########################################################################################################################
## PathSeq Pipeline WDL
########################################################################################################################
##
## Runs the PathSeq pipeline
##
## For further info see the GATK Documentation for the PathSeqPipelineSpark tool:
##   https://software.broadinstitute.org/gatk/documentation/tooldocs/current/org_broadinstitute_hellbender_tools_spark_pathseq_PathSeqPipelineSpark.php
##
########################################################################################################################
##
## Input requirements :
## - Sequencing data in BAM format
## - Host and microbe references files available in the GATK Resource Bundle (available on FTP):
##     https://software.broadinstitute.org/gatk/download/bundle
##
## - Input BAM files must comply with the following requirements:
## - - file must pass validation by ValidateSamFile
## - - all reads must have an RG tag
## - - one or more read groups all belong to a single sample (SM)
##
## Output:
## - BAM file containing microbe-mapped reads and reads of unknown sequence
## - Tab-separated value (.tsv) file of taxonomic abundance scores
## - Picard-style metrics files for the filter and scoring phases of the pipeline
##
########################################################################################################################

task PathSeqFilter {

  # Inputs for this task
  String sample_name
  File input_bam_or_cram

  File kmer_file
  File filter_bwa_image

  # Required if cram is provided
  File? cram_reference_fasta
  File? cram_reference_fasta_index
  File? cram_reference_dict

  Boolean is_host_aligned
  Boolean skip_quality_filters = false
  Boolean skip_pre_bwa_repartition = false
  Boolean filter_duplicates = true
  Boolean gather_metrics = true

  Int filter_bwa_seed_length = 19
  Int host_min_identity = 30
  Int min_clipped_read_length = 60
  Int bam_partition_size = 4000000

  String paired_bam_output_path = "${sample_name}.non_host.paired.bam"
  String unpaired_bam_output_path = "${sample_name}.non_host.unpaired.bam"
  String filter_metrics_output_path = "${sample_name}.pathseq.filter_metrics"

  File? gatk4_jar_override

  # Runtime parameters
  String gatk_docker
  Int mem_gb = 32
  Int preemptible_attempts = 3
  Float additional_disk_gb = 50
  Int cpu = 8
  Boolean use_ssd = false

  Int disk_size = ceil(size(input_bam_or_cram, "GB") + size(kmer_file, "GB") + size(filter_bwa_image, "GB") + additional_disk_gb)

  # Mem is in units of GB but our command and memory runtime values are in MB
  Int machine_mem = mem_gb * 1000
  Int command_mem = ceil((machine_mem - size(filter_bwa_image, "MB")) * 0.8)

  command <<<
    set -euo pipefail
    export GATK_LOCAL_JAR=${default="/root/gatk.jar" gatk4_jar_override}
    touch ${filter_metrics_output_path}
    gatk --java-options "-Xmx${command_mem}m" \
      PathSeqFilterSpark \
      --input ${input_bam_or_cram} \
      --paired-output ${paired_bam_output_path} \
      --unpaired-output ${unpaired_bam_output_path} \
      ~{if defined(cram_reference_fasta) then "--reference ${cram_reference_fasta}" else ""} \
      ~{if gather_metrics then "--filter-metrics ${filter_metrics_output_path}" else ""} \
      --kmer-file ${kmer_file} \
      --filter-bwa-image ${filter_bwa_image} \
      --bam-partition-size ${bam_partition_size} \
      --is-host-aligned ${is_host_aligned} \
      --skip-quality-filters ${skip_quality_filters} \
      --min-clipped-read-length ${min_clipped_read_length} \
      --filter-bwa-seed-length ${filter_bwa_seed_length} \
      --host-min-identity ${host_min_identity} \
      --filter-duplicates ${filter_duplicates} \
      --skip-pre-bwa-repartition ${skip_pre_bwa_repartition}

    if [ ! -f "${paired_bam_output_path}" ]; then
    	echo "File ${paired_bam_output_path} not found, creating empty BAM"
    	printf "@HD     VN:1.5  SO:queryname\n@RG     ID:A    SM:${sample_name}\n" | samtools view -hb - > ${paired_bam_output_path}
    fi

    if [ ! -f "${unpaired_bam_output_path}" ]; then
    	echo "File ${unpaired_bam_output_path} not found, creating empty BAM"
    	printf "@HD     VN:1.5  SO:queryname\n@RG     ID:A    SM:${sample_name}\n" | samtools view -hb - > ${unpaired_bam_output_path}
    fi
  >>>
  runtime {
    docker: gatk_docker
    memory: machine_mem + " MB"
    # Note that the space before SSD and HDD should be included.
    disks: "local-disk " + disk_size + if use_ssd then " SSD" else " HDD"
    preemptible: preemptible_attempts
    cpu: cpu
  }
  output {
    File paired_bam_out = "${paired_bam_output_path}"
    File unpaired_bam_out = "${unpaired_bam_output_path}"
    File filter_metrics = "${sample_name}.pathseq.filter_metrics"
  }
}

task PathSeqAlign {

  # Inputs for this task
  String sample_name
  File input_paired_bam
  File input_unpaired_bam

  File microbe_bwa_image
  File microbe_dict

  String paired_bam_output_path = "${sample_name}.microbe_aligned.paired.bam"
  String unpaired_bam_output_path = "${sample_name}.microbe_aligned.unpaired.bam"

  File? gatk4_jar_override

  # Runtime parameters
  String gatk_docker
  Int mem_gb = 140
  Int preemptible_attempts = 3
  Float additional_disk_gb = 20
  Int cpu = 32
  Boolean use_ssd = false

  Int disk_size = ceil(2*size(input_paired_bam, "GB") + 2*size(input_unpaired_bam, "GB") + size(microbe_bwa_image, "GB") + additional_disk_gb)

  # Mem is in units of GB but our command and memory runtime values are in MB
  Int machine_mem = mem_gb * 1000
  Int command_mem = ceil((machine_mem - size(microbe_bwa_image, "MB")) * 0.8)

  command <<<
    set -e
    export GATK_LOCAL_JAR=${default="/root/gatk.jar" gatk4_jar_override}
    gatk --java-options "-Xmx${command_mem}m" \
      PathSeqBwaSpark \
      --paired-input ${input_paired_bam} \
      --unpaired-input ${input_unpaired_bam} \
      --paired-output ${paired_bam_output_path} \
      --unpaired-output ${unpaired_bam_output_path} \
      --microbe-bwa-image ${microbe_bwa_image} \
      --microbe-dict ${microbe_dict}
  >>>
  runtime {
    docker: gatk_docker
    memory: machine_mem + " MB"
    # Note that the space before SSD and HDD should be included.
    disks: "local-disk " + disk_size + if use_ssd then " SSD" else " HDD"
    preemptible: preemptible_attempts
    cpu: cpu
  }
  output {
    File paired_bam_out = "${paired_bam_output_path}"
    File unpaired_bam_out = "${unpaired_bam_output_path}"
  }
}

task PathSeqScore {

  # Inputs for this task
  String sample_name
  File input_paired_bam
  File input_unpaired_bam

  File taxonomy_file

  Boolean divide_by_genome_length = true
  Float min_score_identity = 0.9
  Float identity_margin = 0.02

  String bam_output_path = "${sample_name}.pathseq.bam"
  String scores_output_path = "${sample_name}.pathseq.tsv"
  String score_metrics_output_path = "${sample_name}.pathseq.score_metrics"

  File? gatk4_jar_override

  # Runtime parameters
  String gatk_docker
  Int mem_gb = 8
  Int preemptible_attempts = 3
  Float additional_disk_gb = 10
  Int cpu = 2
  Boolean use_ssd = false

  Int disk_size = ceil(size(input_paired_bam, "GB") + size(input_unpaired_bam, "GB") + additional_disk_gb)

  # Mem is in units of GB but our command and memory runtime values are in MB
  Int machine_mem = mem_gb * 1000
  Int command_mem = ceil(machine_mem * 0.8)

  command <<<
    set -e
    export GATK_LOCAL_JAR=${default="/root/gatk.jar" gatk4_jar_override}
    gatk --java-options "-Xmx${command_mem}m" \
      PathSeqScoreSpark \
      --paired-input ${input_paired_bam} \
      --unpaired-input ${input_unpaired_bam} \
      --output ${bam_output_path} \
      --scores-output ${scores_output_path} \
      --score-metrics ${score_metrics_output_path} \
      --taxonomy-file ${taxonomy_file} \
      --min-score-identity ${min_score_identity} \
      --identity-margin ${identity_margin} \
      --divide-by-genome-length ${divide_by_genome_length}
  >>>
  runtime {
    docker: gatk_docker
    memory: machine_mem + " MB"
    # Note that the space before SSD and HDD should be included.
    disks: "local-disk " + disk_size + if use_ssd then " SSD" else " HDD"
    preemptible: preemptible_attempts
    cpu: cpu
  }
  output {
    File bam_out = "${sample_name}.pathseq.bam"
    File scores = "${sample_name}.pathseq.tsv"
    File score_metrics = "${sample_name}.pathseq.score_metrics"
  }
}

workflow PathSeqThreeStageWorkflow {

  String sample_name
  File input_bam_or_cram

  File? cram_reference_fasta
  File? cram_reference_fasta_index
  File? cram_reference_dict

  File kmer_file
  File filter_bwa_image
  File microbe_bwa_image
  File microbe_dict
  File taxonomy_file

  Boolean is_host_aligned
  Boolean? filter_duplicates
  Boolean? skip_pre_bwa_repartition
  Int? filter_bam_partition_size
  Boolean? divide_by_genome_length
  Int? filter_bwa_seed_length
  Int? host_min_identity
  Int? min_clipped_read_length
  Float? min_score_identity
  Float? identity_margin

  File? gatk4_jar_override

  # Runtime parameters
  String gatk_docker

  Int? filter_preemptible_attempts
  Int? align_preemptible_attempts
  Int? score_preemptible_attempts

  Int? filter_cpu
  Int? align_cpu
  Int? score_cpu

  Int? filter_mem_gb
  Int? align_mem_gb
  Int? score_mem_gb

  # Optional input to increase all disk sizes in case of outlier sample with strange size behavior
  Int? filter_additional_disk_gb
  Int? align_additional_disk_gb
  Int? score_additional_disk_gb

  call PathSeqFilter {
    input:
      sample_name=sample_name,
      input_bam_or_cram=input_bam_or_cram,
      cram_reference_fasta = cram_reference_fasta,
      cram_reference_fasta_index = cram_reference_fasta_index,
      cram_reference_dict = cram_reference_dict,
      kmer_file=kmer_file,
      filter_bwa_image=filter_bwa_image,
      is_host_aligned=is_host_aligned,
      filter_duplicates=filter_duplicates,
      min_clipped_read_length=min_clipped_read_length,
      bam_partition_size=filter_bam_partition_size,
      host_min_identity=host_min_identity,
      filter_bwa_seed_length=filter_bwa_seed_length,
      skip_pre_bwa_repartition=skip_pre_bwa_repartition,
      gatk4_jar_override=gatk4_jar_override,
      mem_gb=filter_mem_gb,
      gatk_docker=gatk_docker,
      preemptible_attempts=filter_preemptible_attempts,
      additional_disk_gb=filter_additional_disk_gb,
      cpu=filter_cpu
  }

  call PathSeqAlign {
    input:
      sample_name=sample_name,
      input_paired_bam=PathSeqFilter.paired_bam_out,
      input_unpaired_bam=PathSeqFilter.unpaired_bam_out,
      microbe_bwa_image=microbe_bwa_image,
      microbe_dict=microbe_dict,
      gatk4_jar_override=gatk4_jar_override,
      mem_gb=align_mem_gb,
      gatk_docker=gatk_docker,
      preemptible_attempts=align_preemptible_attempts,
      additional_disk_gb=align_additional_disk_gb,
      cpu=align_cpu
  }

  call PathSeqScore {
    input:
      sample_name=sample_name,
      input_paired_bam=PathSeqAlign.paired_bam_out,
      input_unpaired_bam=PathSeqAlign.unpaired_bam_out,
      taxonomy_file=taxonomy_file,
      divide_by_genome_length=divide_by_genome_length,
      min_score_identity=min_score_identity,
      identity_margin=identity_margin,
      gatk4_jar_override=gatk4_jar_override,
      mem_gb=score_mem_gb,
      gatk_docker=gatk_docker,
      preemptible_attempts=score_preemptible_attempts,
      additional_disk_gb=score_additional_disk_gb,
      cpu=score_cpu
  }

  output {
    File final_bam = PathSeqScore.bam_out
    File scores = PathSeqScore.scores
    File filter_metrics = PathSeqFilter.filter_metrics
    File score_metrics = PathSeqScore.score_metrics
  }
}