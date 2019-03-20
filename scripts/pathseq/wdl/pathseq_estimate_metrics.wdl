########################################################################################################################
## PathSeq Estimate Filter Metrics WDL
########################################################################################################################
##
## Estimates filter metrics for PathSeq
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
## - Picard-style metrics files for the filter phase of the pipeline at the downsampled coverage
##
########################################################################################################################

import "pathseq_3_stage_pipeline.wdl" as pathseq_pipeline

# Downsamples BAM to a specified number of reads
task Downsample {
  File input_bam_file
  File input_bam_index_file
  String downsampled_bam_filename

  Int reads_after_downsampling

  Int disk_size_gb
  Int preemptible_tries

  command <<<
    set -e
    set -o pipefail
    NUM_READS=`samtools idxstats ${input_bam_file} | awk '{s+=$3+$4} END {print s}'`
    P_DOWNSAMPLE=`python -c "print ${reads_after_downsampling}/float($NUM_READS)"`
    RESULT=`python -c "print $P_DOWNSAMPLE > 1"`
    if [ "$RESULT" == "True" ]
    then
        P_DOWNSAMPLE="1"
    fi
    java -Xmx2000m -jar /usr/gitc/picard.jar \
      DownsampleSam \
      INPUT=${input_bam_file} \
      OUTPUT=${downsampled_bam_filename} \
      P=$P_DOWNSAMPLE
  >>>
  output {
    File output_bam_file = "${downsampled_bam_filename}"
  }
  runtime {
    preemptible: "${preemptible_tries}"
    docker: "us.gcr.io/broad-gotc-prod/genomes-in-the-cloud:2.2.5-1486412288"
    memory: "3.75 GB"
    cpu: "1"
    disks: "local-disk ${disk_size_gb} HDD"
  }
}

workflow EstimateMetricsPathseq {

  String sample_name
  File input_bam
  File input_bam_index
  Int? reads_after_downsampling

  File kmer_file
  File filter_bwa_image

  Boolean is_host_aligned
  Boolean? filter_duplicates
  Boolean? skip_pre_bwa_repartition
  Int? filter_bwa_seed_length
  Int? host_min_identity
  Int? min_clipped_read_length

  File? gatk4_jar_override

  # Runtime parameters
  String gatk_docker
  Int? preemptible_attempts
  Int? cpu
  Int? filter_bam_partition_size
  Int? filter_mem_gb

  # Optional input to increase all disk sizes in case of outlier sample with strange size behavior
  Int? increase_disk_size

  # Some tasks need wiggle room, and we also need to add a small amount of disk to prevent getting a
  # Cromwell error from asking for 0 disk when the input is less than 1GB.
  # Also Spark requires some temporary storage.
  Int additional_disk = select_first([increase_disk_size, 20])
  Float downsample_disk_space_gb = size(input_bam, "GB") + additional_disk

  # Downsample input bam
  call Downsample {
    input:
      input_bam_file=input_bam,
      input_bam_index_file=input_bam_index,
      downsampled_bam_filename="${sample_name}.downsampled.bam",
      reads_after_downsampling=select_first([reads_after_downsampling, 100000]),
      disk_size_gb=downsample_disk_space_gb,
      preemptible_tries=preemptible_attempts
  }

  Float filter_disk_space_gb = 2*size(Downsample.output_bam_file, "GB") + size(kmer_file, "GB") + size(filter_bwa_image, "GB") + additional_disk
  call pathseq_pipeline.PathSeqThreeStageWorkflow {
    input:
      sample_name=sample_name,
      input_bam=Downsample.output_bam_file,
      kmer_file=kmer_file,
      filter_bwa_image=filter_bwa_image,
      is_host_aligned=is_host_aligned,
      filter_duplicates=filter_duplicates,
      min_clipped_read_length=min_clipped_read_length,
      bam_partition_size=filter_bam_partition_size,
      host_min_identity=host_min_identity,
      filter_bwa_seed_length=filter_bwa_seed_length,
      gatk4_jar_override=gatk4_jar_override,
      mem_gb=filter_mem_gb,
      gatk_docker=gatk_docker,
      preemptible_attempts=preemptible_attempts,
      disk_space_gb=filter_disk_space_gb,
      cpu=cpu
  }

  output {
    File metrics = PathseqFilter.outputFilterMetricsFile
  }
}