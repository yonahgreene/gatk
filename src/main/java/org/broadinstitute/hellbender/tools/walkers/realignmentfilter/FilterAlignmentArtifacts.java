package org.broadinstitute.hellbender.tools.walkers.realignmentfilter;

import htsjdk.samtools.Cigar;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFHeaderLine;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.ArgumentCollection;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.barclay.argparser.ExperimentalFeature;
import org.broadinstitute.barclay.help.DocumentedFeature;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.engine.*;
import org.broadinstitute.hellbender.engine.filters.CountingVariantFilter;
import org.broadinstitute.hellbender.engine.filters.ReadFilter;
import org.broadinstitute.hellbender.engine.filters.VariantFilterLibrary;
import org.broadinstitute.hellbender.tools.walkers.haplotypecaller.graphs.*;
import org.broadinstitute.hellbender.tools.walkers.haplotypecaller.readthreading.MultiDeBruijnVertex;
import org.broadinstitute.hellbender.tools.walkers.haplotypecaller.readthreading.ReadThreadingGraph;
import org.broadinstitute.hellbender.tools.walkers.mutect.Mutect2Engine;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.bwa.BwaMemAlignment;
import org.broadinstitute.hellbender.utils.haplotype.EventMap;
import org.broadinstitute.hellbender.utils.haplotype.Haplotype;
import org.broadinstitute.hellbender.utils.read.CigarUtils;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.smithwaterman.SmithWatermanAligner;
import org.broadinstitute.hellbender.utils.smithwaterman.SmithWatermanJavaAligner;
import org.broadinstitute.hellbender.utils.variant.GATKVCFConstants;
import org.broadinstitute.hellbender.utils.variant.GATKVCFHeaderLines;
import picard.cmdline.programgroups.VariantFilteringProgramGroup;

import java.io.File;
import java.util.*;
import java.util.stream.Collectors;

/**
 * <p>Filter false positive alignment artifacts from a VCF callset.</p>
 *
 * <p>
 *     FilterAlignmentArtifacts identifies alignment artifacts, that is, apparent variants due to reads being mapped to the wrong genomic locus.
 * </p>
 * <p>
 *     Alignment artifacts can occur whenever there is sufficient sequence similarity between two or more regions in the genome
 *     to confuse the alignment algorithm.  This can occur when the aligner for whatever reason overestimate how uniquely a read
 *     maps, thereby assigning it too high of a mapping quality.  It can also occur through no fault of the aligner due to gaps in
 *     the reference, which can also hide the true position to which a read should map.  By using a good alignment algorithm
 *     (the GATK wrapper of BWA-MEM), giving it sensitive settings (which may have been impractically slow for the original
 *     bam alignment) and mapping to the best available reference we can avoid these pitfalls.  The last point is especially important:
 *     one can (and should) use a BWA-MEM index image corresponding to the best reference, regardless of the reference to which
 *     the bam was aligned.
 * </p>
 * <p>
 *     This tool is featured in the Somatic Short Mutation calling Best Practice Workflow.
 *     See <a href="https://software.broadinstitute.org/gatk/documentation/article?id=11136">Tutorial#11136</a> for a
 *     step-by-step description of the workflow and <a href="https://software.broadinstitute.org/gatk/documentation/article?id=11127">Article#11127</a>
 *     for an overview of what traditional somatic calling entails. For the latest pipeline scripts, see the
 *     <a href="https://github.com/broadinstitute/gatk/tree/master/scripts/mutect2_wdl">Mutect2 WDL scripts directory</a>.
 * </p>
 * <p>
 *     The bam input to this tool should be the reassembly bamout produced by HaplotypeCaller or Mutect2 in the process of generating
 *     the input callset.  The original bam will also work but might fail to filter some indels.  The reference passed with the -R argument
 *     must be the reference to which the input bam was realigned.  This does not need to correspond to the reference of the BWA-MEM
 *     index image.  The latter should be derived from the best available reference, for example hg38 in humans as of February 2018.
 * </p>
 *
 * <h3>Usage example</h3>
 * <pre>
 * gatk FilterAlignmentArtifacts \
 *   -R hg19.fasta
 *   -V somatic.vcf.gz \
 *   -I somatic_bamout.bam \
 *   --bwa-mem-index-image hg38.index_image \
 *   -O filtered.vcf.gz
 * </pre>
 *
 */
@CommandLineProgramProperties(
        summary = "Filter alignment artifacts from a vcf callset.",
        oneLineSummary = "Filter alignment artifacts from a vcf callset.",
        programGroup = VariantFilteringProgramGroup.class
)
@DocumentedFeature
@ExperimentalFeature
public class FilterAlignmentArtifacts extends MultiVariantWalkerGroupedOnStart {
    public static final int DEFAULT_DISTANCE_TO_GROUP_VARIANTS = 1000;
    public static final int DEFAULT_REF_PADDING = 100;
    public static final int DEFAULT_MAX_GROUPED_SPAN = 10_000;
    private static final int MIN_UNITIG_LENGTH = 30;
    private static final SmithWatermanAligner ALIGNER = SmithWatermanAligner.getAligner(SmithWatermanAligner.Implementation.FASTEST_AVAILABLE);

    @Argument(fullName = StandardArgumentDefinitions.OUTPUT_LONG_NAME, shortName = StandardArgumentDefinitions.OUTPUT_SHORT_NAME,
            doc="The output filtered VCF file", optional=false)
    private final String outputVcf = null;

    public static final int DEFAULT_INDEL_START_TOLERANCE = 5;
    public static final String INDEL_START_TOLERANCE_LONG_NAME = "indel-start-tolerance";
    @Argument(fullName = INDEL_START_TOLERANCE_LONG_NAME, doc="Max distance between indel start of aligned read in the bam and the variant in the vcf", optional=true)
    private int indelStartTolerance = DEFAULT_INDEL_START_TOLERANCE;

    public static final int DEFAULT_KMER_SIZE = 21;
    public static final String KMER_SIZE_LONG_NAME = "kmer-size";
    @Argument(fullName = KMER_SIZE_LONG_NAME, doc="Kmer size for reassembly", optional=true)
    private int kmerSize = DEFAULT_KMER_SIZE;

    public static final String DONT_SKIP_ALREADY_FILTERED_VARIANTS_LONG_NAME = "dont-skip-filtered-variants";
    @Argument(fullName = DONT_SKIP_ALREADY_FILTERED_VARIANTS_LONG_NAME,
            doc="Try to realign all variants, even ones that have already been filtered.", optional=true)
    private boolean dontSkipFilteredVariants = false;


    @ArgumentCollection
    protected RealignmentArgumentCollection realignmentArgumentCollection = new RealignmentArgumentCollection();

    private VariantContextWriter vcfWriter;
    private RealignmentEngine realignmentEngine;

    @Override
    public List<ReadFilter> getDefaultReadFilters() { return Mutect2Engine.makeStandardMutect2ReadFilters(); }

    @Override
    protected CountingVariantFilter makeVariantFilter() {
        return new CountingVariantFilter(dontSkipFilteredVariants ? VariantFilterLibrary.ALLOW_ALL_VARIANTS : VariantFilterLibrary.PASSES_FILTERS);
    }

    @Override
    public boolean requiresReads() { return true; }

    @Override
    protected int defaultDistanceToGroupVariants() { return DEFAULT_DISTANCE_TO_GROUP_VARIANTS; }

    @Override
    protected int defaultReferenceWindowPadding() { return DEFAULT_REF_PADDING; }

    @Override
    protected int defaultMaxGroupedSpan() {
        return DEFAULT_MAX_GROUPED_SPAN;
    }

    @Override
    public void onTraversalStart() {
        realignmentEngine = new RealignmentEngine(realignmentArgumentCollection);
        vcfWriter = createVCFWriter(new File(outputVcf));

        final VCFHeader inputHeader = getHeaderForVariants();
        final Set<VCFHeaderLine> headerLines = new HashSet<>(inputHeader.getMetaDataInSortedOrder());
        headerLines.add(GATKVCFHeaderLines.getFilterLine(GATKVCFConstants.ALIGNMENT_ARTIFACT_FILTER_NAME));
        headerLines.add(GATKVCFHeaderLines.getInfoLine(GATKVCFConstants.UNITIG_SIZES_KEY));
        headerLines.add(GATKVCFHeaderLines.getInfoLine(GATKVCFConstants.ALIGNMENT_SCORE_DIFFERENCE_KEY));
        headerLines.add(GATKVCFHeaderLines.getInfoLine(GATKVCFConstants.JOINT_ALIGNMENT_COUNT_KEY));
        headerLines.addAll(getDefaultToolVCFHeaderLines());
        final VCFHeader vcfHeader = new VCFHeader(headerLines, inputHeader.getGenotypeSamples());
        vcfWriter.writeHeader(vcfHeader);
    }

    @Override
    public Object onTraversalSuccess() {
        return "SUCCESS";
    }

    @Override
    public void apply(List<VariantContext> variantContexts, ReferenceContext referenceContext, final List<ReadsContext> readsContexts) {


        // for now we do one variant at a time but eventually we will want to combine all reads supporting all variants
        // into a single graph.  This is non-trivial because there may be more than one phasing between variants.
        for (final VariantContext vc : variantContexts) {
            final SeqGraph seqGraph = makeSeqGraph(referenceContext, readsContexts, vc);
            final Optional<KBestHaplotype> bestVariantHaplotype = getBestHaplotype(seqGraph, vc, referenceContext);

            if (!bestVariantHaplotype.isPresent()) {
                vcfWriter.add(vc);
                return;
            }

            final List<byte[]> unitigs = getUnitigs(seqGraph, bestVariantHaplotype.get());

            final VariantContextBuilder vcb = new VariantContextBuilder(vc)
                    .attribute(GATKVCFConstants.UNITIG_SIZES_KEY, unitigs.stream().mapToInt(u -> u.length).toArray());

            final List<List<BwaMemAlignment>> unitigAlignments = unitigs.stream()
                    .map(realignmentEngine::realign).collect(Collectors.toList());

            final List<List<BwaMemAlignment>> jointAlignments = RealignmentEngine.findJointAlignments(unitigAlignments, realignmentArgumentCollection.maxReasonableFragmentLength);
            vcb.attribute(GATKVCFConstants.JOINT_ALIGNMENT_COUNT_KEY, jointAlignments.size());
            if (jointAlignments.size() > 1) {
                jointAlignments.sort(Comparator.comparingInt(FilterAlignmentArtifacts::jointAlignmentScore).reversed());
                final int totalBases = unitigs.stream().mapToInt(unitig -> unitig.length).sum();
                final int scoreDiff = jointAlignmentScore(jointAlignments.get(0)) - jointAlignmentScore(jointAlignments.get(1));
                final int mismatchDiff = totalMismatches(jointAlignments.get(1)) - totalMismatches(jointAlignments.get(0));

                vcb.attribute(GATKVCFConstants.ALIGNMENT_SCORE_DIFFERENCE_KEY, scoreDiff);

                final boolean multimapping = (double) scoreDiff / totalBases < realignmentArgumentCollection.minAlignerScoreDifferencePerBase
                        && (double) mismatchDiff / totalBases < realignmentArgumentCollection.minMismatchDifferencePerBase;

                if (multimapping) {
                    vcb.filter(GATKVCFConstants.ALIGNMENT_ARTIFACT_FILTER_NAME);
                }
            }

            vcfWriter.add(vcb.make());
        }
    }


    private SeqGraph makeSeqGraph(final ReferenceContext referenceContext, final List<ReadsContext> readsContexts, VariantContext vc) {
        final Set<String> variantReadNames = readsContexts.stream().flatMap(Utils::stream)
                .filter(read -> RealignmentEngine.supportsVariant(read, vc, indelStartTolerance))
                .map(GATKRead::getName)
                .collect(Collectors.toSet());

        final List<GATKRead> variantReads = readsContexts.stream().flatMap(Utils::stream)
                .filter(read -> variantReadNames.contains(read.getName()))
                .collect(Collectors.toList());

        final ReadThreadingGraph graph = new ReadThreadingGraph(kmerSize);

        // note: add the reference ot make graph construction easier, but we ultimately discard any assembled sequence
        // supported only by the reference
        graph.addSequence("reference", referenceContext.getBases(), true);
        variantReads.forEach(read -> graph.addSequence(read.getBasesNoCopy(), false));
        graph.buildGraphIfNecessary();

        final ChainPruner<MultiDeBruijnVertex, MultiSampleEdge> pruner = new AdaptiveChainPruner<>(0.001, 1, 5);
        pruner.pruneLowWeightChains(graph);

        final SmithWatermanAligner aligner = SmithWatermanJavaAligner.getInstance();
        graph.recoverDanglingTails(1, 3, false, aligner);
        graph.recoverDanglingHeads(1, 3, false, aligner);
        graph.removePathsNotConnectedToRef();

        return graph.toSequenceGraph();
    }

    private static Optional<KBestHaplotype> getBestHaplotype(SeqGraph seqGraph, final VariantContext vc, final ReferenceContext referenceContext) {
        //seqGraph.zipLinearChains();
        seqGraph.removeSingletonOrphanVertices();
        seqGraph.removeVerticesNotConnectedToRefRegardlessOfEdgeDirection();
        //seqGraph.simplifyGraph();
        seqGraph.removePathsNotConnectedToRef();
        //seqGraph.simplifyGraph();

        final List<KBestHaplotype> bestPaths = new KBestHaplotypeFinder(seqGraph).findBestHaplotypes(10);

        // find best haplotype that supports the variant
        return bestPaths.stream()
                .filter(hap -> {
                    final Haplotype haplotype = hap.haplotype();
                    final Cigar cigar = CigarUtils.calculateCigar(referenceContext.getBases(), haplotype.getBases(), ALIGNER);
                    haplotype.setCigar(cigar);
                    final EventMap eventMap = new EventMap(haplotype, referenceContext.getBases(), referenceContext.getWindow(), "name", 1);
                    haplotype.setEventMap(eventMap);
                    return haplotype.getEventMap().getStartPositions().stream().anyMatch(pos -> Math.abs(pos - vc.getStart()) <= vc.getLengthOnReference());
                })
                .sorted(Comparator.comparingDouble(KBestHaplotype::score).reversed())
                .findFirst();
    }

    private List<byte[]> getUnitigs(final SeqGraph seqGraph, final KBestHaplotype bestVariantHaplotype) {
        final List<StringBuilder> unitigBuilders = new ArrayList<>();
        boolean lastEdgeInUnitig = false;

        for (final BaseEdge edge : bestVariantHaplotype.getEdges()) {
            // if it has support from variant reads, add to current unitig
            final boolean currentEdgeInUnitig = !edge.isRef() || edge.getMultiplicity() > 1;
            if (currentEdgeInUnitig) {
                if (!lastEdgeInUnitig) {
                    unitigBuilders.add(new StringBuilder());
                    unitigBuilders.get(unitigBuilders.size() - 1).append(seqGraph.getEdgeSource(edge).getSequence());
                }
                unitigBuilders.get(unitigBuilders.size() - 1).append(new String(seqGraph.getEdgeTarget(edge).getAdditionalSequence(false)));
            }
            lastEdgeInUnitig = currentEdgeInUnitig;
        }

        return unitigBuilders.stream()
                .map(builder -> builder.toString().getBytes())
                .filter(unitig -> unitig.length > MIN_UNITIG_LENGTH)
                .collect(Collectors.toList());
    }

    private static int jointAlignmentScore(final List<BwaMemAlignment> alignments) {
        return alignments.stream().mapToInt(BwaMemAlignment::getAlignerScore).sum();
    }

    private static int totalMismatches(final List<BwaMemAlignment> alignments) {
        return alignments.stream().mapToInt(BwaMemAlignment::getNMismatches).sum();
    }

    @Override
    public void closeTool() {
        if ( vcfWriter != null ) {
            vcfWriter.close();
        }
    }

}
