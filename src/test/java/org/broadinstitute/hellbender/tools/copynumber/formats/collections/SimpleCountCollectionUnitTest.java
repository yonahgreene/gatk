package org.broadinstitute.hellbender.tools.copynumber.formats.collections;

import com.google.common.collect.Lists;
import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.samtools.SAMTextHeaderCodec;
import htsjdk.samtools.util.BufferedLineReader;
import htsjdk.samtools.util.Log;
import org.apache.commons.math3.linear.Array2DRowRealMatrix;
import org.apache.commons.math3.linear.RealMatrix;
import org.broadinstitute.hellbender.GATKBaseTest;
import org.broadinstitute.hellbender.engine.FeatureDataSource;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.tools.copynumber.formats.metadata.SampleLocatableMetadata;
import org.broadinstitute.hellbender.tools.copynumber.formats.metadata.SimpleSampleLocatableMetadata;
import org.broadinstitute.hellbender.tools.copynumber.formats.records.SimpleCount;
import org.broadinstitute.hellbender.utils.LoggingUtils;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.gcs.BucketUtils;
import org.testng.Assert;
import org.testng.annotations.Test;

import java.io.File;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashSet;
import java.util.List;
import java.util.stream.Collectors;
import java.util.stream.Stream;


public final class SimpleCountCollectionUnitTest extends GATKBaseTest {
    private static final File TEST_SUB_DIR = new File(toolsTestDir, "copynumber/formats/collections");
    private static final File INTEGER_COUNTS_FILE = new File(TEST_SUB_DIR,"simple-count-collection-integer-counts.tsv");
    private static final File INTEGER_COUNTS_HDF5_FILE = new File(TEST_SUB_DIR,"simple-count-collection-integer-counts.hdf5");
    private static final File INTEGER_COUNTS_MISSING_HEADER_FILE = new File(TEST_SUB_DIR,"simple-count-collection-integer-counts-missing-header.tsv");
    private static final File DOUBLE_COUNTS_FILE = new File(TEST_SUB_DIR, "simple-count-collection-double-counts.tsv");

    private static final SampleLocatableMetadata METADATA_EXPECTED = new SimpleSampleLocatableMetadata(
            "test-sample",
            new SAMSequenceDictionary(Collections.singletonList(
                    new SAMSequenceRecord("20", 200000))));

    private static final List<SimpleInterval> INTERVALS_EXPECTED = Arrays.asList(
            new SimpleInterval("20", 1,10000),
            new SimpleInterval("20", 10001,20000),
            new SimpleInterval("20", 20001, 30000),
            new SimpleInterval("20", 30001, 40000),
            new SimpleInterval("20", 40001, 50000),
            new SimpleInterval("20", 50001, 60000),
            new SimpleInterval("20", 60001, 70000),
            new SimpleInterval("20", 70001, 80000),
            new SimpleInterval("20", 80001, 90000),
            new SimpleInterval("20", 90001, 100000),
            new SimpleInterval("20", 100001, 110000),
            new SimpleInterval("20", 110001, 120000),
            new SimpleInterval("20", 120001, 130000),
            new SimpleInterval("20", 130001, 140000),
            new SimpleInterval("20", 140001, 150000),
            new SimpleInterval("20", 150001, 160000));
    private static final RealMatrix READ_COUNTS_EXPECTED = new Array2DRowRealMatrix(
            new double[][]{{0, 0, 0, 0, 0, 0, 94, 210, 22, 21, 24, 84, 247, 181, 27, 72}});

    @Test
    public void testReadIntegerCounts() {
        final SimpleCountCollection scc = SimpleCountCollection.read(INTEGER_COUNTS_FILE);
        final SampleLocatableMetadata metadata = scc.getMetadata();
        final List<SimpleInterval> intervals = scc.getIntervals();
        final RealMatrix readCounts = new Array2DRowRealMatrix(new double[][]{scc.getRecords().stream().mapToDouble(SimpleCount::getCount).toArray()});

        Assert.assertEquals(metadata, METADATA_EXPECTED);
        Assert.assertEquals(intervals, INTERVALS_EXPECTED);
        Assert.assertEquals(readCounts, READ_COUNTS_EXPECTED);
    }

    @Test
    public void testReadIntegerCountsHDF5() {
        final SimpleCountCollection scc = SimpleCountCollection.read(INTEGER_COUNTS_HDF5_FILE);
        final SampleLocatableMetadata metadata = scc.getMetadata();
        final List<SimpleInterval> intervals = scc.getIntervals();
        final RealMatrix readCounts = new Array2DRowRealMatrix(new double[][]{scc.getRecords().stream().mapToDouble(SimpleCount::getCount).toArray()});

        Assert.assertEquals(metadata, METADATA_EXPECTED);
        Assert.assertEquals(intervals, INTERVALS_EXPECTED);
        Assert.assertEquals(readCounts, READ_COUNTS_EXPECTED);
    }

    @Test(expectedExceptions = IllegalArgumentException.class)
    public void testReadIntegerCountsMissingHeader() {
        SimpleCountCollection.read(INTEGER_COUNTS_MISSING_HEADER_FILE);
    }

    @Test(expectedExceptions = UserException.BadInput.class)
    public void testReadDoubleCounts() {
        SimpleCountCollection.read(DOUBLE_COUNTS_FILE);
    }

    @Test
    public void testQuery() {
        LoggingUtils.setLoggingLevel(Log.LogLevel.DEBUG);
        final String localPath = "/home/slee/working/gatk/test_files/test.counts.tsv";
        final String bucketPath = "gs://broad-dsde-methods-slee/test.counts.tsv";

        final File file = new File(localPath);
        final FeatureDataSource<SimpleCount> localSource = new FeatureDataSource<>(localPath);
        final FeatureDataSource<SimpleCount> bucketSource = new FeatureDataSource<>(bucketPath);

        final SimpleInterval interval = new SimpleInterval("1", 1, 500000);

        final SimpleCountCollection counts = SimpleCountCollection.read(file);

        final BufferedLineReader localReader = new BufferedLineReader(BucketUtils.openFile(localPath));
        final SAMFileHeader localHeader = new SAMTextHeaderCodec().decode(localReader, localPath);

        final BufferedLineReader bucketReader = new BufferedLineReader(BucketUtils.openFile(bucketPath));
        final SAMFileHeader bucketHeader = new SAMTextHeaderCodec().decode(bucketReader, bucketPath);

        System.out.println(Stream.of(counts.getMetadata().toHeader().getSAMString(), localHeader.getSAMString(), bucketHeader.getSAMString()).distinct().count() == 1);

        System.out.println(localSource.getHeader().toString());
        System.out.println(bucketSource.getHeader().toString());

        System.out.println(counts.getOverlapDetector().getOverlaps(interval).stream().sorted(counts.getComparator()).collect(Collectors.toList()));
        System.out.println(Lists.newArrayList(localSource.query(interval)));
        System.out.println(Lists.newArrayList(bucketSource.query(interval)));
    }
}