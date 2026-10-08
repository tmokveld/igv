package org.broad.igv.sam.mods;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.samtools.util.CloseableIterator;
import org.broad.igv.prefs.Constants;
import org.broad.igv.prefs.PreferencesManager;
import org.broad.igv.sam.Alignment;
import org.broad.igv.sam.SAMAlignment;
import org.broad.igv.sam.reader.BAMReader;
import org.broad.igv.util.ResourceLocator;
import org.broad.igv.util.TestUtils;
import org.junit.After;
import org.junit.Before;
import org.junit.Test;

import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashMap;
import java.util.Map;

import static junit.framework.Assert.assertEquals;
import static junit.framework.Assert.assertTrue;

public class BaseModificationCountsTest {
    private String previousValidateBaseCount;


    @Before
    public void setup() {
        previousValidateBaseCount = PreferencesManager.getPreferences().get(Constants.BASEMOD_VALIDATE_BASE_COUNT, null);
        PreferencesManager.getPreferences().put(Constants.BASEMOD_VALIDATE_BASE_COUNT, "false");
    }

    @After
    public void restorePreferences() {
        restorePreference(Constants.BASEMOD_VALIDATE_BASE_COUNT, previousValidateBaseCount);
    }

    private void restorePreference(String key, String value) {
        if (value == null) {
            PreferencesManager.getPreferences().remove(key);
        } else {
            PreferencesManager.getPreferences().put(key, value);
        }
    }

    private SAMAlignment alignment(String bases, String cigar, String mm, byte... likelihoods) {
        SAMFileHeader header = new SAMFileHeader();
        header.addSequence(new SAMSequenceRecord("chr1", 1000));
        SAMRecord record = new SAMRecord(header);
        record.setReadName("modifications");
        record.setReferenceName("chr1");
        record.setAlignmentStart(101);
        record.setCigarString(cigar);
        record.setReadString(bases);
        record.setAttribute("MM", mm);
        record.setAttribute("ML", likelihoods);
        return new SAMAlignment(record);
    }

    @Test
    public void zeroAbsentAndUnsignedCallsGiveDifferentCoverage() {
        SAMFileHeader header = new SAMFileHeader();
        header.addSequence(new SAMSequenceRecord("chr1", 1000));
        SAMRecord record = new SAMRecord(header);
        record.setReadName("modifications");
        record.setReferenceName("chr1");
        record.setAlignmentStart(101);
        record.setCigarString("3M");
        record.setReadString("CCC");
        record.setAttribute("MM", "C+mh?,0,1;");
        record.setAttribute("ML", new byte[]{0, 0, (byte) 255, 0});

        BaseModificationCounts counts = new BaseModificationCounts();
        counts.incrementCounts(new SAMAlignment(record));
        BaseModificationKey modified = BaseModificationKey.getKey('C', '+', "m");
        BaseModificationKey unmodified = BaseModificationKey.getKey('C', '+', "NONE_C");
        assertEquals(1, counts.getCount(100, modified, 0, false));
        assertEquals(0, counts.getCount(100, modified, 0.01f, false));
        assertEquals(1, counts.getCount(100, unmodified, 1, true));
        assertEquals(0, counts.getCount(101, modified, 0, false));
        assertEquals(0, counts.getCount(101, unmodified, 0, true));
        assertEquals(1, counts.getCount(102, modified, 1, false));
        assertEquals(255, counts.getLikelihoodSum(102, modified, 1, false));
        assertEquals(1, counts.getCount(102, modified, 1, true));
        assertEquals(0, counts.getCount(102, unmodified, 0, true));
    }

    @Test
    public void keysFollowFirstObservedCallsAndTheirCanonicalBases() {
        BaseModificationCounts counts = new BaseModificationCounts();
        counts.incrementCounts(alignment("CCAC", "4M", "C+m?,2;A+a?,0;C-h?,0;",
                (byte) 190, (byte) 40, (byte) 80));

        BaseModificationKey h = BaseModificationKey.getKey('C', '-', "h");
        BaseModificationKey m = BaseModificationKey.getKey('C', '+', "m");
        BaseModificationKey a = BaseModificationKey.getKey('A', '+', "a");
        BaseModificationKey noneG = BaseModificationKey.getKey('G', '+', "NONE_G");
        BaseModificationKey noneA = BaseModificationKey.getKey('A', '+', "NONE_A");
        BaseModificationKey noneC = BaseModificationKey.getKey('C', '+', "NONE_C");
        assertEquals(Arrays.asList(h, noneG, a, noneA, m, noneC),
                new ArrayList<>(counts.getAllModificationKeys()));
        assertEquals(80, counts.getLikelihoodSum(100, h, 0, false));
        assertEquals(175, counts.getLikelihoodSum(100, noneG, 0.5f, true));
        assertEquals(0, counts.getCount(101, h, 0, false));
        assertEquals(0, counts.getCount(101, noneG, 0, true));
        assertEquals(40, counts.getLikelihoodSum(102, a, 0, false));
        assertEquals(215, counts.getLikelihoodSum(102, noneA, 0.5f, true));
        assertEquals(190, counts.getLikelihoodSum(103, m, 0.5f, true));
        assertEquals(0, counts.getCount(103, noneC, 0, true));
    }

    @Test
    public void collidingModificationNamesKeepTheFirstObservedBaseAndStrand() {
        Map<String, BaseModificationKey> previousKeyCache = BaseModificationKey.keyCache;
        BaseModificationKey.keyCache = Collections.synchronizedMap(new HashMap<>());
        try {
            BaseModificationCounts counts = new BaseModificationCounts();
            // C+'m' and A-'m' have the same existing numeric cache name.
            // Their set order is opposite their first-observed read position.
            counts.incrementCounts(alignment("AC", "2M", "C+m?,0;A-m?,0;",
                    (byte) 200, (byte) 40));
            BaseModificationKey modified = counts.getAllModificationKeys().iterator().next();
            assertEquals('A', modified.getBase());
            assertEquals('-', modified.getStrand());
            assertEquals(1, counts.getCount(100, modified, 0, false));
            assertEquals(40, counts.getLikelihoodSum(100, modified, 0, false));
            assertEquals(1, counts.getCount(101, modified, 0.5f, false));
            assertEquals(200, counts.getLikelihoodSum(101, modified, 0.5f, true));
            BaseModificationKey noneT = BaseModificationKey.getKey('T', '+', "NONE_T");
            BaseModificationKey noneC = BaseModificationKey.getKey('C', '+', "NONE_C");
            assertEquals(Arrays.asList(modified, noneT, noneC), new ArrayList<>(counts.getAllModificationKeys()));
            assertEquals(215, counts.getLikelihoodSum(100, noneT, 0.5f, true));
            assertEquals(0, counts.getCount(101, noneC, 0, true));
        } finally {
            BaseModificationKey.keyCache = previousKeyCache;
        }
    }

    @Test
    public void duplicateSetOccurrencesContributeToResidualAndFirstWinnerWinsTies() {
        BaseModificationCounts counts = new BaseModificationCounts();
        SAMAlignment read = alignment("C", "1M", "C+mmh?,0;", (byte) 20, (byte) 70, (byte) 70);
        counts.incrementCounts(read);
        counts.incrementCounts(read);

        BaseModificationKey m = BaseModificationKey.getKey('C', '+', "m");
        BaseModificationKey h = BaseModificationKey.getKey('C', '+', "h");
        BaseModificationKey noneC = BaseModificationKey.getKey('C', '+', "NONE_C");
        assertEquals(Arrays.asList(m, h, noneC), new ArrayList<>(counts.getAllModificationKeys()));
        assertEquals(2, counts.getCount(100, m, 69.75f / 255, false));
        assertEquals(0, counts.getCount(100, m, 70.25f / 255, false));
        assertEquals(140, counts.getLikelihoodSum(100, m, 0, false));
        assertEquals(2, counts.getCount(100, m, 0, true));
        assertEquals(140, counts.getLikelihoodSum(100, m, 0, true));
        assertEquals(0, counts.getCount(100, h, 0, false));
        assertEquals(0, counts.getCount(100, noneC, 0, true));

        BaseModificationCounts residualTie = new BaseModificationCounts();
        residualTie.incrementCounts(alignment("C", "1M", "C+mh?,0;", (byte) 85, (byte) 85));
        assertEquals(85, residualTie.getLikelihoodSum(100, m, 0, true));
        assertEquals(0, residualTie.getCount(100, h, 0, true));
        assertEquals(0, residualTie.getCount(100, noneC, 0, true));
    }

    @Test
    public void softClipsInsertionsAndSplitBlockReadOffsetsDoNotShiftCalls() {
        String previousSoftClips = PreferencesManager.getPreferences().get(Constants.SAM_SHOW_SOFT_CLIPPED, null);
        PreferencesManager.getPreferences().put(Constants.SAM_SHOW_SOFT_CLIPPED, true);
        try {
            BaseModificationCounts counts = new BaseModificationCounts();
            counts.incrementCounts(alignment("CCCCCC", "1S1M1I1M2N1M1S",
                    "C+f?,0,4;C+m?,1,1,0;C+h?,2;",
                    (byte) 255, (byte) 255, (byte) 200, (byte) 130, (byte) 0, (byte) 255));

            BaseModificationKey m = BaseModificationKey.getKey('C', '+', "m");
            BaseModificationKey noneC = BaseModificationKey.getKey('C', '+', "NONE_C");
            BaseModificationKey f = BaseModificationKey.getKey('C', '+', "f");
            BaseModificationKey h = BaseModificationKey.getKey('C', '+', "h");
            assertEquals(Arrays.asList(m, noneC), new ArrayList<>(counts.getAllModificationKeys()));
            assertEquals(200, counts.getLikelihoodSum(100, m, 0, false));
            assertEquals(130, counts.getLikelihoodSum(101, m, 0.5f, false));
            assertEquals(0, counts.getCount(102, m, 0, false));
            assertEquals(0, counts.getCount(103, m, 0, false));
            assertEquals(1, counts.getCount(104, m, 0, false));
            assertEquals(0, counts.getCount(104, m, 0.01f, false));
            assertEquals(255, counts.getLikelihoodSum(104, noneC, 1, true));
            assertEquals(0, counts.getCount(99, f, 0, false));
            assertEquals(0, counts.getCount(105, f, 0, false));
            assertEquals(0, counts.getCount(101, h, 0, false));
        } finally {
            restorePreference(Constants.SAM_SHOW_SOFT_CLIPPED, previousSoftClips);
        }
    }

    @Test
    public void incrementCounts() throws IOException {

        String bamfile = TestUtils.DATA_DIR + "bam/chr20_mod_call_sample.bam";
        String chr = "20";
        int start = 13846149;
        int end = 13846294;

        BAMReader bamreader = new BAMReader(new ResourceLocator(bamfile), true);
        CloseableIterator<SAMAlignment> bamiter = bamreader.query(chr, start, end, false);
        int readCount = 0;

        BaseModificationCounts counts = new BaseModificationCounts();
        while (bamiter.hasNext()) {
            Alignment alignment = bamiter.next();
            counts.incrementCounts(alignment);
            readCount++;
        }
        assertTrue("No data retrieved:  " + readCount, readCount > 0);

        int[] expectedPositions = {13846181, 13846182, 13846227, 13846228, 13846232, 13846233, 13846234};
        int[] expectedCounts =    {1,        1,        1,        2,        1,        1,        1       };

        BaseModificationKey key =  BaseModificationKey.getKey('C', '+', "m");

        boolean includeNoMods = false;
        for(int i=0; i<expectedPositions.length; i++) {
            int c = counts.getCount(expectedPositions[i] - 1, key, 0, includeNoMods);
            assertEquals("Unexpected count at position " + expectedPositions[i], expectedCounts[i], c);
        }
    }

    @Test
    public void incrementCounts2() throws IOException {

        String bamfile = "https://www.dropbox.com/s/q32hk7tvsejryjt/HG002_chr11_119076212_119102218_2.bam";
        String indexFile = "https://www.dropbox.com/s/ax1ljny7ja5fdcu/HG002_chr11_119076212_119102218_2.bam.bai";
        String chr = "chr11";
        int start = 119094722;
        int end = 119094724;
        boolean includeNoMods = true;

        ResourceLocator locator = new ResourceLocator(bamfile);
        locator.setIndexPath(indexFile);
        BAMReader bamreader = new BAMReader(locator, true);
        CloseableIterator<SAMAlignment> bamiter = bamreader.query(chr, start, end, false);
        int readCount = 0;

        BaseModificationCounts counts = new BaseModificationCounts();
        while (bamiter.hasNext()) {
            Alignment alignment = bamiter.next();
            counts.incrementCounts(alignment);
            readCount++;
        }
        assertTrue("No data retrieved:  " + readCount, readCount > 0);

        BaseModificationKey cmKey = BaseModificationKey.getKey('C', '+', "m");
        int aboveThreshold = counts.getCount(119094723, cmKey, 0.5f, includeNoMods);
        assertEquals("Counts above threshold", 3, aboveThreshold);

        BaseModificationKey noModKey = BaseModificationKey.getKey('C', '+', "NONE_C");
        int belowThreshold = counts.getCount( 119094723, noModKey, 0.5f, includeNoMods);
        assertEquals("Counts below threshold", 12, belowThreshold);
    }
}