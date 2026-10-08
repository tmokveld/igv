/*
 * The MIT License (MIT)
 *
 * Copyright (c) 2007-2015 Broad Institute
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 *
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
 * THE SOFTWARE.
 */

package org.broad.igv.sam;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMSequenceRecord;

import org.broad.igv.AbstractHeadlessTest;
import org.broad.igv.sam.reader.AlignmentReader;
import org.broad.igv.sam.reader.AlignmentReaderFactory;
import org.broad.igv.util.TestUtils;
import org.junit.Test;
import java.util.Arrays;
import java.util.Iterator;
import static org.junit.Assert.*;

/**
 * @author jacob
 * @date 2013-Oct-31
 */
public class DenseAlignmentCountsTest extends AbstractHeadlessTest {

    @Test
    public void countsBothStrandsQualitiesAndDynamicBasesWithinSmallInterval() {
        DenseAlignmentCounts counts = new DenseAlignmentCounts(100, 104, null);
        counts.incCounts(alignment(99, "nAr=tN", new byte[]{10, 20, 30, 40, 50, 60}, false));
        counts.incCounts(alignment(100, "aR=t", new byte[]{7, 8, 9, 10}, true));
        counts.incCounts(alignment(101, "r", new byte[]{11}, false));

        byte[] bases = {'a', 'r', '=', 't'};
        int[] totals = {2, 3, 2, 2};
        int[] positiveCounts = {1, 2, 1, 1};
        int[] qualitySums = {27, 49, 49, 60};
        for (int i = 0; i < bases.length; i++) {
            int position = 100 + i;
            assertEquals(totals[i], counts.getTotalCount(position));
            assertEquals(totals[i], counts.getCount(position, bases[i]));
            assertEquals(positiveCounts[i], counts.getPosCount(position, bases[i]));
            assertEquals(1, counts.getNegCount(position, bases[i]));
            assertEquals(positiveCounts[i], counts.getTotalPositiveCount(position));
            assertEquals(1, counts.getTotalNegativeCount(position));
            assertEquals(qualitySums[i], counts.getQuality(position, bases[i]));
            assertEquals(qualitySums[i], counts.getTotalQuality(position));
        }
        assertTrue(counts.getBases().contains((byte) 'R'));
        assertTrue(counts.getBases().contains((byte) '='));
        assertFalse(counts.getBases().contains((byte) 'r'));
        assertEquals(0, counts.getTotalCount(99));
        assertEquals(0, counts.getTotalCount(104));
        assertEquals(0, counts.getCount(99, (byte) 'n'));
        assertEquals(0, counts.getCount(104, (byte) 'N'));
        assertEquals(0, counts.getTotalQuality(99));
        assertEquals(0, counts.getTotalQuality(104));
        assertEquals(3, counts.getMaxCount(100, 104));
        assertEquals(3, counts.getMaxCount(0, 200));
    }

    @Test
    public void maxCoverageIncludesLastBaseOfLargeIntervalWithoutCountingOutsideBases() {
        int start = 1000;
        int end = 2027;
        DenseAlignmentCounts counts = new DenseAlignmentCounts(start, end, null);
        String spanningBases = "A".repeat(end - start + 2);
        byte[] qualities = new byte[spanningBases.length()];
        Arrays.fill(qualities, (byte) 25);
        counts.incCounts(alignment(start - 1, spanningBases, qualities, false));
        counts.incCounts(alignment(start - 1, spanningBases, qualities, true));
        counts.incCounts(alignment(end - 1, "a", new byte[]{30}, false));

        assertEquals(2, counts.getTotalCount(start));
        assertEquals(2, counts.getTotalCount(end - 2));
        assertEquals(3, counts.getTotalCount(end - 1));
        assertEquals(2, counts.getPosCount(end - 1, (byte) 'A'));
        assertEquals(1, counts.getNegCount(end - 1, (byte) 'a'));
        assertEquals(80, counts.getQuality(end - 1, (byte) 'a'));
        assertEquals(80, counts.getTotalQuality(end - 1));
        assertEquals(0, counts.getTotalCount(start - 1));
        assertEquals(0, counts.getTotalCount(end));
        assertEquals(0, counts.getTotalQuality(start - 1));
        assertEquals(0, counts.getTotalQuality(end));
        assertEquals(3, counts.getMaxCount(start, end));
        assertEquals(3, counts.getMaxCount(start - 100, end + 100));
        assertEquals(3, counts.getMaxCount(end - 1, end));
    }

    @Test
    public void testEqBase() throws Exception {
        String chr = "chr1";
        int start = 59300;
        int end = 59400;
        String path = TestUtils.DATA_DIR + "bam/squeeze.sam";
        DenseAlignmentCounts daCounts = new DenseAlignmentCounts(start, end, null);

        AlignmentReader reader = AlignmentReaderFactory.getReader(path, false);
        Iterator<Alignment> iter = reader.iterator();
        while (iter.hasNext()) {
            Alignment al = iter.next();
            daCounts.incCounts(al);
        }
        int count = daCounts.getCount(59303, (byte) '=');
        assertEquals(3, count);

    }

    private Alignment alignment(int start, String bases, byte[] qualities, boolean negativeStrand) {
        SAMFileHeader header = new SAMFileHeader();
        header.addSequence(new SAMSequenceRecord("chr1", 10000));
        SAMRecord record = new SAMRecord(header);
        record.setReadName("coverage");
        record.setReferenceName("chr1");
        record.setAlignmentStart(start + 1);
        record.setCigarString(bases.length() + "M");
        record.setReadString(bases);
        record.setBaseQualities(qualities);
        record.setReadNegativeStrandFlag(negativeStrand);
        return new SAMAlignment(record);
    }
}
