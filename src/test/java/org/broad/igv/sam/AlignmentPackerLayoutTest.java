package org.broad.igv.sam;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMSequenceRecord;
import org.broad.igv.Globals;
import org.broad.igv.feature.Chromosome;
import org.broad.igv.feature.genome.Genome;
import org.broad.igv.feature.genome.GenomeManager;
import org.broad.igv.prefs.IGVPreferences;
import org.broad.igv.prefs.PreferencesManager;
import org.broad.igv.track.Track;
import org.broad.igv.ui.panel.ReferenceFrame;
import org.junit.After;
import org.junit.Before;
import org.junit.Test;

import java.util.ArrayList;
import java.util.List;

import static org.broad.igv.prefs.Constants.SAM_SHOW_SOFT_CLIPPED;
import static org.junit.Assert.assertEquals;

/** Exact public packing layouts, with in-memory SAM records and no external data. */
public class AlignmentPackerLayoutTest {
    private static final String CHR = "chr1";
    private static final int CHROMOSOME_LENGTH = 20_000_000;

    private Genome previousGenome;
    private boolean previousHeadless;
    private boolean previousBatch;
    private IGVPreferences preferences;
    private String previousSoftClipped;
    private SAMFileHeader header;

    @Before
    public void setUp() {
        previousHeadless = Globals.isHeadless();
        previousBatch = Globals.isBatch();
        Globals.setHeadless(true);
        Globals.setBatch(false);
        previousGenome = GenomeManager.getInstance().getCurrentGenome();
        GenomeManager.getInstance().setCurrentGenomeForTest(new Genome("packing-layout-test",
                List.of(new Chromosome(0, CHR, CHROMOSOME_LENGTH))));
        preferences = PreferencesManager.getPreferences();
        previousSoftClipped = preferences.get(SAM_SHOW_SOFT_CLIPPED, null);
        preferences.put(SAM_SHOW_SOFT_CLIPPED, "true");
        header = new SAMFileHeader();
        header.addSequence(new SAMSequenceRecord(CHR, CHROMOSOME_LENGTH));
        header.setSortOrder(SAMFileHeader.SortOrder.coordinate);
    }

    @After
    public void tearDown() {
        if (previousSoftClipped == null) {
            // Remove rather than put(null), so an inherited/default boolean is
            // resolved again instead of leaving a cached false value behind.
            preferences.remove(SAM_SHOW_SOFT_CLIPPED);
        } else {
            preferences.put(SAM_SHOW_SOFT_CLIPPED, previousSoftClipped);
        }
        GenomeManager.getInstance().setCurrentGenomeForTest(previousGenome);
        Globals.setBatch(previousBatch);
        Globals.setHeadless(previousHeadless);
    }

    @Test
    public void densePackingPreservesQueueOrderAndRestartsAcrossEmptiedBucketsAndGaps() {
        List<Alignment> reads = List.of(
                read("tie-a", 1000, "61M"),
                read("tie-b", 1000, "61M"),
                read("tie-c", 1000, "61M"),
                read("tie-d", 1000, "61M"),
                read("short", 1000, "10M"),
                read("one-base-too-close", 1062, "1M"),
                read("bit-63", 1063, "1M"),
                read("bit-64", 1064, "1M"),
                read("bit-128", 1128, "1M"),
                read("spacing-boundary", 1131, "2M"),
                read("far", 2_001_000, "1M"),
                read("far-overlap", 2_001_001, "10M"));

        // Bucket positions are relative to 1000. Searches hit bit 63 exactly, start
        // between bits 64 and 128, cross a two-million-base gap, and finish beyond
        // the last occupied position but inside the range. Later rows revisit empty
        // buckets. Equal-end PriorityQueue entries deliberately are NOT stable:
        // with the shorter fifth entry, their removal order is a, b, d, c.
        // The long reads end at 1061, so start 1062 is too close for all four;
        // it can only follow the shorter read on the last row.
        assertLayout(pack(reads, options(), Track.DisplayMode.EXPANDED, null), List.of(""), List.of(List.of(
                List.of("tie-a", "bit-63", "bit-128", "spacing-boundary", "far"),
                List.of("tie-b", "bit-64", "far-overlap"),
                List.of("tie-d"),
                List.of("tie-c"),
                List.of("short", "one-base-too-close"))));
    }

    @Test
    public void softClippedStartsRetainFirstBucketClampAndSpacing() {
        List<Alignment> reads = List.of(
                read("zero", 0, "2M"),
                read("negative-clip", 1, "5S4M"),
                read("too-close", 6, "1M"),
                read("boundary", 7, "1M"),
                read("plain", 60, "2M"),
                read("positive-clip", 64, "5S4M"),
                read("tail", 70, "1M"));

        // Input is sorted by SAM alignment start, not the displayed clipped start.
        // The second read starts at -4 and is clamped into bucket zero, where its
        // later end wins over "zero". The other clipped start is 59, before "plain".
        assertEquals(-4, reads.get(1).getStart());
        assertEquals(59, reads.get(5).getStart());
        assertLayout(pack(reads, options(), Track.DisplayMode.EXPANDED, null), List.of(""), List.of(List.of(
                List.of("negative-clip", "boundary", "positive-clip", "tail"),
                List.of("zero", "too-close", "plain"))));
    }

    @Test
    public void groupedPackingPreservesGroupOrderAndEveryRow() {
        List<Alignment> reads = List.of(
                taggedRead("alpha-long", 100, "20M", "Alpha"),
                taggedRead("alpha-short", 100, "5M", "Alpha"),
                taggedRead("zeta-first", 101, "1M", "zeta"),
                taggedRead("zeta-next", 104, "1M", "zeta"),
                read("missing-tag", 107, "1M"),
                taggedRead("alpha-next", 122, "2M", "Alpha"));
        AlignmentTrack.RenderOptions options = options();
        options.setGroupByOption(AlignmentTrack.GroupOption.TAG);
        options.setGroupByTag("ZG");

        assertLayout(pack(reads, options, Track.DisplayMode.EXPANDED, null), List.of("Alpha", "zeta", ""), List.of(
                List.of(List.of("alpha-long", "alpha-next"), List.of("alpha-short")),
                List.of(List.of("zeta-first", "zeta-next")),
                List.of(List.of("missing-tag"))));
    }

    @Test
    public void pairedPackingPreservesBothMembersAndSecondaryReads() {
        SAMAlignment firstP = pairedRead("p", "p/1", 100, "10M", 150, true);
        SAMAlignment firstQ = pairedRead("q", "q/1", 100, "5M", 130, true);
        SAMRecord secondary = record("p", "p/secondary", 102, "3M");
        secondary.setSecondaryAlignment(true);
        secondary.setReadPairedFlag(true);
        secondary.setMateReferenceName(CHR);
        secondary.setMateAlignmentStart(151);
        List<Alignment> reads = List.of(firstP, firstQ, new SAMAlignment(secondary),
                pairedRead("q", "q/2", 130, "5M", 100, false),
                read("between-pairs", 137, "2M"),
                pairedRead("p", "p/2", 150, "10M", 100, false),
                read("tail", 162, "2M"));
        AlignmentTrack.RenderOptions options = options();
        options.setViewPairs(true);

        // Member-specific labels distinguish the two real SAM records sharing each
        // read name, and distinguish the primary pair from its secondary alignment.
        assertLayout(pack(reads, options, Track.DisplayMode.EXPANDED, null), List.of(""), List.of(List.of(
                List.of("pair(p/1,p/2)", "tail"),
                List.of("pair(q/1,q/2)", "between-pairs"),
                List.of("p/secondary"))));
        assertLayout(pack(reads, options, Track.DisplayMode.FULL, frame(0, 1000)), List.of(""), List.of(List.of(
                List.of("pair(p/1,p/2)"),
                List.of("pair(q/1,q/2)"),
                List.of("p/secondary"),
                List.of("between-pairs"),
                List.of("tail"))));
    }

    @Test
    public void layoutsMatchImmediatelyBelowAndAtSparseThreshold() {
        // The span is measured from the first displayed start to the maximum end:
        // 9,999,999 uses dense buckets; exactly 10,000,000 uses sparse buckets.
        for (int span : new int[]{9_999_999, 10_000_000}) {
            List<Alignment> reads = List.of(
                    read("long", 1000, "4M"),
                    read("short", 1000, "1M"),
                    read("next", 1006, "1M"),
                    read("far", 1000 + span - 1, "1M"));
            assertLayout(pack(reads, options(), Track.DisplayMode.EXPANDED, null), List.of(""), List.of(List.of(
                    List.of("long", "next", "far"),
                    List.of("short"))));
        }
    }

    @Test
    public void fullPackingRetainsInputOrderAndViewportOverlapBoundaries() {
        SAMRecord unmapped = record("unmapped", "unmapped", 150, "1M");
        unmapped.setReadUnmappedFlag(true);
        List<Alignment> reads = List.of(
                read("ends-at-left-edge", 90, "10M"),
                read("overlaps-left-edge", 95, "10M"),
                read("same-start-short", 100, "1M"),
                read("same-start-long", 100, "10M"),
                new SAMAlignment(unmapped),
                read("overlaps-right-edge", 199, "2M"),
                read("starts-at-right-edge", 200, "1M"));

        assertLayout(pack(reads, options(), Track.DisplayMode.FULL, frame(100, 200)), List.of(""), List.of(List.of(
                List.of("overlaps-left-edge"),
                List.of("same-start-short"),
                List.of("same-start-long"),
                List.of("overlaps-right-edge"))));
    }

    private AlignmentTrack.RenderOptions options() {
        AlignmentTrack.RenderOptions options = new AlignmentTrack.RenderOptions(null);
        options.setGroupByOption(AlignmentTrack.GroupOption.NONE);
        options.setLinkedReads(false);
        options.setViewPairs(false);
        return options;
    }

    private PackedAlignments pack(List<Alignment> reads, AlignmentTrack.RenderOptions options,
                                  Track.DisplayMode mode, ReferenceFrame frame) {
        AlignmentInterval interval = new AlignmentInterval(CHR, 0, CHROMOSOME_LENGTH, reads,
                null, null, null, frame);
        return new AlignmentPacker().packAlignments(interval, options, frame, mode);
    }

    private SAMAlignment read(String name, int start, String cigar) {
        return new SAMAlignment(record(name, name, start, cigar));
    }

    private SAMAlignment taggedRead(String name, int start, String cigar, String group) {
        SAMRecord record = record(name, name, start, cigar);
        record.setAttribute("ZG", group);
        return new SAMAlignment(record);
    }

    private SAMAlignment pairedRead(String name, String member, int start, String cigar, int mateStart,
                                    boolean first) {
        SAMRecord record = record(name, member, start, cigar);
        record.setReadPairedFlag(true);
        record.setFirstOfPairFlag(first);
        record.setSecondOfPairFlag(!first);
        record.setReadNegativeStrandFlag(!first);
        record.setMateNegativeStrandFlag(first);
        record.setMateReferenceName(CHR);
        record.setMateAlignmentStart(mateStart + 1);
        return new SAMAlignment(record);
    }

    /** Starts are zero-based IGV coordinates; SAMRecord itself uses one-based starts. */
    private SAMRecord record(String name, String label, int start, String cigar) {
        SAMRecord record = new SAMRecord(header);
        record.setReadName(name);
        record.setReferenceName(CHR);
        record.setAlignmentStart(start + 1);
        record.setCigarString(cigar);
        record.setMappingQuality(60);
        record.setReadString("A".repeat(record.getCigar().getReadLength()));
        record.setBaseQualityString("I".repeat(record.getReadLength()));
        record.setAttribute("ZI", label);
        return record;
    }

    private ReferenceFrame frame(int start, int end) {
        ReferenceFrame frame = new ReferenceFrame("packing-layout-test");
        frame.setBounds(0, 100);
        frame.jumpTo(CHR, start, end);
        return frame;
    }

    private void assertLayout(PackedAlignments packed, List<String> groups,
                              List<List<List<String>>> expectedRows) {
        assertEquals("Ordered groups", groups, new ArrayList<>(packed.keySet()));
        List<List<List<String>>> actualRows = new ArrayList<>();
        for (List<Row> rows : packed.values()) {
            List<List<String>> groupRows = new ArrayList<>();
            for (Row row : rows) {
                List<String> members = new ArrayList<>();
                for (Alignment alignment : row.alignments) {
                    members.add(label(alignment));
                }
                groupRows.add(members);
            }
            actualRows.add(groupRows);
        }
        assertEquals("Complete ordered rows and members", expectedRows, actualRows);
    }

    private String label(Alignment alignment) {
        if (alignment instanceof PairedAlignment) {
            PairedAlignment pair = (PairedAlignment) alignment;
            return "pair(" + label(pair.getFirstAlignment()) + "," + label(pair.getSecondAlignment()) + ")";
        }
        return (String) alignment.getAttribute("ZI");
    }
}
