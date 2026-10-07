package org.broad.igv.sam;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMFileWriter;
import htsjdk.samtools.SAMFileWriterFactory;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMSequenceRecord;
import org.broad.igv.Globals;
import org.broad.igv.feature.genome.Genome;
import org.broad.igv.feature.genome.load.GenomeConfig;
import org.broad.igv.feature.genome.GenomeManager;
import org.broad.igv.feature.genome.InMemorySequence;
import org.broad.igv.prefs.IGVPreferences;
import org.broad.igv.prefs.PreferencesManager;
import org.broad.igv.renderer.DataRange;
import org.broad.igv.renderer.SequenceRenderer;
import org.broad.igv.track.RenderContext;
import org.broad.igv.track.Track;
import org.broad.igv.ui.panel.FrameManager;
import org.broad.igv.ui.panel.ReferenceFrame;
import org.broad.igv.util.ResourceLocator;
import org.junit.After;
import org.junit.Before;
import org.junit.Rule;
import org.junit.Test;
import org.junit.rules.TemporaryFolder;

import javax.swing.JPanel;
import java.awt.AlphaComposite;
import java.awt.Color;
import java.awt.Graphics2D;
import java.awt.Rectangle;
import java.awt.RenderingHints;
import java.awt.image.BufferedImage;
import java.io.File;
import java.lang.reflect.Field;
import java.nio.file.Files;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;

import static org.broad.igv.prefs.Constants.*;
import static org.junit.Assert.*;

/** Pixel contracts at the public render seam; no assertions about raster-cache internals. */
public class CoverageTrackRenderTest {
    private static final int WIDTH = 240;
    private static final int VIEW_HEIGHT = 140;
    private static final String CHR = "chrCoverageRaster";
    private static final Rectangle DEFAULT_RECT = new Rectangle(0, 0, WIDTH, 60);
    private static final Color COVERAGE_COLOR = new Color(93, 117, 137);

    @Rule
    public TemporaryFolder temporaryFolder = new TemporaryFolder();

    private boolean previousHeadless;
    private boolean previousBatch;
    private Genome previousGenome;
    private List<ReferenceFrame> previousFrames;
    private Map<Character, Color> previousPalette;
    private IGVPreferences preferences;
    private final Map<String, String> previousPreferences = new LinkedHashMap<>();
    private final List<AlignmentTrack> tracks = new ArrayList<>();
    private Genome genome;
    private File bam;
    private ReferenceFrame frame;
    private AlignmentTrack track;
    // This is only global fixture isolation, not renderer/count instrumentation.
    private Field knownSnpsField;
    private Object previousKnownSnps;

    @Before
    public void setUp() throws Exception {
        previousHeadless = Globals.isHeadless();
        previousBatch = Globals.isBatch();
        previousGenome = GenomeManager.getInstance().getCurrentGenome();
        previousFrames = FrameManager.getFrames();
        previousPalette = SequenceRenderer.nucleotideColors;
        Globals.setHeadless(true);
        Globals.setBatch(false);
        preferences = PreferencesManager.getPreferences();
        for (String key : List.of(SAM_DOWNSAMPLE_READS, SAM_REDUCED_MEMORY_MODE,
                SAM_SHOW_REF_SEQ, SAM_SHOW_SOFT_CLIPPED, SAM_FILTER_ALIGNMENTS,
                SAM_QUALITY_THRESHOLD, SAM_ALIGNMENT_SCORE_THRESHOLD, SAM_DISPLAY_PAIRED,
                SAM_MAX_VISIBLE_RANGE, SAM_ALLELE_THRESHOLD, SAM_ALLELE_USE_QUALITY,
                ENABLE_ANTIALISING, KNOWN_SNPS)) {
            previousPreferences.put(key, preferences.get(key, null));
        }
        preferences.put(SAM_DOWNSAMPLE_READS, "false");
        preferences.put(SAM_REDUCED_MEMORY_MODE, "false");
        preferences.put(SAM_SHOW_REF_SEQ, "false");
        preferences.put(SAM_SHOW_SOFT_CLIPPED, "false");
        preferences.put(SAM_FILTER_ALIGNMENTS, "false");
        preferences.put(SAM_QUALITY_THRESHOLD, "0");
        preferences.put(SAM_ALIGNMENT_SCORE_THRESHOLD, "0");
        preferences.put(SAM_DISPLAY_PAIRED, "false");
        preferences.put(SAM_MAX_VISIBLE_RANGE, "1000");
        preferences.put(SAM_ALLELE_THRESHOLD, "0.2");
        preferences.put(SAM_ALLELE_USE_QUALITY, "false");
        preferences.put(ENABLE_ANTIALISING, "false");
        preferences.put(KNOWN_SNPS, null);
        knownSnpsField = BaseAlignmentCounts.class.getDeclaredField("knownSnps");
        knownSnpsField.setAccessible(true);
        previousKnownSnps = knownSnpsField.get(null);
        knownSnpsField.set(null, null);
        SequenceRenderer.nucleotideColors = palette();

        byte[] reference = new byte[10000];
        Arrays.fill(reference, (byte) 'A');
        GenomeConfig config = new GenomeConfig();
        config.setId("coverage-raster-test");
        config.setName("coverage-raster-test");
        config.setSequence(new InMemorySequence(CHR, reference));
        genome = new Genome(config);
        GenomeManager.getInstance().setCurrentGenomeForTest(genome);
        bam = writeBam();
        frame = new ReferenceFrame("coverage-raster-test");
        frame.setBounds(0, WIDTH);
        frame.jumpTo(CHR, 960, 1140);
        FrameManager.setFrames(List.of(frame));
        track = newTrack();

        AlignmentCounts counts = interval(track).getCounts();
        assertTrue("Tiny fixture must retain per-base counts", counts.hasBaseCounts());
        assertEquals("All deterministic fixture reads must contribute", 16, counts.getTotalCount(1040));
        assertEquals("Low-quality allele spike must be present", 4, counts.getCount(1040, (byte) 'c'));
        assertEquals("High-quality allele spike must be present", 8, counts.getCount(1065, (byte) 'c'));
    }

    @After
    public void tearDown() throws Exception {
        try {
            for (AlignmentTrack loaded : new ArrayList<>(tracks)) closeTrack(loaded);
        } finally {
            if (knownSnpsField != null) knownSnpsField.set(null, previousKnownSnps);
            if (preferences != null) previousPreferences.forEach(preferences::put);
            SequenceRenderer.nucleotideColors = previousPalette;
            GenomeManager.getInstance().setCurrentGenomeForTest(previousGenome);
            if (previousFrames != null) FrameManager.setFrames(previousFrames);
            Globals.setHeadless(previousHeadless);
            Globals.setBatch(previousBatch);
        }
    }

    @Test
    public void fullAndPartialPaintsMatchAtOneAndTwoTimesIncludingBottomBorder() throws Exception {
        for (int deviceScale : new int[]{1, 2}) {
            for (Rectangle rect : List.of(DEFAULT_RECT, new Rectangle(13, 27, WIDTH, 60))) {
                BufferedImage baseline = fresh(rect, deviceScale, 0, GraphicsCase.NORMAL, false, false);
                assertPixels("Full paint", baseline, paint(track, rect, deviceScale, 0, null));
                assertPixels("Warm full paint", baseline, paint(track, rect, deviceScale, 0, null));
                assertEquals("Bottom border really exists outside the height-60 data area", Color.GRAY.getRGB(),
                        baseline.getRGB(100 * deviceScale, (rect.y + 60) * deviceScale));
                for (int offset : new int[]{0, 20, 40, 59, 60, 61}) {
                    int scroll = rect.y + offset;
                    BufferedImage full = fresh(rect, deviceScale, scroll, GraphicsCase.NORMAL, false, false);
                    assertPixels("Scrolled warm paint at " + offset, full,
                            paint(track, rect, deviceScale, scroll, null));
                    for (Rectangle damage : List.of(
                            new Rectangle(rect.x, scroll, rect.width, 1),
                            new Rectangle(rect.x + 33, scroll + 2, 107, 19),
                            new Rectangle(rect.x + 80, scroll, 1, VIEW_HEIGHT),
                            new Rectangle(rect.x, rect.y + 60, rect.width, 1),
                            new Rectangle(rect.x, scroll + 100, rect.width, 20))) {
                        assertDamage(full, paint(track, rect, deviceScale, scroll, damage),
                                viewport(rect, scroll, false), damage, deviceScale);
                    }
                    if (offset == 60) {
                        assertEquals("Y60 exposure must retain the live bottom border", Color.GRAY.getRGB(),
                                full.getRGB(100 * deviceScale, 0));
                    } else if (offset == 61) {
                        assertWhite("Y61 is below all coverage and border pixels", full);
                    }
                }
            }
        }
    }

    @Test
    public void partialInitialPaintDoesNotPoisonLaterFullOrUncoveredExposure() throws Exception {
        for (int deviceScale : new int[]{1, 2}) {
            for (Rectangle rect : List.of(DEFAULT_RECT, new Rectangle(13, 27, WIDTH, 60))) {
                AlignmentTrack initiallyClipped = newTrack();
                try {
                    BufferedImage full = fresh(rect, deviceScale, 0, GraphicsCase.NORMAL, false, false);
                    Rectangle firstDamage = new Rectangle(rect.x + 61, rect.y + 48, 1, 2);
                    assertDamage(full, paint(initiallyClipped, rect, deviceScale, 0, firstDamage),
                            viewport(rect, 0, false), firstDamage, deviceScale);
                    assertPixels("Full exposure after a one-pixel initial clip", full,
                            paint(initiallyClipped, rect, deviceScale, 0, null));
                    Rectangle uncovered = new Rectangle(rect.x + 100, rect.y, 70, 61);
                    assertDamage(full, paint(initiallyClipped, rect, deviceScale, 0, uncovered),
                            viewport(rect, 0, false), uncovered, deviceScale);
                    BufferedImage scrolled = fresh(rect, deviceScale, rect.y + 40,
                            GraphicsCase.NORMAL, false, false);
                    assertPixels("Vertical exposure reuses complete data, not the initial damage", scrolled,
                            paint(initiallyClipped, rect, deviceScale, rect.y + 40, null));
                } finally {
                    closeTrack(initiallyClipped);
                }
            }
        }
    }

    @Test
    public void disjointInitialDamageDoesNotSuppressTheNextVisiblePaint() throws Exception {
        for (int deviceScale : new int[]{1, 2}) {
            AlignmentTrack initiallyHidden = newTrack();
            try {
                Rectangle damage = new Rectangle(0, 90, WIDTH, 30);
                assertWhite("Read-only damage must not paint coverage", paint(initiallyHidden,
                        DEFAULT_RECT, deviceScale, 0, damage));
                assertPixels("Coverage subsequently exposed after disjoint damage",
                        fresh(DEFAULT_RECT, deviceScale, 0, GraphicsCase.NORMAL, false, false),
                        paint(initiallyHidden, DEFAULT_RECT, deviceScale, 0, null));
            } finally {
                closeTrack(initiallyHidden);
            }
        }
    }

    @Test
    public void horizontalPlacementPreservesLeftEdgeRoundingAndVerticalRelocation() throws Exception {
        frame.setOrigin(960.5);
        assertEquals("Left-edge regression needs fractional genomic-to-pixel mapping", 0.75, frame.getScale(), 0);
        AlignmentInterval original = interval(track);
        Rectangle shiftedX = new Rectangle(10, 0, WIDTH, 60);
        Rectangle shiftedXY = new Rectangle(10, 27, WIDTH, 60);
        for (int deviceScale : new int[]{1, 2}) {
            BufferedImage atZero = fresh(DEFAULT_RECT, deviceScale, 0, GraphicsCase.NORMAL, false, false);
            assertPixels("Initial X=0 left-edge paint", atZero, paint(track, DEFAULT_RECT, deviceScale, 0, null));
            BufferedImage atTen = fresh(shiftedX, deviceScale, 0, GraphicsCase.NORMAL, false, false);
            assertDifferent("Integer truncation must genuinely affect the left-edge allele", atZero, atTen);
            assertEquals("X=0 retains the just-before-origin C spike", SequenceRenderer.nucleotideColors.get('c').getRGB(),
                    atZero.getRGB(0, 57 * deviceScale));
            assertEquals("X=10 clips that allele before the viewport", COVERAGE_COLOR.getRGB(),
                    atTen.getRGB(0, 57 * deviceScale));
            assertPixels("Same frame and width at X=10", atTen, paint(track, shiftedX, deviceScale, 0, null));
            assertPixels("Vertical relocation after X=10", fresh(shiftedXY, deviceScale, 0,
                    GraphicsCase.NORMAL, false, false), paint(track, shiftedXY, deviceScale, 0, null));
            assertPixels("Returning to X=0 restores original rounding", atZero,
                    paint(track, DEFAULT_RECT, deviceScale, 0, null));
            assertSame("Placement changes must not replace loaded data", original, interval(track));
            assertEquals("Placement changes keep the same frame origin", 960.5, frame.getOrigin(), 0);
        }
    }

    @Test
    public void originScaleAndDimensionsInvalidateWarmPixels() throws Exception {
        for (int deviceScale : new int[]{1, 2}) {
            frame.setBounds(0, WIDTH);
            frame.jumpTo(CHR, 960, 1140);
            BufferedImage before = paint(track, DEFAULT_RECT, deviceScale, 0, null);
            frame.setOrigin(970);
            assertChangedAndFresh("Origin", before, DEFAULT_RECT, deviceScale);

            before = paint(track, DEFAULT_RECT, deviceScale, 0, null);
            double previousScale = frame.getScale();
            frame.jumpTo(CHR, 970, 1210);
            assertNotEquals("Scale mutation must take effect", previousScale, frame.getScale(), 0);
            assertChangedAndFresh("Genomic scale", before, DEFAULT_RECT, deviceScale);

            before = paint(track, DEFAULT_RECT, deviceScale, 0, null);
            Rectangle narrow = new Rectangle(0, 0, 180, 60);
            // Keep the canvas/viewport fixed so a stale wider raster is observable.
            BufferedImage expected = fresh(narrow, deviceScale, 0, GraphicsCase.NORMAL, false, true);
            BufferedImage actual = paint(track, narrow, deviceScale, 0, null,
                    GraphicsCase.NORMAL, false, true);
            assertPixels("Narrower width with wider viewport", expected, actual);
            assertNotEquals("Wider renderer really paints bars in the removed columns", Color.WHITE.getRGB(),
                    before.getRGB(220 * deviceScale, 59 * deviceScale));
            assertEquals("No stale bars beyond the narrower renderer boundary", Color.WHITE.getRGB(),
                    actual.getRGB(220 * deviceScale, 59 * deviceScale));

            Rectangle tall = new Rectangle(0, 0, WIDTH, 80);
            assertChangedAndFresh("Track height", before, tall, deviceScale);

            // Matching viewport geometry exercises the ordinary cache-eligible width case as well.
            frame.setBounds(0, 180);
            expected = fresh(narrow, deviceScale, 0, GraphicsCase.NORMAL, false, false);
            assertPixels("New matching frame/rectangle width", expected,
                    paint(track, narrow, deviceScale, 0, null));
            assertPixels("Warm matching frame/rectangle width", expected,
                    paint(track, narrow, deviceScale, 0, null));
        }
    }

    @Test
    public void rangeAndLogMutationsInvalidateWarmPixels() throws Exception {
        for (int deviceScale : new int[]{1, 2}) {
            coverage(track).setDataRange(new DataRange(0, 0, 24));
            BufferedImage before = paint(track, DEFAULT_RECT, deviceScale, 0, null);
            coverage(track).setDataRange(new DataRange(0, 0, 32));
            assertChangedAndFresh("Numeric maximum", before, DEFAULT_RECT, deviceScale);
            before = paint(track, DEFAULT_RECT, deviceScale, 0, null);
            coverage(track).getDataRange().setType(DataRange.Type.LOG);
            assertChangedAndFresh("In-place log-range mutation", before, DEFAULT_RECT, deviceScale);
        }
    }

    @Test
    public void coverageColorAndActualLowercasePaletteValuesInvalidateWarmPixels() throws Exception {
        for (int deviceScale : new int[]{1, 2}) {
            coverage(track).setColor(COVERAGE_COLOR);
            SequenceRenderer.nucleotideColors = palette();
            BufferedImage before = paint(track, DEFAULT_RECT, deviceScale, 0, null);
            coverage(track).setColor(new Color(173, 59, 97));
            assertChangedAndFresh("Coverage color", before, DEFAULT_RECT, deviceScale);
            before = paint(track, DEFAULT_RECT, deviceScale, 0, null);
            // Only lower-case c changes: that is the value CoverageTrack actually consumes.
            SequenceRenderer.nucleotideColors.put('c', new Color(237, 31, 153));
            assertChangedAndFresh("In-place lower-case nucleotide palette", before, DEFAULT_RECT, deviceScale);
            before = paint(track, DEFAULT_RECT, deviceScale, 0, null);
            Map<Character, Color> replacement = new HashMap<>(SequenceRenderer.nucleotideColors);
            replacement.put('c', new Color(11, 197, 83));
            SequenceRenderer.nucleotideColors = replacement;
            assertChangedAndFresh("Replacement custom nucleotide palette", before, DEFAULT_RECT, deviceScale);
        }
    }

    @Test
    public void alleleThresholdAndQualityWeightingChangeRealSpikePixels() throws Exception {
        for (int deviceScale : new int[]{1, 2}) {
            coverage(track).setSnpThreshold(0.2f);
            preferences.put(SAM_ALLELE_USE_QUALITY, "false");
            BufferedImage before = paint(track, DEFAULT_RECT, deviceScale, 0, null);
            assertSpikeColor(before, 1040, deviceScale, SequenceRenderer.nucleotideColors.get('c'));
            coverage(track).setSnpThreshold(0.8f);
            BufferedImage suppressed = assertChangedAndFresh("Allele threshold", before, DEFAULT_RECT, deviceScale);
            assertSpikeColor(suppressed, 1040, deviceScale, COVERAGE_COLOR);

            coverage(track).setSnpThreshold(0.2f);
            before = paint(track, DEFAULT_RECT, deviceScale, 0, null);
            preferences.put(SAM_ALLELE_USE_QUALITY, "true");
            suppressed = assertChangedAndFresh("Allele quality weighting", before, DEFAULT_RECT, deviceScale);
            assertSpikeColor(suppressed, 1040, deviceScale, COVERAGE_COLOR);
            assertSpikeColor(suppressed, 1065, deviceScale, SequenceRenderer.nucleotideColors.get('c'));
        }
    }

    @Test
    public void actualReloadReplacesIntervalAndCountsAtAnUnchangedLocus() throws Exception {
        BufferedImage before = paint(track, DEFAULT_RECT, 2, 0, null);
        AlignmentInterval original = interval(track);
        double origin = frame.getOrigin();
        double scale = frame.getScale();
        // Four low-MAPQ reads carry the low-quality C spike. A normal reload filters them.
        preferences.put(SAM_QUALITY_THRESHOLD, "30");
        track.getDataManager().clear();
        track.load(frame);
        AlignmentInterval reloaded = interval(track);
        assertNotSame("A real fixture load must replace the published interval", original, reloaded);
        assertNotSame("A real fixture load must replace the finalized counts", original.getCounts(), reloaded.getCounts());
        assertEquals("Locus must not hide interval-identity invalidation", origin, frame.getOrigin(), 0);
        assertEquals("Mapping remains unchanged during the reload", scale, frame.getScale(), 0);
        assertEquals("Reload must change actual coverage", 12, reloaded.getCounts().getTotalCount(1040));
        assertEquals("Filtered reads remove the low-quality allele", 0, reloaded.getCounts().getCount(1040, (byte) 'c'));
        assertChangedAndFresh("Published interval/count replacement", before, DEFAULT_RECT, 2);
        assertPixels("Replacement remains correct when warm", fresh(DEFAULT_RECT, 2, 0,
                GraphicsCase.NORMAL, false, false), paint(track, DEFAULT_RECT, 2, 0, null));
    }

    @Test
    public void independentlyLoadedTracksCannotReuseEachOthersPixels() throws Exception {
        AlignmentTrack other = newTrack();
        try {
            coverage(other).setColor(new Color(191, 83, 29));
            coverage(other).setSnpThreshold(0.8f);
            for (int deviceScale : new int[]{1, 2}) {
                BufferedImage first = paint(track, DEFAULT_RECT, deviceScale, 0, null);
                BufferedImage second = paint(other, DEFAULT_RECT, deviceScale, 0, null);
                assertDifferent("Separate tracks must genuinely paint differently", first, second);
                assertPixels("First track after other track paints", first,
                        paint(track, DEFAULT_RECT, deviceScale, 0, null));
                assertPixels("Other track after first track paints", second,
                        paint(other, DEFAULT_RECT, deviceScale, 0, null));
                coverage(other).setColor(new Color(43, 149, 181));
                paint(other, DEFAULT_RECT, deviceScale, 0, null);
                assertPixels("Changing another track must not damage this track", first,
                        paint(track, DEFAULT_RECT, deviceScale, 0, null));
            }
        } finally {
            closeTrack(other);
        }
    }

    @Test
    public void unsupportedGraphicsMatchFreshRenderingAfterAnOrdinaryWarmPaint() throws Exception {
        for (int deviceScale : new int[]{1, 2}) {
            for (GraphicsCase graphicsCase : List.of(GraphicsCase.ANTIALIASED,
                    GraphicsCase.FRACTIONAL_TRANSLATION, GraphicsCase.FRACTIONAL_SCALE,
                    GraphicsCase.TRANSLUCENT_COMPOSITE, GraphicsCase.NULL_CLIP)) {
                paint(track, DEFAULT_RECT, deviceScale, 0, null);
                assertPixels("Unsupported " + graphicsCase,
                        fresh(DEFAULT_RECT, deviceScale, 0, graphicsCase, false, false),
                        paint(track, DEFAULT_RECT, deviceScale, 0, null, graphicsCase, false, false));
                assertPixels("Ordinary rendering after " + graphicsCase,
                        fresh(DEFAULT_RECT, deviceScale, 0, GraphicsCase.NORMAL, false, false),
                        paint(track, DEFAULT_RECT, deviceScale, 0, null));
            }
        }
    }

    @Test
    public void translucentColorsAndPaletteMatchFreshRenderingAndReturnToOpaqueCorrectly() throws Exception {
        for (int deviceScale : new int[]{1, 2}) {
            coverage(track).setColor(COVERAGE_COLOR);
            SequenceRenderer.nucleotideColors = palette();
            BufferedImage before = paint(track, DEFAULT_RECT, deviceScale, 0, null);
            coverage(track).setColor(new Color(93, 117, 137, 93));
            assertChangedAndFresh("Translucent coverage color", before, DEFAULT_RECT, deviceScale);
            coverage(track).setColor(COVERAGE_COLOR);
            assertPixels("Opaque coverage after translucent fallback",
                    fresh(DEFAULT_RECT, deviceScale, 0, GraphicsCase.NORMAL, false, false),
                    paint(track, DEFAULT_RECT, deviceScale, 0, null));
            before = paint(track, DEFAULT_RECT, deviceScale, 0, null);
            SequenceRenderer.nucleotideColors.put('c', new Color(37, 71, 211, 87));
            assertChangedAndFresh("Translucent allele palette", before, DEFAULT_RECT, deviceScale);
            SequenceRenderer.nucleotideColors = palette();
            assertPixels("Opaque palette after translucent fallback",
                    fresh(DEFAULT_RECT, deviceScale, 0, GraphicsCase.NORMAL, false, false),
                    paint(track, DEFAULT_RECT, deviceScale, 0, null));
        }
    }

    @Test
    public void specializedColorModesAndMultiframeMatchFreshRendering() throws Exception {
        for (int deviceScale : new int[]{1, 2}) {
            for (AlignmentTrack.ColorOption colorOption : List.of(AlignmentTrack.ColorOption.READ_STRAND,
                    AlignmentTrack.ColorOption.BISULFITE, AlignmentTrack.ColorOption.BASE_MODIFICATION,
                    AlignmentTrack.ColorOption.BASE_MODIFICATION_2COLOR)) {
                track.getRenderOptions().setColorOption(AlignmentTrack.ColorOption.NONE);
                paint(track, DEFAULT_RECT, deviceScale, 0, null);
                track.getRenderOptions().setColorOption(colorOption);
                assertPixels("Color-mode fallback " + colorOption,
                        fresh(DEFAULT_RECT, deviceScale, 0, GraphicsCase.NORMAL, false, false),
                        paint(track, DEFAULT_RECT, deviceScale, 0, null));
            }
            track.getRenderOptions().setColorOption(AlignmentTrack.ColorOption.NONE);
            BufferedImage single = paint(track, DEFAULT_RECT, deviceScale, 0, null);
            BufferedImage multiple = paint(track, DEFAULT_RECT, deviceScale, 0, null,
                    GraphicsCase.NORMAL, true, false);
            assertDifferent("Multiframe must really suppress the live scale label", single, multiple);
            assertPixels("Multiframe fallback", fresh(DEFAULT_RECT, deviceScale, 0,
                    GraphicsCase.NORMAL, true, false), multiple);
            assertPixels("Scale label returns on the next ordinary paint", single,
                    paint(track, DEFAULT_RECT, deviceScale, 0, null));
        }
    }

    @Test
    public void widerViewportPreservesLegacyBarOverhangInsteadOfCroppingToTrackBounds() throws Exception {
        Rectangle placed = new Rectangle(13, 27, WIDTH, 60);
        for (int deviceScale : new int[]{1, 2}) {
            paint(track, placed, deviceScale, 0, null);
            BufferedImage expected = fresh(placed, deviceScale, 0, GraphicsCase.NORMAL, false, true);
            BufferedImage actual = paint(track, placed, deviceScale, 0, null,
                    GraphicsCase.NORMAL, false, true);
            assertPixels("Wider viewport requires uncropped fallback", expected, actual);
            assertNotEquals("Fixture really paints a left-side data overhang", Color.WHITE.getRGB(),
                    expected.getRGB((placed.x - 1) * deviceScale, (placed.y + 59) * deviceScale));
            assertNotEquals("Fixture really paints a right-side data overhang", Color.WHITE.getRGB(),
                    expected.getRGB((placed.x + placed.width) * deviceScale, (placed.y + 59) * deviceScale));
            Rectangle damage = new Rectangle(placed.x - 2, placed.y + 55, placed.width + 3, 6);
            assertDamage(expected, paint(track, placed, deviceScale, 0, damage,
                    GraphicsCase.NORMAL, false, true), viewport(placed, 0, true), damage, deviceScale);
        }
    }

    @Test
    public void knownSnpLoadedThroughPreferenceAndConstructorSuppressesOnlyItsSpike() throws Exception {
        BufferedImage before = paint(track, DEFAULT_RECT, 2, 0, null);
        assertSpikeColor(before, 1040, 2, SequenceRenderer.nucleotideColors.get('c'));
        AlignmentInterval original = interval(track);
        File snps = temporaryFolder.newFile("known-snps.txt");
        Files.writeString(snps.toPath(), CHR + "\t1041\n");
        preferences.put(KNOWN_SNPS, snps.getAbsolutePath());
        // Another real load initializes the global filter without changing this track's counts or mapping.
        AlignmentTrack filterLoader = newTrack();
        closeTrack(filterLoader);
        assertSame("Global filtering must invalidate even an unchanged loaded interval", original, interval(track));
        BufferedImage filtered = assertChangedAndFresh("Known-SNP filtering", before, DEFAULT_RECT, 2);
        assertSpikeColor(filtered, 1040, 2, COVERAGE_COLOR);
        assertSpikeColor(filtered, 1065, 2, SequenceRenderer.nucleotideColors.get('c'));
        // Clearing the preference does not unload the global table in legacy IGV.
        preferences.put(KNOWN_SNPS, null);
        assertPixels("Loaded global filter remains effective after preference removal", filtered,
                paint(track, DEFAULT_RECT, 2, 0, null));
    }

    private File writeBam() throws Exception {
        SAMFileHeader header = new SAMFileHeader();
        header.addSequence(new SAMSequenceRecord(CHR, 10000));
        header.setSortOrder(SAMFileHeader.SortOrder.coordinate);
        File file = temporaryFolder.newFile("coverage-raster.bam");
        try (SAMFileWriter writer = new SAMFileWriterFactory().setCreateIndex(true)
                .makeBAMWriter(header, true, file)) {
            for (int i = 0; i < 2; i++) {
                SAMRecord read = record(header, "span-" + i, 940, 280, 60);
                byte[] bases = read.getReadBases();
                // A lone allele just before a fractional origin exposes X-dependent integer truncation.
                bases[20] = 'C';
                read.setReadBases(bases);
                writer.addAlignment(read);
            }
            for (int i = 0; i < 10; i++) {
                SAMRecord read = record(header, "main-" + i, 1000, 100, i < 4 ? 20 : 60);
                byte[] bases = read.getReadBases();
                byte[] qualities = read.getBaseQualities();
                if (i < 4) {
                    bases[40] = 'C';
                    qualities[40] = 5;
                }
                if (i < 8) bases[65] = 'C';
                read.setReadBases(bases);
                read.setBaseQualities(qualities);
                writer.addAlignment(read);
            }
            for (int i = 0; i < 4; i++) writer.addAlignment(record(header, "late-" + i, 1030, 110, 60));
        }
        return file;
    }

    private SAMRecord record(SAMFileHeader header, String name, int start, int length, int mappingQuality) {
        SAMRecord record = new SAMRecord(header);
        record.setReadName(name);
        record.setReferenceName(CHR);
        record.setAlignmentStart(start + 1);
        record.setCigarString(length + "M");
        record.setMappingQuality(mappingQuality);
        record.setReadString("A".repeat(length));
        byte[] qualities = new byte[length];
        Arrays.fill(qualities, (byte) 35);
        record.setBaseQualities(qualities);
        return record;
    }

    private AlignmentTrack newTrack() throws Exception {
        ResourceLocator locator = new ResourceLocator(bam.getAbsolutePath());
        AlignmentTrack created = new AlignmentTrack(locator, new AlignmentDataManager(locator, genome), genome);
        tracks.add(created);
        created.setDisplayMode(Track.DisplayMode.EXPANDED);
        AlignmentTrack.RenderOptions options = created.getRenderOptions();
        options.setColorOption(AlignmentTrack.ColorOption.NONE);
        options.setShadeAlignmentsOption(AlignmentTrack.ShadeAlignmentsOption.NONE);
        options.setGroupByOption(AlignmentTrack.GroupOption.NONE);
        options.setDuplicatesOption(AlignmentTrack.DuplicatesOption.SHOW);
        options.setShowMismatches(true);
        options.setShowAllBases(false);
        options.setQuickConsensusMode(false);
        options.setShadeBasesOption(false);
        coverage(created).setHeight(60);
        coverage(created).setColor(COVERAGE_COLOR);
        coverage(created).setDataRange(new DataRange(0, 0, 24));
        coverage(created).setSnpThreshold(0.2f);
        created.load(frame);
        assertNotNull("Indexed temporary BAM must load through the real manager", interval(created));
        return created;
    }

    private void closeTrack(AlignmentTrack loaded) {
        loaded.getCoverageTrack().unload();
        loaded.getSpliceJunctionTrack().unload();
        loaded.unload();
        tracks.remove(loaded);
    }

    private CoverageTrack coverage(AlignmentTrack loaded) {
        return loaded.getCoverageTrack();
    }

    private AlignmentInterval interval(AlignmentTrack loaded) {
        return loaded.getDataManager().getLoadedInterval(frame);
    }

    private BufferedImage fresh(Rectangle rect, int scale, int scroll, GraphicsCase graphicsCase,
                                boolean multiframe, boolean wideViewport) throws Exception {
        AlignmentTrack reference = newTrack();
        try {
            CoverageTrack source = coverage(track);
            coverage(reference).setColor(source.getColor());
            coverage(reference).setSnpThreshold(source.getSnpThreshold());
            DataRange range = source.getDataRange();
            coverage(reference).setDataRange(new DataRange(range.getMinimum(), range.getBaseline(),
                    range.getMaximum(), range.isDrawBaseline(), range.isLog()));
            reference.getRenderOptions().setColorOption(track.getRenderOptions().getColorOption());
            assertNotSame("Reference must have independently loaded counts", interval(track).getCounts(),
                    interval(reference).getCounts());
            return paint(reference, rect, scale, scroll, null, graphicsCase, multiframe, wideViewport);
        } finally {
            closeTrack(reference);
        }
    }

    private BufferedImage assertChangedAndFresh(String reason, BufferedImage before, Rectangle rect, int scale)
            throws Exception {
        BufferedImage expected = fresh(rect, scale, 0, GraphicsCase.NORMAL, false, false);
        assertDifferent(reason + " must affect real pixels", before, expected);
        BufferedImage actual = paint(track, rect, scale, 0, null);
        assertPixels(reason + " must match a newly initialized track", expected, actual);
        assertPixels(reason + " warm repaint", expected, paint(track, rect, scale, 0, null));
        return actual;
    }

    private Rectangle viewport(Rectangle rect, int scroll, boolean wide) {
        return new Rectangle(wide ? 0 : rect.x, scroll, wide ? 280 : rect.width, VIEW_HEIGHT);
    }

    private BufferedImage paint(AlignmentTrack loaded, Rectangle rect, int scale, int scroll, Rectangle damage) {
        return paint(loaded, rect, scale, scroll, damage, GraphicsCase.NORMAL, false, false);
    }

    private BufferedImage paint(AlignmentTrack loaded, Rectangle rect, int scale, int scroll, Rectangle damage,
                                GraphicsCase graphicsCase, boolean multiframe, boolean wideViewport) {
        Rectangle viewport = viewport(rect, scroll, wideViewport);
        BufferedImage image = new BufferedImage(viewport.width * scale, VIEW_HEIGHT * scale, BufferedImage.TYPE_INT_ARGB);
        Graphics2D graphics = image.createGraphics();
        graphics.setColor(Color.WHITE);
        graphics.fillRect(0, 0, image.getWidth(), image.getHeight());
        graphics.setRenderingHint(RenderingHints.KEY_ANTIALIASING, RenderingHints.VALUE_ANTIALIAS_OFF);
        graphics.setRenderingHint(RenderingHints.KEY_TEXT_ANTIALIASING, RenderingHints.VALUE_TEXT_ANTIALIAS_OFF);
        graphics.scale(scale, scale);
        graphics.translate(-viewport.x, -scroll);
        if (graphicsCase == GraphicsCase.ANTIALIASED) {
            graphics.setRenderingHint(RenderingHints.KEY_ANTIALIASING, RenderingHints.VALUE_ANTIALIAS_ON);
        } else if (graphicsCase == GraphicsCase.FRACTIONAL_TRANSLATION) {
            graphics.translate(0.25, 0.25);
        } else if (graphicsCase == GraphicsCase.FRACTIONAL_SCALE) {
            graphics.scale(1.25, 1.25);
        } else if (graphicsCase == GraphicsCase.TRANSLUCENT_COMPOSITE) {
            graphics.setComposite(AlphaComposite.getInstance(AlphaComposite.SRC_OVER, 0.5f));
        }
        graphics.setClip(graphicsCase == GraphicsCase.NULL_CLIP ? null :
                damage == null ? viewport : viewport.intersection(damage));
        JPanel panel = new JPanel();
        panel.setBackground(Color.WHITE);
        RenderContext context = new RenderContext(panel, graphics, frame, viewport);
        context.multiframe = multiframe;
        try {
            coverage(loaded).render(context, new Rectangle(rect));
            assertEquals("Damage must not replace the visible viewport", viewport, context.getVisibleRect());
            return image;
        } finally {
            context.dispose();
            graphics.dispose();
        }
    }

    private void assertSpikeColor(BufferedImage image, int position, int scale, Color expected) {
        int x = (int) ((position - frame.getOrigin()) / frame.getScale());
        // The C segment lies above the A segment; y=24 intersects both fixture C spikes.
        assertEquals("Allele pixel at genomic position " + position, expected.getRGB(),
                image.getRGB(x * scale, 24 * scale));
    }

    private void assertDamage(BufferedImage full, BufferedImage partial, Rectangle viewport, Rectangle damage, int scale) {
        Rectangle clipped = viewport.intersection(damage);
        Rectangle deviceDamage = new Rectangle((clipped.x - viewport.x) * scale,
                (clipped.y - viewport.y) * scale, clipped.width * scale, clipped.height * scale);
        assertEquals(full.getWidth(), partial.getWidth());
        assertEquals(full.getHeight(), partial.getHeight());
        for (int y = 0; y < full.getHeight(); y++) {
            for (int x = 0; x < full.getWidth(); x++) {
                assertEquals("Damage pixel " + x + "," + y + " at " + scale + "x in " + damage,
                        deviceDamage.contains(x, y) ? full.getRGB(x, y) : Color.WHITE.getRGB(), partial.getRGB(x, y));
            }
        }
    }

    private void assertPixels(String reason, BufferedImage expected, BufferedImage actual) {
        assertEquals(reason + " image width", expected.getWidth(), actual.getWidth());
        assertEquals(reason + " image height", expected.getHeight(), actual.getHeight());
        int[] expectedPixels = expected.getRGB(0, 0, expected.getWidth(), expected.getHeight(), null, 0, expected.getWidth());
        int[] actualPixels = actual.getRGB(0, 0, actual.getWidth(), actual.getHeight(), null, 0, actual.getWidth());
        assertArrayEquals(reason, expectedPixels, actualPixels);
    }

    private void assertDifferent(String reason, BufferedImage first, BufferedImage second) {
        int[] firstPixels = first.getRGB(0, 0, first.getWidth(), first.getHeight(), null, 0, first.getWidth());
        int[] secondPixels = second.getRGB(0, 0, second.getWidth(), second.getHeight(), null, 0, second.getWidth());
        assertFalse(reason, Arrays.equals(firstPixels, secondPixels));
    }

    private void assertWhite(String reason, BufferedImage image) {
        int[] pixels = image.getRGB(0, 0, image.getWidth(), image.getHeight(), null, 0, image.getWidth());
        for (int pixel : pixels) assertEquals(reason, Color.WHITE.getRGB(), pixel);
    }

    private Map<Character, Color> palette() {
        Map<Character, Color> colors = new HashMap<>();
        colors.put('a', new Color(19, 173, 83));
        colors.put('c', new Color(37, 71, 211));
        colors.put('g', new Color(221, 149, 23));
        colors.put('t', new Color(203, 41, 109));
        colors.put('n', Color.GRAY);
        for (char base : new char[]{'a', 'c', 'g', 't', 'n'}) colors.put(Character.toUpperCase(base), colors.get(base));
        return colors;
    }

    private enum GraphicsCase {
        NORMAL, ANTIALIASED, FRACTIONAL_TRANSLATION, FRACTIONAL_SCALE, TRANSLUCENT_COMPOSITE, NULL_CLIP
    }
}
