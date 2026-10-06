package org.broad.igv.sam;

import htsjdk.samtools.*;
import org.broad.igv.Globals;
import org.broad.igv.feature.genome.Genome;
import org.broad.igv.feature.genome.GenomeManager;
import org.broad.igv.feature.genome.InMemorySequence;
import org.broad.igv.feature.genome.load.GenomeConfig;
import org.broad.igv.prefs.IGVPreferences;
import org.broad.igv.prefs.PreferencesManager;
import org.broad.igv.renderer.SequenceRenderer;
import org.broad.igv.track.RenderContext;
import org.broad.igv.track.Track;
import org.broad.igv.ui.panel.ReferenceFrame;
import org.broad.igv.util.ResourceLocator;
import org.junit.*;
import org.junit.rules.TemporaryFolder;

import java.awt.*;
import java.awt.geom.Rectangle2D;
import java.awt.image.BufferedImage;
import java.io.File;
import java.util.Arrays;
import java.util.HashMap;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;

import static org.broad.igv.prefs.Constants.*;
import static org.junit.Assert.*;

/** Headless paints against the old quality formula and individual Java2D rectangles. */
public class AlignmentRendererQualityColorTest {
    private static final int WIDTH = 1600;
    private static final int START = 1000;
    private static final Color BACKGROUND = new Color(245, 239, 227);
    private static final Rectangle ROW = new Rectangle(0, 2, WIDTH, 6);
    private static final String BASES = "ACGTN?";

    @Rule
    public TemporaryFolder temporaryFolder = new TemporaryFolder();

    private Genome previousGenome;
    private boolean previousHeadless;
    private IGVPreferences preferences;
    private final Map<String, String> previousPreferences = new LinkedHashMap<>();
    private HashMap<Character, Color> previousSamPalette;
    private Map<Character, Color> previousSequencePalette;
    private AlignmentTrack track;
    private AlignmentTrack.RenderOptions options;
    private AlignmentRenderer renderer;
    private SAMFileHeader header;

    @Before
    public void setUp() throws Exception {
        previousHeadless = Globals.isHeadless();
        Globals.setHeadless(true);
        previousGenome = GenomeManager.getInstance().getCurrentGenome();
        byte[] reference = new byte[10000];
        Arrays.fill(reference, (byte) 'A');
        GenomeConfig config = new GenomeConfig();
        config.setId("quality-color-raster-test");
        config.setName("quality-color-raster-test");
        config.setSequence(new InMemorySequence("chr16", reference));
        Genome genome = new Genome(config);
        GenomeManager.getInstance().setCurrentGenomeForTest(genome);

        preferences = PreferencesManager.getPreferences();
        for (String key : List.of(SAM_SHOW_SOFT_CLIPPED, SAM_SHOW_CENTER_LINE, SAM_FLAG_CLIPPING,
                SAM_FLAG_UNMAPPED_PAIR, SAM_FLAG_ZERO_QUALITY, SAM_FLAG_LARGE_INDELS,
                SAM_BASE_QUALITY_MIN, SAM_BASE_QUALITY_MAX,
                SAM_COLOR_A, SAM_COLOR_C, SAM_COLOR_G, SAM_COLOR_T, SAM_COLOR_N)) {
            previousPreferences.put(key, preferences.get(key, null));
        }
        preferences.put(SAM_SHOW_SOFT_CLIPPED, "true");
        preferences.put(SAM_SHOW_CENTER_LINE, "false");
        preferences.put(SAM_FLAG_CLIPPING, "false");
        preferences.put(SAM_FLAG_UNMAPPED_PAIR, "false");
        preferences.put(SAM_FLAG_ZERO_QUALITY, "false");
        preferences.put(SAM_FLAG_LARGE_INDELS, "false");
        thresholds(5, 20);
        previousSamPalette = AlignmentRenderer.nucleotideColors;
        previousSequencePalette = SequenceRenderer.nucleotideColors;

        header = new SAMFileHeader();
        header.addSequence(new SAMSequenceRecord("chr16", reference.length));
        header.setSortOrder(SAMFileHeader.SortOrder.coordinate);
        File bam = temporaryFolder.newFile("quality-colors.bam");
        try (SAMFileWriter writer = new SAMFileWriterFactory().setCreateIndex(true).makeBAMWriter(header, true, bam)) {
            writer.addAlignment(record("fixture", "20M", "A".repeat(20), new byte[20]));
        }
        ResourceLocator locator = new ResourceLocator(bam.getAbsolutePath());
        track = new AlignmentTrack(locator, new AlignmentDataManager(locator, genome), genome);
        track.setDisplayMode(Track.DisplayMode.EXPANDED);
        track.setColor(new Color(92, 117, 133, 173));
        options = track.getRenderOptions();
        options.setColorOption(AlignmentTrack.ColorOption.NONE);
        options.setShadeAlignmentsOption(AlignmentTrack.ShadeAlignmentsOption.NONE);
        options.setShowMismatches(true);
        options.setShowAllBases(true);
        options.setQuickConsensusMode(false);
        options.setShadeBasesOption(true);
        options.setHideSmallIndels(false);
        options.setDuplicatesOption(AlignmentTrack.DuplicatesOption.SHOW);
        renderer = new AlignmentRenderer(track);
        synchronizePalettePreferences();
        AlignmentRenderer.nucleotideColors = customPalette();
        SequenceRenderer.nucleotideColors = insertionPalette();
    }

    @After
    public void tearDown() {
        try {
            if (preferences != null) previousPreferences.forEach(preferences::put);
            // Refresh the renderer's preference-string snapshot before restoring a
            // user's direct palette map; later fixtures must not lose custom edits.
            if (renderer != null) synchronizePalettePreferences();
        } finally {
            if (previousSamPalette != null) AlignmentRenderer.nucleotideColors = previousSamPalette;
            SequenceRenderer.nucleotideColors = previousSequencePalette;
            if (track != null) {
                track.getCoverageTrack().unload();
                track.getSpliceJunctionTrack().unload();
                track.unload();
            }
            GenomeManager.getInstance().setCurrentGenomeForTest(previousGenome);
            Globals.setHeadless(previousHeadless);
        }
    }

    @Test
    public void customTranslucentPaletteAndUnknownFallbackRenderEverySignedQualityWithShadingOffAndOn() {
        Alignment read = qualityMatrix(false);
        for (boolean shade : new boolean[]{false, true, false, true}) {
            options.setShadeBasesOption(shade);
            for (double scale : new double[]{1, 2}) {
                assertRow(read, customPalette(), shade, 5, 20, scale, scale == 1 ? 0 : 2);
            }
        }
        // Soft clips use the same quality lookup but have a distinct batching path.
        options.setShowAllBases(false);
        assertRow(qualityMatrix(true), customPalette(), true, 5, 20, 2, 2);
    }

    @Test
    public void sameRendererObservesPaletteReplacementInPlaceEditsAndLiveThresholds() {
        Alignment read = qualityMatrix(false);
        assertRow(read, customPalette(), true, 5, 20, 1, 0);
        HashMap<Character, Color> replacement = customPalette();
        replacement.put('C', new Color(229, 17, 101, 67));
        AlignmentRenderer.nucleotideColors = replacement;
        assertRow(read, new HashMap<>(replacement), true, 5, 20, 1, 0);
        replacement.put('C', new Color(11, 193, 157, 129));
        replacement.put('A', new Color(11, 193, 157, 0));
        replacement.remove('N');
        assertRow(read, new HashMap<>(replacement), true, 5, 20, 1, 0);
        for (int[] range : new int[][]{{0, 60}, {-20, -5}, {10, 10}, {20, 5}, {-300, 300},
                {Integer.MIN_VALUE, Integer.MAX_VALUE}, {5, 20}}) {
            thresholds(range[0], range[1]);
            assertRow(read, new HashMap<>(replacement), true, range[0], range[1], 1, 0);
        }
        options.setShadeBasesOption(false);
        assertRow(read, new HashMap<>(replacement), false, 5, 20, 1, 0);
        options.setShadeBasesOption(true);
        assertRow(read, new HashMap<>(replacement), true, 5, 20, 1, 0);
    }

    @Test
    public void changedSamPreferenceStringsRefreshNextRowWithoutReplacingUnchangedCustomEdits() {
        Alignment read = new SAMAlignment(record("preferences", "256M", "Cc".repeat(128), signedQualities(256)));
        HashMap<Character, Color> expected = customPalette();
        expected.put('c', expected.get('C'));
        AlignmentRenderer.nucleotideColors = new HashMap<>(expected);
        assertRow(read, expected, true, 5, 20, 1, 0);
        preferences.put(SAM_COLOR_A, "19,173,83");
        preferences.put(SAM_COLOR_C, "227,31,113");
        preferences.put(SAM_COLOR_G, "221,149,23");
        preferences.put(SAM_COLOR_T, "203,41,109");
        preferences.put(SAM_COLOR_N, "109,53,139");
        expected.put('C', new Color(227, 31, 113));
        expected.put('c', new Color(227, 31, 113));
        assertRow(read, expected, true, 5, 20, 1, 0);
        AlignmentRenderer.nucleotideColors.put('C', new Color(17, 199, 71, 49));
        AlignmentRenderer.nucleotideColors.put('c', new Color(17, 199, 71, 121));
        expected.put('C', new Color(17, 199, 71, 49));
        expected.put('c', new Color(17, 199, 71, 121));
        assertRow(read, expected, true, 5, 20, 1, 0);
        preferences.put(SAM_COLOR_C, "41,67,239");
        expected.put('C', new Color(41, 67, 239));
        expected.put('c', new Color(41, 67, 239));
        assertRow(read, expected, true, 5, 20, 1, 0);
    }

    @Test
    public void bisulfiteAndGraySpecializedCallersRetainQualityShading() {
        Alignment read = new SAMAlignment(record("specialized", "256M", "C".repeat(256), signedQualities(256)));
        HashMap<Character, Color> expected = new HashMap<>();
        for (AlignmentTrack.ColorOption mode : new AlignmentTrack.ColorOption[]{
                AlignmentTrack.ColorOption.BISULFITE, AlignmentTrack.ColorOption.BASE_MODIFICATION,
                AlignmentTrack.ColorOption.SMRT_SUBREAD_IPD}) {
            options.setColorOption(mode);
            // C against reference A is a methylated-color CHARACTER, with no
            // context-dependent rectangle expansion and no modification tags.
            expected.put('C', mode == AlignmentTrack.ColorOption.BISULFITE ? Color.RED : Color.GRAY);
            for (boolean shade : new boolean[]{false, true}) {
                options.setShadeBasesOption(shade);
                assertRow(read, expected, shade, 5, 20, 1, 0);
            }
            thresholds(0, 60);
            assertRow(read, expected, true, 0, 60, 1, 0);
            thresholds(5, 20);
        }
    }

    @Test
    public void expandedInsertionsUseDistinctSequencePaletteAndSnapshotEachInvocation() {
        String inserted = matrixBases();
        byte[] qualities = new byte[inserted.length() + 2];
        System.arraycopy(signedQualities(inserted.length()), 0, qualities, 1, inserted.length());
        Alignment read = new SAMAlignment(record("insertion", "1M" + inserted.length() + "I1M",
                "A" + inserted + "A", qualities));
        InsertionMarker marker = new InsertionMarker(START + 1, inserted.length());
        Map<Character, Color> expected = insertionPalette();
        for (boolean shade : new boolean[]{false, true}) {
            options.setShadeBasesOption(shade);
            assertInsertion(read, marker, expected, shade, 5, 20);
        }
        // SAM palette replacement must not leak into expanded insertions.
        AlignmentRenderer.nucleotideColors = new HashMap<>();
        assertInsertion(read, marker, expected, true, 5, 20);
        Map<Character, Color> replacement = insertionPalette();
        replacement.put('G', new Color(241, 7, 113, 53));
        SequenceRenderer.nucleotideColors = replacement;
        assertInsertion(read, marker, new HashMap<>(replacement), true, 5, 20);
        replacement.put('G', new Color(43, 181, 71, 121));
        replacement.remove('N');
        assertInsertion(read, marker, new HashMap<>(replacement), true, 5, 20);
        thresholds(0, 60);
        assertInsertion(read, marker, new HashMap<>(replacement), true, 0, 60);
        options.setShadeBasesOption(false);
        assertInsertion(read, marker, new HashMap<>(replacement), false, 0, 60);
        options.setColorOption(AlignmentTrack.ColorOption.BASE_MODIFICATION);
        options.setShadeBasesOption(true);
        HashMap<Character, Color> gray = new HashMap<>();
        for (char base : BASES.toCharArray()) gray.put(base, Color.GRAY);
        assertInsertion(read, marker, gray, true, 0, 60);
    }

    private void assertRow(Alignment read, Map<Character, Color> palette, boolean shade,
                           int min, int max, double scale, int tolerance) {
        try (Raster actual = new Raster(START, scale); Raster expected = new Raster(START, scale)) {
            renderer.renderAlignments(List.of(read), null, actual.context, ROW, options);
            boolean showAll = options.isShowAllBases();
            boolean showMismatches = options.isShowMismatches();
            try {
                options.setShowAllBases(false);
                options.setShowMismatches(false);
                renderer.renderAlignments(List.of(read), null, expected.context, ROW, options);
            } finally {
                options.setShowAllBases(showAll);
                options.setShowMismatches(showMismatches);
            }
            Graphics2D paint = (Graphics2D) expected.graphics.create();
            try {
                paint.setComposite(AlphaComposite.getInstance(AlphaComposite.SRC_OVER, 0.75f));
                for (AlignmentBlock block : read.getAlignmentBlocks()) {
                    if (!showAll && !block.isSoftClip() && options.getColorOption() == AlignmentTrack.ColorOption.NONE) continue;
                    for (int i = 0; i < block.getLength(); i++) {
                        char base = (char) block.getBases().getByte(i);
                        Color color = palette.getOrDefault(base, Color.BLACK);
                        if (shade) color = BaseRendererQualityColorTest.legacyColor(color, block.getQuality(i), min, max);
                        int x = (int) ((block.getStart() + i - START) / scale);
                        paint.setColor(color);
                        paint.fillRect(x, ROW.y, (int) Math.max(1, 1 / scale), ROW.height - 2);
                    }
                }
            } finally {
                paint.dispose();
            }
            assertPixels(expected.image, actual.image, tolerance);
        }
    }

    private void assertInsertion(Alignment read, InsertionMarker marker, Map<Character, Color> palette,
                                 boolean shade, int min, int max) {
        try (Raster actual = new Raster(marker.position, 1); Raster expected = new Raster(marker.position, 1)) {
            BaseRenderer.drawExpandedInsertions(marker, List.of(read), actual.context, ROW, true, options);
            AlignmentBlock insertion = read.getInsertionAt(marker.position);
            assertNotNull("Fixture must contain an expanded insertion", insertion);
            Graphics2D paint = (Graphics2D) expected.graphics.create();
            try {
                paint.setRenderingHint(RenderingHints.KEY_ANTIALIASING, RenderingHints.VALUE_ANTIALIAS_ON);
                for (int i = 0; i < insertion.getBasesLength(); i++) {
                    char base = (char) insertion.getBases().getByte(i);
                    Color color = palette.getOrDefault(base, Color.BLACK);
                    if (shade) color = BaseRendererQualityColorTest.legacyColor(color, insertion.getQuality(i), min, max);
                    paint.setColor(color);
                    paint.fill(new Rectangle2D.Double(i, ROW.y, 1, ROW.height));
                }
            } finally {
                paint.dispose();
            }
            assertPixels(expected.image, actual.image, 0);
        }
    }

    private void synchronizePalettePreferences() {
        try (Raster raster = new Raster(START, 1)) {
            renderer.renderAlignments(List.of(), null, raster.context, ROW, options);
        }
    }

    private void thresholds(int min, int max) {
        preferences.put(SAM_BASE_QUALITY_MIN, Integer.toString(min));
        preferences.put(SAM_BASE_QUALITY_MAX, Integer.toString(max));
    }

    private Alignment qualityMatrix(boolean softClip) {
        String bases = matrixBases();
        SAMRecord record = record("matrix", softClip ? bases.length() + "S1M" : bases.length() + "M",
                softClip ? bases + "A" : bases, signedQualities(bases.length() + (softClip ? 1 : 0)));
        // SAMAlignment extends the alignment start leftward by the clip length.
        if (softClip) record.setAlignmentStart(START + bases.length() + 1);
        return new SAMAlignment(record);
    }

    private static String matrixBases() {
        StringBuilder bases = new StringBuilder();
        for (char base : BASES.toCharArray()) bases.append(String.valueOf(base).repeat(256));
        return bases.toString();
    }

    private static byte[] signedQualities(int length) {
        byte[] qualities = new byte[length];
        for (int i = 0; i < length; i++) qualities[i] = (byte) i;
        return qualities;
    }

    private SAMRecord record(String name, String cigar, String bases, byte[] qualities) {
        SAMRecord record = new SAMRecord(header);
        record.setReadName(name);
        record.setReferenceName("chr16");
        record.setAlignmentStart(START + 1);
        record.setMappingQuality(60);
        record.setCigarString(cigar);
        record.setReadString(bases);
        record.setBaseQualities(qualities);
        return record;
    }

    private static HashMap<Character, Color> customPalette() {
        HashMap<Character, Color> palette = new HashMap<>();
        palette.put('A', new Color(19, 173, 83, 43));
        palette.put('C', new Color(37, 71, 211, 0));
        palette.put('G', new Color(37, 71, 211, 91));
        palette.put('T', new Color(203, 41, 109, 203));
        palette.put('N', new Color(109, 53, 139));
        return palette;
    }

    private static HashMap<Character, Color> insertionPalette() {
        HashMap<Character, Color> palette = customPalette();
        palette.put('A', new Color(211, 67, 31, 71));
        palette.put('C', new Color(113, 199, 47, 0));
        palette.put('G', new Color(113, 199, 47, 139));
        return palette;
    }

    private static void assertPixels(BufferedImage expected, BufferedImage actual, int tolerance) {
        for (int y = 0; y < expected.getHeight(); y++) {
            for (int x = 0; x < expected.getWidth(); x++) {
                int reference = expected.getRGB(x, y);
                int observed = actual.getRGB(x, y);
                String pixel = "Pixel " + x + "," + y;
                assertEquals(pixel + " alpha", reference >>> 24, observed >>> 24);
                assertEquals(pixel + " coverage", reference != BACKGROUND.getRGB(), observed != BACKGROUND.getRGB());
                if (tolerance == 0 || reference == BACKGROUND.getRGB()) {
                    assertEquals(pixel, reference, observed);
                } else {
                    for (int shift : new int[]{16, 8, 0}) {
                        assertTrue(pixel + " channel " + shift,
                                Math.abs(((reference >>> shift) & 255) - ((observed >>> shift) & 255)) <= tolerance);
                    }
                }
            }
        }
    }

    private static final class Raster implements AutoCloseable {
        final BufferedImage image = new BufferedImage(WIDTH, 12, BufferedImage.TYPE_INT_ARGB_PRE);
        final Graphics2D graphics = image.createGraphics();
        final RenderContext context;

        Raster(double origin, double scale) {
            graphics.setColor(BACKGROUND);
            graphics.fillRect(0, 0, WIDTH, image.getHeight());
            graphics.setRenderingHint(RenderingHints.KEY_ANTIALIASING, RenderingHints.VALUE_ANTIALIAS_OFF);
            graphics.clip(new Rectangle(0, 0, WIDTH, image.getHeight()));
            ReferenceFrame frame = new ReferenceFrame("quality-color-raster-test") {
                @Override
                public double getScale() {
                    return scale;
                }
            };
            frame.setBounds(0, WIDTH);
            frame.jumpTo("chr16", (int) origin, (int) (origin + WIDTH * scale));
            frame.setOrigin(origin);
            context = new RenderContext(null, graphics, frame, new Rectangle(0, 0, WIDTH, image.getHeight()));
        }

        @Override
        public void close() {
            context.dispose();
            graphics.dispose();
        }
    }
}
