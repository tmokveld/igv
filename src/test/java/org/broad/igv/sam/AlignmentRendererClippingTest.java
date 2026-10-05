package org.broad.igv.sam;

import htsjdk.samtools.*;
import org.apache.batik.dom.GenericDOMImplementation;
import org.apache.batik.svggen.SVGGraphics2D;
import org.broad.igv.Globals;
import org.broad.igv.feature.Chromosome;
import org.broad.igv.feature.genome.Genome;
import org.broad.igv.feature.genome.GenomeManager;
import org.broad.igv.prefs.IGVPreferences;
import org.broad.igv.prefs.PreferencesManager;
import org.broad.igv.track.RenderContext;
import org.broad.igv.track.Track;
import org.broad.igv.ui.panel.ReferenceFrame;
import org.broad.igv.util.ResourceLocator;
import org.junit.*;
import org.junit.rules.TemporaryFolder;
import org.w3c.dom.Document;
import org.w3c.dom.Element;
import org.w3c.dom.NodeList;

import java.awt.*;
import java.awt.image.BufferedImage;
import java.io.File;
import java.util.Arrays;
import java.util.HashMap;
import java.util.List;
import java.util.LinkedHashMap;
import java.util.Map;

import static org.broad.igv.prefs.Constants.*;
import static org.junit.Assert.*;

public class AlignmentRendererClippingTest {
    @Rule
    public TemporaryFolder temporaryFolder = new TemporaryFolder();

    private Genome previousGenome;
    private boolean previousHeadless;
    private IGVPreferences preferences;
    private final Map<String, String> previousPreferences = new LinkedHashMap<>();
    private HashMap<Character, Color> previousNucleotideColors;
    private Map<String, Color> previousShadedColors;
    private AlignmentTrack track;
    private SAMFileHeader header;

    @Before
    public void setUp() throws Exception {
        previousHeadless = Globals.isHeadless();
        Globals.setHeadless(true);
        previousGenome = GenomeManager.getInstance().getCurrentGenome();
        Genome genome = new Genome("clipping-test", List.of(new Chromosome(0, "chr16", 50000)));
        GenomeManager.getInstance().setCurrentGenomeForTest(genome);
        preferences = PreferencesManager.getPreferences();
        for (String key : List.of(SAM_SHOW_SOFT_CLIPPED, SAM_SHOW_CENTER_LINE,
                SAM_BASE_QUALITY_MIN, SAM_BASE_QUALITY_MAX)) {
            previousPreferences.put(key, preferences.get(key, null));
        }
        preferences.put(SAM_SHOW_SOFT_CLIPPED, "true");
        preferences.put(SAM_SHOW_CENTER_LINE, "false");
        preferences.put(SAM_BASE_QUALITY_MIN, "5");
        preferences.put(SAM_BASE_QUALITY_MAX, "20");
        previousNucleotideColors = AlignmentRenderer.nucleotideColors;
        AlignmentRenderer.nucleotideColors = new HashMap<>(previousNucleotideColors);
        AlignmentRenderer.nucleotideColors.put('A', Color.GREEN);
        AlignmentRenderer.nucleotideColors.put('C', Color.BLUE);
        AlignmentRenderer.nucleotideColors.put('G', Color.ORANGE);
        AlignmentRenderer.nucleotideColors.put('T', Color.RED);
        synchronized (BaseRenderer.shadedColorCache) {
            previousShadedColors = new HashMap<>(BaseRenderer.shadedColorCache);
            BaseRenderer.shadedColorCache.clear();
        }

        header = new SAMFileHeader();
        header.addSequence(new SAMSequenceRecord("chr16", 50000));
        header.setSortOrder(SAMFileHeader.SortOrder.coordinate);
        File bam = temporaryFolder.newFile("clipped.bam");
        try (SAMFileWriter writer = new SAMFileWriterFactory().setCreateIndex(true).makeBAMWriter(header, true, bam)) {
            writer.addAlignment(record("5000S20M5000S", "A".repeat(10020)));
        }
        ResourceLocator locator = new ResourceLocator(bam.getAbsolutePath());
        track = new AlignmentTrack(locator, new AlignmentDataManager(locator, genome), genome);
        track.setDisplayMode(Track.DisplayMode.EXPANDED);
        track.setColor(Color.GRAY);
        AlignmentTrack.RenderOptions options = track.getRenderOptions();
        options.setColorOption(AlignmentTrack.ColorOption.NONE);
        options.setShadeAlignmentsOption(AlignmentTrack.ShadeAlignmentsOption.NONE);
        options.setShowMismatches(true);
        options.setShowAllBases(false);
        options.setQuickConsensusMode(false);
        options.setShadeBasesOption(false);
        options.setHideSmallIndels(false);
    }

    @After
    public void tearDown() {
        if (track != null) {
            track.getCoverageTrack().unload();
            track.getSpliceJunctionTrack().unload();
            track.unload();
        }
        if (preferences != null) previousPreferences.forEach(preferences::put);
        if (previousNucleotideColors != null) AlignmentRenderer.nucleotideColors = previousNucleotideColors;
        if (previousShadedColors != null) {
            synchronized (BaseRenderer.shadedColorCache) {
                BaseRenderer.shadedColorCache.clear();
                BaseRenderer.shadedColorCache.putAll(previousShadedColors);
            }
        }
        GenomeManager.getInstance().setCurrentGenomeForTest(previousGenome);
        Globals.setHeadless(previousHeadless);
    }

    @Test
    public void zoomedOutHomopolymerClipsMergeIntoVectorRuns() {
        SAMAlignment alignment = new SAMAlignment(record("5000S20M5000S", "A".repeat(10020)));
        Document document = GenericDOMImplementation.getDOMImplementation()
                .createDocument("http://www.w3.org/2000/svg", "svg", null);
        SVGGraphics2D graphics = new SVGGraphics2D(document);
        try {
            render(alignment, graphics, 0, 20000, 400);
            Element root = graphics.getRoot();
            NodeList rectangles = root.getElementsByTagName("rect");
            assertTrue("Homopolymer columns should be batched into vector runs",
                    rectangles.getLength() < 20);
            assertTrue("Both clip ends must remain visible", rectangles.getLength() >= 2);
            assertEquals("SVG must not embed a raster strip", 0, root.getElementsByTagName("image").getLength());
        } finally {
            graphics.dispose();
        }
    }

    @Test
    public void bothClippedEndsStayVisibleWithDetailedColorsWhenZoomedIn() {
        SAMAlignment alignment = new SAMAlignment(record("40S20M40S",
                "ACGT".repeat(10) + "A".repeat(20) + "TGCA".repeat(10)));
        BufferedImage overview = image(alignment, 9920, 10080, 80);
        BufferedImage detailed = image(alignment, 9920, 10080, 160);
        BufferedImage reference = originalSoftClipImage(alignment, 9920, 10080, 80, 1);
        for (int x = 22; x < 39; x++) {
            assertChannelsWithinTwo("Left clip column " + x, reference.getRGB(x, 5), overview.getRGB(x, 5));
        }
        for (int x = 52; x < 68; x++) {
            assertChannelsWithinTwo("Right clip column " + x, reference.getRGB(x, 5), overview.getRGB(x, 5));
        }
        assertNotEquals(detailed.getRGB(48, 5), detailed.getRGB(49, 5));
        assertNotEquals(detailed.getRGB(49, 5), detailed.getRGB(50, 5));
        assertNotEquals(detailed.getRGB(50, 5), detailed.getRGB(51, 5));
    }

    @Test
    public void overviewBlendsDrawableBasesWithoutPaintingTrailingEquals() {
        SAMRecord record = record("40S20M", "A=C=".repeat(10) + "A".repeat(20));
        byte[] qualities = record.getBaseQualities();
        for (int i = 2; i < 40; i += 4) qualities[i] = 0;
        record.setBaseQualities(qualities);
        SAMAlignment alignment = new SAMAlignment(record);
        track.getRenderOptions().setShadeBasesOption(true);
        BufferedImage overview = image(alignment, 9920, 10080, 40);
        BufferedImage reference = originalSoftClipImage(alignment, 9920, 10080, 40, 1);
        for (int x = 12; x < 19; x++) {
            assertChannelsWithinTwo("Drawable base quality, column " + x,
                    reference.getRGB(x, 5), overview.getRGB(x, 5));
        }
        track.getRenderOptions().setShadeBasesOption(false);
        assertNotEquals(image(alignment, 9920, 10080, 40).getRGB(18, 5), overview.getRGB(18, 5));
    }

    @Test
    public void equalsOnlyColumnsRespectShowAllBases() {
        SAMAlignment alignment = new SAMAlignment(record("40S20M", "=".repeat(40) + "A".repeat(20)));
        BufferedImage mismatchesOnly = image(alignment, 9920, 10080, 40);
        track.getRenderOptions().setShowMismatches(false);
        assertEquals(image(alignment, 9920, 10080, 40).getRGB(18, 5), mismatchesOnly.getRGB(18, 5));
        track.getRenderOptions().setShowAllBases(true);
        BufferedImage allBases = image(alignment, 9920, 10080, 40);
        assertNotEquals(mismatchesOnly.getRGB(18, 5), allBases.getRGB(18, 5));
        assertChannelsWithinTwo("Equals bases use the ordinary fallback color",
                originalSoftClipImage(alignment, 9920, 10080, 40, 1).getRGB(18, 5), allBases.getRGB(18, 5));
    }

    @Test
    public void partialColumnsAtViewportEdgesRetainClippedBases() {
        SAMAlignment alignment = new SAMAlignment(record("80S20M", "C".repeat(80) + "A".repeat(20)));
        BufferedImage overview = image(alignment, 9933, 9976, 20);
        BufferedImage reference = originalSoftClipImage(alignment, 9933, 9976, 20, 1);
        for (int x = 0; x < 20; x++) {
            assertChannelsWithinTwo("Clipped column " + x, reference.getRGB(x, 5), overview.getRGB(x, 5));
        }
    }


    @Test
    public void scaledRasterClipsKeepSharpColumnsAndOrderedColors() {
        SAMAlignment alignment = new SAMAlignment(record("40S20M40S",
                "AACCGGTT".repeat(5) + "A".repeat(20) + "TTGGCCAA".repeat(5)));
        for (int scale : new int[]{1, 2}) {
            BufferedImage overview = image(alignment, 9920, 10080, 80, scale);
            BufferedImage reference = originalSoftClipImage(alignment, 9920, 10080, 80, scale);
            for (int x = 22; x < 68; x++) {
                if (x >= 39 && x < 52) continue; // Aligned block and end indicators.
                int expected = reference.getRGB(scale * x, 5 * scale);
                for (int dx = 0; dx < scale; dx++) {
                    for (int y = 2 * scale; y < 8 * scale; y++) {
                        assertChannelsWithinTwo("Scale " + scale + ", column " + x + ", row " + y,
                                expected, overview.getRGB(scale * x + dx, y));
                    }
                }
            }
            assertEquals(Color.WHITE.getRGB(), overview.getRGB(5 * scale, 5 * scale));
        }
    }

    @Test
    public void shadedClipsBlendMixedColorsAndQualitiesInReadOrderAtBothScales() {
        SAMRecord record = record("48S20M", "TAC=TAG=TAT=".repeat(4) + "A".repeat(20));
        byte[] qualities = record.getBaseQualities();
        for (int i = 2; i < 48; i += 4) {
            qualities[i - 2] = 0;
            qualities[i - 1] = 10;
            qualities[i] = (byte) ((i / 4) % 3 == 0 ? 0 : (i / 4) % 3 == 1 ? 10 : 30);
            qualities[i + 1] = 40; // A trailing '=' contributes no layer.
        }
        record.setBaseQualities(qualities);
        SAMAlignment alignment = new SAMAlignment(record);
        track.getRenderOptions().setShadeBasesOption(true);
        BufferedImage detailed = image(alignment, 9920, 10080, 160);
        for (int scale : new int[]{1, 2}) {
            BufferedImage overview = image(alignment, 9920, 10080, 40, scale);
            BufferedImage reference = originalSoftClipImage(alignment, 9920, 10080, 40, scale);
            for (int x = 10; x < 19; x++) {
                for (int dx = 0; dx < scale; dx++) {
                    String message = "Scale " + scale + ", shaded column " + x;
                    int expected = reference.getRGB(scale * x + dx, 5 * scale);
                    int actual = overview.getRGB(scale * x + dx, 5 * scale);
                    assertChannelsWithinTwo(message, expected, actual);
                }
            }
        }
        track.getRenderOptions().setShadeBasesOption(false);
        BufferedImage unshaded = image(alignment, 9920, 10080, 160);
        assertNotEquals("Low quality remains translucent", unshaded.getRGB(46, 5), detailed.getRGB(46, 5));
        assertNotEquals("Intermediate quality remains translucent", unshaded.getRGB(50, 5), detailed.getRGB(50, 5));
        assertEquals("High quality retains the opaque base color", unshaded.getRGB(54, 5), detailed.getRGB(54, 5));
    }

    @Test
    public void lowQualityClipsAccumulateEveryDrawableBaseAtBothScales() {
        SAMRecord record = record("80S20M", "C".repeat(80) + "A".repeat(20));
        byte[] qualities = record.getBaseQualities();
        Arrays.fill(qualities, 0, 80, (byte) 0);
        record.setBaseQualities(qualities);
        SAMAlignment alignment = new SAMAlignment(record);
        track.getRenderOptions().setShadeBasesOption(true);
        int singleBaseColor = image(alignment, 9920, 10080, 160).getRGB(32, 5);
        for (int scale : new int[]{1, 2}) {
            BufferedImage actual = image(alignment, 9920, 10080, 40, scale);
            BufferedImage reference = originalSoftClipImage(alignment, 9920, 10080, 40, scale);
            assertNotEquals("Several faint bases must not collapse to one faint base",
                    singleBaseColor, reference.getRGB(8 * scale, 5 * scale));
            for (int x = 2; x < 19; x++) {
                for (int dx = 0; dx < scale; dx++) {
                    assertChannelsWithinTwo("Low quality, scale " + scale + ", column " + x,
                            reference.getRGB(x * scale + dx, 5 * scale),
                            actual.getRGB(x * scale + dx, 5 * scale));
                }
            }
        }
    }

    @Test
    public void equalsOnlyHolesBetweenColoredColumnsStayUnpainted() {
        SAMAlignment alignment = new SAMAlignment(record("48S20M",
                "AA==CC==GG==TT==".repeat(3) + "A".repeat(20)));
        BufferedImage overview = image(alignment, 9920, 10080, 80);
        BufferedImage reference = originalSoftClipImage(alignment, 9920, 10080, 80, 1);
        SAMAlignment blankClip = new SAMAlignment(record("48S20M", "=".repeat(48) + "A".repeat(20)));
        BufferedImage blank = image(blankClip, 9920, 10080, 80);
        for (int x = 18; x < 38; x++) {
            assertChannelsWithinTwo("Column " + x, reference.getRGB(x, 5), overview.getRGB(x, 5));
            if (x % 2 == 1) {
                assertEquals("Equals-only hole " + x, blank.getRGB(x, 5), overview.getRGB(x, 5));
            } else {
                assertNotEquals("Colored column " + x, blank.getRGB(x, 5), overview.getRGB(x, 5));
            }
        }
    }

    @Test
    public void reusedScratchDoesNotLeakColorsAcrossBlocksReadsOrWidths() {
        List<Alignment> alignments = List.of(
                new SAMAlignment(record("80S20M80S", "A".repeat(80) + "A".repeat(20) + "C".repeat(80))),
                new SAMAlignment(record("20S20M60S", "==GG".repeat(5) + "A".repeat(20) + "==TT".repeat(15))),
                new SAMAlignment(record("40S20M20S", "CC==".repeat(10) + "A".repeat(20) + "GG==".repeat(5))),
                new SAMAlignment(record("100S20M100S", "T".repeat(100) + "A".repeat(20) + "G".repeat(100))),
                new SAMAlignment(record("20S20M60S", "==GG".repeat(5) + "A".repeat(20) + "==TT".repeat(15))));
        BufferedImage overview = new BufferedImage(120, 12 * alignments.size(), BufferedImage.TYPE_INT_RGB);
        Graphics2D graphics = overview.createGraphics();
        ReferenceFrame frame = frame(9920, 10160, 120);
        RenderContext context = null;
        try {
            graphics.setColor(Color.WHITE);
            graphics.fillRect(0, 0, overview.getWidth(), overview.getHeight());
            context = new RenderContext(null, graphics, frame, new Rectangle(0, 0, 120, overview.getHeight()));
            AlignmentRenderer renderer = new AlignmentRenderer(track);
            for (int row = 0; row < alignments.size(); row++) {
                renderer.renderAlignments(List.of(alignments.get(row)), null, context,
                        new Rectangle(0, 12 * row, 120, 12), track.getRenderOptions());
            }
        } finally {
            if (context != null) context.dispose();
            graphics.dispose();
        }
        int[][] clipRanges = {{2, 38, 52, 88}, {32, 38, 52, 78}, {22, 38, 52, 58},
                {2, 38, 52, 98}, {32, 38, 52, 78}};
        for (int row = 0; row < alignments.size(); row++) {
            BufferedImage reference = originalSoftClipImage(alignments.get(row), 9920, 10160, 120, 1);
            int[] ranges = clipRanges[row];
            for (int range = 0; range < ranges.length; range += 2) {
                for (int x = ranges[range]; x < ranges[range + 1]; x++) {
                    assertChannelsWithinTwo("Read " + row + ", column " + x,
                            reference.getRGB(x, 5), overview.getRGB(x, 12 * row + 5));
                }
            }
        }
        for (int row : new int[]{1, 2, 4}) {
            assertEquals("No old wide strip beyond read " + row, Color.WHITE.getRGB(),
                    overview.getRGB(110, 12 * row + 5));
        }
    }

    @Test
    public void fractionalOriginAndViewportClipPreserveBothEdgeColumns() {
        SAMAlignment alignment = new SAMAlignment(record("80S20M", "C".repeat(80) + "A".repeat(20)));
        for (int scale : new int[]{1, 2}) {
            BufferedImage image = new BufferedImage(20 * scale, 12 * scale, BufferedImage.TYPE_INT_RGB);
            Graphics2D graphics = image.createGraphics();
            BufferedImage reference;
            RenderContext context = null;
            try {
                graphics.setColor(Color.WHITE);
                graphics.fillRect(0, 0, image.getWidth(), image.getHeight());
                graphics.scale(scale, scale);
                graphics.setClip(new Rectangle(3, 2, 14, 7));
                ReferenceFrame frame = frame(9933, 9976, 20);
                frame.setOrigin(9933.5);
                context = new RenderContext(null, graphics, frame, new Rectangle(3, 0, 14, 12));
                new AlignmentRenderer(track).renderAlignments(List.of(alignment), null, context,
                        new Rectangle(3, 0, 14, 12), track.getRenderOptions());
                reference = originalSoftClipImage(alignment, 9933.5, 9976, 20, scale);
            } finally {
                if (context != null) context.dispose();
                graphics.dispose();
            }
            for (int x = 0; x < image.getWidth(); x++) {
                assertChannelsWithinTwo("Scale " + scale + ", viewport column " + x,
                        x >= 3 * scale && x < 17 * scale ? reference.getRGB(x, 5 * scale) : Color.WHITE.getRGB(),
                        image.getRGB(x, 5 * scale));
                assertEquals("Above viewport clip", Color.WHITE.getRGB(), image.getRGB(x, scale));
                assertEquals("Below viewport clip", Color.WHITE.getRGB(), image.getRGB(x, 10 * scale));
            }
        }
    }

    @Test
    public void alignedBasesStillCompositeEveryBaseWithTheOriginalGraphicsComposite() {
        track.getRenderOptions().setShowAllBases(true);
        SAMAlignment alignment = new SAMAlignment(record("80M", "AC".repeat(40)));
        BufferedImage overview = image(alignment, 9920, 10080, 80);
        BufferedImage reference = new BufferedImage(1, 1, BufferedImage.TYPE_INT_RGB);
        Graphics2D graphics = reference.createGraphics();
        try {
            graphics.setColor(Color.WHITE);
            graphics.fillRect(0, 0, 1, 1);
            graphics.setComposite(AlphaComposite.getInstance(AlphaComposite.SRC_OVER, 0.75f));
            // Aligned bases retain the original A-then-C overdraw.
            for (Color color : List.of(track.getColor(), AlignmentRenderer.nucleotideColors.get('A'),
                    AlignmentRenderer.nucleotideColors.get('C'))) {
                graphics.setColor(color);
                graphics.fillRect(0, 0, 1, 1);
            }
        } finally {
            graphics.dispose();
        }
        for (int x = 43; x < 75; x++) {
            assertEquals("Aligned column " + x, reference.getRGB(0, 0), overview.getRGB(x, 5));
        }
    }

    @Test
    public void specializedModesDoNotPaintOrdinaryNucleotideSoftClips() {
        SAMAlignment alignment = new SAMAlignment(record("40S20M40S",
                "ACGT".repeat(10) + "A".repeat(20) + "TGCA".repeat(10)));
        for (AlignmentTrack.ColorOption mode : List.of(AlignmentTrack.ColorOption.BISULFITE,
                AlignmentTrack.ColorOption.NOMESEQ)) {
            track.getRenderOptions().setColorOption(mode);
            track.getRenderOptions().setShowMismatches(false);
            track.getRenderOptions().setShowAllBases(false);
            BufferedImage background = image(alignment, 9920, 10080, 80);
            track.getRenderOptions().setShowMismatches(true);
            track.getRenderOptions().setShowAllBases(true);
            BufferedImage specialized = image(alignment, 9920, 10080, 80);
            for (int x = 22; x < 68; x++) {
                if (x >= 39 && x < 52) continue;
                assertEquals(mode + " soft-clip column " + x,
                        background.getRGB(x, 5), specialized.getRGB(x, 5));
            }
        }
    }

    @Test
    public void stripAppliesExtraAlphaToEveryLayerAndAccumulatesLowAlphaAtBothScales() {
        Color[][] columns = blendingColumns();
        Color background = new Color(37, 79, 113);
        for (float extraAlpha : new float[]{0.5f, 0.75f}) {
            for (int scale : new int[]{1, 2}) {
                BufferedImage actual = stripImage(columns, background, extraAlpha, scale, true);
                BufferedImage reference = stripImage(columns, background, extraAlpha, scale, false);
                assertNotEquals("Read order changes mixed-color output",
                        reference.getRGB(scale, 5 * scale), reference.getRGB(3 * scale, 5 * scale));
                assertNotEquals("Low-alpha layers accumulate visibly", background.getRGB(),
                        reference.getRGB(5 * scale, 5 * scale));
                for (int y = 0; y < actual.getHeight(); y++) {
                    for (int x = 0; x < actual.getWidth(); x++) {
                        assertChannelsWithinTwo("Extra alpha " + extraAlpha + ", scale " + scale +
                                ", pixel " + x + "," + y, reference.getRGB(x, y), actual.getRGB(x, y));
                    }
                }
            }
        }
    }

    @Test
    public void vectorStripColorsCompositeLikeOrderedRasterLayersWithoutEmbeddedImages() {
        Color[][] columns = blendingColumns();
        Color background = new Color(37, 79, 113);
        for (int scale : new int[]{1, 2}) {
            Document document = GenericDOMImplementation.getDOMImplementation()
                    .createDocument("http://www.w3.org/2000/svg", "svg", null);
            SVGGraphics2D graphics = new SVGGraphics2D(document);
            try {
                graphics.scale(scale, scale);
                graphics.setComposite(AlphaComposite.getInstance(AlphaComposite.SRC_OVER, 0.5f));
                BaseRenderer.ColorStrip strip = new BaseRenderer.ColorStrip();
                strip.reset(graphics, 0, columns.length);
                blendColumns(strip, columns);
                strip.draw(2, 6);
                Element root = graphics.getRoot();
                assertEquals("SVG must not embed a raster strip", 0, root.getElementsByTagName("image").getLength());
                BufferedImage raster = stripImage(columns, background, 0.5f, scale, true);
                BufferedImage reference = stripImage(columns, background, 0.5f, scale, false);
                for (int x = 0; x < columns.length; x++) {
                    int svgColor = vectorPixel(root, x, 5, background);
                    assertChannelsWithinTwo("SVG column " + x + ", scale " + scale,
                            reference.getRGB(x * scale, 5 * scale), svgColor);
                    assertChannelsWithinTwo("SVG/raster column " + x + ", scale " + scale,
                            raster.getRGB(x * scale, 5 * scale), svgColor);
                    assertEquals("Vector gaps remain untouched", columns[x].length != 0,
                            vectorRectangleCovers(root, x, 5));
                }
            } finally {
                graphics.dispose();
            }
        }
    }

    private static Color[][] blendingColumns() {
        Color red = new Color(255, 0, 0, 51);
        Color green = new Color(0, 255, 0, 128);
        Color blue = new Color(0, 0, 255, 51);
        Color[] lowAlpha = new Color[16];
        Arrays.fill(lowAlpha, new Color(0, 0, 255, 16));
        return new Color[][]{{}, {red, green, blue}, {red, green, blue},
                {blue, green, red}, {}, lowAlpha, {}, {}};
    }

    private static void blendColumns(BaseRenderer.ColorStrip strip, Color[][] columns) {
        for (int x = 0; x < columns.length; x++) {
            for (Color color : columns[x]) strip.blendColor(x, color);
        }
    }

    private static BufferedImage stripImage(Color[][] columns, Color background, float extraAlpha,
                                            int scale, boolean batched) {
        BufferedImage image = new BufferedImage(columns.length * scale, 12 * scale, BufferedImage.TYPE_INT_RGB);
        Graphics2D graphics = image.createGraphics();
        try {
            graphics.setColor(background);
            graphics.fillRect(0, 0, image.getWidth(), image.getHeight());
            graphics.scale(scale, scale);
            graphics.setComposite(AlphaComposite.getInstance(AlphaComposite.SRC_OVER, extraAlpha));
            if (batched) {
                BaseRenderer.ColorStrip strip = new BaseRenderer.ColorStrip();
                strip.reset(graphics, 0, columns.length);
                blendColumns(strip, columns);
                strip.draw(2, 6);
            } else {
                for (int x = 0; x < columns.length; x++) {
                    for (Color color : columns[x]) {
                        BaseRenderer.drawBase(graphics, color, 'N', x, 2, 1, 6, false, null);
                    }
                }
            }
        } finally {
            graphics.dispose();
        }
        return image;
    }

    private static boolean vectorRectangleCovers(Element root, int x, int y) {
        NodeList rectangles = root.getElementsByTagName("rect");
        for (int i = 0; i < rectangles.getLength(); i++) {
            if (rectangleContains((Element) rectangles.item(i), x, y)) return true;
        }
        return false;
    }

    private static boolean rectangleContains(Element rectangle, int x, int y) {
        double left = Double.parseDouble(rectangle.getAttribute("x"));
        double top = Double.parseDouble(rectangle.getAttribute("y"));
        double width = Double.parseDouble(rectangle.getAttribute("width"));
        double height = Double.parseDouble(rectangle.getAttribute("height"));
        return x >= left && x < left + width && y >= top && y < top + height;
    }

    // SVGGraphics2D exports presentation attributes inherited from enclosing groups.
    // Replay rectangle paints through Java2D rather than reproducing the strip's blending formula.
    // Coordinates are logical/user-space, so the same samples apply at either SVG scale.
    private static int vectorPixel(Element root, int x, int y, Color background) {
        BufferedImage image = new BufferedImage(1, 1, BufferedImage.TYPE_INT_RGB);
        Graphics2D graphics = image.createGraphics();
        try {
            graphics.setColor(background);
            graphics.fillRect(0, 0, 1, 1);
            NodeList rectangles = root.getElementsByTagName("rect");
            for (int i = 0; i < rectangles.getLength(); i++) {
                Element rectangle = (Element) rectangles.item(i);
                if (!rectangleContains(rectangle, x, y)) continue;
                String fill = inheritedAttribute(rectangle, "fill", "black");
                assertTrue("Expected an SVG RGB fill, got " + fill, fill.startsWith("rgb("));
                String[] channels = fill.substring(4, fill.length() - 1).split(",");
                Color color = new Color(Integer.parseInt(channels[0].trim()),
                        Integer.parseInt(channels[1].trim()), Integer.parseInt(channels[2].trim()));
                float alpha = Float.parseFloat(inheritedAttribute(rectangle, "fill-opacity", "1"));
                for (org.w3c.dom.Node node = rectangle; node instanceof Element; node = node.getParentNode()) {
                    String opacity = ((Element) node).getAttribute("opacity");
                    if (!opacity.isEmpty()) alpha *= Float.parseFloat(opacity);
                }
                graphics.setComposite(AlphaComposite.getInstance(AlphaComposite.SRC_OVER, alpha));
                graphics.setColor(color);
                graphics.fillRect(0, 0, 1, 1);
            }
        } finally {
            graphics.dispose();
        }
        return image.getRGB(0, 0);
    }

    private static String inheritedAttribute(Element element, String name, String fallback) {
        for (org.w3c.dom.Node node = element; node instanceof Element; node = node.getParentNode()) {
            String value = ((Element) node).getAttribute(name);
            if (!value.isEmpty()) return value;
        }
        return fallback;
    }

    private BufferedImage originalSoftClipImage(Alignment alignment, double start, int end, int width, int scale) {
        AlignmentTrack.RenderOptions options = track.getRenderOptions();
        boolean showAllBases = options.isShowAllBases();
        boolean showMismatches = options.isShowMismatches();
        BufferedImage reference = new BufferedImage(width * scale, 12 * scale, BufferedImage.TYPE_INT_RGB);
        Graphics2D graphics = reference.createGraphics();
        RenderContext context = null;
        ReferenceFrame frame = frame((int) start, end, width);
        frame.setOrigin(start);
        try {
            graphics.setColor(Color.WHITE);
            graphics.fillRect(0, 0, reference.getWidth(), reference.getHeight());
            graphics.scale(scale, scale);
            context = new RenderContext(null, graphics, frame, new Rectangle(0, 0, width, 12));
            options.setShowAllBases(false);
            options.setShowMismatches(false);
            new AlignmentRenderer(track).renderAlignments(List.of(alignment), null, context,
                    new Rectangle(0, 0, width, 12), options);
            options.setShowAllBases(showAllBases);
            options.setShowMismatches(showMismatches);
            graphics.setComposite(AlphaComposite.getInstance(AlphaComposite.SRC_OVER, 0.75f));
            for (AlignmentBlock block : alignment.getAlignmentBlocks()) {
                if (!block.isSoftClip()) continue;
                int first = Math.max((int) Math.floor(frame.getOrigin()), block.getStart());
                int last = Math.min((int) Math.ceil(frame.getEnd()), block.getEnd());
                for (int loc = first; loc < last; loc++) {
                    int offset = loc - block.getStart();
                    char base = (char) block.getBases().getByte(offset);
                    if (base == '=' && !showAllBases) continue;
                    Color color = AlignmentRenderer.nucleotideColors.getOrDefault(base, Color.BLACK);
                    if (options.getShadeBasesOption()) {
                        color = BaseRenderer.getShadedColor(color, block.getQuality(offset),
                                options.getBaseQualityMin(), options.getBaseQualityMax());
                    }
                    int x = (int) ((loc - frame.getOrigin()) / frame.getScale());
                    BaseRenderer.drawBase(graphics, color, base, x, 0, 1, 10, false, null);
                }
            }
        } finally {
            options.setShowAllBases(showAllBases);
            options.setShowMismatches(showMismatches);
            if (context != null) context.dispose();
            graphics.dispose();
        }
        return reference;
    }


    private static void assertChannelsWithinTwo(String message, int expected, int actual) {
        assertEquals(message + " alpha", expected >>> 24, actual >>> 24);
        for (int shift : new int[]{16, 8, 0}) {
            int expectedChannel = (expected >>> shift) & 255;
            int actualChannel = (actual >>> shift) & 255;
            assertTrue(message + " channel " + shift + ": expected " + expectedChannel + ", got " + actualChannel,
                    Math.abs(expectedChannel - actualChannel) <= 2);
        }
    }

    private BufferedImage image(Alignment alignment, int start, int end, int width) {
        return image(alignment, start, end, width, 1);
    }

    private BufferedImage image(Alignment alignment, int start, int end, int width, int scale) {
        BufferedImage image = new BufferedImage(width * scale, 12 * scale, BufferedImage.TYPE_INT_RGB);
        Graphics2D graphics = image.createGraphics();
        try {
            graphics.setColor(Color.WHITE);
            graphics.fillRect(0, 0, image.getWidth(), image.getHeight());
            graphics.scale(scale, scale);
            // A caller's smoothing hint must not blur categorical base columns.
            graphics.setRenderingHint(RenderingHints.KEY_INTERPOLATION, RenderingHints.VALUE_INTERPOLATION_BICUBIC);
            render(alignment, graphics, start, end, width);
        } finally {
            graphics.dispose();
        }
        return image;
    }

    private SAMRecord record(String cigar, String bases) {
        SAMRecord record = new SAMRecord(header);
        record.setReadName("clipped");
        record.setReferenceName("chr16");
        record.setAlignmentStart(10001);
        record.setMappingQuality(60);
        record.setCigarString(cigar);
        record.setReadString(bases);
        byte[] qualities = new byte[bases.length()];
        Arrays.fill(qualities, (byte) 30);
        record.setBaseQualities(qualities);
        return record;
    }

    private void render(Alignment alignment, Graphics2D graphics, int start, int end, int width) {
        ReferenceFrame frame = frame(start, end, width);
        RenderContext context = new RenderContext(null, graphics, frame, new Rectangle(0, 0, width, 12));
        try {
            new AlignmentRenderer(track).renderAlignments(List.of(alignment), null, context,
                    new Rectangle(0, 0, width, 12), track.getRenderOptions());
        } finally {
            context.dispose();
        }
    }

    private ReferenceFrame frame(int start, int end, int width) {
        ReferenceFrame frame = new ReferenceFrame("clipping-test");
        frame.setBounds(0, width);
        frame.jumpTo("chr16", start, end);
        return frame;
    }
}
