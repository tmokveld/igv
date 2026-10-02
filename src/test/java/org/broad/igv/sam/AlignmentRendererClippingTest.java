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
            // Two homogeneous clips need two runs, not one rectangle per column.
            assertEquals("Unmerged homopolymer runs", 2, rectangles.getLength());
            assertVectorRectangle((Element) rectangles.item(0), 100, 100);
            assertVectorRectangle((Element) rectangles.item(1), 200, 101);
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
        // Compare interior columns to their final base at one base/pixel. Check
        // both clip ends and every nucleotide, not the alignment's outline/arrow.
        for (int x = 22; x < 39; x++) {
            assertEquals(detailed.getRGB(2 * x + 1, 5), overview.getRGB(x, 5));
        }
        for (int x = 52; x < 68; x++) {
            assertEquals(detailed.getRGB(2 * x + 1, 5), overview.getRGB(x, 5));
        }
        assertNotEquals(detailed.getRGB(48, 5), detailed.getRGB(49, 5));
        assertNotEquals(detailed.getRGB(49, 5), detailed.getRGB(50, 5));
        assertNotEquals(detailed.getRGB(50, 5), detailed.getRGB(51, 5));
    }

    @Test
    public void overviewUsesQualityOfLastDrawableBaseNotTrailingEquals() {
        SAMRecord record = record("40S20M", "A=C=".repeat(10) + "A".repeat(20));
        byte[] qualities = record.getBaseQualities();
        for (int i = 2; i < 40; i += 4) qualities[i] = 0;
        record.setBaseQualities(qualities);
        SAMAlignment alignment = new SAMAlignment(record);
        track.getRenderOptions().setShadeBasesOption(true);
        BufferedImage overview = image(alignment, 9920, 10080, 40);
        BufferedImage detailed = image(alignment, 9920, 10080, 160);
        assertChannelsWithinTwo("Last drawable base quality", detailed.getRGB(74, 5), overview.getRGB(18, 5));
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
        assertEquals(image(alignment, 9920, 10080, 160).getRGB(75, 5), allBases.getRGB(18, 5));
    }

    @Test
    public void partialColumnsAtViewportEdgesRetainClippedBases() {
        SAMAlignment alignment = new SAMAlignment(record("80S20M", "C".repeat(80) + "A".repeat(20)));
        BufferedImage overview = image(alignment, 9933, 9976, 20);
        BufferedImage detailed = image(alignment, 9920, 10080, 160);
        int clippedColor = detailed.getRGB(55, 5);
        for (int x = 0; x < 20; x++) {
            assertEquals("Clipped column " + x, clippedColor, overview.getRGB(x, 5));
        }
    }

    @Test
    public void scaledSvgMergesColoredRunsWithoutBridgingEqualsOnlyGaps() {
        SAMAlignment alignment = new SAMAlignment(record("80S20M", "AAAA====".repeat(10) + "A".repeat(20)));
        Document document = GenericDOMImplementation.getDOMImplementation()
                .createDocument("http://www.w3.org/2000/svg", "svg", null);
        SVGGraphics2D graphics = new SVGGraphics2D(document);
        try {
            graphics.scale(2, 2);
            render(alignment, graphics, 9920, 10080, 80);
            Element root = graphics.getRoot();
            assertEquals("SVG must remain vector graphics", 0, root.getElementsByTagName("image").getLength());
            NodeList rectangles = root.getElementsByTagName("rect");
            assertEquals("Ten colored runs separated by untouched holes", 10, rectangles.getLength());
            for (int run = 0; run < rectangles.getLength(); run++) {
                assertVectorRectangle((Element) rectangles.item(run), 4 * run, 2);
            }
        } finally {
            graphics.dispose();
        }
    }

    @Test
    public void scaledRasterClipsKeepSharpColumnsAndDetailedColors() {
        SAMAlignment alignment = new SAMAlignment(record("40S20M40S",
                "AACCGGTT".repeat(5) + "A".repeat(20) + "TTGGCCAA".repeat(5)));
        BufferedImage detailed = image(alignment, 9920, 10080, 160);
        for (int scale : new int[]{1, 2}) {
            BufferedImage overview = image(alignment, 9920, 10080, 80, scale);
            for (int x = 22; x < 68; x++) {
                if (x >= 39 && x < 52) continue; // Aligned block and end indicators.
                int expected = detailed.getRGB(2 * x + 1, 5);
                for (int dx = 0; dx < scale; dx++) {
                    for (int y = 2 * scale; y < 8 * scale; y++) {
                        assertEquals("Scale " + scale + ", column " + x + ", row " + y,
                                expected, overview.getRGB(scale * x + dx, y));
                    }
                }
            }
            assertEquals(Color.WHITE.getRGB(), overview.getRGB(5 * scale, 5 * scale));
        }
    }

    @Test
    public void shadedClipsSelectTheLastDrawableBaseAndItsQualityAtBothScales() {
        SAMRecord record = record("48S20M", "TAC=TAG=TAT=".repeat(4) + "A".repeat(20));
        byte[] qualities = record.getBaseQualities();
        for (int i = 2; i < 48; i += 4) {
            qualities[i] = (byte) ((i / 4) % 3 == 0 ? 0 : (i / 4) % 3 == 1 ? 10 : 30);
            qualities[i + 1] = 40; // A trailing '=' must not replace the drawable sample.
        }
        record.setBaseQualities(qualities);
        SAMAlignment alignment = new SAMAlignment(record);
        track.getRenderOptions().setShadeBasesOption(true);
        BufferedImage detailed = image(alignment, 9920, 10080, 160);
        for (int scale : new int[]{1, 2}) {
            BufferedImage overview = image(alignment, 9920, 10080, 40, scale);
            for (int x = 10; x < 19; x++) {
                for (int dx = 0; dx < scale; dx++) {
                    String message = "Scale " + scale + ", shaded column " + x;
                    int expected = detailed.getRGB(4 * x + 2, 5);
                    int actual = overview.getRGB(scale * x + dx, 5 * scale);
                    if ((x - 8) % 3 == 2) {
                        assertEquals(message + " opaque sample", expected, actual);
                    } else {
                        assertChannelsWithinTwo(message, expected, actual);
                    }
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
    public void equalsOnlyHolesBetweenColoredColumnsStayUnpainted() {
        SAMAlignment alignment = new SAMAlignment(record("48S20M",
                "AA==CC==GG==TT==".repeat(3) + "A".repeat(20)));
        BufferedImage overview = image(alignment, 9920, 10080, 80);
        BufferedImage detailed = image(alignment, 9920, 10080, 160);
        SAMAlignment blankClip = new SAMAlignment(record("48S20M", "=".repeat(48) + "A".repeat(20)));
        BufferedImage blank = image(blankClip, 9920, 10080, 80);
        for (int x = 18; x < 38; x++) {
            assertEquals("Column " + x, detailed.getRGB(2 * x + 1, 5), overview.getRGB(x, 5));
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
            BufferedImage detailed = image(alignments.get(row), 9920, 10160, 240);
            int[] ranges = clipRanges[row];
            for (int range = 0; range < ranges.length; range += 2) {
                for (int x = ranges[range]; x < ranges[range + 1]; x++) {
                    assertEquals("Read " + row + ", column " + x,
                            detailed.getRGB(2 * x + 1, 5), overview.getRGB(x, 12 * row + 5));
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
        int expectedColor = image(alignment, 9920, 10080, 160).getRGB(55, 5);
        for (int scale : new int[]{1, 2}) {
            BufferedImage image = new BufferedImage(20 * scale, 12 * scale, BufferedImage.TYPE_INT_RGB);
            Graphics2D graphics = image.createGraphics();
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
            } finally {
                if (context != null) context.dispose();
                graphics.dispose();
            }
            for (int x = 0; x < image.getWidth(); x++) {
                assertEquals("Scale " + scale + ", viewport column " + x,
                        x >= 3 * scale && x < 17 * scale ? expectedColor : Color.WHITE.getRGB(),
                        image.getRGB(x, 5 * scale));
                assertEquals("Above viewport clip", Color.WHITE.getRGB(), image.getRGB(x, scale));
                assertEquals("Below viewport clip", Color.WHITE.getRGB(), image.getRGB(x, 10 * scale));
            }
        }
    }

    @Test
    public void alignedBasesStillCompositeEveryBaseInsteadOfUsingTheSoftClipPolicy() {
        track.getRenderOptions().setShowAllBases(true);
        SAMAlignment alignment = new SAMAlignment(record("80M", "AC".repeat(40)));
        BufferedImage overview = image(alignment, 9920, 10080, 80);
        BufferedImage reference = new BufferedImage(1, 1, BufferedImage.TYPE_INT_RGB);
        Graphics2D graphics = reference.createGraphics();
        try {
            graphics.setColor(Color.WHITE);
            graphics.fillRect(0, 0, 1, 1);
            graphics.setComposite(AlphaComposite.getInstance(AlphaComposite.SRC_OVER, 0.75f));
            // Aligned bases retain the original A-then-C overdraw, unlike soft clips.
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

    private static void assertVectorRectangle(Element rectangle, int x, int width) {
        assertEquals("Run x", x, Double.parseDouble(rectangle.getAttribute("x")), 0.0);
        assertEquals("Run width", width, Double.parseDouble(rectangle.getAttribute("width")), 0.0);
        assertEquals("Run y", 0.0, Double.parseDouble(rectangle.getAttribute("y")), 0.0);
        assertEquals("Run height", 10.0, Double.parseDouble(rectangle.getAttribute("height")), 0.0);
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
