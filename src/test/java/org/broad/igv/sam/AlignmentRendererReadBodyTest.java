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
import org.broad.igv.ui.color.ColorUtilities;
import org.broad.igv.util.ResourceLocator;
import org.junit.*;
import org.junit.rules.TemporaryFolder;
import org.w3c.dom.Document;
import org.w3c.dom.Element;
import org.w3c.dom.Node;
import org.w3c.dom.NodeList;

import java.awt.*;
import java.awt.image.BufferedImage;
import java.io.File;
import java.util.Arrays;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;

import static org.broad.igv.prefs.Constants.*;
import static org.junit.Assert.*;

public class AlignmentRendererReadBodyTest {
    @Rule
    public TemporaryFolder temporaryFolder = new TemporaryFolder();

    private Genome previousGenome;
    private boolean previousHeadless;
    private IGVPreferences preferences;
    private final Map<String, String> previousPreferences = new LinkedHashMap<>();
    private AlignmentTrack track;
    private SAMFileHeader header;

    @Before
    public void setUp() throws Exception {
        previousHeadless = Globals.isHeadless();
        Globals.setHeadless(true);
        previousGenome = GenomeManager.getInstance().getCurrentGenome();
        Genome genome = new Genome("read-body-test", List.of(new Chromosome(0, "chr16", 50000)));
        GenomeManager.getInstance().setCurrentGenomeForTest(genome);
        preferences = PreferencesManager.getPreferences();
        for (String key : List.of(SAM_SHOW_SOFT_CLIPPED, SAM_SHOW_CENTER_LINE, SAM_FLAG_CLIPPING,
                SAM_CLIPPING_THRESHOLD, SAM_FLAG_UNMAPPED_PAIR, SAM_FLAG_ZERO_QUALITY, SAM_FLAG_LARGE_INDELS)) {
            previousPreferences.put(key, preferences.get(key, null));
        }
        preferences.put(SAM_SHOW_SOFT_CLIPPED, "false");
        preferences.put(SAM_SHOW_CENTER_LINE, "false");
        preferences.put(SAM_FLAG_CLIPPING, "false");
        preferences.put(SAM_CLIPPING_THRESHOLD, "0");
        preferences.put(SAM_FLAG_UNMAPPED_PAIR, "true");
        preferences.put(SAM_FLAG_ZERO_QUALITY, "true");
        preferences.put(SAM_FLAG_LARGE_INDELS, "false");

        header = new SAMFileHeader();
        header.addSequence(new SAMSequenceRecord("chr16", 50000));
        header.setSortOrder(SAMFileHeader.SortOrder.coordinate);
        File bam = temporaryFolder.newFile("body.bam");
        try (SAMFileWriter writer = new SAMFileWriterFactory().setCreateIndex(true).makeBAMWriter(header, true, bam)) {
            writer.addAlignment(record("40M", false));
        }
        ResourceLocator locator = new ResourceLocator(bam.getAbsolutePath());
        track = new AlignmentTrack(locator, new AlignmentDataManager(locator, genome), genome);
        track.setDisplayMode(Track.DisplayMode.EXPANDED);
        track.setColor(Color.GRAY);
        AlignmentTrack.RenderOptions options = track.getRenderOptions();
        options.setColorOption(AlignmentTrack.ColorOption.NONE);
        options.setShadeAlignmentsOption(AlignmentTrack.ShadeAlignmentsOption.NONE);
        options.setShowMismatches(false);
        options.setShowAllBases(false);
        options.setQuickConsensusMode(false);
        options.setShadeBasesOption(false);
        options.setHideSmallIndels(false);
        options.setDuplicatesOption(AlignmentTrack.DuplicatesOption.SHOW);
    }

    @After
    public void tearDown() {
        if (track != null) {
            track.getCoverageTrack().unload();
            track.getSpliceJunctionTrack().unload();
            track.unload();
        }
        if (preferences != null) previousPreferences.forEach(preferences::put);
        GenomeManager.getInstance().setCurrentGenomeForTest(previousGenome);
        Globals.setHeadless(previousHeadless);
    }

    @Test
    public void rectangularBodiesMatchLegacyPolygonsAtViewportEdgesAndDeviceScales() {
        // Short rows, clipped strand tips, and tips rounded down to zero are all rectangles.
        double[][] windows = {{9980, 10140}, {9999.75, 10015.75}, {10015.25, 10065.25},
                {10001.25, 10039.25}, {9900, 10900}};
        for (boolean negative : new boolean[]{false, true}) {
            for (int height : new int[]{4, 6, 10}) {
                for (double[] window : windows) {
                    for (double deviceScale : new double[]{1, 1.25, 2}) {
                        assertLegacy(List.of(new SAMAlignment(record("40M", negative))), window[0], window[1],
                                100, height, deviceScale);
                    }
                }
            }
        }
        track.setDisplayMode(Track.DisplayMode.SQUISHED);
        assertLegacy(List.of(new SAMAlignment(record("40M", true))), 9980, 10140, 100, 6, 2);
    }

    @Test
    public void terminalArrowsAndLaterBlocksRetainStrandAndLeftmostBookkeeping() {
        for (boolean negative : new boolean[]{false, true}) {
            SAMAlignment alignment = new SAMAlignment(record("80M", negative));
            for (double[] window : new double[][]{{9990, 10090}, {10020, 10120}, {9970, 10070}}) {
                assertLegacy(List.of(alignment), window[0], window[1], 100, 10, 1);
            }
            // The first block at x=0 is rectangular. It must still consume `leftmost`,
            // otherwise the negative-strand second block gains an incorrect arrow.
            assertLegacy(List.of(new SAMAlignment(record("20M5D30M", negative))), 10000, 10100, 100, 10, 2);
            // First block rounds to zero width; filling it must not paint a pixel or
            // leave the next block incorrectly marked as the first painted block.
            assertLegacy(List.of(new SAMAlignment(record("1M2D40M", negative))), 10000.4, 10300.4, 100, 10, 1);
        }
    }

    @Test
    public void repeatedBodyAndAllBaseOverlaysKeepOrderedAlpha() {
        SAMAlignment first = new SAMAlignment(record("40M", false));
        SAMRecord secondRecord = record("40M", true);
        secondRecord.setReadName("second");
        SAMAlignment second = new SAMAlignment(secondRecord);
        List<Alignment> alignments = List.of(first, second, first);
        assertLegacy(alignments, 9980, 10140, 80, 6, 2);
        track.getRenderOptions().setShowAllBases(true);
        assertLegacy(alignments, 9980, 10140, 80, 6, 1);
    }

    @Test
    public void selectionQualityMateOutlinesTextureAndClippingKeepLegacyBoundaries() {
        for (boolean negative : new boolean[]{false, true}) {
            for (int height : new int[]{6, 10}) {
                SAMAlignment ordinary = new SAMAlignment(record("40M", negative));
                track.getSelectedReadNames().put(ordinary.getReadName(), Color.BLUE);
                assertLegacy(List.of(ordinary), 9980, 10080, 100, height, 2);
                track.getSelectedReadNames().clear();

                SAMRecord zeroQuality = record("40M", negative);
                zeroQuality.setMappingQuality(0);
                assertLegacy(List.of(new SAMAlignment(zeroQuality)), 9980, 10080, 100, height, 1);

                SAMRecord unmappedMate = record("40M", negative);
                unmappedMate.setReadPairedFlag(true);
                unmappedMate.setMateUnmappedFlag(true);
                assertLegacy(List.of(new SAMAlignment(unmappedMate)), 9980, 10080, 100, height, 1);

                SAMRecord duplicate = record("40M", negative);
                duplicate.setDuplicateReadFlag(true);
                track.getRenderOptions().setDuplicatesOption(AlignmentTrack.DuplicatesOption.TEXTURE);
                assertLegacy(List.of(new SAMAlignment(duplicate)), 9980.25, 10080.25, 100, height, 2);
                track.getRenderOptions().setDuplicatesOption(AlignmentTrack.DuplicatesOption.SHOW);

                preferences.put(SAM_FLAG_CLIPPING, "true");
                assertLegacy(List.of(new SAMAlignment(record("8H40M9H", negative))), 9980, 10080, 100, height, 2);
                preferences.put(SAM_FLAG_CLIPPING, "false");
            }
        }
    }

    @Test
    public void svgKeepsVectorBodyCoverageAlphaAndRealStrandTips() {
        for (boolean negative : new boolean[]{false, true}) {
            for (int height : new int[]{6, 10}) {
                Document document = GenericDOMImplementation.getDOMImplementation()
                        .createDocument("http://www.w3.org/2000/svg", "svg", null);
                SVGGraphics2D graphics = new SVGGraphics2D(document);
                try {
                    graphics.scale(2, 2);
                    SAMAlignment alignment = new SAMAlignment(record("40M", negative));
                    render(List.of(alignment, alignment), graphics, frame(9980, 10080, 100), height);
                    Element root = graphics.getRoot();
                    assertEquals("Bodies and arrows must remain vectors", 0, root.getElementsByTagName("image").getLength());
                    assertEquals("Repeated fills must retain 0.75 alpha", compositeGray(2), vectorPixel(root, 30, 3));
                    assertEquals(Color.WHITE.getRGB(), vectorPixel(root, 30, height + 3));
                    if (height == 10) {
                        double tipX = negative ? 18 : 62;
                        double otherSide = negative ? 62 : 18;
                        assertEquals("Strand tip remains filled", compositeGray(2), vectorPixel(root, tipX, 7));
                        assertEquals("Opposite edge is rectangular", Color.WHITE.getRGB(), vectorPixel(root, otherSide, 7));
                    } else {
                        assertEquals(Color.WHITE.getRGB(), vectorPixel(root, 18, 3));
                        assertEquals(Color.WHITE.getRGB(), vectorPixel(root, 62, 3));
                    }
                } finally {
                    graphics.dispose();
                }
            }
        }
    }

    private SAMRecord record(String cigar, boolean negative) {
        SAMRecord record = new SAMRecord(header);
        record.setReadName("body");
        record.setReferenceName("chr16");
        record.setAlignmentStart(10001);
        record.setMappingQuality(60);
        record.setCigarString(cigar);
        int length = record.getCigar().getReadLength();
        record.setReadString("ACGT".repeat((length + 3) / 4).substring(0, length));
        byte[] qualities = new byte[length];
        Arrays.fill(qualities, (byte) 30);
        record.setBaseQualities(qualities);
        record.setReadNegativeStrandFlag(negative);
        return record;
    }

    private ReferenceFrame frame(double start, double end, int width) {
        ReferenceFrame frame = new ReferenceFrame("read-body-test") {
            @Override
            public double getScale() {
                return (end - start) / width;
            }
        };
        frame.setBounds(0, width);
        frame.jumpTo("chr16", (int) start, (int) end);
        frame.setOrigin(start);
        return frame;
    }

    private void render(List<Alignment> alignments, Graphics2D graphics, ReferenceFrame frame, int bodyHeight) {
        Rectangle row = new Rectangle(0, 2, frame.getWidthInPixels(),
                bodyHeight + (track.getDisplayMode() == Track.DisplayMode.SQUISHED ? 0 : 2));
        RenderContext context = new RenderContext(null, graphics, frame, row);
        try {
            new AlignmentRenderer(track).renderAlignments(alignments, null, context, row, track.getRenderOptions());
        } finally {
            context.dispose();
        }
    }

    private void assertLegacy(List<Alignment> alignments, double start, double end, int width, int height, double scale) {
        BufferedImage actual = image(width, height, scale);
        BufferedImage expected = image(width, height, scale);
        ReferenceFrame frame = frame(start, end, width);
        Graphics2D actualGraphics = graphics(actual, scale);
        Graphics2D expectedGraphics = graphics(expected, scale);
        try {
            render(alignments, actualGraphics, frame, height);
            for (Alignment alignment : alignments) legacyPolygon(expectedGraphics, alignment, frame, width, height);
        } finally {
            actualGraphics.dispose();
            expectedGraphics.dispose();
        }
        for (int y = 0; y < actual.getHeight(); y++) {
            for (int x = 0; x < actual.getWidth(); x++) {
                assertEquals("Pixel " + x + "," + y + "; origin=" + start + ", height=" + height + ", scale=" + scale,
                        expected.getRGB(x, y), actual.getRGB(x, y));
            }
        }
    }

    private static BufferedImage image(int width, int height, double scale) {
        return new BufferedImage((int) Math.ceil(width * scale), (int) Math.ceil((height + 6) * scale), BufferedImage.TYPE_INT_RGB);
    }

    private static Graphics2D graphics(BufferedImage image, double scale) {
        Graphics2D graphics = image.createGraphics();
        graphics.setColor(Color.WHITE);
        graphics.fillRect(0, 0, image.getWidth(), image.getHeight());
        graphics.setColor(Color.BLACK);
        graphics.scale(scale, scale);
        return graphics;
    }

    // Independent pre-optimization oracle: every block is the original six-point
    // polygon, including collapsed tips/zero-width blocks. No renderer shape helper
    // or optimized body drawing is used to construct the expected image.
    private void legacyPolygon(Graphics2D graphics, Alignment alignment, ReferenceFrame frame, int width, int height) {
        double origin = frame.getOrigin();
        double scale = frame.getScale();
        final int y = 2;
        if (alignment.getGaps() != null) {
            for (Gap gap : alignment.getGaps()) {
                if (gap.getStart() + gap.getnBases() <= origin || gap.getStart() >= frame.getEnd()) continue;
                Graphics2D gapGraphics = (Graphics2D) graphics.create();
                try {
                    gapGraphics.setColor(Color.BLACK);
                    if (height > 5) gapGraphics.setStroke(new BasicStroke(2));
                    gapGraphics.drawLine((int) ((Math.max(origin, gap.getStart()) - origin) / scale), y + height / 2,
                            (int) ((Math.min(Math.ceil(frame.getEnd()), gap.getStart() + gap.getnBases()) - origin) / scale), y + height / 2);
                } finally {
                    gapGraphics.dispose();
                }
            }
        }
        Color color = track.getColor();
        if (alignment.getMappingQuality() == 0) {
            color = ColorUtilities.getCompositeColor(Color.WHITE, color, 0.15f);
        }
        Color outline = null;
        float outlineWidth = 1;
        if (track.getSelectedReadNames().containsKey(alignment.getReadName())) {
            color = track.getSelectedReadNames().get(alignment.getReadName());
            if (color == null) color = Color.BLUE;
            outline = Color.BLACK;
            outlineWidth = 2;
        } else if (alignment.isPaired() && !alignment.getMate().isMapped()) {
            outline = Color.RED;
        } else if (alignment.getMappingQuality() == 0) {
            outline = new Color(185, 185, 185);
        }
        boolean texture = alignment.isDuplicate() && track.getRenderOptions().getDuplicatesOption() == AlignmentTrack.DuplicatesOption.TEXTURE;
        AlignmentBlock[] blocks = alignment.getAlignmentBlocks();
        double pixelLength = alignment.getLengthOnReference() / scale;
        int arrow = pixelLength == 0 ? 0 : (int) Math.min(Math.min(5, height / 2), pixelLength / 6);
        boolean first = true;
        for (int i = 0; i < blocks.length; i++) {
            AlignmentBlock block = blocks[i];
            int start = (int) Math.round((block.getStart() - origin) / scale);
            int end = (int) Math.round((block.getEnd() - origin) / scale);
            if (end < 0) continue;
            if (start > width) break;
            boolean last = i == blocks.length - 1;
            int spacing = (int) (AlignmentPacker.MIN_ALIGNMENT_SPACING / scale);
            if (spacing < arrow) arrow = Math.max(0, arrow - spacing);
            boolean pointed = height > 6 && (first && start > 0 || last && end < width);
            start = Math.max(0, start);
            end = Math.min(width, end);
            int leftTip = start - (first && alignment.isNegativeStrand() && pointed ? arrow : 0);
            int rightTip = end + (last && !alignment.isNegativeStrand() && pointed ? arrow : 0);
            Polygon body = new Polygon(new int[]{leftTip, start, end, rightTip, end, start},
                    new int[]{y + height / 2, y, y, y + height / 2, y + height, y + height}, 6);
            Graphics2D fill = (Graphics2D) graphics.create();
            try {
                fill.setComposite(AlphaComposite.getInstance(AlphaComposite.SRC_OVER, 0.75f));
                fill.setColor(color);
                if (texture) fill.setPaint(legacyTexture(color));
                fill.fill(body);
            } finally {
                fill.dispose();
            }
            Graphics2D decoration = (Graphics2D) graphics.create();
            try {
                if (outline != null) {
                    decoration.setColor(outline);
                    decoration.setStroke(new BasicStroke(outlineWidth));
                    decoration.draw(body);
                }
                decoration.setColor(new Color(255, 20, 147));
                decoration.setStroke(new BasicStroke(height > 5 ? 1.2f : 1));
                if (preferences.getAsBoolean(SAM_FLAG_CLIPPING)) {
                    if (first && alignment.getClippingCounts().getLeft() > 0) {
                        decoration.drawLine(leftTip, y + height / 2, start, y + height);
                        decoration.drawLine(start, y - 1, leftTip, y + height / 2);
                    }
                    if (last && alignment.getClippingCounts().getRight() > 0) {
                        decoration.drawLine(end, y + height, rightTip, y + height / 2);
                        decoration.drawLine(rightTip, y + height / 2, end, y - 1);
                    }
                }
            } finally {
                decoration.dispose();
            }
            first = false;
        }
        if (track.getRenderOptions().isShowAllBases()) {
            Graphics2D bases = (Graphics2D) graphics.create();
            try {
                bases.setComposite(AlphaComposite.getInstance(AlphaComposite.SRC_OVER, 0.75f));
                for (AlignmentBlock block : blocks) {
                    for (int loc = Math.max((int) Math.floor(origin), block.getStart());
                         loc < Math.min(Math.ceil(frame.getEnd()), block.getEnd()); loc++) {
                        char base = (char) block.getBases().getByte(loc - block.getStart());
                        bases.setColor(AlignmentRenderer.nucleotideColors.getOrDefault(base, Color.BLACK));
                        bases.fillRect((int) ((loc - origin) / scale), y, (int) Math.max(1, 1 / scale), height);
                    }
                }
            } finally {
                bases.dispose();
            }
        }
    }

    private static TexturePaint legacyTexture(Color color) {
        BufferedImage image = new BufferedImage(5, 5, BufferedImage.TYPE_INT_RGB);
        Graphics2D graphics = image.createGraphics();
        try {
            graphics.setColor(color);
            graphics.fillRect(0, 0, 5, 5);
            graphics.setColor(color.darker());
            graphics.drawLine(0, 2, 0, 3);
            graphics.drawLine(1, 3, 1, 4);
            graphics.drawLine(2, 4, 2, 4);
            graphics.drawLine(2, 0, 2, 0);
            graphics.drawLine(3, 0, 3, 1);
            graphics.drawLine(4, 1, 4, 2);
        } finally {
            graphics.dispose();
        }
        return new TexturePaint(image, new Rectangle(0, 0, 5, 5));
    }

    private static int compositeGray(int layers) {
        BufferedImage image = image(1, 1, 1);
        Graphics2D graphics = graphics(image, 1);
        try {
            graphics.setComposite(AlphaComposite.getInstance(AlphaComposite.SRC_OVER, 0.75f));
            graphics.setColor(Color.GRAY);
            for (int i = 0; i < layers; i++) graphics.fillRect(0, 0, 1, 1);
        } finally {
            graphics.dispose();
        }
        return image.getRGB(0, 0);
    }

    private static int vectorPixel(Element root, double x, double y) {
        BufferedImage image = image(1, 1, 1);
        Graphics2D graphics = graphics(image, 1);
        try {
            NodeList elements = root.getElementsByTagName("*");
            for (int i = 0; i < elements.getLength(); i++) {
                Element element = (Element) elements.item(i);
                Shape shape;
                if (element.getTagName().equals("rect")) {
                    shape = new java.awt.geom.Rectangle2D.Double(number(element, "x"), number(element, "y"),
                            number(element, "width"), number(element, "height"));
                } else if (element.getTagName().equals("polygon")) {
                    String[] points = element.getAttribute("points").trim().split("[ ,]+");
                    java.awt.geom.Path2D path = new java.awt.geom.Path2D.Double();
                    path.moveTo(Double.parseDouble(points[0]), Double.parseDouble(points[1]));
                    for (int p = 2; p < points.length; p += 2) path.lineTo(Double.parseDouble(points[p]), Double.parseDouble(points[p + 1]));
                    path.closePath();
                    shape = path;
                } else continue;
                if (!shape.contains(x, y)) continue;
                String fill = inherited(element, "fill", "black");
                Color color = ColorUtilities.stringToColorNoDefault(fill);
                assertNotNull("Expected a valid SVG body fill: " + fill, color);
                graphics.setColor(color);
                float alpha = Float.parseFloat(inherited(element, "fill-opacity", "1"));
                for (Node node = element; node instanceof Element; node = node.getParentNode()) {
                    String opacity = ((Element) node).getAttribute("opacity");
                    if (!opacity.isEmpty()) alpha *= Float.parseFloat(opacity);
                }
                graphics.setComposite(AlphaComposite.getInstance(AlphaComposite.SRC_OVER, alpha));
                graphics.fillRect(0, 0, 1, 1);
            }
        } finally {
            graphics.dispose();
        }
        return image.getRGB(0, 0);
    }

    private static double number(Element element, String name) {
        String value = element.getAttribute(name);
        return value.isEmpty() ? 0 : Double.parseDouble(value);
    }

    private static String inherited(Element element, String name, String fallback) {
        for (Node node = element; node instanceof Element; node = node.getParentNode()) {
            String value = ((Element) node).getAttribute(name);
            if (!value.isEmpty()) return value;
        }
        return fallback;
    }
}
