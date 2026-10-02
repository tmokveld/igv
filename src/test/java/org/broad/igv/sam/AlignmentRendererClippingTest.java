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

import java.awt.*;
import java.awt.image.BufferedImage;
import java.io.File;
import java.util.Arrays;
import java.util.List;

import static org.broad.igv.prefs.Constants.SAM_SHOW_SOFT_CLIPPED;
import static org.junit.Assert.*;

public class AlignmentRendererClippingTest {
    @Rule
    public TemporaryFolder temporaryFolder = new TemporaryFolder();

    private Genome previousGenome;
    private boolean previousHeadless;
    private IGVPreferences preferences;
    private String previousSoftClips;
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
        previousSoftClips = preferences.get(SAM_SHOW_SOFT_CLIPPED, null);
        preferences.put(SAM_SHOW_SOFT_CLIPPED, "true");

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
        if (preferences != null) preferences.put(SAM_SHOW_SOFT_CLIPPED, previousSoftClips);
        GenomeManager.getInstance().setCurrentGenomeForTest(previousGenome);
        Globals.setHeadless(previousHeadless);
    }

    @Test
    public void zoomedOutClipsPaintAtMostOncePerPixel() {
        SAMAlignment alignment = new SAMAlignment(record("5000S20M5000S", "A".repeat(10020)));
        Document document = GenericDOMImplementation.getDOMImplementation()
                .createDocument("http://www.w3.org/2000/svg", "svg", null);
        SVGGraphics2D graphics = new SVGGraphics2D(document);
        try {
            render(alignment, graphics, 0, 20000, 400);
            int rectangles = graphics.getRoot().getElementsByTagName("rect").getLength();
            // Each 5000-base clip covers at most 101 columns at 50 bases/pixel.
            assertTrue("Unbounded per-base overdraw: " + rectangles, rectangles <= 202);
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
        assertEquals(detailed.getRGB(74, 5), overview.getRGB(18, 5));
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

    private BufferedImage image(Alignment alignment, int start, int end, int width) {
        BufferedImage image = new BufferedImage(width, 12, BufferedImage.TYPE_INT_RGB);
        Graphics2D graphics = image.createGraphics();
        try {
            graphics.setColor(Color.WHITE);
            graphics.fillRect(0, 0, width, 12);
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
        ReferenceFrame frame = new ReferenceFrame("clipping-test");
        frame.setBounds(0, width);
        frame.jumpTo("chr16", start, end);
        RenderContext context = new RenderContext(null, graphics, frame, new Rectangle(0, 0, width, 12));
        try {
            new AlignmentRenderer(track).renderAlignments(List.of(alignment), null, context,
                    new Rectangle(0, 0, width, 12), track.getRenderOptions());
        } finally {
            context.dispose();
        }
    }
}
