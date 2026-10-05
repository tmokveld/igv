package org.broad.igv.sam;

import htsjdk.samtools.*;
import org.broad.igv.Globals;
import org.broad.igv.feature.Chromosome;
import org.broad.igv.feature.genome.Genome;
import org.broad.igv.feature.genome.GenomeManager;
import org.broad.igv.prefs.IGVPreferences;
import org.broad.igv.prefs.PreferencesManager;
import org.broad.igv.renderer.SequenceRenderer;
import org.broad.igv.track.RenderContext;
import org.broad.igv.track.Track;
import org.broad.igv.ui.panel.ReferenceFrame;
import org.broad.igv.util.ResourceLocator;
import org.junit.*;
import org.junit.rules.TemporaryFolder;

import javax.swing.*;
import java.awt.*;
import java.awt.image.BufferedImage;
import java.io.File;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;

import static org.broad.igv.prefs.Constants.*;
import static org.junit.Assert.*;

public class AlignmentTrackDamageTest {
    private static final int WIDTH = 240;
    private static final int VIEW_HEIGHT = 120;

    @Rule
    public TemporaryFolder temporaryFolder = new TemporaryFolder();

    private Genome previousGenome;
    private boolean previousHeadless;
    private Map<Character, Color> previousSequenceColors;
    private IGVPreferences preferences;
    private final Map<String, String> previousPreferences = new LinkedHashMap<>();
    private AlignmentTrack track;
    private ReferenceFrame frame;

    @Before
    public void setUp() throws Exception {
        previousHeadless = Globals.isHeadless();
        Globals.setHeadless(true);
        previousSequenceColors = SequenceRenderer.nucleotideColors;
        SequenceRenderer.setNucleotideColors();
        previousGenome = GenomeManager.getInstance().getCurrentGenome();
        Genome genome = new Genome("damage-test", List.of(
                new Chromosome(0, "chr1", 50000), new Chromosome(1, "chr2", 50000)));
        GenomeManager.getInstance().setCurrentGenomeForTest(genome);
        preferences = PreferencesManager.getPreferences();
        for (String key : List.of(SAM_DOWNSAMPLE_READS, SAM_SHOW_REF_SEQ, SAM_SHOW_SOFT_CLIPPED,
                SAM_SHOW_CENTER_LINE, SAM_SHOW_GROUP_SEPARATOR, SAM_FLAG_CLIPPING,
                SAM_CLIPPING_THRESHOLD, SAM_SHOW_CONNECTED_CHR_NAME, SAM_FILTER_ALIGNMENTS,
                SAM_QUALITY_THRESHOLD, SAM_DISPLAY_PAIRED)) {
            previousPreferences.put(key, preferences.get(key, null));
        }
        preferences.put(SAM_DOWNSAMPLE_READS, "false");
        preferences.put(SAM_SHOW_REF_SEQ, "false");
        preferences.put(SAM_SHOW_SOFT_CLIPPED, "true");
        preferences.put(SAM_SHOW_CENTER_LINE, "true");
        preferences.put(SAM_SHOW_GROUP_SEPARATOR, "true");
        preferences.put(SAM_FLAG_CLIPPING, "true");
        preferences.put(SAM_CLIPPING_THRESHOLD, "0");
        preferences.put(SAM_SHOW_CONNECTED_CHR_NAME, "true");
        preferences.put(SAM_FILTER_ALIGNMENTS, "false");
        preferences.put(SAM_QUALITY_THRESHOLD, "0");
        preferences.put(SAM_DISPLAY_PAIRED, "false");

        SAMFileHeader header = new SAMFileHeader();
        header.addSequence(new SAMSequenceRecord("chr1", 50000));
        header.addSequence(new SAMSequenceRecord("chr2", 50000));
        header.setSortOrder(SAMFileHeader.SortOrder.coordinate);
        for (int group = 0; group < 3; group++) {
            SAMReadGroupRecord readGroup = new SAMReadGroupRecord("group-" + group);
            readGroup.setSample("sample-" + group);
            header.addReadGroup(readGroup);
        }
        File bam = temporaryFolder.newFile("damage.bam");
        try (SAMFileWriter writer = new SAMFileWriterFactory().setCreateIndex(true)
                .makeBAMWriter(header, true, bam)) {
            for (int group = 0; group < 3; group++) {
                for (int read = 0; read < 8; read++) {
                    SAMRecord record = new SAMRecord(header);
                    record.setReadName("read-" + group + "-" + read);
                    record.setReferenceName("chr1");
                    record.setAlignmentStart(1001);
                    record.setMappingQuality(60);
                    record.setCigarString("10S35M5I25M20D40M10S");
                    record.setReadString("ACGT".repeat(31) + "A");
                    byte[] qualities = new byte[125];
                    Arrays.fill(qualities, (byte) 30);
                    record.setBaseQualities(qualities);
                    record.setReadNegativeStrandFlag(read % 2 == 1);
                    record.setAttribute("RG", "group-" + group);
                    record.setAttribute("SA", "chr2,2001,+,115S10M,60,0;");
                    writer.addAlignment(record);
                }
            }
        }
        ResourceLocator locator = new ResourceLocator(bam.getAbsolutePath());
        track = new AlignmentTrack(locator, new AlignmentDataManager(locator, genome), genome);
        track.setDisplayMode(Track.DisplayMode.EXPANDED);
        track.setColor(Color.GRAY);
        AlignmentTrack.RenderOptions options = track.getRenderOptions();
        options.setColorOption(AlignmentTrack.ColorOption.NONE);
        options.setShadeAlignmentsOption(AlignmentTrack.ShadeAlignmentsOption.NONE);
        options.setGroupByOption(AlignmentTrack.GroupOption.NONE);
        options.setShowMismatches(true);
        options.setShowAllBases(false);
        options.setQuickConsensusMode(false);
        options.setShadeBasesOption(false);
        options.setHideSmallIndels(false);
        options.setIndelQualColoring(true);
        frame = new ReferenceFrame("damage-test");
        frame.setBounds(0, WIDTH);
        frame.jumpTo("chr1", 960, 1140);
        track.load(frame);
        assertNotNull("Temporary BAM must be loaded", track.getDataManager().getLoadedInterval(frame));
    }

    @After
    public void tearDown() {
        if (track != null) {
            track.getCoverageTrack().unload();
            track.getSpliceJunctionTrack().unload();
            track.unload();
        }
        if (preferences != null) previousPreferences.forEach(preferences::put);
        SequenceRenderer.nucleotideColors = previousSequenceColors;
        GenomeManager.getInstance().setCurrentGenomeForTest(previousGenome);
        Globals.setHeadless(previousHeadless);
    }

    @Test
    public void narrowDamageMatchesFullRenderAcrossRowsAndScrolledViewports() {
        for (int scale : new int[]{1, 2}) {
            for (int scroll : new int[]{0, 41, 127}) {
                BufferedImage full = render(scale, scroll, null, false, false);
                for (Rectangle damage : List.of(new Rectangle(0, scroll, WIDTH, 1),
                        new Rectangle(29, scroll + 29, 178, 16),
                        new Rectangle(0, scroll + 63, WIDTH, 1),
                        new Rectangle(0, scroll + VIEW_HEIGHT - 1, WIDTH, 1))) {
                    assertDamageMatches(full, render(scale, scroll, damage, false, false), damage, scale, scroll);
                }
            }
        }
    }

    @Test
    public void squishedAndGroupedDamageKeepsViewportLayoutAndGroupPositions() {
        for (Track.DisplayMode mode : List.of(Track.DisplayMode.EXPANDED, Track.DisplayMode.SQUISHED)) {
            track.setDisplayMode(mode);
            track.getRenderOptions().setGroupByOption(AlignmentTrack.GroupOption.SAMPLE);
            track.getDataManager().packAlignments(track.getRenderOptions(), mode);
            for (int scale : new int[]{1, 2}) {
                int scroll = mode == Track.DisplayMode.SQUISHED ? 7 : 111;
                BufferedImage full = render(scale, scroll, null, false, false);
                List<Row> rows = rows();
                Map<Row, double[]> geometry = new LinkedHashMap<>();
                for (Row row : rows) {
                    if (row.y <= scroll + VIEW_HEIGHT && row.y + row.h > scroll) {
                        geometry.put(row, new double[]{row.y, row.h});
                        row.y = -1000;
                        row.h = -1000;
                    }
                }
                Rectangle damage = new Rectangle(0, scroll + 51, WIDTH, 3);
                assertDamageMatches(full, render(scale, scroll, damage, false, false), damage, scale, scroll);
                for (Row row : geometry.keySet()) {
                    assertEquals("Row position must survive damage culling", geometry.get(row)[0], row.y, 0);
                    assertEquals("Squished height comes from viewport, not damage", geometry.get(row)[1], row.h, 0);
                    int hitY = (int) (row.y + row.h / 2);
                    if (hitY >= scroll && hitY < scroll + VIEW_HEIGHT) {
                        assertSame("Undamaged visible rows remain selectable", row.alignments.get(0),
                                track.getAlignmentAt(1020, hitY, frame));
                    }
                }
                // Group labels and dividers below the damaged row retain their full-group offsets.
                Rectangle bottom = new Rectangle(0, scroll + VIEW_HEIGHT - 18, WIDTH, 18);
                assertDamageMatches(full, render(scale, scroll, bottom, false, false), bottom, scale, scroll);
            }
        }
    }


    @Test
    public void selectedOutlineSpillingAboveRowIsRestoredByPartialPaint() {
        for (Row row : rows()) {
            track.getSelectedReadNames().put(row.alignments.get(0).getReadName(), Color.BLUE);
        }
        for (int scale : new int[]{1, 2}) {
            BufferedImage full = render(scale, 0, null, false, false);
            int spillY = (int) rows().get(5).y - 1;
            assertNotEquals("Selection stroke must really spill above the row", Color.WHITE.getRGB(),
                    full.getRGB(80 * scale, spillY * scale));
            Rectangle damage = new Rectangle(70, spillY, 100, 1);
            assertDamageMatches(full, render(scale, 0, damage, false, false), damage, scale, 0);
        }
    }

    @Test
    public void nullGraphicsClipRetainsFullVisibleRendering() {
        track.getRenderOptions().setGroupByOption(AlignmentTrack.GroupOption.SAMPLE);
        track.getDataManager().packAlignments(track.getRenderOptions(), track.getDisplayMode());
        for (int scale : new int[]{1, 2}) {
            int scroll = 53;
            Rectangle viewport = new Rectangle(0, scroll, WIDTH, VIEW_HEIGHT);
            assertDamageMatches(render(scale, scroll, null, false, false),
                    render(scale, scroll, null, true, false), viewport, scale, scroll);
        }
    }

    @Test
    public void expandedInsertionDamageMatchesFullInsertionRendering() {
        for (Track.DisplayMode mode : List.of(Track.DisplayMode.EXPANDED, Track.DisplayMode.SQUISHED)) {
            track.setDisplayMode(mode);
            track.getRenderOptions().setGroupByOption(AlignmentTrack.GroupOption.SAMPLE);
            track.getDataManager().packAlignments(track.getRenderOptions(), mode);
            for (int scale : new int[]{1, 2}) {
                int scroll = 17;
                BufferedImage full = render(scale, scroll, null, false, true);
                for (Rectangle damage : List.of(new Rectangle(45, scroll, 150, 1),
                        new Rectangle(45, scroll + 42, 150, 7),
                        new Rectangle(45, scroll + VIEW_HEIGHT - 1, 150, 1))) {
                    assertDamageMatches(full, render(scale, scroll, damage, false, true), damage, scale, scroll);
                }
            }
        }
    }

    @Test
    public void flowQualityBarsFarAboveTheirRowsSurviveDamage() throws Exception {
        SAMFileHeader header = new SAMFileHeader();
        header.addSequence(new SAMSequenceRecord("chr1", 50000));
        header.setSortOrder(SAMFileHeader.SortOrder.coordinate);
        SAMReadGroupRecord readGroup = new SAMReadGroupRecord("flow");
        readGroup.setPlatform("ULTIMA");
        readGroup.setSample("flow");
        header.addReadGroup(readGroup);
        File bam = temporaryFolder.newFile("flow-damage.bam");
        try (SAMFileWriter writer = new SAMFileWriterFactory().setCreateIndex(true)
                .makeBAMWriter(header, true, bam)) {
            for (int i = 0; i < 8; i++) {
                SAMRecord record = new SAMRecord(header);
                record.setReadName("flow-" + i);
                record.setReferenceName("chr1");
                record.setAlignmentStart(1001);
                record.setMappingQuality(60);
                record.setCigarString("30M1I30M");
                record.setReadString("A".repeat(30) + "C" + "A".repeat(30));
                byte[] qualities = new byte[61];
                Arrays.fill(qualities, (byte) 93);
                record.setBaseQualities(qualities);
                record.setAttribute("RG", "flow");
                record.setAttribute("tp", new byte[61]);
                writer.addAlignment(record);
            }
        }
        track.getCoverageTrack().unload();
        track.getSpliceJunctionTrack().unload();
        track.unload();
        Genome genome = GenomeManager.getInstance().getCurrentGenome();
        ResourceLocator locator = new ResourceLocator(bam.getAbsolutePath());
        track = new AlignmentTrack(locator, new AlignmentDataManager(locator, genome), genome);
        track.setDisplayMode(Track.DisplayMode.EXPANDED);
        AlignmentTrack.RenderOptions options = track.getRenderOptions();
        options.setColorOption(AlignmentTrack.ColorOption.NONE);
        options.setGroupByOption(AlignmentTrack.GroupOption.NONE);
        options.setShowMismatches(false);
        options.setShowAllBases(false);
        options.setHideSmallIndels(false);
        options.setIndelQualUsesMin(true);
        track.load(frame);
        for (int scale : new int[]{1, 2}) {
            options.setIndelQualColoring(false);
            BufferedImage ordinary = render(scale, 0, null, false, false);
            options.setIndelQualColoring(true);
            BufferedImage full = render(scale, 0, null, false, false);
            int x = (int) ((1030 - frame.getOrigin()) / frame.getScale());
            Row lastRow = rows().get(7);
            int spillY = (int) (lastRow.y + 11 * ((42 - 93.0) / 42)) - 1;
            assertTrue("Quality decoration must be outside normal row padding", lastRow.y - spillY > 6);
            assertNotEquals("Quality coloring must really paint the spilled pixel",
                    ordinary.getRGB(x * scale, spillY * scale), full.getRGB(x * scale, spillY * scale));
            Rectangle damage = new Rectangle(x - 2, spillY, 6, 1);
            assertDamageMatches(full, render(scale, 0, damage, false, false), damage, scale, 0);
        }
    }

    private PackedAlignments groups() {
        AlignmentDataManager manager = track.getDataManager();
        return manager.getGroups(manager.getLoadedInterval(frame), track.getRenderOptions());
    }

    private List<Row> rows() {
        List<Row> result = new ArrayList<>();
        groups().values().forEach(result::addAll);
        return result;
    }

    private BufferedImage render(int scale, int scroll, Rectangle damage, boolean nullClip, boolean insertion) {
        BufferedImage image = new BufferedImage(WIDTH * scale, VIEW_HEIGHT * scale, BufferedImage.TYPE_INT_ARGB);
        Graphics2D graphics = image.createGraphics();
        graphics.setColor(Color.WHITE);
        graphics.fillRect(0, 0, image.getWidth(), image.getHeight());
        graphics.setColor(Color.BLACK);
        graphics.scale(scale, scale);
        graphics.translate(0, -scroll);
        Rectangle viewport = new Rectangle(0, scroll, WIDTH, VIEW_HEIGHT);
        graphics.setClip(nullClip ? null : damage == null ? viewport : viewport.intersection(damage));
        JPanel panel = new JPanel();
        panel.setBackground(Color.WHITE);
        RenderContext context = new RenderContext(panel, graphics, frame, viewport);
        try {
            if (insertion) {
                track.renderExpandedInsertion(new InsertionMarker(1035, 5), context, new Rectangle(0, 0, WIDTH, 400));
            } else {
                track.render(context, new Rectangle(0, 0, WIDTH, 400));
            }
            assertEquals("Painting must not replace the viewport with damage", viewport, context.getVisibleRect());
            return image;
        } finally {
            context.dispose();
            graphics.dispose();
        }
    }

    private void assertDamageMatches(BufferedImage full, BufferedImage partial, Rectangle damage, int scale, int scroll) {
        Rectangle deviceDamage = new Rectangle(damage.x * scale, (damage.y - scroll) * scale,
                damage.width * scale, damage.height * scale);
        for (int y = 0; y < full.getHeight(); y++) {
            for (int x = 0; x < full.getWidth(); x++) {
                if (deviceDamage.contains(x, y)) {
                    assertEquals("Damaged pixel " + x + "," + y + " at scale " + scale + ", scroll " + scroll,
                            full.getRGB(x, y), partial.getRGB(x, y));
                } else {
                    assertEquals("Painting must stay inside graphics damage", Color.WHITE.getRGB(), partial.getRGB(x, y));
                }
            }
        }
    }
}
