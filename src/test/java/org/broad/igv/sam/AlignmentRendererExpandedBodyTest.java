package org.broad.igv.sam;

import htsjdk.samtools.*;
import org.broad.igv.Globals;
import org.broad.igv.feature.Strand;
import org.broad.igv.feature.genome.Genome;
import org.broad.igv.feature.genome.GenomeManager;
import org.broad.igv.feature.genome.InMemorySequence;
import org.broad.igv.feature.genome.load.GenomeConfig;
import org.broad.igv.prefs.IGVPreferences;
import org.broad.igv.prefs.PreferencesManager;
import org.broad.igv.track.RenderContext;
import org.broad.igv.track.Track;
import org.broad.igv.ui.color.ColorUtilities;
import org.broad.igv.ui.panel.ReferenceFrame;
import org.broad.igv.util.ResourceLocator;
import org.junit.*;
import org.junit.rules.TemporaryFolder;

import java.awt.*;
import java.awt.geom.Ellipse2D;
import java.awt.image.BufferedImage;
import java.io.File;
import java.util.Arrays;
import java.util.HashMap;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;
import java.util.function.Consumer;

import static org.broad.igv.prefs.Constants.*;
import static org.junit.Assert.*;

/** Body/gap/decoration/base oracle never invokes the candidate renderer. */
public class AlignmentRendererExpandedBodyTest {
    private static final int WIDTH = 100;
    private static final Color BACKGROUND = new Color(231, 211, 179, 193);
    private static final Color BODY = new Color(47, 113, 149, 117);
    private static final float BODY_ALPHA = 0.75f;

    @Rule
    public TemporaryFolder temporaryFolder = new TemporaryFolder();

    private Genome previousGenome;
    private boolean previousHeadless;
    private IGVPreferences preferences;
    private final Map<String, String> previousPreferences = new LinkedHashMap<>();
    private HashMap<Character, Color> previousNucleotideColors;
    private SAMFileHeader header;
    private AlignmentTrack track;

    @Before
    public void setUp() throws Exception {
        previousHeadless = Globals.isHeadless();
        Globals.setHeadless(true);
        previousGenome = GenomeManager.getInstance().getCurrentGenome();
        byte[] sequence = new byte[5000];
        Arrays.fill(sequence, (byte) 'A');
        GenomeConfig config = new GenomeConfig();
        config.setId("expanded-body-raster-test");
        config.setName("expanded-body-raster-test");
        config.setSequence(new InMemorySequence("chr16", sequence));
        Genome genome = new Genome(config);
        GenomeManager.getInstance().setCurrentGenomeForTest(genome);
        preferences = PreferencesManager.getPreferences();
        for (String key : List.of(SAM_SHOW_SOFT_CLIPPED, SAM_SHOW_CENTER_LINE, SAM_FLAG_CLIPPING,
                SAM_CLIPPING_THRESHOLD, SAM_FLAG_UNMAPPED_PAIR, SAM_FLAG_ZERO_QUALITY, SAM_FLAG_LARGE_INDELS,
                SAM_BASE_QUALITY_MIN, SAM_BASE_QUALITY_MAX)) {
            previousPreferences.put(key, preferences.get(key, null));
        }
        preferences.put(SAM_SHOW_SOFT_CLIPPED, "true");
        preferences.put(SAM_SHOW_CENTER_LINE, "false");
        preferences.put(SAM_FLAG_CLIPPING, "false");
        preferences.put(SAM_CLIPPING_THRESHOLD, "0");
        preferences.put(SAM_FLAG_UNMAPPED_PAIR, "true");
        preferences.put(SAM_FLAG_ZERO_QUALITY, "true");
        preferences.put(SAM_FLAG_LARGE_INDELS, "false");
        preferences.put(SAM_BASE_QUALITY_MIN, "5");
        preferences.put(SAM_BASE_QUALITY_MAX, "20");
        previousNucleotideColors = AlignmentRenderer.nucleotideColors;
        AlignmentRenderer.nucleotideColors = new HashMap<>(previousNucleotideColors);
        AlignmentRenderer.nucleotideColors.put('A', new Color(17, 179, 61));
        AlignmentRenderer.nucleotideColors.put('C', new Color(211, 43, 137));
        AlignmentRenderer.nucleotideColors.put('G', new Color(37, 61, 211));
        AlignmentRenderer.nucleotideColors.put('T', new Color(211, 149, 23));
        header = new SAMFileHeader();
        header.addSequence(new SAMSequenceRecord("chr16", 5000));
        header.setSortOrder(SAMFileHeader.SortOrder.coordinate);
        File bam = temporaryFolder.newFile("expanded-body.bam");
        try (SAMFileWriter writer = new SAMFileWriterFactory().setCreateIndex(true).makeBAMWriter(header, true, bam)) {
            writer.addAlignment(record("fixture", 1000, "40M", false));
        }
        ResourceLocator locator = new ResourceLocator(bam.getAbsolutePath());
        track = new AlignmentTrack(locator, new AlignmentDataManager(locator, genome), genome);
        track.setDisplayMode(Track.DisplayMode.EXPANDED);
        track.setColor(BODY);
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
        if (previousNucleotideColors != null) AlignmentRenderer.nucleotideColors = previousNucleotideColors;
        GenomeManager.getInstance().setCurrentGenomeForTest(previousGenome);
        Globals.setHeadless(previousHeadless);
    }

    @Test
    public void expandedAdjacentOverlappingRevisitedCollapsedReversedAndTerminalBodiesAreExact() {
        Alignment read = blocks("rectangles", false,
                block(1000, 20), block(1020, 10), block(1018, 12), block(1002, 14),
                block(1040, 0), block(1041, 1), block(1050, -3), block(1100, 14));
        for (int type : new int[]{BufferedImage.TYPE_INT_ARGB, BufferedImage.TYPE_INT_ARGB_PRE}) {
            for (Color background : new Color[]{BACKGROUND, new Color(0, 0, 0, 0), new Color(239, 227, 193)}) {
                for (int alpha : new int[]{1, 51, 117, 173, 254, 255}) {
                    track.setColor(new Color(BODY.getRed(), BODY.getGreen(), BODY.getBlue(), alpha));
                    for (double bpp : new double[]{2.5, 8, 31.75}) {
                        for (int deviceScale : new int[]{1, 2}) {
                            assertReplay(List.of(read, read, read), frame(998.25, bpp), row(6), type, background,
                                    g -> g.scale(deviceScale, deviceScale));
                        }
                    }
                }
            }
        }
    }

    @Test
    public void partialAndCurvedDamageIntegerTranslationsAndRowDeviceEdgesAreExact() {
        Alignment read = blocks("damage", false, block(985, 40), block(1000, 40), block(1020, 24),
                block(1001, 8), block(1260, 80));
        for (int scale : new int[]{1, 2}) {
            for (int tx : new int[]{-3, 0, 7}) {
                for (Shape damage : List.of(new Rectangle(0, 2, 1, 8), new Rectangle(9, 3, 24, 4),
                        new Rectangle(96, 1, 4, 9), new Ellipse2D.Double(5.25, 1.25, 36.5, 9.5))) {
                    assertReplay(List.of(read, read), frame(999.75, 3.25), row(6), BufferedImage.TYPE_INT_ARGB_PRE,
                            BACKGROUND, g -> {
                                g.translate(tx, 3);
                                g.scale(scale, scale);
                                g.clip(damage);
                            });
                }
            }
        }
        // Body left clamp is logical zero, not row.x; retain this original geometry.
        assertReplay(List.of(read), frame(999.75, 3.25), new Rectangle(7, 2, 93, 8),
                BufferedImage.TYPE_INT_ARGB, BACKGROUND, g -> g.scale(2, 2));
    }

    @Test
    public void gapsThenAllBodiesThenBasesRemainOrderedAcrossReads() {
        SAMRecord source = record("gap-base", 1000, "24M2I10D24M", false);
        byte[] bases = source.getReadBases();
        bases[3] = 'C';
        bases[29] = 'G';
        source.setReadBases(bases);
        Alignment first = new SAMAlignment(source);
        Alignment second = blocks("overlap-base", false, block(1006, 44), block(1006, 44), block(1020, 12));
        track.getRenderOptions().setShowMismatches(true);
        track.getRenderOptions().setShadeBasesOption(true);
        for (int scale : new int[]{1, 2}) {
            for (int type : new int[]{BufferedImage.TYPE_INT_ARGB, BufferedImage.TYPE_INT_ARGB_PRE}) {
                for (List<Alignment> reads : List.of(List.of(first, second, first), List.of(second, first, second))) {
                    assertReplay(reads, frame(985.25, 4), row(6), type, BACKGROUND, g -> g.scale(scale, scale));
                }
            }
        }
    }

    @Test
    public void pendingRectanglesFlushBeforeMissingBaseSoftClipAndMismatchGraphicsThenResume() {
        AlignmentBlockImpl softClip = block(1012, 18, false, 'S');
        softClip.setSoftClipped(true);
        Alignment read = blocks("alternate-graphics", false,
                block(1000, 44), block(1004, 32), softClip, block(1008, 28),
                block(1016, 20, false, 'X'), block(1002, 46));
        for (int scale : new int[]{1, 2}) {
            assertReplay(List.of(read, read), frame(985.25, 4), row(6), BufferedImage.TYPE_INT_ARGB_PRE,
                    BACKGROUND, g -> g.scale(scale, scale));
        }
    }

    @Test
    public void pendingRectanglesFlushBeforeOverlappingTerminalPolygonAndClippingDecorations() {
        for (boolean negative : new boolean[]{false, true}) {
            track.getRenderOptions().setShowMismatches(true);
            Alignment read = blocks("arrows", negative, mismatchBlock(1000, 44, 5), block(1004, 40), block(1001, 48));
            for (int scale : new int[]{1, 2}) {
                assertReplay(List.of(read, read), frame(970.25, 4), row(10), BufferedImage.TYPE_INT_ARGB_PRE,
                        BACKGROUND, g -> g.scale(scale, scale));
            }
            preferences.put(SAM_FLAG_CLIPPING, "true");
            SAMRecord source = record("clipped-ends", 1000, "8H44M10D44M9H", negative);
            Alignment clipped = new SAMAlignment(source) {
                @Override
                public AlignmentBlock[] getAlignmentBlocks() {
                    return new AlignmentBlock[]{mismatchBlock(1000, 44, 5), block(1004, 40), block(1001, 48)};
                }
            };
            for (int scale : new int[]{1, 2}) {
                assertReplay(List.of(clipped, clipped), frame(970.25, 4), row(10), BufferedImage.TYPE_INT_ARGB,
                        BACKGROUND, g -> g.scale(scale, scale));
            }
            preferences.put(SAM_FLAG_CLIPPING, "false");
        }
    }

    @Test
    public void selectedQualityMateOutlinesAndTexturedDuplicatesRetainOriginalPaints() {
        for (boolean negative : new boolean[]{false, true}) {
            Alignment selected = blocks("selected", negative, block(1000, 30), block(1004, 32), block(1001, 34));
            track.getSelectedReadNames().put(selected.getReadName(), new Color(43, 179, 71, 113));
            assertReplay(List.of(selected, selected), frame(985.25, 4), row(6), BufferedImage.TYPE_INT_ARGB_PRE,
                    BACKGROUND, g -> g.scale(2, 2));
            track.getSelectedReadNames().clear();
            SAMRecord zero = record("zero-quality", 1000, "30M4D30M", negative);
            zero.setMappingQuality(0);
            assertReplay(List.of(new SAMAlignment(zero)), frame(985.25, 4), row(6), BufferedImage.TYPE_INT_ARGB_PRE,
                    BACKGROUND, g -> g.scale(2, 2));
            SAMRecord mate = record("unmapped-mate", 1000, "30M4D30M", negative);
            mate.setReadPairedFlag(true);
            mate.setMateUnmappedFlag(true);
            assertReplay(List.of(new SAMAlignment(mate)), frame(985.25, 4), row(10), BufferedImage.TYPE_INT_ARGB,
                    BACKGROUND, g -> { });
            SAMRecord duplicate = record("textured", 1000, "30M4D30M", negative);
            duplicate.setDuplicateReadFlag(true);
            track.getRenderOptions().setDuplicatesOption(AlignmentTrack.DuplicatesOption.TEXTURE);
            assertReplay(List.of(new SAMAlignment(duplicate), new SAMAlignment(duplicate)), frame(985.25, 4), row(6),
                    BufferedImage.TYPE_INT_ARGB_PRE, BACKGROUND, g -> g.scale(2, 2));
            track.getRenderOptions().setDuplicatesOption(AlignmentTrack.DuplicatesOption.SHOW);
        }
    }

    @Test
    public void expandedAndSquishedDetailedAndNonordinaryModesKeepIndependentGeometry() {
        Alignment read = blocks("fallback", false, block(1000, 20), block(1004, 16), block(1030, 20));
        for (Track.DisplayMode mode : List.of(Track.DisplayMode.EXPANDED, Track.DisplayMode.SQUISHED)) {
            track.setDisplayMode(mode);
            for (double bpp : new double[]{0.5, 1, 4, 100}) {
                assertReplay(List.of(read, read), frame(995.25, bpp), row(6), BufferedImage.TYPE_INT_ARGB_PRE,
                        BACKGROUND, g -> g.scale(2, 2));
            }
        }
        track.setDisplayMode(Track.DisplayMode.SQUISHED);
        // Single-layer thin strokes have an exact primitive oracle; repeated
        // translucent thin bodies retain the existing rounding-bounded tests.
        Color previousBodyColor = track.getColor();
        track.setColor(new Color(BODY.getRed(), BODY.getGreen(), BODY.getBlue()));
        Alignment thin = blocks("thin-fallback", false, block(1000, 16), block(1080, 16), block(1160, 16));
        for (int scale : new int[]{1, 2}) {
            assertReplay(List.of(thin), frame(995.25, 4), row(1), BufferedImage.TYPE_INT_ARGB_PRE,
                    BACKGROUND, g -> g.scale(scale, scale));
        }
        track.setColor(previousBodyColor);
        track.setDisplayMode(Track.DisplayMode.EXPANDED);
        for (AlignmentTrack.ColorOption option : List.of(AlignmentTrack.ColorOption.READ_STRAND,
                AlignmentTrack.ColorOption.BISULFITE, AlignmentTrack.ColorOption.NOMESEQ,
                AlignmentTrack.ColorOption.BASE_MODIFICATION, AlignmentTrack.ColorOption.BASE_MODIFICATION_2COLOR,
                AlignmentTrack.ColorOption.SMRT_CCS_FWD_IPD)) {
            track.getRenderOptions().setColorOption(option);
            assertReplay(List.of(read, read), frame(995.25, 4), row(6), BufferedImage.TYPE_INT_ARGB_PRE,
                    BACKGROUND, g -> g.scale(2, 2));
        }
        track.getRenderOptions().setColorOption(AlignmentTrack.ColorOption.NONE);
        track.getRenderOptions().setQuickConsensusMode(true);
        assertReplay(List.of(read, read), frame(995.25, 4), row(6), BufferedImage.TYPE_INT_ARGB_PRE, BACKGROUND, g -> { });
        track.getRenderOptions().setQuickConsensusMode(false);
        track.getRenderOptions().setHideSmallIndels(true);
        assertReplay(List.of(new SAMAlignment(record("hidden-indels", 1000, "12M1D12M1D12M", false))),
                frame(995.25, 4), row(6), BufferedImage.TYPE_INT_ARGB_PRE, BACKGROUND, g -> g.scale(2, 2));
        track.getRenderOptions().setHideSmallIndels(false);
        track.getRenderOptions().setShowAllBases(true);
        assertReplay(List.of(read, read), frame(995.25, 4), row(6), BufferedImage.TYPE_INT_ARGB_PRE,
                BACKGROUND, g -> g.scale(2, 2));
    }

    @Test
    public void unsupportedTransformsAntialiasingAndNullDamageKeepOriginalRectangles() {
        Alignment read = blocks("graphics-fallback", false, block(1000, 44), block(1004, 36), block(1001, 48));
        for (Consumer<Graphics2D> configure : List.<Consumer<Graphics2D>>of(
                g -> g.scale(1.25, 1.25), g -> g.scale(1, 2), g -> g.translate(0.25, 0.75),
                g -> g.shear(0.1, 0), g -> g.rotate(0.1),
                g -> g.setRenderingHint(RenderingHints.KEY_ANTIALIASING, RenderingHints.VALUE_ANTIALIAS_ON),
                g -> g.setClip(null))) {
            assertReplay(List.of(read, read), frame(985.25, 4), row(6), BufferedImage.TYPE_INT_ARGB_PRE,
                    BACKGROUND, configure);
        }
    }

    @Test
    public void consecutiveWideNarrowEmptyAndCopiedContextPaintsCannotStealPendingWork() {
        Alignment read = blocks("reset", false, block(1000, 40), block(1004, 32), block(1001, 38));
        ReferenceFrame frame = frame(985.25, 4);
        try (Raster actual = new Raster(frame, BufferedImage.TYPE_INT_ARGB_PRE, BACKGROUND, g -> g.scale(2, 2));
             Raster expected = new Raster(frame, BufferedImage.TYPE_INT_ARGB_PRE, BACKGROUND, g -> g.scale(2, 2))) {
            for (Rectangle row : List.of(row(6), new Rectangle(0, 14, 13, 5), new Rectangle(0, 22, 100, 8))) {
                render(List.of(read, read), actual.context, row);
                replay(List.of(read, read), expected.graphics, frame, row);
            }
            Graphics2D pending = actual.context.getGraphics2D("pending-parent-body");
            pending.setColor(new Color(211, 43, 137, 179));
            pending.setComposite(AlphaComposite.getInstance(AlphaComposite.SRC_OVER, 0.37f));
            ExpandedReadBody parentBody = actual.context.getExpandedReadBody();
            Rectangle pendingRow = new Rectangle(0, 22, 100, 8);
            assertTrue(parentBody.reset(pending, pendingRow, 6));
            assertTrue(parentBody.fill(2, 24));
            RenderContext child = new RenderContext(actual.context);
            try {
                render(List.of(read), child, pendingRow);
                replay(List.of(read), expected.graphics, frame, pendingRow);
                assertTrue(parentBody.fill(8, 17));
                parentBody.draw();
                Graphics2D oracle = (Graphics2D) expected.graphics.create();
                try {
                    oracle.setColor(pending.getColor());
                    oracle.setComposite(pending.getComposite());
                    oracle.fillRect(2, pendingRow.y, 22, 6);
                    oracle.fillRect(8, pendingRow.y, 9, 6);
                } finally {
                    oracle.dispose();
                }
            } finally {
                child.dispose();
                child.getGraphics().dispose();
            }
            actual.context.clearGraphicsCache();
            actual.graphics.setClip(new Rectangle(0, 0, 0, 0));
            render(List.of(read), actual.context, row(6));
            ExpandedReadBodyTest.assertExactPixels(expected.image, actual.image);
        }
    }

    private Rectangle row(int height) {
        return new Rectangle(0, 2, WIDTH, height + (track.getDisplayMode() == Track.DisplayMode.SQUISHED ? 0 : 2));
    }

    private void assertReplay(List<Alignment> reads, ReferenceFrame frame, Rectangle row, int type,
                              Color background, Consumer<Graphics2D> configure) {
        try (Raster actual = new Raster(frame, type, background, configure);
             Raster expected = new Raster(frame, type, background, configure)) {
            render(reads, actual.context, row);
            replay(reads, expected.graphics, frame, row);
            ExpandedReadBodyTest.assertExactPixels(expected.image, actual.image);
        }
    }

    private void render(List<Alignment> reads, RenderContext context, Rectangle row) {
        new AlignmentRenderer(track).renderAlignments(reads, null, context, row, track.getRenderOptions());
    }

    // Independent pre-batching body traversal: preserve Math.round endpoints,
    // inclusive gap strokes, per-block original fillRect/polygon, decorations,
    // and only then the read's individual base paints. No candidate code is an oracle.
    private void replay(List<Alignment> reads, Graphics2D source, ReferenceFrame frame, Rectangle row) {
        AlignmentTrack.RenderOptions options = track.getRenderOptions();
        double origin = frame.getOrigin();
        double bpp = frame.getScale();
        double viewEnd = origin + frame.getWidthInPixels() * bpp;
        boolean margin = track.getDisplayMode() != Track.DisplayMode.SQUISHED;
        int height = Math.max(1, row.height - (margin ? 2 : 0));
        int lastPixel = -1;
        for (Alignment read : reads) {
            double pixelStart = (read.getStart() - origin) / bpp;
            double pixelEnd = (read.getEnd() - origin) / bpp;
            if (pixelEnd < row.x || pixelStart > row.getMaxX()) continue;
            Color color = bodyColor(read);
            if (pixelEnd - pixelStart < 2) {
                if (pixelEnd <= lastPixel && color == track.getColor()) continue;
                Graphics2D tiny = (Graphics2D) source.create();
                try {
                    tiny.setColor(color);
                    tiny.setComposite(AlphaComposite.getInstance(AlphaComposite.SRC_OVER, BODY_ALPHA));
                    int tinyHeight = Math.max(1, row.height - 2);
                    int width = Math.max(1, (int) (pixelEnd - pixelStart));
                    tiny.fillRect((int) pixelStart, row.y + (row.height - tinyHeight) / 2, width, tinyHeight);
                    lastPixel = (int) pixelStart + width;
                } finally {
                    tiny.dispose();
                }
                continue;
            }
            replayGaps(source, read, origin, viewEnd, bpp, row, height);
            Color outline = null;
            float outlineWidth = 1;
            if (track.getSelectedReadNames().containsKey(read.getReadName())) {
                color = track.getSelectedReadNames().get(read.getReadName());
                if (color == null) color = Color.BLUE;
                outline = source.getColor();
                outlineWidth = 2;
            } else if (options.isFlagUnmappedPairs() && read.isPaired() && !read.getMate().isMapped()) {
                outline = Color.RED;
            } else if (options.isFlagZeroQualityAlignments() && read.getMappingQuality() == 0) {
                outline = new Color(185, 185, 185);
            }
            boolean textured = read.isDuplicate() && options.getDuplicatesOption() == AlignmentTrack.DuplicatesOption.TEXTURE;
            boolean flagClipping = preferences.getAsBoolean(SAM_FLAG_CLIPPING);
            int threshold = preferences.getAsInt(SAM_CLIPPING_THRESHOLD);
            boolean leftClipped = flagClipping && read.getClippingCounts().getLeft() > threshold;
            boolean rightClipped = flagClipping && read.getClippingCounts().getRight() > threshold;
            AlignmentBlock[] blocks = read.getAlignmentBlocks();
            double pixelLength = read.getLengthOnReference() / bpp;
            int arrow = pixelLength == 0 ? 0 : (int) Math.min(Math.min(5, height / 2), pixelLength / 6);
            int chromStart = blocks[0].getStart();
            boolean first = true;
            for (int index = 0; index < blocks.length; index++) {
                AlignmentBlock block = blocks[index];
                int chromEnd = block.getStart() + block.getLength();
                int start = (int) Math.round((chromStart - origin) / bpp);
                int end = (int) Math.round((chromEnd - origin) / bpp);
                if (end < 0) {
                    if (index + 1 < blocks.length) chromStart = blocks[index + 1].getStart();
                    continue;
                }
                if (start > row.getMaxX()) break;
                boolean last = index + 1 == blocks.length;
                if (!last) {
                    if (options.isHideSmallIndels() && blocks[index + 1].getStart() - chromEnd < options.getSmallIndelThreshold()) continue;
                    chromStart = blocks[index + 1].getStart();
                }
                Graphics2D fill = (Graphics2D) source.create();
                try {
                    fill.setColor(color);
                    fill.setComposite(AlphaComposite.getInstance(AlphaComposite.SRC_OVER, BODY_ALPHA));
                    if (height == 1) {
                        fill.drawLine(start, row.y, end, row.y);
                    } else {
                        int spacing = (int) (AlignmentPacker.MIN_ALIGNMENT_SPACING / bpp);
                        if (spacing < arrow) arrow = Math.max(0, arrow - spacing);
                        boolean pointed = height > 6 && (first && start > 0 || last && end < row.getMaxX());
                        start = Math.max(0, start);
                        end = Math.min(row.x + row.width, end);
                        int leftTip = first && read.isNegativeStrand() && pointed ? arrow : 0;
                        int rightTip = last && !read.isNegativeStrand() && pointed ? arrow : 0;
                        boolean leftDecoration = first && leftClipped;
                        boolean rightDecoration = last && rightClipped;
                        if (!block.hasBases() && (block.isSoftClip() || block.getCigarOperator() == 'X')) {
                            fill.setColor(block.isSoftClip() ? Color.BLACK : Color.RED);
                            fill.setComposite(source.getComposite());
                        }
                        if (leftTip == 0 && rightTip == 0 && outline == null && !textured && !leftDecoration && !rightDecoration && end >= start) {
                            fill.fillRect(start, row.y, end - start, height);
                        } else {
                            Polygon polygon = new Polygon(new int[]{start - leftTip, start, end, end + rightTip, end, start},
                                    new int[]{row.y + height / 2, row.y, row.y, row.y + height / 2, row.y + height, row.y + height}, 6);
                            if (textured) fill.setPaint(texture(fill.getColor()));
                            fill.fill(polygon);
                            Graphics2D decoration = (Graphics2D) source.create();
                            try {
                                if (outline != null) {
                                    decoration.setColor(outline);
                                    decoration.setStroke(new BasicStroke(outlineWidth));
                                    decoration.draw(polygon);
                                }
                                decoration.setColor(new Color(255, 20, 147));
                                if (height > 5) decoration.setStroke(new BasicStroke(1.2f));
                                if (leftDecoration) {
                                    decoration.drawLine(start - leftTip, row.y + height / 2, start, row.y + height);
                                    decoration.drawLine(start, row.y - 1, start - leftTip, row.y + height / 2);
                                }
                                if (rightDecoration) {
                                    decoration.drawLine(end, row.y + height, end + rightTip, row.y + height / 2);
                                    decoration.drawLine(end + rightTip, row.y + height / 2, end, row.y - 1);
                                }
                            } finally {
                                decoration.dispose();
                            }
                        }
                    }
                } finally {
                    fill.dispose();
                }
                first = false;
            }
            replayBases(source, read, origin, viewEnd, bpp, row, height);
            replayInsertions(source, read, origin, bpp, row, height);
        }
    }

    private void replayGaps(Graphics2D source, Alignment read, double origin, double viewEnd, double bpp, Rectangle row, int height) {
        if (read.getGaps() == null) return;
        for (Gap gap : read.getGaps()) {
            int end = gap.getStart() + gap.getnBases();
            if (end <= origin) continue;
            if (gap.getStart() >= Math.ceil(viewEnd)) break;
            if (track.getRenderOptions().isHideSmallIndels() && gap.getnBases() < track.getRenderOptions().getSmallIndelThreshold()) continue;
            Graphics2D graphics = (Graphics2D) source.create();
            try {
                if (gap.getType() == SAMAlignment.UNKNOWN) graphics.setColor(new Color(0, 150, 0));
                else if (gap.getType() == SAMAlignment.SKIPPED_REGION) graphics.setColor(new Color(150, 184, 200));
                else {
                    graphics.setColor(height > 5 ? source.getColor() : AlignmentRenderer.deletionColor);
                    if (height > 5) graphics.setStroke(new BasicStroke(2));
                }
                graphics.drawLine((int) ((Math.max(origin, gap.getStart()) - origin) / bpp), row.y + height / 2,
                        (int) ((Math.min(Math.ceil(viewEnd), end) - origin) / bpp), row.y + height / 2);
            } finally {
                graphics.dispose();
            }
        }
    }

    private void replayBases(Graphics2D source, Alignment read, double origin, double viewEnd, double bpp, Rectangle row, int height) {
        AlignmentTrack.RenderOptions options = track.getRenderOptions();
        if (bpp >= 100 || !options.isShowMismatches() && !options.isShowAllBases()) return;
        Graphics2D graphics = (Graphics2D) source.create();
        try {
            graphics.setComposite(AlphaComposite.getInstance(AlphaComposite.SRC_OVER, BODY_ALPHA));
            for (AlignmentBlock block : read.getAlignmentBlocks()) {
                if (!block.hasBases() || block.getLength() <= 0) continue;
                for (int loc = Math.max((int) Math.floor(origin), block.getStart());
                     loc < Math.min(Math.ceil(viewEnd), block.getEnd()); loc++) {
                    char base = (char) block.getBases().getByte(loc - block.getStart());
                    boolean mismatch = block.isSoftClip() ? base != '=' : base != 'A' && base != '=';
                    if (!options.isShowAllBases() && !mismatch) continue;
                    int x = (int) ((loc - origin) / bpp);
                    int width = (int) Math.max(1, 1 / bpp);
                    if (x > row.getMaxX()) break;
                    if (x + width < row.x) continue;
                    Color color = AlignmentRenderer.nucleotideColors.getOrDefault(base, Color.BLACK);
                    if (options.getShadeBasesOption()) {
                        byte quality = block.getQuality(loc - block.getStart());
                        // Sparse high-quality markers avoid the pre-existing strip's
                        // low-alpha flattening differences, without any tolerance.
                        assertTrue("Fixture mismatch quality is above max=20", quality >= 20);
                    }
                    graphics.setColor(color);
                    graphics.fillRect(x, row.y, width, height);
                }
            }
        } finally {
            graphics.dispose();
        }
    }

    private void replayInsertions(Graphics2D source, Alignment read, double origin, double bpp, Rectangle row, int height) {
        if (read.getInsertions() == null) return;
        for (AlignmentBlock insertion : read.getInsertions()) {
            int x = (int) ((insertion.getStart() - origin) / bpp);
            if (x > row.getMaxX()) break;
            if (x < row.x) continue;
            if (track.getRenderOptions().isHideSmallIndels() && insertion.getBasesLength() < track.getRenderOptions().getSmallIndelThreshold()) continue;
            Graphics2D graphics = (Graphics2D) source.create();
            try {
                graphics.setColor(AlignmentRenderer.purple);
                int y = row.y + (row.height - height) / 2 - (track.getDisplayMode() == Track.DisplayMode.SQUISHED ? 0 : 1);
                int wing = height > 10 ? 2 : height > 5 ? 1 : 0;
                graphics.fillRect(x, y, 2, height);
                graphics.fillRect(x - wing, y, 2 + 2 * wing, 2);
                graphics.fillRect(x - wing, y + height - 2, 2 + 2 * wing, 2);
            } finally {
                graphics.dispose();
            }
        }
    }

    private Color bodyColor(Alignment read) {
        Color color;
        switch (track.getRenderOptions().getColorOption()) {
            case READ_STRAND:
                color = read.isNegativeStrand() ? new Color(150, 150, 230) : new Color(230, 150, 150);
                break;
            case BISULFITE:
            case BASE_MODIFICATION:
            case BASE_MODIFICATION_2COLOR:
            case SMRT_CCS_FWD_IPD:
                color = read.getFirstOfPairStrand() == Strand.POSITIVE ? new Color(195, 195, 195) : new Color(195, 210, 195);
                break;
            case NOMESEQ:
                color = new Color(195, 195, 195);
                break;
            default:
                color = track.getColor();
        }
        if (read.getMappingQuality() == 0 && track.getRenderOptions().isFlagZeroQualityAlignments()) {
            color = ColorUtilities.getCompositeColor(Color.WHITE, color, 0.15f);
        }
        return color;
    }

    private Alignment blocks(String name, boolean negative, AlignmentBlock... blocks) {
        return new SAMAlignment(record(name, 1000, "400M", negative)) {
            @Override
            public AlignmentBlock[] getAlignmentBlocks() {
                return blocks;
            }
        };
    }

    private static AlignmentBlockImpl mismatchBlock(int start, int length, int offset) {
        byte[] sequence = new byte[length];
        byte[] qualities = new byte[length];
        Arrays.fill(sequence, (byte) 'A');
        Arrays.fill(qualities, (byte) 30);
        sequence[offset] = 'C';
        return new AlignmentBlockImpl(start, sequence, qualities);
    }

    private static AlignmentBlockImpl block(int start, int length) {
        return block(start, length, true, 'M');
    }

    private static AlignmentBlockImpl block(int start, int length, boolean bases, char operator) {
        byte[] sequence = new byte[bases ? Math.max(0, length) : 0];
        byte[] qualities = new byte[sequence.length];
        Arrays.fill(sequence, (byte) 'A');
        Arrays.fill(qualities, (byte) 30);
        return new AlignmentBlockImpl(start, sequence, qualities, 0, length, operator);
    }

    private SAMRecord record(String name, int start, String cigar, boolean negative) {
        SAMRecord record = new SAMRecord(header);
        record.setReadName(name);
        record.setReferenceName("chr16");
        record.setAlignmentStart(start + 1);
        record.setMappingQuality(60);
        record.setCigarString(cigar);
        int length = record.getCigar().getReadLength();
        record.setReadString("A".repeat(length));
        byte[] qualities = new byte[length];
        Arrays.fill(qualities, (byte) 30);
        record.setBaseQualities(qualities);
        record.setReadNegativeStrandFlag(negative);
        return record;
    }

    private static ReferenceFrame frame(double origin, double bpp) {
        ReferenceFrame frame = new ReferenceFrame("expanded-body-raster-test") {
            @Override
            public double getScale() {
                return bpp;
            }
        };
        frame.setBounds(0, WIDTH);
        frame.jumpTo("chr16", (int) origin, (int) (origin + WIDTH * bpp));
        frame.setOrigin(origin);
        return frame;
    }

    private static TexturePaint texture(Color color) {
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

    private static final class Raster implements AutoCloseable {
        final BufferedImage image;
        final Graphics2D graphics;
        final RenderContext context;

        Raster(ReferenceFrame frame, int type, Color background, Consumer<Graphics2D> configure) {
            image = new BufferedImage(220, 90, type);
            graphics = image.createGraphics();
            graphics.setComposite(AlphaComposite.Src);
            graphics.setColor(background);
            graphics.fillRect(0, 0, 220, 90);
            graphics.setComposite(AlphaComposite.SrcOver);
            graphics.setColor(Color.BLACK);
            graphics.setRenderingHint(RenderingHints.KEY_ANTIALIASING, RenderingHints.VALUE_ANTIALIAS_OFF);
            graphics.setClip(new Rectangle(0, 0, 220, 90));
            configure.accept(graphics);
            context = new RenderContext(null, graphics, frame, new Rectangle(0, 0, WIDTH, 40));
        }

        @Override
        public void close() {
            context.dispose();
            graphics.dispose();
        }
    }
}
