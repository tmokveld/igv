package org.broad.igv.sam;

import htsjdk.samtools.*;
import org.broad.igv.Globals;
import org.broad.igv.feature.Chromosome;
import org.broad.igv.feature.genome.Genome;
import org.broad.igv.feature.genome.GenomeManager;
import org.broad.igv.feature.genome.InMemorySequence;
import org.broad.igv.feature.genome.load.GenomeConfig;
import org.broad.igv.prefs.IGVPreferences;
import org.broad.igv.prefs.PreferencesManager;
import org.broad.igv.track.RenderContext;
import org.broad.igv.track.Track;
import org.broad.igv.ui.panel.ReferenceFrame;
import org.broad.igv.util.ResourceLocator;
import org.junit.*;
import org.junit.rules.TemporaryFolder;

import java.awt.*;
import java.awt.geom.Ellipse2D;
import java.awt.image.BufferedImage;
import java.io.File;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.function.Consumer;

import static org.broad.igv.prefs.Constants.*;
import static org.junit.Assert.*;

/** Pixel-level comparisons against ordered Java2D paints, not batching internals. */
public class AlignmentRendererBatchingTest {
    private static final Color BACKGROUND = new Color(245, 239, 227);
    private static final int WIDTH = 100;
    private static final byte[] QUALITIES = {0, 4, 5, 10, 15, 19, 20, 30};

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
        byte[] sequence = new byte[5000];
        Arrays.fill(sequence, (byte) 'A');
        // A worked ambiguity fixture: G/R, C/R, T/Y, =/A, C/A, ?/missing, G/N, N/A.
        byte[] ambiguousReference = {'R', 'R', 'Y', 'A', 'A', 0, 'N', 'A'};
        System.arraycopy(ambiguousReference, 0, sequence, 1200, ambiguousReference.length);
        GenomeConfig config = new GenomeConfig();
        config.setId("batching-raster-test");
        config.setName("batching-raster-test");
        config.setSequence(new InMemorySequence("chr16", sequence));
        Genome genome = new Genome(config);
        GenomeManager.getInstance().setCurrentGenomeForTest(genome);

        preferences = PreferencesManager.getPreferences();
        for (String key : List.of(SAM_SHOW_SOFT_CLIPPED, SAM_SHOW_CENTER_LINE, SAM_FLAG_CLIPPING,
                SAM_FLAG_UNMAPPED_PAIR, SAM_FLAG_ZERO_QUALITY, SAM_FLAG_LARGE_INDELS,
                SAM_BASE_QUALITY_MIN, SAM_BASE_QUALITY_MAX)) {
            previousPreferences.put(key, preferences.get(key, null));
        }
        preferences.put(SAM_SHOW_SOFT_CLIPPED, "true");
        preferences.put(SAM_SHOW_CENTER_LINE, "false");
        preferences.put(SAM_FLAG_CLIPPING, "false");
        preferences.put(SAM_FLAG_UNMAPPED_PAIR, "false");
        preferences.put(SAM_FLAG_ZERO_QUALITY, "false");
        preferences.put(SAM_FLAG_LARGE_INDELS, "false");
        preferences.put(SAM_BASE_QUALITY_MIN, "5");
        preferences.put(SAM_BASE_QUALITY_MAX, "20");
        previousNucleotideColors = AlignmentRenderer.nucleotideColors;
        AlignmentRenderer.nucleotideColors = new HashMap<>(previousNucleotideColors);
        AlignmentRenderer.nucleotideColors.put('A', new Color(19, 173, 83));
        AlignmentRenderer.nucleotideColors.put('C', new Color(37, 71, 211));
        AlignmentRenderer.nucleotideColors.put('G', new Color(221, 149, 23));
        AlignmentRenderer.nucleotideColors.put('T', new Color(203, 41, 109));
        AlignmentRenderer.nucleotideColors.put('N', new Color(109, 53, 139));

        header = new SAMFileHeader();
        header.addSequence(new SAMSequenceRecord("chr16", 5000));
        header.setSortOrder(SAMFileHeader.SortOrder.coordinate);
        File bam = temporaryFolder.newFile("batching.bam");
        try (SAMFileWriter writer = new SAMFileWriterFactory().setCreateIndex(true).makeBAMWriter(header, true, bam)) {
            writer.addAlignment(record("fixture", 1000, "40M", "A".repeat(40)));
        }
        ResourceLocator locator = new ResourceLocator(bam.getAbsolutePath());
        track = new AlignmentTrack(locator, new AlignmentDataManager(locator, genome), genome);
        track.setDisplayMode(Track.DisplayMode.EXPANDED);
        track.setColor(new Color(92, 117, 133, 173));
        AlignmentTrack.RenderOptions options = track.getRenderOptions();
        options.setColorOption(AlignmentTrack.ColorOption.NONE);
        options.setShadeAlignmentsOption(AlignmentTrack.ShadeAlignmentsOption.NONE);
        options.setShowMismatches(true);
        options.setShowAllBases(false);
        options.setQuickConsensusMode(false);
        options.setShadeBasesOption(true);
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
    public void overviewCombinesSmallBlocksAndBothClipsInReadOrderWithQualityThresholds() {
        String cigar = "24S" + "1M1D".repeat(40) + "1M24S";
        Alignment first = new SAMAlignment(record("first", 1000, cigar, pattern("CTG=A", 89)));
        Alignment second = new SAMAlignment(record("second", 1002, cigar, pattern("GTC=A", 89)));
        Rectangle row = new Rectangle(0, 2, WIDTH, 6);
        for (double basesPerPixel : new double[]{2.5, 8, 31.75}) {
            ReferenceFrame frame = frame(970.25, basesPerPixel);
            for (double deviceScale : new double[]{1, 2}) {
                Consumer<Graphics2D> transform = g -> g.scale(deviceScale, deviceScale);
                for (List<Alignment> reads : List.of(List.of(first, second, first), List.of(second, first, second))) {
                    assertReplay(reads, frame, row, transform, null, false, 2);
                }
            }
        }
        ReferenceFrame frame = frame(970.25, 2.5);
        try (Raster forward = new Raster(frame, g -> { }); Raster reverse = new Raster(frame, g -> { })) {
            render(List.of(first, second), forward.context, row);
            render(List.of(second, first), reverse.context, row);
            assertTrue("Overlapping read order must remain visible", differentPixels(forward.image, reverse.image));
        }
    }

    @Test
    public void revisitedColumnsRetainEarlierSoftClipAndMismatchContributions() {
        SAMRecord record = record("revisited", 1000, "12S4M4D4M4D4M12S", pattern("CTG=", 36));
        Alignment alignment = new SAMAlignment(record) {
            @Override
            public AlignmentBlock[] getAlignmentBlocks() {
                AlignmentBlock[] blocks = super.getAlignmentBlocks();
                // Keep the same genomic spans, but revisit columns after painting the right clip.
                return new AlignmentBlock[]{blocks[0], blocks[4], blocks[2], blocks[1], blocks[3]};
            }
        };
        for (double scale : new double[]{1, 2}) {
            assertReplay(List.of(alignment), frame(985.25, 8), new Rectangle(0, 2, WIDTH, 6),
                    g -> g.scale(scale, scale), null, false, 2);
        }
    }

    @Test
    public void ambiguousMissingAndEqualsBasesKeepMismatchSemanticsAtOverviewAndDetailedScales() {
        Alignment alignment = new SAMAlignment(record("ambiguity", 1200, "6S8M6S", "N?A=CTGCT=C?GNT=NC?A"));
        Set<Integer> mismatches = Set.of(1201, 1204, 1207);
        for (double basesPerPixel : new double[]{0.5, 1, 4, 8}) {
            for (double scale : new double[]{1, 2}) {
                assertReplay(List.of(alignment), frame(1185.25, basesPerPixel), new Rectangle(0, 2, WIDTH, 6),
                        g -> g.scale(scale, scale), mismatches, false, basesPerPixel <= 1 ? 0 : 2);
            }
        }
        SAMRecord missingQualities = record("missing-qualities", 1200, "6S8M6S", "N?A=CTGCT=C?GNT=NC?A");
        missingQualities.setBaseQualities(SAMRecord.NULL_QUALS);
        assertReplay(List.of(new SAMAlignment(missingQualities)), frame(1185.25, 4),
                new Rectangle(0, 2, WIDTH, 6), g -> { }, mismatches, false, 2);
        SAMRecord missingBases = record("missing-bases", 1200, "20M", "A".repeat(20));
        missingBases.setReadBases(SAMRecord.NULL_SEQUENCE);
        missingBases.setBaseQualities(SAMRecord.NULL_QUALS);
        assertReplay(List.of(new SAMAlignment(missingBases)), frame(1185.25, 4),
                new Rectangle(0, 2, WIDTH, 6), g -> { }, Set.of(), false, 0);
    }

    @Test
    public void missingReferenceSkipsAlignedBasesButStillPaintsSoftClips() {
        GenomeManager.getInstance().setCurrentGenomeForTest(new Genome("no-batching-reference",
                List.of(new Chromosome(0, "chr16", 5000))));
        Alignment alignment = new SAMAlignment(record("no-reference", 1000, "12S40M12S", pattern("CTG=", 64)));
        assertReplay(List.of(alignment), frame(975.25, 4), new Rectangle(0, 2, WIDTH, 6),
                g -> g.scale(2, 2), Set.of(), false, 2);
    }

    @Test
    public void partialDamageAndCurvedClipsMatchFullRenderingExactly() {
        Alignment alignment = new SAMAlignment(record("damage", 1000,
                "24S" + "1M1D".repeat(20) + "1M24S", pattern("GCT=A", 69)));
        for (double scale : new double[]{1, 2}) {
            for (Shape damage : List.of(new Rectangle(4, 3, 21, 3), new Ellipse2D.Double(3.5, 1, 24, 8))) {
                for (double translation : new double[]{0, 0.25}) {
                    Consumer<Graphics2D> transform = g -> {
                        g.translate(translation, translation);
                        g.scale(scale, scale);
                    };
                    ReferenceFrame frame = frame(979.75, 3.25);
                    Rectangle row = new Rectangle(0, 2, WIDTH, 6);
                    try (Raster full = new Raster(frame, transform);
                         Raster actual = new Raster(frame, g -> { transform.accept(g); g.clip(damage); });
                         Raster expected = new Raster(frame, g -> { transform.accept(g); g.clip(damage); })) {
                        render(List.of(alignment, alignment), full.context, row);
                        render(List.of(alignment, alignment), actual.context, row);
                        // Copy the full raster through the same device damage mask.
                        // This checks clipping independently of compositing-rounding
                        // differences already present in the legacy soft-clip strip.
                        expected.graphics.setTransform(new java.awt.geom.AffineTransform());
                        expected.graphics.setComposite(AlphaComposite.Src);
                        expected.graphics.drawImage(full.image, 0, 0, null);
                        assertPixels(expected.image, actual.image, 0);
                    }
                }
            }
        }
    }

    @Test
    public void repeatedSelectedThinBodiesAccumulateEveryInclusiveLineBeforeTheirBases() {
        track.setDisplayMode(Track.DisplayMode.SQUISHED);
        Alignment alignment = new SAMAlignment(record("selected", 1000,
                "1M1D".repeat(20) + "1M", pattern("CTG=", 21)));
        track.getSelectedReadNames().put(alignment.getReadName(), new Color(47, 129, 89, 113));
        Rectangle row = new Rectangle(0, 2, WIDTH, 1);
        for (boolean bases : new boolean[]{false, true}) {
            track.getRenderOptions().setShowMismatches(bases);
            for (double scaleX : new double[]{1, 1.25, 2}) {
                for (double scaleY : new double[]{1, 1.25, 2}) {
                    for (Shape damage : List.of(new Rectangle(2, 1, 30, 3), new Ellipse2D.Double(2, 0, 30, 5))) {
                        assertReplay(List.of(alignment, alignment, alignment), frame(976.25, 3.25), row,
                                g -> { g.translate(0.25, 0.75); g.scale(scaleX, scaleY); g.clip(damage); },
                                null, true, 2);
                    }
                }
            }
        }
        track.getRenderOptions().setShowMismatches(false);
        for (double scale : new double[]{1, 1.25, 2}) {
            assertReplay(List.of(alignment, alignment), frame(976.25, 3.25),
                    new Rectangle(8, 2, 14, 1), g -> g.scale(scale, scale), null, true, 2);
        }
        Alignment continuous = new SAMAlignment(record("selected", 1000, "40M", "A".repeat(40)));
        try (Raster once = new Raster(frame(976.25, 3.25), g -> { });
             Raster repeated = new Raster(frame(976.25, 3.25), g -> { })) {
            render(List.of(continuous), once.context, row);
            render(List.of(continuous, continuous, continuous), repeated.context, row);
            assertNotEquals("Repeated faint bodies are counts, not a coverage union",
                    once.image.getRGB(9, 2), repeated.image.getRGB(9, 2));
        }
    }

    @Test
    public void detailedRenderingRemainsExactAndShowAllBasesKeepsFallbackSemantics() {
        Alignment alignment = new SAMAlignment(record("detailed", 1000, "8S12M8S", pattern("CTG=A", 28)));
        for (double basesPerPixel : new double[]{0.5, 1}) {
            assertReplay(List.of(alignment, alignment), frame(988.25, basesPerPixel), new Rectangle(0, 2, WIDTH, 6),
                    g -> g.scale(2, 2), null, false, 0);
        }
        track.getRenderOptions().setShowAllBases(true);
        assertReplay(List.of(alignment, alignment), frame(988.25, 3.25), new Rectangle(0, 2, WIDTH, 6),
                g -> { }, null, false, 2);
        Alignment alignedOnly = new SAMAlignment(record("all-aligned", 1000, "28M", pattern("CTG=A", 28)));
        assertReplay(List.of(alignedOnly, alignedOnly), frame(988.25, 3.25), new Rectangle(0, 2, WIDTH, 6),
                g -> g.scale(2, 2), null, false, 0);
    }

    @Test
    public void reusedContextResetsScratchBetweenWideNarrowAndEmptyOverlays() {
        Alignment wide = new SAMAlignment(record("wide", 1000, "60S20M60S", pattern("CTG=", 140)));
        Alignment narrow = new SAMAlignment(record("narrow", 1000, "8S20M8S", pattern("GTC=", 36)));
        Alignment empty = new SAMAlignment(record("empty", 1000, "60S20M60S", "=".repeat(140)));
        ReferenceFrame frame = frame(930.25, 2.5);
        try (Raster actual = new Raster(frame, g -> { }); Raster expected = new Raster(frame, g -> { })) {
            Alignment[] reads = {wide, narrow, empty};
            for (int i = 0; i < reads.length; i++) {
                Rectangle row = new Rectangle(0, 2 + i * 10, i == 1 ? 48 : WIDTH, 6);
                render(List.of(reads[i]), actual.context, row);
                replay(List.of(reads[i]), expected, row, null, false);
            }
            assertPixels(expected.image, actual.image, 2);
        }
    }

    @Test
    public void copiedContextsKeepInterleavedPendingStripsIndependentAcrossRendererPaints() {
        ReferenceFrame frame = frame(975, 4);
        Alignment alignment = new SAMAlignment(record("copy", 1000, "12S40M12S", pattern("CTG=", 64)));
        try (Raster actual = new Raster(frame, g -> { }); Raster expected = new Raster(frame, g -> { })) {
            Graphics2D parentGraphics = actual.context.getGraphics2D("pending-parent");
            parentGraphics.setComposite(AlphaComposite.getInstance(AlphaComposite.SRC_OVER, 0.75f));
            BaseRenderer.ColorStrip parentStrip = actual.context.getBaseColorStrip();
            parentStrip.reset(parentGraphics, 2, 4);
            parentStrip.blendColor(2, new Color(203, 41, 109, 51));
            RenderContext copy = new RenderContext(actual.context);
            try {
                // Rendering through the copy must not reset the parent's unfinished overlay.
                render(List.of(alignment), copy, new Rectangle(0, 2, WIDTH, 6));
                replay(List.of(alignment), expected, new Rectangle(0, 2, WIDTH, 6), null, false);
                Graphics2D childGraphics = copy.getGraphics2D("pending-child");
                childGraphics.setComposite(parentGraphics.getComposite());
                BaseRenderer.ColorStrip childStrip = copy.getBaseColorStrip();
                childStrip.reset(childGraphics, 12, 3);
                childStrip.blendColor(12, new Color(37, 71, 211, 102));
                parentStrip.blendColor(4, new Color(19, 173, 83, 178));
                parentStrip.blendColor(2, new Color(221, 149, 23, 102));
                childStrip.blendColor(14, new Color(203, 41, 109, 229));
                childStrip.blendColor(12, new Color(19, 173, 83, 51));
                parentStrip.draw(16, 4);
                childStrip.draw(24, 4);
                try (GraphicsPaint paint = new GraphicsPaint(expected.graphics)) {
                    paint.fill(2, 16, new Color(203, 41, 109, 51));
                    paint.fill(4, 16, new Color(19, 173, 83, 178));
                    paint.fill(2, 16, new Color(221, 149, 23, 102));
                    paint.fill(12, 24, new Color(37, 71, 211, 102));
                    paint.fill(14, 24, new Color(203, 41, 109, 229));
                    paint.fill(12, 24, new Color(19, 173, 83, 51));
                }
                assertPixels(expected.image, actual.image, 2);
            } finally {
                copy.dispose();
                copy.getGraphics().dispose();
            }
        }
    }

    @Test
    public void referenceTileHitsAndMissesPreserveOrdinaryAndDetailedPixelsAcrossTheBoundary() throws Exception {
        assertReferenceTileRendering(List.of(AlignmentTrack.ColorOption.NONE));
    }

    @Test
    public void referenceTileHitsAndMissesPreserveWholeBlockSpecializedContext() throws Exception {
        assertReferenceTileRendering(List.of(AlignmentTrack.ColorOption.BISULFITE,
                AlignmentTrack.ColorOption.NOMESEQ, AlignmentTrack.ColorOption.BASE_MODIFICATION,
                AlignmentTrack.ColorOption.BASE_MODIFICATION_2COLOR, AlignmentTrack.ColorOption.SMRT_CCS_FWD_IPD));
    }

    private void assertReferenceTileRendering(List<AlignmentTrack.ColorOption> colorOptions) throws Exception {
        byte[] sequence = new byte[2_000_000];
        Arrays.fill(sequence, (byte) 'A');
        byte[] context = pattern("ACGCGTNRYS", 80).getBytes(java.nio.charset.StandardCharsets.US_ASCII);
        System.arraycopy(context, 0, sequence, 999_960, context.length);
        InMemorySequence directReference = new InMemorySequence("chr16", sequence);
        GenomeConfig config = new GenomeConfig();
        config.setId("reference-tile-raster-test");
        config.setName("reference-tile-raster-test");
        config.setSequence(directReference);
        Genome cachedGenome = new Genome(config);
        // Independent API seam: the expected renderer bypasses SequenceWrapper entirely.
        Genome directGenome = new Genome(config) {
            @Override
            public byte[] getSequence(String chr, int start, int end) {
                return directReference.getSequence(getCanonicalChrName(chr), start, Math.min(end, sequence.length));
            }
        };
        header = new SAMFileHeader();
        header.addSequence(new SAMSequenceRecord("chr16", sequence.length));
        List<Alignment> reads = new ArrayList<>();
        for (int start : new int[]{999_980, 999_994, 1_000_003}) {
            SAMRecord record = record("tile-" + start, start, "6S12M6S", pattern("TCG=AC", 24));
            record.setAttribute("MM", "C+m,0,0,0;");
            record.setAttribute("ML", new byte[]{(byte) 255, (byte) 128, (byte) 64});
            record.setAttribute("fi", new short[24]);
            reads.add(new SAMAlignment(record));
        }
        for (Track.DisplayMode mode : List.of(Track.DisplayMode.EXPANDED, Track.DisplayMode.SQUISHED)) {
            track.setDisplayMode(mode);
            Rectangle row = new Rectangle(0, 2, WIDTH, mode == Track.DisplayMode.SQUISHED ? 1 : 6);
            for (AlignmentTrack.ColorOption colorOption : colorOptions) {
                track.getRenderOptions().setColorOption(colorOption);
                for (double scale : new double[]{0.1, 0.5, 1, 4}) {
                    ReferenceFrame frame = frame(999_972.25, scale);
                    for (int paint = 0; paint < 2; paint++) {
                        try (Raster expected = new Raster(frame, g -> g.scale(2, 2));
                             Raster actual = new Raster(frame, g -> g.scale(2, 2))) {
                            GenomeManager.getInstance().setCurrentGenomeForTest(directGenome);
                            render(reads, expected.context, row);
                            GenomeManager.getInstance().setCurrentGenomeForTest(cachedGenome);
                            render(reads, actual.context, row);
                            assertPixels(expected.image, actual.image, 0);
                        }
                    }
                }
            }
        }
    }

    private void assertReplay(List<Alignment> reads, ReferenceFrame frame, Rectangle row,
                              Consumer<Graphics2D> configuration, Set<Integer> mismatches,
                              boolean manualThinBodies, int tolerance) {
        try (Raster actual = new Raster(frame, configuration); Raster expected = new Raster(frame, configuration)) {
            render(reads, actual.context, row);
            replay(reads, expected, row, mismatches, manualThinBodies);
            assertPixels(expected.image, actual.image, tolerance);
        }
    }

    // Reuse actual body/gap/decorations rendering with bases disabled. The independent
    // oracle is a sequence of individual fillRect paints, applied AFTER each read's body.
    // Thin bodies instead replay drawLine, retaining Java2D's device-space square caps.
    private void replay(List<Alignment> reads, Raster raster, Rectangle row,
                        Set<Integer> mismatches, boolean manualThinBodies) {
        AlignmentTrack.RenderOptions options = track.getRenderOptions();
        boolean showMismatches = options.isShowMismatches();
        boolean showAll = options.isShowAllBases();
        try {
            for (Alignment read : reads) {
                options.setShowMismatches(false);
                options.setShowAllBases(false);
                if (manualThinBodies) {
                    try (GraphicsPaint paint = new GraphicsPaint(raster.graphics)) {
                        paint.graphics.setColor(track.getSelectedReadNames().getOrDefault(read.getReadName(), track.getColor()));
                        // Deletion connectors retain their independent opaque paint.
                        if (read.getGaps() != null) {
                            try (GraphicsPaint gaps = new GraphicsPaint(raster.graphics)) {
                                gaps.graphics.setComposite(AlphaComposite.SrcOver);
                                gaps.graphics.setColor(AlignmentRenderer.deletionColor);
                                for (Gap gap : read.getGaps()) {
                                    int left = (int) ((Math.max(raster.context.getOrigin(), gap.getStart()) - raster.context.getOrigin()) / raster.context.getScale());
                                    int right = (int) ((Math.min(Math.ceil(raster.context.getEndLocation()), gap.getStart() + gap.getnBases()) - raster.context.getOrigin()) / raster.context.getScale());
                                    gaps.graphics.drawLine(left, row.y, right, row.y);
                                }
                            }
                        }
                        for (AlignmentBlock block : read.getAlignmentBlocks()) {
                            int left = (int) Math.round((block.getStart() - raster.context.getOrigin()) / raster.context.getScale());
                            int right = (int) Math.round((block.getEnd() - raster.context.getOrigin()) / raster.context.getScale());
                            if (right >= 0 && left <= row.getMaxX()) paint.graphics.drawLine(left, row.y, right, row.y);
                        }
                    }
                } else {
                    render(List.of(read), raster.context, row);
                }
                if (!showMismatches && !showAll) continue;
                int height = row.height - (track.getDisplayMode() == Track.DisplayMode.SQUISHED ? 0 : 2);
                int width = (int) Math.max(1, 1 / raster.context.getScale());
                try (GraphicsPaint paint = new GraphicsPaint(raster.graphics)) {
                    for (AlignmentBlock block : read.getAlignmentBlocks()) {
                        if (!block.hasBases()) continue;
                        int first = Math.max(block.getStart(), (int) Math.floor(raster.context.getOrigin()));
                        int end = Math.min(block.getEnd(), (int) Math.ceil(raster.context.getEndLocation()));
                        for (int loc = first; loc < end; loc++) {
                            char base = (char) block.getBases().getByte(loc - block.getStart());
                            boolean mismatch = block.isSoftClip() ? base != '=' :
                                    mismatches == null ? base != 'A' && base != '=' : mismatches.contains(loc);
                            if (!showAll && !mismatch) continue;
                            int x = (int) ((loc - raster.context.getOrigin()) / raster.context.getScale());
                            if (x > row.getMaxX() || x + width < row.x) continue;
                            Color color = AlignmentRenderer.nucleotideColors.getOrDefault(base, Color.BLACK);
                            if (options.getShadeBasesOption()) color = referenceQualityColor(color, block.getQuality(loc - block.getStart()));
                            paint.graphics.setColor(color);
                            paint.graphics.fillRect(x, row.y, width, height);
                        }
                    }
                }
            }
        } finally {
            options.setShowMismatches(showMismatches);
            options.setShowAllBases(showAll);
        }
    }

    // Explicit expectations for min=5/max=20 and the qualities used by this fixture.
    private static Color referenceQualityColor(Color color, byte quality) {
        if (quality >= 20) return color;
        int alpha;
        switch (quality) {
            case 0: case 4: case 5: alpha = 51; break;
            case 10: alpha = 102; break;
            case 15: alpha = 178; break;
            case 19: alpha = 229; break;
            default: throw new AssertionError("Unexpected fixture quality: " + quality);
        }
        return new Color(color.getRed(), color.getGreen(), color.getBlue(), alpha);
    }

    private void render(List<Alignment> reads, RenderContext context, Rectangle row) {
        new AlignmentRenderer(track).renderAlignments(reads, null, context, row, track.getRenderOptions());
    }

    private SAMRecord record(String name, int start, String cigar, String bases) {
        SAMRecord record = new SAMRecord(header);
        record.setReadName(name);
        record.setReferenceName("chr16");
        record.setAlignmentStart(start + 1);
        record.setMappingQuality(60);
        record.setCigarString(cigar);
        record.setReadString(bases);
        byte[] qualities = new byte[bases.length()];
        for (int i = 0; i < qualities.length; i++) qualities[i] = QUALITIES[i % QUALITIES.length];
        record.setBaseQualities(qualities);
        return record;
    }

    private static String pattern(String motif, int length) {
        return motif.repeat((length + motif.length() - 1) / motif.length()).substring(0, length);
    }

    private ReferenceFrame frame(double origin, double scale) {
        ReferenceFrame frame = new ReferenceFrame("batching-raster-test") {
            @Override
            public double getScale() {
                return scale;
            }
        };
        frame.setBounds(0, WIDTH);
        frame.jumpTo("chr16", (int) origin, (int) (origin + WIDTH * scale));
        frame.setOrigin(origin);
        return frame;
    }

    private static boolean differentPixels(BufferedImage first, BufferedImage second) {
        for (int y = 0; y < first.getHeight(); y++) {
            for (int x = 0; x < first.getWidth(); x++) {
                if (first.getRGB(x, y) != second.getRGB(x, y)) return true;
            }
        }
        return false;
    }

    private static void assertPixels(BufferedImage expected, BufferedImage actual, int tolerance) {
        for (int y = 0; y < expected.getHeight(); y++) {
            for (int x = 0; x < expected.getWidth(); x++) {
                int reference = expected.getRGB(x, y);
                int observed = actual.getRGB(x, y);
                String pixel = "Pixel " + x + "," + y;
                assertEquals(pixel + " alpha", reference >>> 24, observed >>> 24);
                assertEquals(pixel + " paint coverage", reference != BACKGROUND.getRGB(), observed != BACKGROUND.getRGB());
                if (tolerance == 0 || reference == BACKGROUND.getRGB()) {
                    assertEquals(pixel, reference, observed);
                } else {
                    for (int shift : new int[]{16, 8, 0}) {
                        int difference = Math.abs(((reference >>> shift) & 255) - ((observed >>> shift) & 255));
                        assertTrue(pixel + " channel " + shift + ": expected " + new Color(reference, true)
                                + ", got " + new Color(observed, true), difference <= tolerance);
                    }
                }
            }
        }
    }

    private static final class GraphicsPaint implements AutoCloseable {
        final Graphics2D graphics;

        GraphicsPaint(Graphics2D source) {
            graphics = (Graphics2D) source.create();
            graphics.setComposite(AlphaComposite.getInstance(AlphaComposite.SRC_OVER, 0.75f));
        }

        void fill(int x, int y, Color color) {
            graphics.setColor(color);
            graphics.fillRect(x, y, 1, 4);
        }

        @Override
        public void close() {
            graphics.dispose();
        }
    }

    private static final class Raster implements AutoCloseable {
        final BufferedImage image = new BufferedImage(205, 70, BufferedImage.TYPE_INT_ARGB_PRE);
        final Graphics2D graphics = image.createGraphics();
        final RenderContext context;

        Raster(ReferenceFrame frame, Consumer<Graphics2D> configuration) {
            graphics.setColor(BACKGROUND);
            graphics.fillRect(0, 0, image.getWidth(), image.getHeight());
            graphics.setRenderingHint(RenderingHints.KEY_ANTIALIASING, RenderingHints.VALUE_ANTIALIAS_OFF);
            graphics.setRenderingHint(RenderingHints.KEY_INTERPOLATION, RenderingHints.VALUE_INTERPOLATION_BICUBIC);
            configuration.accept(graphics);
            // Swing provides an explicit viewport/damage clip; thin batching needs it.
            graphics.clip(new Rectangle(0, 0, WIDTH, 32));
            context = new RenderContext(null, graphics, frame, new Rectangle(0, 0, WIDTH, 32));
        }

        @Override
        public void close() {
            context.dispose();
            graphics.dispose();
        }
    }
}
