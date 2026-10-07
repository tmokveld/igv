/*
 * The MIT License (MIT)
 *
 * Copyright (c) 2007-2015 Broad Institute
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 *
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
 * THE SOFTWARE.
 */

package org.broad.igv.sam;

import com.google.common.primitives.Ints;
import org.broad.igv.Globals;
import org.broad.igv.feature.Strand;
import org.broad.igv.feature.genome.Genome;
import org.broad.igv.feature.genome.GenomeManager;
import org.broad.igv.logging.LogManager;
import org.broad.igv.logging.Logger;
import org.broad.igv.prefs.IGVPreferences;
import org.broad.igv.prefs.PreferencesManager;
import org.broad.igv.renderer.GraphicUtils;
import org.broad.igv.sam.AlignmentTrack.ColorOption;
import org.broad.igv.sam.BisulfiteBaseInfo.DisplayStatus;
import org.broad.igv.sam.mods.BaseModificationRenderer;
import org.broad.igv.sam.smrt.SMRTKineticsRenderer;
import org.broad.igv.track.RenderContext;
import org.broad.igv.track.Track;
import org.broad.igv.ui.FontManager;
import org.broad.igv.ui.color.*;
import org.broad.igv.ultima.render.ColorByTagValueList;
import org.broad.igv.ultima.render.FlowIndelRendering;
import org.broad.igv.util.ChromosomeColors;

import java.awt.*;
import java.awt.geom.AffineTransform;
import java.awt.geom.Rectangle2D;
import java.awt.image.BufferedImage;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Objects;

import static org.broad.igv.prefs.Constants.*;

/**
 * @author jrobinso
 */
public class AlignmentRenderer {

    private static final Logger log = LogManager.getLogger(AlignmentRenderer.class);


    private static final Color negStrandColor = new Color(150, 150, 230);
    private static final Color posStrandColor = new Color(230, 150, 150);

    private static final Color firstOfPairColor = new Color(150, 150, 230);
    private static final Color secondOfPairColor = new Color(230, 150, 150);
    private static final Color firstAndSecondofPairColor = new Color(150, 150, 0);
    private static final Color neitherForOrSecondOfPair = new Color(198, 106, 245, 255);

    private static final Color RL_COLOR = new Color(0, 150, 0);
    private static final Color RR_COLOR = new Color(20, 50, 200);
    private static final Color LL_COLOR = new Color(0, 150, 150);
    private static final Color smallISizeColor = new Color(0, 0, 150);
    private static final Color largeISizeColor = new Color(200, 0, 0);
    private static final Color OUTLINE_COLOR = new Color(185, 185, 185);

    //chimeric read connecting on the same contig colors
    private static final Color INVERSION_COLOR = new Color(200, 0, 0);
    private static final Color NON_INVERSION_COLOR = new Color(230, 150, 150);

    // Clipping colors
    private static final Color clippedColor = new Color(255, 20, 147);

    // Indel colors
    public static Color purple = new Color(118, 24, 220);
    public static Color deletionColor = Color.black;
    private static Color skippedColor = new Color(150, 184, 200);
    private static Color unknownGapColor = new Color(0, 150, 0);

    // Bisulfite colors
    private static final Color bisulfiteColorFw1 = new Color(195, 195, 195);
    private static final Color bisulfiteColorRev1 = new Color(195, 210, 195);
    private static final Color nomeseqColor = new Color(195, 195, 195);

    private static Map<String, AlignmentTrack.OrientationType> frOrientationTypes;
    private static Map<String, AlignmentTrack.OrientationType> f1f2OrientationTypes;
    private static Map<String, AlignmentTrack.OrientationType> f2f1OrientationTypes;
    private static Map<String, AlignmentTrack.OrientationType> rfOrientationTypes;
    private static Map<AlignmentTrack.OrientationType, Color> typeToColorMap;
    private static final Map<Color, TexturePaint> duplicateTextures = new HashMap<>();

    public static final HSLColorTable tenXColorTable1 = new HSLColorTable(30);
    public static final HSLColorTable tenXColorTable2 = new HSLColorTable(270);
    public static final GreyscaleColorTable tenXColorTable3 = new GreyscaleColorTable();

    public static final Color GROUP_DIVIDER_COLOR = new Color(200, 200, 200);
    // A "dummy" reference for soft-clipped reads.
    private static final byte[] softClippedReference = new byte[1000];

    private static ColorTable readGroupColors;
    private static ColorTable sampleColors;
    private static ColorTable movieColors;
    private static ColorTable zmwColors;
    private static Map<String, ColorTable> tagValueColors;
    private static ColorTable defaultTagColors;
    public static volatile HashMap<Character, Color> nucleotideColors;
    private static String nucleotideA, nucleotideC, nucleotideT, nucleotideG, nucleotideN;

    final private static ColorByTagValueList colorByTagValueList = new ColorByTagValueList();
    final private static FlowIndelRendering flowIndelRendering = new FlowIndelRendering();

    private static void initializeTagTypes() {
        // pre-seed from orientation colors

        // fr Orientations (e.g. Illumina paired-end libraries)
        frOrientationTypes = new HashMap<>();
        //LR
        frOrientationTypes.put("F1R2", AlignmentTrack.OrientationType.LR);
        frOrientationTypes.put("F2R1", AlignmentTrack.OrientationType.LR);
        frOrientationTypes.put("F R ", AlignmentTrack.OrientationType.LR);
        frOrientationTypes.put("FR", AlignmentTrack.OrientationType.LR);
        //LL
        frOrientationTypes.put("F1F2", AlignmentTrack.OrientationType.LL);
        frOrientationTypes.put("F2F1", AlignmentTrack.OrientationType.LL);
        frOrientationTypes.put("F F ", AlignmentTrack.OrientationType.LL);
        frOrientationTypes.put("FF", AlignmentTrack.OrientationType.LL);
        //RR
        frOrientationTypes.put("R1R2", AlignmentTrack.OrientationType.RR);
        frOrientationTypes.put("R2R1", AlignmentTrack.OrientationType.RR);
        frOrientationTypes.put("R R ", AlignmentTrack.OrientationType.RR);
        frOrientationTypes.put("RR", AlignmentTrack.OrientationType.RR);
        //RL
        frOrientationTypes.put("R1F2", AlignmentTrack.OrientationType.RL);
        frOrientationTypes.put("R2F1", AlignmentTrack.OrientationType.RL);
        frOrientationTypes.put("R F ", AlignmentTrack.OrientationType.RL);
        frOrientationTypes.put("RF", AlignmentTrack.OrientationType.RL);

        // rf orienation  (e.g. Illumina mate-pair libraries)
        rfOrientationTypes = new HashMap<>();
        //LR
        rfOrientationTypes.put("R1F2", AlignmentTrack.OrientationType.LR);
        rfOrientationTypes.put("R2F1", AlignmentTrack.OrientationType.LR);
        //rfOrientationTypes.put("R F ", AlignmentTrack.OrientationType.LR);
        //rfOrientationTypes.put("RF", AlignmentTrack.OrientationType.LR);
        //LL
        rfOrientationTypes.put("R1R2", AlignmentTrack.OrientationType.LL);
        rfOrientationTypes.put("R2R1", AlignmentTrack.OrientationType.LL);
        rfOrientationTypes.put("R R ", AlignmentTrack.OrientationType.LL);
        rfOrientationTypes.put("RR ", AlignmentTrack.OrientationType.LL);

        rfOrientationTypes.put("F1F2", AlignmentTrack.OrientationType.RR);
        rfOrientationTypes.put("F2F1", AlignmentTrack.OrientationType.RR);
        rfOrientationTypes.put("F F ", AlignmentTrack.OrientationType.RR);
        rfOrientationTypes.put("FF", AlignmentTrack.OrientationType.RR);
        //RL
        rfOrientationTypes.put("F1R2", AlignmentTrack.OrientationType.RL);
        rfOrientationTypes.put("F2R1", AlignmentTrack.OrientationType.RL);
        rfOrientationTypes.put("F R ", AlignmentTrack.OrientationType.RL);
        rfOrientationTypes.put("FR", AlignmentTrack.OrientationType.RL);

        // f1f2 orienation  (e.g. SOLID libraries, AlignmentTrack.OrientationType.second read appears first on + strand (leftmost))
        f2f1OrientationTypes = new HashMap<>();
        //LR
        f2f1OrientationTypes.put("F2F1", AlignmentTrack.OrientationType.LR);
        f2f1OrientationTypes.put("R1R2", AlignmentTrack.OrientationType.LR);

        //LL
        f2f1OrientationTypes.put("F2R1", AlignmentTrack.OrientationType.LL);
        f2f1OrientationTypes.put("R1F2", AlignmentTrack.OrientationType.LL);

        //RR
        f2f1OrientationTypes.put("R2F1", AlignmentTrack.OrientationType.RR);
        f2f1OrientationTypes.put("F1R2", AlignmentTrack.OrientationType.RR);

        //RL
        f2f1OrientationTypes.put("R2R1", AlignmentTrack.OrientationType.RL);
        f2f1OrientationTypes.put("F1F2", AlignmentTrack.OrientationType.RL);

        // f1f2 orientation  (e.g. SOLID libraries, AlignmentTrack.OrientationType.actually is this one even possible?)
        f1f2OrientationTypes = new HashMap<>();
        //LR
        f1f2OrientationTypes.put("F1F2", AlignmentTrack.OrientationType.LR);
        f1f2OrientationTypes.put("R2R1", AlignmentTrack.OrientationType.LR);
        //LL
        f1f2OrientationTypes.put("F1R2", AlignmentTrack.OrientationType.LL);
        f1f2OrientationTypes.put("R2F1", AlignmentTrack.OrientationType.LL);
        //RR
        f1f2OrientationTypes.put("R1F2", AlignmentTrack.OrientationType.RR);
        f1f2OrientationTypes.put("F2R1", AlignmentTrack.OrientationType.RR);
        //RL
        f1f2OrientationTypes.put("R1R2", AlignmentTrack.OrientationType.RL);
        f1f2OrientationTypes.put("F2F1", AlignmentTrack.OrientationType.RL);
    }

    private static void initializeTagColors() {
        ColorPalette palette = ColorUtilities.getPalette("Pastel 1");  // TODO let user choose
        readGroupColors = new PaletteColorTable(palette);
        sampleColors = new PaletteColorTable(palette);
        movieColors = new PaletteColorTable(palette);
        zmwColors = new PaletteColorTable(palette);
        defaultTagColors = new PaletteColorTable(palette);
        tagValueColors = new HashMap<>();

        typeToColorMap = new HashMap<>(5);
        typeToColorMap.put(AlignmentTrack.OrientationType.LL, LL_COLOR);
        typeToColorMap.put(AlignmentTrack.OrientationType.RL, RL_COLOR);
        typeToColorMap.put(AlignmentTrack.OrientationType.RR, RR_COLOR);
    }

    private static synchronized HashMap<Character, Color> getNucleotideColors() {
        IGVPreferences prefs = PreferencesManager.getPreferences();
        String aPreference = prefs.get(SAM_COLOR_A);
        String cPreference = prefs.get(SAM_COLOR_C);
        String tPreference = prefs.get(SAM_COLOR_T);
        String gPreference = prefs.get(SAM_COLOR_G);
        String nPreference = prefs.get(SAM_COLOR_N);
        if (nucleotideColors != null &&
                Objects.equals(nucleotideA, aPreference) && Objects.equals(nucleotideC, cPreference) &&
                Objects.equals(nucleotideT, tPreference) && Objects.equals(nucleotideG, gPreference) &&
                Objects.equals(nucleotideN, nPreference)) return nucleotideColors;

        Color a = ColorUtilities.stringToColor(aPreference, Color.green);
        Color c = ColorUtilities.stringToColor(cPreference, Color.blue);
        Color t = ColorUtilities.stringToColor(tPreference, Color.red);
        Color g = ColorUtilities.stringToColor(gPreference, new Color(209, 113, 5));
        Color n = ColorUtilities.stringToColor(nPreference, new Color(64, 64, 64));
        HashMap<Character, Color> colors = new HashMap<>();
        colors.put('A', a);
        colors.put('a', a);
        colors.put('C', c);
        colors.put('c', c);
        colors.put('T', t);
        colors.put('t', t);
        colors.put('G', g);
        colors.put('g', g);
        colors.put('N', n);
        colors.put('n', n);
        colors.put('-', Color.lightGray);
        nucleotideA = aPreference;
        nucleotideC = cPreference;
        nucleotideT = tPreference;
        nucleotideG = gPreference;
        nucleotideN = nPreference;
        // Publish only a fully populated palette. Existing paints retain their snapshot.
        nucleotideColors = colors;
        return colors;

    }

    static {
        initializeTagTypes();
        getNucleotideColors();
        initializeTagColors();
    }


    AlignmentTrack track;
    private volatile BaseRenderer.PaletteColors baseColors;

    public AlignmentRenderer(AlignmentTrack track) {
        this.track = track;
    }

    private void initializeGraphics(RenderContext context) {

        Font font = FontManager.getFont(10);

        Graphics2D g1 = context.getGraphics2D("ALIGNMENT");
        float alpha = 0.75f;
        int type = AlphaComposite.SRC_OVER;
        Composite alignmentAlphaComposite = AlphaComposite.getInstance(type, alpha);
        g1.setComposite(alignmentAlphaComposite);
        g1.setFont(font);

        Graphics2D g2 = context.getGraphics2D("LINK_LINE");
        alpha = 0.3f;
        alignmentAlphaComposite = AlphaComposite.getInstance(type, alpha);
        g2.setComposite(alignmentAlphaComposite);

        Graphics2D g3 = context.getGraphics2D("THICK_STROKE");
        g3.setStroke(new BasicStroke(2.0f));

        Graphics2D g4 = context.getGraphics2D("STRAND");
        g4.setColor(Color.DARK_GRAY);

        Graphics2D g5 = context.getGraphics2D("SOFT_CLIP");
        g5.setColor(Color.BLACK);

        Graphics2D g6 = context.getGraphics2D("MISMATCH");
        g6.setColor(Color.RED);
    }

    /**
     * Render a row of alignments in the given rectangle.
     */
    public void renderAlignments(List<Alignment> alignments,
                                 AlignmentCounts alignmentCounts,
                                 RenderContext context,
                                 Rectangle rowRect,
                                 AlignmentTrack.RenderOptions renderOptions
    ) {

        initializeGraphics(context);
        HashMap<Character, Color> palette = getNucleotideColors();
        boolean shadeBases = renderOptions.getShadeBasesOption();
        int minQ = renderOptions.getBaseQualityMin();
        int maxQ = renderOptions.getBaseQualityMax();
        BaseRenderer.PaletteColors colors = baseColors;
        if (colors == null || !colors.matches(palette, shadeBases, minQ, maxQ)) {
            colors = new BaseRenderer.PaletteColors(palette, shadeBases, minQ, maxQ);
            baseColors = colors;
        }

        final boolean leaveMargin = this.track.getDisplayMode() != Track.DisplayMode.SQUISHED;
        ThinReadBody thinReadBody = Math.max(1, rowRect.height - (leaveMargin ? 2 : 0)) == 1
                ? ThinReadBody.create(context, rowRect, overviewBatching(context, renderOptions)) : null;

        double origin = context.getOrigin();
        double locScale = context.getScale();

        Color defaultColor = track.getColor();

        List<Alignment> trimmedAlignments;

        if (renderOptions.isHideTailSbx()) {
            trimmedAlignments = new ArrayList<Alignment>();
            for (Alignment al : alignments) {
                if (al instanceof SAMAlignment) {
                    Alignment trimmed = al.trimSimplexTails();
                    if (trimmed.getAlignmentStart() < trimmed.getAlignmentEnd()) {
                        trimmedAlignments.add(trimmed);
                    }
                }
            }

        } else {
            trimmedAlignments = alignments;
        }

        if ((trimmedAlignments != null) && (trimmedAlignments.size() > 0)) {

            int lastPixelDrawn = -1;

            for (Alignment alignment : trimmedAlignments) {
                // Compute the start and end of the alignment in pixels
                double pixelStart = ((alignment.getStart() - origin) / locScale);
                double pixelEnd = ((alignment.getEnd() - origin) / locScale);

                // If any part of the feature fits in the track rectangle draw it
                if (pixelEnd < rowRect.x || pixelStart > rowRect.getMaxX()) {
                    continue;
                }


                // If the alignment is 3 pixels or less,  draw alignment as a single block,
                // further detail would not be seen and just add to drawing overhead
                // Does the change for Bisulfite kill some machines?
                double pixelWidth = pixelEnd - pixelStart;
                Color alignmentColor = getAlignmentColor(alignment, track);
                final ColorOption colorOption = renderOptions.getColorOption();
                if ((pixelWidth < 2) &&
                        !((AlignmentTrack.isBisulfiteColorType(colorOption) ||
                                colorOption.isBaseMod() ||
                                colorOption.isSMRTKinetics()) &&
                                (pixelWidth >= 1))) {
                    // Optimization for really zoomed out views.  If this alignment occupies screen space already taken,
                    // and it is the default color, skip drawing.
                    if (pixelEnd <= lastPixelDrawn && alignmentColor == defaultColor) {
                        continue;
                    }
                    Graphics2D g = context.getGraphics2D("ALIGNMENT");
                    g.setColor(alignmentColor);
                    int w = Math.max(1, (int) (pixelWidth));
                    int h = (int) Math.max(1, rowRect.getHeight() - 2);
                    int y = (int) (rowRect.getY() + (rowRect.getHeight() - h) / 2);
                    g.fillRect((int) pixelStart, y, w, h);
                    lastPixelDrawn = (int) pixelStart + w;
                } else if (alignment instanceof PairedAlignment) {
                    drawPairedAlignment((PairedAlignment) alignment, rowRect, context, renderOptions, leaveMargin, alignmentCounts, thinReadBody, colors);
                } else if (alignment instanceof LinkedAlignment) {
                    drawLinkedAlignment((LinkedAlignment) alignment, rowRect, context, renderOptions, leaveMargin, alignmentCounts, thinReadBody, colors);
                } else {
                    drawAlignment(alignment, rowRect, context, alignmentColor, renderOptions, leaveMargin, alignmentCounts, false, thinReadBody, colors);
                }
            }

            // Optionally draw a border around the center base
            IGVPreferences prefs = this.track.getPreferences();
            boolean showCenterLine = prefs.getAsBoolean(SAM_SHOW_CENTER_LINE);
            final int bottom = rowRect.y + rowRect.height;
            if (showCenterLine) {
                // Calculate center lines
                double center = (int) (context.getReferenceFrame().getCenter() - origin);
                int centerLeftP = (int) (center / locScale);
                int centerRightP = (int) ((center + 1) / locScale);
                //float transparency = Math.max(0.5f, (float) Math.round(10 * (1 - .75 * locScale)) / 10);
                Graphics2D g = context.getGraphics();
                g.setColor(Color.black);
                GraphicUtils.drawDottedDashLine(g, centerLeftP, rowRect.y, centerLeftP, bottom);
                if ((centerRightP - centerLeftP > 2)) {
                    GraphicUtils.drawDottedDashLine(g, centerRightP, rowRect.y, centerRightP, bottom);
                }
            }
        }
    }

    private void drawLinkedAlignment(LinkedAlignment alignment, Rectangle rowRect, RenderContext context,
                                     AlignmentTrack.RenderOptions renderOptions, boolean leaveMargin,
                                     AlignmentCounts alignmentCounts, ThinReadBody thinReadBody, BaseRenderer.PaletteColors colors) {

        double origin = context.getOrigin();
        double locScale = context.getScale();
        Color alignmentColor = getAlignmentColor(alignment, track);
        List<Alignment> barcodedAlignments = alignment.alignments;

        if (barcodedAlignments.size() > 0) {
            boolean mixedStrand = alignment.getStrand() == Strand.NONE;
            Alignment firstAlignment = barcodedAlignments.get(0);
            if (barcodedAlignments.size() > 1) {
                Graphics2D gline = context.getGraphics2D("LINK_LINE");
                gline.setColor(alignmentColor);
                int startX = (int) ((firstAlignment.getEnd() - origin) / locScale);
                int endX = (int) ((barcodedAlignments.get(barcodedAlignments.size() - 1).getStart() - origin) / locScale);
                int h = (int) Math.max(1, rowRect.getHeight() - (leaveMargin ? 2 : 0));
                int y = (int) (rowRect.getY());
                startX = Math.max(rowRect.x, startX);
                endX = Math.min(rowRect.x + rowRect.width, endX);
                gline.drawLine(startX, y + h / 2, endX, y + h / 2);
            }
            for (int i = 0; i < barcodedAlignments.size(); i++) {
                boolean overlapped = false;
                Alignment al = barcodedAlignments.get(i);
                if (al.isNegativeStrand()) {
                    if (mixedStrand) alignmentColor = negStrandColor;
                    overlapped = (i > 0 && barcodedAlignments.get(i - 1).getAlignmentEnd() > al.getAlignmentStart());
                } else {
                    if (mixedStrand) alignmentColor = posStrandColor;
                    overlapped = i < barcodedAlignments.size() - 1 && al.getAlignmentEnd() > barcodedAlignments.get(i + 1).getAlignmentStart();
                }
                drawAlignment(al, rowRect, context, alignmentColor, renderOptions, leaveMargin, alignmentCounts, overlapped, thinReadBody, colors);
            }
        }
    }


    /**
     * Method for drawing alignments without "blocks" (e.g. DotAlignedAlignment)
     */
    private void drawSimpleAlignment(Alignment alignment,
                                     Rectangle rect,
                                     Graphics2D g,
                                     RenderContext context,
                                     boolean flagUnmappedPair) {

        double origin = context.getOrigin();
        double locScale = context.getScale();
        int x = (int) ((alignment.getStart() - origin) / locScale);
        int length = alignment.getEnd() - alignment.getStart();
        int w = (int) Math.ceil(length / locScale);
        int h = (int) Math.max(1, rect.getHeight() - 2);
        int y = (int) (rect.getY() + (rect.getHeight() - h) / 2);
        int arrowLength = Math.min(5, w / 6);
        int d = Math.max(0, (int) (arrowLength + 2 - AlignmentPacker.MIN_ALIGNMENT_SPACING / context.getScale()));
        if (alignment.isNegativeStrand()) x += d;
        if (!alignment.isNegativeStrand()) w -= d;

        int[] xPoly = null;
        int[] yPoly = {y, y, y + h / 2, y + h, y + h};

        // Don't draw off edge of clipping rect
        if (x < rect.x && (x + w) > (rect.x + rect.width)) {
            x = rect.x;
            w = rect.width;
            arrowLength = 0;
        } else if (x < rect.x) {
            int delta = rect.x - x;
            x = rect.x;
            w -= delta;
            if (alignment.isNegativeStrand()) {
                arrowLength = 0;
            }
        } else if ((x + w) > (rect.x + rect.width)) {
            w -= ((x + w) - (rect.x + rect.width));
            if (!alignment.isNegativeStrand()) {
                arrowLength = 0;
            }
        }


        if (alignment.isNegativeStrand()) {
            xPoly = new int[]{x + w, x, x - arrowLength, x, x + w};
        } else {
            xPoly = new int[]{x, x + w, x + w + arrowLength, x + w, x};
        }
        g.fillPolygon(xPoly, yPoly, xPoly.length);

        if (flagUnmappedPair && alignment.isPaired() && !alignment.getMate().isMapped()) {
            g.setColor(Color.red);
            g.drawPolygon(xPoly, yPoly, xPoly.length);
        }
    }

    /**
     * Draw a pair of alignments as a single "template".
     */
    private void drawPairedAlignment(
            PairedAlignment pair,
            Rectangle rowRect,
            RenderContext context,
            AlignmentTrack.RenderOptions renderOptions,
            boolean leaveMargin,
            AlignmentCounts alignmentCounts,
            ThinReadBody thinReadBody, BaseRenderer.PaletteColors colors) {

        double locScale = context.getScale();

        Color alignmentColor1 = getAlignmentColor(pair.firstAlignment, track);
        Color alignmentColor2 = null;

        boolean overlapped = pair.secondAlignment != null && (pair.firstAlignment.getChr().equals(pair.secondAlignment.getChr())) &&
                pair.firstAlignment.getAlignmentEnd() > pair.secondAlignment.getAlignmentStart();

        Graphics2D g = context.getGraphics2D("ALIGNMENT");
        g.setColor(alignmentColor1);

        drawAlignment(pair.firstAlignment, rowRect, context, alignmentColor1, renderOptions, leaveMargin, alignmentCounts, overlapped, thinReadBody, colors);

        //If the paired alignment is in memory, we draw it.
        //However, we get the coordinates from the first alignment
        if (pair.secondAlignment != null) {
            if (alignmentColor2 == null) {
                alignmentColor2 = getAlignmentColor(pair.secondAlignment, track);
            }
            g.setColor(alignmentColor2);

            drawAlignment(pair.secondAlignment, rowRect, context, alignmentColor2, renderOptions, leaveMargin, alignmentCounts, overlapped, thinReadBody, colors);
        } else {
            return;
        }

        Color lineColor = track.getColor();
        if (alignmentColor1.equals(alignmentColor2) || pair.secondAlignment == null) {
            lineColor = alignmentColor1;
        }
        g.setColor(lineColor);

        double origin = context.getOrigin();
        int startX = (int) ((pair.firstAlignment.getEnd() - origin) / locScale);
        int endX = (int) ((pair.firstAlignment.getMate().getStart() - origin) / locScale);
        int h = (int) Math.max(1, rowRect.getHeight() - (leaveMargin ? 2 : 0));
        int y = (int) (rowRect.getY());
        startX = Math.max(rowRect.x, startX);
        endX = Math.min(rowRect.x + rowRect.width, endX);
        g.drawLine(startX, y + h / 2, endX, y + h / 2);
    }

    private static boolean overviewBatching(RenderContext context, AlignmentTrack.RenderOptions options) {
        return context.getScale() > 1 && context.getScale() < 100 &&
                options.getColorOption() == ColorOption.NONE && !options.isShowAllBases() &&
                !options.isQuickConsensusMode() && !options.isHideSmallIndels();
    }

    private static boolean rasterBaseStrip(Graphics2D graphics) {
        if (!graphics.getClass().getName().equals("sun.java2d.SunGraphics2D") ||
                !(graphics.getComposite() instanceof AlphaComposite composite) ||
                composite.getRule() != AlphaComposite.SRC_OVER ||
                graphics.getRenderingHint(RenderingHints.KEY_ANTIALIASING) == RenderingHints.VALUE_ANTIALIAS_ON) {
            return false;
        }
        AffineTransform transform = graphics.getTransform();
        // Integer device boundaries make image columns equivalent to fillRect.
        // Fractional scales/translations and vector/custom graphics keep their
        // existing per-base path (including per-block soft-clip batching).
        return transform.getShearX() == 0 && transform.getShearY() == 0 &&
                (transform.getScaleX() == 1 || transform.getScaleX() == 2) &&
                (transform.getScaleY() == 1 || transform.getScaleY() == 2) &&
                transform.getTranslateX() == Math.rint(transform.getTranslateX()) &&
                transform.getTranslateY() == Math.rint(transform.getTranslateY());
    }

    /**
     * Raster-only replacement for a normalized, square-capped one-pixel stroke.
     * Work in device pixels: a logical fillRect is not equivalent at 2x or 1.25x.
     * The cached graphics retains the original device clip. At overview zoom,
     * coverage counts retain repeated SRC_OVER contributions from overlapping caps.
     */
    private static final class ThinReadBody {
        private static final BasicStroke DEFAULT_STROKE = new BasicStroke();
        private static final double DEVICE_LIMIT = 1 << 20;
        private final Graphics2D graphics;
        private final double scaleX;
        private final double translateX;
        private final int leftInset;
        private final int top;
        private final int height;
        private final int startX;
        private final int width;
        private final int[] coverage;
        private final BaseRenderer.ColorStrip colorStrip;

        private ThinReadBody(Graphics2D graphics, AffineTransform transform, Rectangle row,
                             RenderContext context, boolean batch) {
            this.graphics = graphics;
            scaleX = transform.getScaleX();
            translateX = transform.getTranslateX();
            leftInset = scaleX == 2 ? -1 : 0;
            double scaleY = transform.getScaleY();
            int centerY = (int) Math.floor(row.y * scaleY + transform.getTranslateY() + 0.25);
            top = centerY + (scaleY == 2 ? -1 : 0);
            height = scaleY == 2 ? 2 : 1;
            // Bound scratch by the device damage clip, not the layout rectangle:
            // normalized stroke caps may paint outside that rectangle. With no
            // explicit clip, retain the existing fill path rather than guess the
            // raster's bounds.
            Rectangle clip = graphics.getClipBounds();
            startX = clip == null ? 0 : clip.x;
            width = clip == null ? 0 : clip.width;
            batch = batch && width > 0;
            coverage = batch ? context.getReadBodyCoverage(width) : null;
            colorStrip = batch ? context.getBaseColorStrip() : null;
        }

        private static ThinReadBody create(RenderContext context, Rectangle row, boolean batch) {
            Graphics2D source = context.getGraphics2D("ALIGNMENT");
            // Do not change SVG/vector semantics or assume a custom Graphics2D rasterizer.
            if (!source.getClass().getName().equals("sun.java2d.SunGraphics2D") ||
                    !DEFAULT_STROKE.equals(source.getStroke()) ||
                    !(source.getComposite() instanceof AlphaComposite composite) ||
                    composite.getRule() != AlphaComposite.SRC_OVER ||
                    source.getRenderingHint(RenderingHints.KEY_ANTIALIASING) == RenderingHints.VALUE_ANTIALIAS_ON ||
                    source.getRenderingHint(RenderingHints.KEY_STROKE_CONTROL) == RenderingHints.VALUE_STROKE_PURE) {
                return null;
            }
            AffineTransform transform = source.getTransform();
            if (transform.getShearX() != 0 || transform.getShearY() != 0 ||
                    !supportedScale(transform.getScaleX()) || !supportedScale(transform.getScaleY()) ||
                    !quarterPixel(transform.getTranslateX()) || !quarterPixel(transform.getTranslateY()) ||
                    Math.abs(row.y * transform.getScaleY() + transform.getTranslateY()) > DEVICE_LIMIT ||
                    Math.abs(row.x * transform.getScaleX() + transform.getTranslateX()) > DEVICE_LIMIT ||
                    Math.abs(row.getMaxX() * transform.getScaleX() + transform.getTranslateX()) > DEVICE_LIMIT) {
                return null;
            }
            Graphics2D graphics = context.getGraphics2D("THIN_READ_BODY");
            graphics.setTransform(transform);
            graphics.setClip(source.getClip());
            graphics.setTransform(new AffineTransform());
            graphics.setComposite(source.getComposite());
            graphics.setRenderingHint(RenderingHints.KEY_ANTIALIASING, RenderingHints.VALUE_ANTIALIAS_OFF);
            return new ThinReadBody(graphics, transform, row, context, batch);
        }

        private static boolean supportedScale(double scale) {
            return scale == 1 || scale == 1.25 || scale == 2;
        }

        private static boolean quarterPixel(double coordinate) {
            return Math.abs(coordinate) <= DEVICE_LIMIT && coordinate * 4 == Math.rint(coordinate * 4);
        }

        private void reset(Color color) {
            graphics.setColor(color);
            if (coverage != null) Arrays.fill(coverage, 0, width, 0);
        }

        private void draw() {
            if (coverage == null) return;
            colorStrip.reset(graphics, startX, width);
            Color color = graphics.getColor();
            for (int i = 0; i < width; i++) {
                for (int layer = 0; layer < coverage[i]; layer++) {
                    colorStrip.blendColor(startX + i, color);
                }
            }
            colorStrip.draw(top, height);
        }

        private boolean fill(int start, int end) {
            double first = Math.min(start, end) * scaleX + translateX;
            double last = Math.max(start, end) * scaleX + translateX;
            if (Math.abs(first) > DEVICE_LIMIT || Math.abs(last) > DEVICE_LIMIT) return false;
            // Java2D normalizes endpoints to floor(deviceCoordinate + .25) + .25.
            // Pixel-center coverage of the square caps includes both endpoints,
            // even when they coincide or arrive in reverse order.
            int left = (int) Math.floor(first + 0.25) + leftInset;
            int right = (int) Math.floor(last + 0.25) + 1;
            if (coverage == null) {
                graphics.fillRect(left, top, right - left, height);
            } else {
                int firstColumn = Math.max(0, left - startX);
                int lastColumn = Math.min(width, right - startX);
                for (int i = firstColumn; i < lastColumn; i++) coverage[i]++;
            }
            return true;
        }
    }

    // Keep expanded scratch state out of the shared thin-body/base traversal.
    private void prepareExpandedReadBody(RenderContext context, Rectangle row, int height,
                                         Graphics2D graphics, Graphics2D outline,
                                         AlignmentTrack.RenderOptions options) {
        boolean eligible = track.getDisplayMode() == Track.DisplayMode.EXPANDED &&
                outline == null && overviewBatching(context, options);
        ExpandedReadBody body = eligible ? context.getExpandedReadBody() : context.getExpandedReadBodyIfPresent();
        if (body != null) body.reset(eligible ? graphics : null, row, height);
    }

    private static void flushExpandedReadBody(RenderContext context) {
        ExpandedReadBody body = context.getExpandedReadBodyIfPresent();
        if (body != null) body.draw();
    }

    /**
     * Draw a (possibly gapped) alignment
     * <p>
     * NOTE: This is a large method, but every attempt to break it up results in methods with very long argument lists.
     */
    private void drawAlignment(
            Alignment alignment,
            Rectangle rowRect,
            RenderContext context,
            Color alignmentColor,
            AlignmentTrack.RenderOptions renderOptions,
            boolean leaveMargin,
            AlignmentCounts alignmentCounts,
            boolean overlapped,
            ThinReadBody thinReadBody, BaseRenderer.PaletteColors colors) {

        AlignmentBlock[] blocks = alignment.getAlignmentBlocks();

        Graphics2D gAlignment = context.getGraphics2D("ALIGNMENT");
        gAlignment.setColor(alignmentColor);

        // No blocks.  Note: SAM/BAM alignments always have at least 1 block
        if (blocks == null || blocks.length == 0) {
            drawSimpleAlignment(alignment, rowRect, gAlignment, context, renderOptions.isFlagUnmappedPairs());
            return;
        }

        IGVPreferences prefs = this.track.getPreferences();
        boolean flagLargeIndels = prefs.getAsBoolean(SAM_FLAG_LARGE_INDELS);
        int largeInsertionsThreshold = prefs.getAsInt(SAM_LARGE_INDELS_THRESHOLD);
        boolean hideSmallIndelsBP = renderOptions.isHideSmallIndels();
        int indelThresholdBP = renderOptions.getSmallIndelThreshold();
        boolean quickConsensus = renderOptions.isQuickConsensusMode();
        final float snpThreshold = prefs.getAsFloat(SAM_ALLELE_THRESHOLD);

        // Scale and position of the alignment rendering.
        double locScale = context.getScale();
        int h = (int) Math.max(1, rowRect.getHeight() - (leaveMargin ? 2 : 0));
        int y = (int) (rowRect.getY());


        // Get a graphics context for drawing clipping indicators.
        Graphics2D clippedGraphics = context.getGraphic2DForColor(clippedColor);
        if (h > 5) {
            clippedGraphics.setStroke(new BasicStroke(1.2f));
        }

        // Get a graphics context for drawing strand indicators.
        Graphics2D strandGraphics = context.getGraphics2D("STRAND");

        // Define a graphics context for indel labels.
        Graphics2D largeIndelGraphics = context.getGraphics2D("INDEL_LABEL");
        largeIndelGraphics.setFont(FontManager.getFont(Font.BOLD, h - 2));

        // Get a graphics context for drawing individual basepairs.
        Graphics2D bpGraphics = context.getGraphics2D("BASE");
        int dX = (int) Math.max(1, (1.0 / locScale));
        if (dX >= 8) {
            Font f = FontManager.getFont(Font.BOLD, Math.min(dX, h));
            bpGraphics.setFont(f);
        }

        /* Clipping */
        boolean flagClipping = prefs.getAsBoolean(SAM_FLAG_CLIPPING);
        int clippingThreshold = prefs.getAsInt(SAM_CLIPPING_THRESHOLD);
        ClippingCounts clipping = alignment.getClippingCounts();
        boolean leftClipped = flagClipping && (clipping.getLeft() > clippingThreshold);
        boolean rightClipped = flagClipping && (clipping.getRight() > clippingThreshold);

        double bpStart = context.getOrigin();
        double bpEnd = Math.ceil(context.getEndLocation());

        // Draw gaps (deletions and junctions)
        List<Gap> gaps = alignment.getGaps();
        if (gaps != null) {
            for (Gap gap : gaps) {
                int gapStart = gap.getStart();
                int gapWidth = gap.getnBases();
                int gapEnd = gapStart + gapWidth;

                int gapPxStart = (int) ((Math.max(bpStart, gapStart) - bpStart) / locScale);
                int gapPxEnd = (int) ((Math.min(bpEnd, gapEnd) - bpStart) / locScale);

                if (gapEnd <= bpStart) { // gap ends before the visible context
                    continue; // move to next gap
                } else if (gapStart >= bpEnd) { // gap starts after the visible context
                    break; // done examining gaps
                }

                // Draw the gap if it is sufficiently large
                boolean drawGap = (!hideSmallIndelsBP || gapWidth >= indelThresholdBP);

                //Check deletion consistency  -- TODO this needs work, disabled for now
                //if (drawGap && quickConsensus) {
                //    drawGap = alignmentCounts.isConsensusDeletion(gapChromStart, gap.getnBases(), snpThreshold);
                //}
                if (!drawGap) {
                    continue;
                }

                // Draw the gap line.
                Graphics2D gapGraphics = context.getGraphics2D("GAP");
                if (gap.getType() == SAMAlignment.UNKNOWN) {
                    gapGraphics.setColor(unknownGapColor);
                } else if (gap.getType() == SAMAlignment.SKIPPED_REGION) {
                    gapGraphics.setColor(skippedColor);
                } else {
                    gapGraphics.setColor(deletionColor);
                    if (h > 5) {
                        gapGraphics = context.getGraphics2D("THICK_STROKE");
                    }
                }

                gapGraphics.drawLine(gapPxStart, y + h / 2, gapPxEnd, y + h / 2);

                // Label the size of the deletion if it is "large" and the label fits.
                if (flagLargeIndels && (gap.getType() == SAMAlignment.DELETION) && gapWidth > largeInsertionsThreshold) {
                    drawLargeIndelLabel(largeIndelGraphics,
                            false,
                            Globals.DECIMAL_FORMAT.format(gapWidth),
                            ((gapPxStart + gapPxEnd) / 2),
                            y,
                            h,
                            gapPxEnd - gapPxStart - 2,
                            context.translateX,
                            null,
                            alignment,
                            context);
                }

                // gap extensions
                if (flowIndelRendering.handlesAlignment(alignment) && flowIndelRendering.handlesGap(gap)) {
                    flowIndelRendering.renderDeletionGap(alignment, gap, y, h, gapPxStart, gapPxEnd - gapPxStart, context, renderOptions);
                }
            }
        }

        // Draw blocks
        // Get a graphics context for outlining alignment blocks.
        Graphics2D outlineGraphics = null;
        final HashMap<String, Color> selectedReadNames = this.track.getSelectedReadNames();
        final String readName = alignment.getReadName();
        if (selectedReadNames.containsKey(readName)) {
            Color c = selectedReadNames.get(readName);
            c = (c == null) ? Color.blue : c;
            outlineGraphics = context.getGraphics2D("THICK_STROKE");
            gAlignment.setColor(c);
        } else if (renderOptions.isFlagUnmappedPairs() && alignment.isPaired() && !alignment.getMate().isMapped()) {
            outlineGraphics = context.getGraphics2D("OUTLINE");
            outlineGraphics.setColor(Color.red);
        } else if (alignment.getMappingQuality() == 0 && renderOptions.isFlagZeroQualityAlignments()) {
            outlineGraphics = context.getGraphic2DForColor(OUTLINE_COLOR);
        }
        if (thinReadBody != null) thinReadBody.reset(gAlignment.getColor());
        if (h > 1) prepareExpandedReadBody(context, rowRect, h, gAlignment, outlineGraphics, renderOptions);


        // Compute arrow width from total length of alignment on reference
        double pixelLengthOnReference = alignment.getLengthOnReference() / locScale;
        int arrowPxWidth = pixelLengthOnReference == 0 ? 0 : (int) Math.min(Math.min(5, h / 2), pixelLengthOnReference / 6);

        int blockChromStart = blocks[0].getStart();
        boolean leftmost = true;
        for (int blockIx = 0; blockIx < blocks.length; blockIx++) {
            AlignmentBlock block = blocks[blockIx];

            int blockChromEnd = block.getStart() + block.getLength();
            int blockPxStart = (int) Math.round((blockChromStart - bpStart) / locScale);
            int blockPxEnd = (int) Math.round((blockChromEnd - bpStart) / locScale);

            // Check if block is in visible rectangle
            if (blockPxEnd < 0) {
                if (blockIx + 1 < blocks.length) {
                    blockChromStart = blocks[blockIx + 1].getStart(); // start position for the next block
                }
                continue;
            } else if (blockPxStart > rowRect.x + rowRect.width) {
                break;
            }

            boolean rightmost = blockIx + 1 == blocks.length;

            if (!rightmost) { // consider waiting to draw the block unless it is rightmost
                if (hideSmallIndelsBP && (blocks[blockIx + 1].getStart() - blockChromEnd) < indelThresholdBP) {
                    continue; // small indel between this block and the next; wait to draw
                } else {
                    blockChromStart = blocks[blockIx + 1].getStart(); // start position for the next block
                }
            }

            if (h == 1) {
                if (thinReadBody == null || !thinReadBody.fill(blockPxStart, blockPxEnd)) {
                    // Flush earlier contributions before an unsupported endpoint's
                    // original line paint, then resume with empty coverage.
                    if (thinReadBody != null) {
                        thinReadBody.draw();
                        thinReadBody.reset(gAlignment.getColor());
                    }
                    gAlignment.drawLine(blockPxStart, y, blockPxEnd, y);
                }
            } else {
                if (!overlapped) {
                    int pixelGap = (int) (AlignmentPacker.MIN_ALIGNMENT_SPACING / locScale);
                    if (pixelGap < arrowPxWidth) {
                        arrowPxWidth = Math.max(0, arrowPxWidth - pixelGap);
                    }
                }

                boolean drawArrow = h > 6 && (leftmost && blockPxStart > 0 || rightmost && blockPxEnd < rowRect.x + rowRect.width);
                blockPxStart = Math.max(0, blockPxStart);
                blockPxEnd = Math.min(rowRect.x + rowRect.width, blockPxEnd);

                boolean drawLeftClip = leftmost && leftClipped;
                boolean drawRightClip = rightmost && rightClipped;
                boolean textured = renderOptions.getDuplicatesOption() == AlignmentTrack.DuplicatesOption.TEXTURE && alignment.isDuplicate();
                int leftArrowWidth = leftmost && alignment.isNegativeStrand() && drawArrow ? arrowPxWidth : 0;
                int rightArrowWidth = rightmost && !alignment.isNegativeStrand() && drawArrow ? arrowPxWidth : 0;

                Graphics2D g = gAlignment;
                if (!block.hasBases()) {
                    if (block.isSoftClip()) g = context.getGraphics2D("SOFT_CLIP");
                    else if (block.getCigarOperator() == 'X') g = context.getGraphics2D("MISMATCH");
                }

                // A terminal block can still be rectangular: its strand's tip may be
                // clipped, on another block, or rounded down to zero pixels.
                if (leftArrowWidth == 0 && rightArrowWidth == 0 && outlineGraphics == null &&
                        !textured && !drawLeftClip && !drawRightClip && blockPxEnd >= blockPxStart) {
                    ExpandedReadBody expandedReadBody = context.getExpandedReadBodyIfPresent();
                    if (expandedReadBody == null || g != gAlignment || !expandedReadBody.fill(blockPxStart, blockPxEnd)) {
                        flushExpandedReadBody(context);
                        g.fillRect(blockPxStart, y, blockPxEnd - blockPxStart, h);
                    }
                } else {
                    // Preserve traversal order when a decoration interrupts rectangles.
                    flushExpandedReadBody(context);
                    // Keep the original geometry for arrows, outlines, textures and clipping decorations.
                    int[] xPoly = {blockPxStart - leftArrowWidth, blockPxStart, blockPxEnd,
                            blockPxEnd + rightArrowWidth, blockPxEnd, blockPxStart};
                    int[] yPoly = {y + h / 2, y, y, y + h / 2, y + h, y + h};
                    Shape blockShape = new Polygon(xPoly, yPoly, xPoly.length);
                    if (textured) {
                        final Graphics2D tg = (Graphics2D) g.create();
                        tg.setPaint(getDuplicateTexture(tg.getColor()));
                        tg.fill(blockShape);
                        tg.dispose();
                    } else {
                        g.fill(blockShape);
                    }

                    if (outlineGraphics != null) {
                        outlineGraphics.draw(blockShape);
                    }

                    if (drawLeftClip || drawRightClip) {
                        final SupplementaryAlignment.SupplementaryNeighbors supplementaryRenderingInfo =
                                SupplementaryAlignment.getAdjacentSupplementaryReads(alignment);
                        drawClippedEnds(clippedGraphics, xPoly, yPoly, drawLeftClip, drawRightClip, supplementaryRenderingInfo);
                    }
                }

            }
            leftmost = false;
        }
        if (thinReadBody != null) thinReadBody.draw();
        if (h > 1) flushExpandedReadBody(context);

        // Draw bases for an alignment block.  The bases are "overlaid" on the block with a transparency value (alpha)
        // that is proportional to the base quality score, or flow signal deviation, whichever is selected.

        final ColorOption colorOption = renderOptions.getColorOption();
        if (locScale < 100) {

            boolean showAllBases = renderOptions.isShowAllBases() &&
                    !(colorOption == ColorOption.BISULFITE || colorOption == ColorOption.NOMESEQ); // Disable showAllBases in bisulfite mode

            if (renderOptions.isShowMismatches() || showAllBases) {
                BaseRenderer.ColorStrip alignmentStrip = null;
                if (overviewBatching(context, renderOptions) && rasterBaseStrip(gAlignment)) {
                    Graphics2D stripGraphics = context.getGraphics2D("BASE_COLOR_STRIP");
                    stripGraphics.setTransform(gAlignment.getTransform());
                    stripGraphics.setClip(gAlignment.getClip());
                    stripGraphics.setComposite(gAlignment.getComposite());
                    alignmentStrip = context.getBaseColorStrip();
                    // Preserve partial edge columns and the existing inclusive
                    // right-edge test, independently of the graphics damage clip.
                    alignmentStrip.reset(stripGraphics, rowRect.x - 1, rowRect.width + 2);
                }


                for (AlignmentBlock block : alignment.getAlignmentBlocks()) {

                    boolean haveBases = (block.hasBases() && block.getLength() > 0);
                    if (!haveBases) {
                        continue;   // nothing to draw
                    }

                    String chr = context.getChr();
                    final int start = block.getStart();
                    final int end = block.getEnd();
//                    if (end <= bpStart) { // block ends before the visible context
//                        continue; // move to next block
//                    } else if (start >= bpEnd) { // block starts after the visible context
//                        break; // done examining blocks
//                    }

                    boolean isSoftClip = block.isSoftClip();

                    // Get the reference sequence
                    Genome genome = GenomeManager.getInstance().getCurrentGenome();
                    final byte[] reference = isSoftClip ? softClippedReference : genome.getSequence(chr, start, end);

                    if (reference == null && !showAllBases) {
                        continue;   // No reference to compare to.
                    }

                    // Compute bounds
                    int pY = (int) rowRect.getY();
                    int dY = (int) rowRect.getHeight();
                    dX = (int) Math.max(1, (1.0 / locScale));

                    BisulfiteBaseInfo bisinfo = null;
                    boolean nomeseqMode = (colorOption.equals(AlignmentTrack.ColorOption.NOMESEQ));
                    boolean bisulfiteMode = AlignmentTrack.isBisulfiteColorType(colorOption);
                    if (nomeseqMode) {
                        bisinfo = new BisulfiteBaseInfoNOMeseq(reference, alignment, block, renderOptions.bisulfiteContext);
                    } else if (bisulfiteMode) {
                        bisinfo = new BisulfiteBaseInfo(reference, alignment, block, renderOptions.bisulfiteContext);
                    }

                    ByteSubarray blockBases = block.getBases();
                    final int s = (int) Math.max(Math.floor(bpStart), start);
                    final int e = (int) Math.min(Math.ceil(bpEnd), end);
                    final boolean batchSoftClips = isSoftClip && locScale > 1 && !bisulfiteMode &&
                            gAlignment.getComposite() instanceof AlphaComposite composite &&
                            composite.getRule() == AlphaComposite.SRC_OVER;
                    BaseRenderer.ColorStrip colorStrip = alignmentStrip;
                    if (colorStrip == null && batchSoftClips && s < e) {
                        int firstPixel = Math.max(rowRect.x - 1, (int) ((s - bpStart) / locScale));
                        int lastPixel = Math.min((int) rowRect.getMaxX(), (int) ((e - 1.0 - bpStart) / locScale));
                        if (firstPixel > lastPixel) continue;
                        Graphics2D stripGraphics = context.getGraphics2D("BASE_COLOR_STRIP");
                        stripGraphics.setComposite(gAlignment.getComposite());
                        colorStrip = context.getBaseColorStrip();
                        colorStrip.reset(stripGraphics, firstPixel, lastPixel - firstPixel + 1);
                    }
                    for (int loc = s; loc < e; loc++) {

                        int idx = loc - start;

                        boolean misMatch = AlignmentUtils.isMisMatch(reference, blockBases, isSoftClip, idx);

                        // if base is modified paint rectangle
                        if (showAllBases || (!bisulfiteMode && misMatch) ||
                                (bisulfiteMode && (DisplayStatus.NOTHING != bisinfo.getDisplayStatus(idx)))) {
                            char c = (char) blockBases.getByte(idx);

                            double bisulfiteXaxisShift = (bisulfiteMode) ? bisinfo.getXaxisShift(idx) : 0;
                            int pX = (int) (((double) loc + bisulfiteXaxisShift - bpStart) / locScale);

                            // Don't draw out of clipping rect
                            if (pX > rowRect.getMaxX()) {
                                break;
                            } else if (pX + dX < rowRect.getX()) {
                                continue;
                            }
                            Color color = null;
                            if (bisulfiteMode) {
                                color = bisinfo.getDisplayColor(idx);
                            } else if (colorOption.isBaseMod() ||
                                    colorOption.isSMRTKinetics()) {
                                color = Color.GRAY;
                            } else {
                                byte quality = colors.shade ? block.getQuality(idx) : 0;
                                color = colors.getColor((byte) c, quality);
                            }
                            if (color == null) {
                                color = Color.black;
                            }

                            if (colors.shade && (bisulfiteMode || colorOption.isBaseMod() || colorOption.isSMRTKinetics())) {
                                byte qual = block.getQuality(idx);
                                color = BaseRenderer.getShadedColor(color, qual, colors.minQ, colors.maxQ);
                            }

                            BisulfiteBaseInfo.DisplayStatus bisstatus = (bisinfo == null) ? null : bisinfo.getDisplayStatus(idx);

                            final boolean showBase = showAllBases ||
                                    isSoftClip ||
                                    bisulfiteMode ||
                                    // In "quick consensus" mode, only show mismatches at positions with a consistent alternative basepair.
                                    (!quickConsensus || alignmentCounts.isConsensusMismatch(loc, reference[idx], chr, snpThreshold));
                            if (showBase) {
                                if (colorStrip != null) {
                                    colorStrip.blendColor(pX, color);
                                } else {
                                    BaseRenderer.drawBase(gAlignment, color, c, pX, pY, dX, dY - (leaveMargin ? 2 : 0), bisulfiteMode, bisstatus);
                                }
                            }
                        }
                    }
                    if (colorStrip != null && alignmentStrip == null) {
                        colorStrip.draw(pY, dY - (leaveMargin ? 2 : 0));
                    }
                }
                if (alignmentStrip != null) {
                    alignmentStrip.draw(rowRect.y, h);
                }
            }
        }

        // Base modification
        if (colorOption.isBaseMod()) {
            BaseModificationRenderer.drawModifications(alignment, bpStart, locScale, rowRect, context.getGraphics(), renderOptions);
        }

        // Kinetic data
        if (colorOption.isSMRTKinetics()) {
            SMRTKineticsRenderer.drawSmrtKinetics(alignment, bpStart, locScale, rowRect, context.getGraphics(), colorOption);
        }

        // DRAW Insertions
        AlignmentBlock[] insertions = alignment.getInsertions();
        if (insertions != null) {


            for (AlignmentBlock aBlock : insertions) {

                if (aBlock.getStart() == context.expandedInsertionPosition) continue;   // Skip, will be drawn expanded

                int x = (int) ((aBlock.getStart() - bpStart) / locScale);
                int bpWidth = aBlock.getBasesLength();
                double pxWidthExact = ((double) bpWidth) / locScale;
                h = (int) Math.max(1, rowRect.getHeight() - (leaveMargin ? 2 : 0));
                y = (int) (rowRect.getY() + (rowRect.getHeight() - h) / 2) - (leaveMargin ? 1 : 0);

                // Don't draw out of clipping rect
                if (x > rowRect.getMaxX()) {
                    break;
                } else if (x < rowRect.getX()) {
                    continue;
                }

                if ((!hideSmallIndelsBP || bpWidth >= indelThresholdBP)) {
                    // && (!quickConsensus || alignmentCounts.isConsensusInsertion(aBlock.getStart(), snpThreshold))) {
                    if (flagLargeIndels && bpWidth > largeInsertionsThreshold) {
                        drawLargeIndelLabel(context.getGraphics2D("INDEL_LABEL"),
                                true,
                                Globals.DECIMAL_FORMAT.format(bpWidth),
                                x - 1,
                                y,
                                h,
                                (int) pxWidthExact,
                                context.translateX,
                                aBlock,
                                alignment,
                                context);
                    } else {
                        int pxWing = (h > 10 ? 2 : (h > 5) ? 1 : 0);
                        Graphics2D ig = context.getGraphics();
                        Color color = purple;
                        boolean sbxFullDiscordantInsertionRecolor = renderOptions.isIndelQualSbx() &&
                                SbxUtils.isFullDiscordantInsertion(aBlock);
                        if (sbxFullDiscordantInsertionRecolor) {
                            color = SbxUtils.fullDiscordantInsertionColor;
                        } else {
                            boolean sbxDiscordantInsertionRecolor = renderOptions.isIndelQualSbx() &&
                                    SbxUtils.isDiscordantInsertion(alignment, aBlock);
                            color = sbxDiscordantInsertionRecolor ? SbxUtils.discordantInsertionColor : purple;
                        }
                        ig.setColor(color);
                        if (flowIndelRendering.handlesAlignment(alignment) && flowIndelRendering.handlesBlock(aBlock)) {
                            flowIndelRendering.renderSmallInsertion(alignment, aBlock, context, h, x, y, renderOptions);
                        } else {
                            ig.fillRect(x, y, 2, h);
                            ig.fillRect(x - pxWing, y, 2 + 2 * pxWing, 2);
                            ig.fillRect(x - pxWing, y + h - 2, 2 + 2 * pxWing, 2);
                        }

                        aBlock.setPixelRange(context.translateX + x - pxWing, context.translateX + x + 2 + pxWing);
                    }
                }
            }
        }

        // draw simplex tails
        if (renderOptions.isTailQualSbx() && !renderOptions.isHideTailSbx()) {
            int simplexTailLeftEnd = 0;
            int alignmentLength = alignment.getEnd() - alignment.getStart();
            for (; simplexTailLeftEnd < alignmentLength && (alignment.getPhred(alignment.getStart() + simplexTailLeftEnd) <= SAMAlignment.SBX_LOW_BASEQ_TAIL_THRESHOLD || alignment.getBase(alignment.getStart() + simplexTailLeftEnd) == 0); simplexTailLeftEnd++) {

            }
            int simplexTailRightStart = alignmentLength - 1;
            for (; simplexTailRightStart >= 0 && (alignment.getPhred(alignment.getStart() + simplexTailRightStart) <= SAMAlignment.SBX_LOW_BASEQ_TAIL_THRESHOLD || alignment.getBase(alignment.getStart() + simplexTailRightStart) == 0); simplexTailRightStart--)
                ;
            int[] tailStarts = new int[]{};
            int[] tailEnds = new int[]{};
            if (simplexTailLeftEnd > 0 && simplexTailRightStart < alignmentLength - 1 && simplexTailRightStart > 0) {
                tailStarts = new int[]{alignment.getStart(), alignment.getStart() + simplexTailRightStart + 1};
                tailEnds = new int[]{simplexTailLeftEnd, alignmentLength};
            } else if (simplexTailLeftEnd != 0) {
                tailStarts = new int[]{alignment.getStart()};
                tailEnds = new int[]{simplexTailLeftEnd};
            } else if (simplexTailRightStart < alignmentLength - 1) {
                tailStarts = new int[]{alignment.getStart() + simplexTailRightStart + 1};
                tailEnds = new int[]{alignmentLength};
            }
            for (int i = 0; i < tailStarts.length; i++) {
                int blockChromEnd1 = tailStarts[i];
                int blockChromEnd2 = alignment.getStart() + tailEnds[i];
                int blockPxStart = (int) Math.round((blockChromEnd1 - bpStart) / locScale);
                int blockPxEnd = (int) Math.round((blockChromEnd2 - bpStart) / locScale);
                gAlignment.setColor(new Color(0, 0, 0, 15));
                gAlignment.fillRect(blockPxStart, y, blockPxEnd - blockPxStart, h);

                Color tailColor = new Color(150, 50, 200, 20);

                Graphics2D ig = context.getGraphics();
                ig.setColor(tailColor);
                ig.fillRect(blockPxStart, y, blockPxEnd - blockPxStart, h);

                Graphics2D foo = context.getGraphics2D("INDEL_LABEL");
                foo.setColor(tailColor);
                foo.fillRect(blockPxStart, y, blockPxEnd - blockPxStart, h);
            }
        }

    }

    /**
     * get a texture to apply to duplicate reads, caches the created textures according to their color
     * @param baseColor the color to render the read
     * @return a texture matching the base color with shading
     */
    private TexturePaint getDuplicateTexture(Color baseColor) {
        return duplicateTextures.computeIfAbsent(baseColor, color -> {
            BufferedImage texture = new BufferedImage(5, 5, BufferedImage.TYPE_INT_RGB);
            Graphics2D big = texture.createGraphics();
            //Render into the BufferedImage graphics to create the texture
            big.setColor(baseColor);
            big.fillRect(0, 0, 5, 5);

            final Color darker = baseColor.darker();
            big.setColor(darker);
            // hash pattern, probably better to save it as a bitmap somewhere
            big.drawLine(0, 2, 0, 3);
            big.drawLine(1, 3, 1, 4);
            big.drawLine(2, 4, 2, 4);
            big.drawLine(2, 0, 2, 0);
            big.drawLine(3, 0, 3, 1);
            big.drawLine(4, 1, 4, 2);
            big.dispose();
            // Create a texture paint from the buffered image
            Rectangle r = new Rectangle(0, 0, 5, 5);
            return new TexturePaint(texture, r);
        });
    }

    private static void drawClippedEnds(final Graphics2D g, final int[] xPoly, final int[] yPoly,
                                        final boolean drawLeftClip, final boolean drawRightClip,
                                        final SupplementaryAlignment.SupplementaryNeighbors sri) {
        /*
                5       4
             0 <|=======|> 3
                1       2
         */
        final int xLeftPoint = xPoly[0];
        final int xLeft = xPoly[1];
        final int xRightPoint = xPoly[3];
        final int xRight = xPoly[2];
        final int yMiddle = yPoly[0];
        final int yBottom = yPoly[5];
        final int yTop = yPoly[1];
        final Color savedColor = g.getColor();
        final int arrowWidth = 3;
        String leftContigLabel = null;
        String rightContigLabel = null;
        try {
            //left side
            if (drawLeftClip) {
                if (sri != null && sri.previous() != null) {
                    final SupplementaryAlignment previous = sri.previous();
                    if (previous.contigsMatch(sri.alignment())) {
                        g.setColor(previous.getStrand() == sri.alignment().getReadStrand() ? NON_INVERSION_COLOR : INVERSION_COLOR);
                    } else {
                        if (previous.getContig() != null) {
                            g.setColor(ChromosomeColors.getColor(previous.getContig()));
                            leftContigLabel = previous.getContig();
                        } else {
                            log.warn("previous is missing contig: " + previous);
                        }
                    }
                    final Polygon thickLeftArrow = new Polygon(
                            new int[]{xLeftPoint, xLeft, xLeft + arrowWidth, xLeft + arrowWidth, xLeft},
                            new int[]{yMiddle, yBottom, yBottom, yTop, yTop},
                            5);
                    g.drawPolygon(thickLeftArrow);
                    g.fillPolygon(thickLeftArrow);
                    g.setColor(savedColor);
                }
                g.drawLine(xLeftPoint, yMiddle, xLeft, yBottom);
                g.drawLine(xLeft, yTop - 1, xLeftPoint, yMiddle);
            }
            //right side
            if (drawRightClip) {
                if (sri != null && sri.next() != null) {
                    final SupplementaryAlignment next = sri.next();
                    if (next.contigsMatch(sri.alignment())) {
                        //TODO Set to scaled color by contig position
                        g.setColor(next.getStrand() == sri.alignment().getReadStrand() ? NON_INVERSION_COLOR : INVERSION_COLOR);
                    } else {
                        g.setColor(ChromosomeColors.getColor(next.getContig()));
                        rightContigLabel = next.getContig();
                    }
                    Polygon thickRightArrow = new Polygon(new int[]{xRightPoint, xRight, xRight - arrowWidth, xRight - arrowWidth, xRight},
                            new int[]{yMiddle, yBottom, yBottom, yTop, yTop}, 5);
                    g.drawPolygon(thickRightArrow);
                    g.fillPolygon(thickRightArrow);
                    g.setColor(savedColor);
                }
                g.drawLine(xRight, yBottom, xRightPoint, yMiddle);
                g.drawLine(xRightPoint, yMiddle, xRight, yTop - 1);
            }
            //contig names
            if (PreferencesManager.getPreferences().getAsBoolean(SAM_SHOW_CONNECTED_CHR_NAME)) {
                if (leftContigLabel != null || rightContigLabel != null) {
                    final int height = yBottom - yTop;
                    int textHeight = g.getFontMetrics().getHeight();
                    if (textHeight <= height + 4) {
                        int leftTextWidth = leftContigLabel != null ? g.getFontMetrics().stringWidth(leftContigLabel) : 0;
                        int rightContigWidth = rightContigLabel != null ? g.getFontMetrics().stringWidth(rightContigLabel) : 0;

                        final int width = xRight - xLeft;
                        int totalTextWidth = leftTextWidth + rightContigWidth + g.getFontMetrics().stringWidth(" ");
                        final int distanceFromArrow = 4;
                        if (totalTextWidth <= width + 2 * distanceFromArrow) {
                            if (leftContigLabel != null) {
                                g.setColor(ChromosomeColors.getColor(leftContigLabel));
                                g.drawString(leftContigLabel, xLeft + distanceFromArrow, yBottom - 1);
                            }
                            if (rightContigLabel != null) {
                                g.setColor(ChromosomeColors.getColor(rightContigLabel));
                                g.drawString(rightContigLabel, xRight - (rightContigWidth + distanceFromArrow), yBottom - 1);
                            }
                        }
                    }
                }
            }
        } finally {
            g.setColor(savedColor);
        }
    }

    private void drawLargeIndelLabel(Graphics2D g, boolean isInsertion, String labelText, int pxCenter,
                                     int pxTop, int pxH, int pxWmax, int translateX, AlignmentBlock insertionBlock, Alignment alignment, RenderContext context) {

        final int pxPad = 2;   // text padding in the label
        final int pxWing = (pxH > 10 ? 2 : 1);  // width of the cursor "wing"
        final int minTextHeight = 8; // min height to draw text

        // Calculate the width required to draw the label
        Rectangle2D textBounds = g.getFontMetrics().getStringBounds(labelText, g);
        int pxTextW = 2 * pxPad + (int) textBounds.getWidth();
        boolean doesTextFit = (pxH >= minTextHeight) && (pxTextW < pxWmax);

        if (!doesTextFit && !isInsertion) {
            return;
        } // only label deletions when the text fits

        // Calculate the pixel bounds of the label
        int pxW = (int) Math.max(2, Math.min(pxTextW, pxWmax)),
                pxLeft = pxCenter - (int) Math.ceil(pxW / 2),
                pxRight = pxLeft + pxW;

        // Draw the label
        Color alternateColor = purple;
        boolean sbxFullDiscordantInsertionRecolor = track.renderOptions.isIndelQualSbx() &&
                SbxUtils.isFullDiscordantInsertion(insertionBlock);
        if (sbxFullDiscordantInsertionRecolor) {
            alternateColor = SbxUtils.fullDiscordantInsertionColor;
        } else {
            boolean sbxDiscordantInsertionRecolor = track.renderOptions.isIndelQualSbx() &&
                    SbxUtils.isDiscordantInsertion(alignment, insertionBlock);
            alternateColor = sbxDiscordantInsertionRecolor ? SbxUtils.discordantInsertionColor : purple;
        }

        g.setColor(isInsertion ? alternateColor : Color.white);
        g.fillRect(pxLeft, pxTop, pxRight - pxLeft, pxH);

        if (isInsertion && pxH > 5) {
            if (flowIndelRendering.handlesAlignment(alignment) && flowIndelRendering.handlesBlock(insertionBlock)) {
                flowIndelRendering.renderSmallInsertionWings(alignment, insertionBlock, context, pxH, pxTop, pxRight, pxLeft, track.renderOptions);
            } else {
                g.fillRect(pxLeft - pxWing, pxTop, pxRight - pxLeft + 2 * pxWing, 2);
                g.fillRect(pxLeft - pxWing, pxTop + pxH - 2, pxRight - pxLeft + 2 * pxWing, 2);
            }
        } // draw "wings" For insertions

        if (doesTextFit) {
            g.setColor(isInsertion ? Color.white : alternateColor);
            GraphicUtils.drawCenteredText(labelText, pxLeft, pxTop, pxW, pxH, g);
        } // draw the text if it fits

        if (insertionBlock != null) {
            insertionBlock.setPixelRange(context.translateX + pxLeft, context.translateX + pxRight);
        }
    }

    private Color getAlignmentColor(Alignment alignment, AlignmentTrack track) {

        // Set color used to draw the feature.  Highlight features that intersect the
        // center line.  Also restorePersistentState row "score" if alignment intersects center line

        Color defaultColor = track.getColor();
        Color c = defaultColor;
        AlignmentTrack.RenderOptions renderOptions = track.getRenderOptions();
        ColorOption colorOption = renderOptions.getColorOption();
        String readNameParts[];

        switch (colorOption) {

            case YC_TAG:

                Color ycColor = alignment.getYcColor();
                if (ycColor != null) {
                    c = ycColor;
                }
                break;

            case BISULFITE:
            case BASE_MODIFICATION:
            case BASE_MODIFICATION_2COLOR:
            case SMRT_SUBREAD_IPD:
            case SMRT_SUBREAD_PW:
            case SMRT_CCS_FWD_IPD:
            case SMRT_CCS_FWD_PW:
            case SMRT_CCS_REV_IPD:
            case SMRT_CCS_REV_PW:
                // Just a simple forward/reverse strand color scheme that won't clash with the
                // methylation rectangles.
                c = (alignment.getFirstOfPairStrand() == Strand.POSITIVE) ? bisulfiteColorFw1 : bisulfiteColorRev1;

//                if (alignment.isNegativeStrand()) {
//                    c = (alignment.isSecondOfPair()) ? bisulfiteColorRev2 : bisulfiteColorRev1;
//                } else {
//                    c = (alignment.isSecondOfPair()) ? bisulfiteColorFw2 : bisulfiteColorFw1;
//                }
                break;
            case NOMESEQ:
                c = nomeseqColor;
                break;

            case UNEXPECTED_PAIR:
            case PAIR_ORIENTATION:
                c = getOrientationColor(alignment, getPEStats(alignment, renderOptions));
                if (c != defaultColor || colorOption == ColorOption.PAIR_ORIENTATION) {
                    break;
                }
            case INSERT_SIZE:
                boolean isPairedAlignment = alignment instanceof PairedAlignment;
                ReadMate mate = alignment.getMate();
                if ((alignment.isPaired() && mate != null && mate.isMapped()) || isPairedAlignment) {
//                   boolean sameChr = isPairedAlignment || alignment.getMate().getChr().equals(alignment.getChr());
                    String mateChr = mate == null ? null : mate.getChr();
                    boolean sameChr = isPairedAlignment || (mateChr != null && mateChr.equals(alignment.getChr()));
                    if (sameChr) {
                        int readDistance = Math.abs(alignment.getInferredInsertSize());
                        if (readDistance != 0) {

                            int minThreshold = renderOptions.getMinInsertSize();
                            int maxThreshold = renderOptions.getMaxInsertSize();
                            PEStats peStats = getPEStats(alignment, renderOptions);
                            if (renderOptions.isComputeIsizes() && peStats != null) {
                                minThreshold = peStats.getMinThreshold();
                                maxThreshold = peStats.getMaxThreshold();
                            }

                            if (readDistance < minThreshold) {
                                c = smallISizeColor;
                            } else if (readDistance > maxThreshold) {
                                c = largeISizeColor;
                            }
                        }

                        //return renderOptions.insertSizeColorScale.getColor(readDistance);
                    } else {
                        c = ChromosomeColors.getColor(alignment.getMate().getChr());
                        if (c == null) {
                            c = Color.black;
                        }
                    }
                }


                break;
            case READ_STRAND:
                if (alignment.isNegativeStrand()) {
                    c = negStrandColor;
                } else {
                    c = posStrandColor;
                }
                break;
            case FIRST_OF_PAIR_STRAND:
                final Strand fragmentStrand = alignment.getFirstOfPairStrand();
                if (fragmentStrand == Strand.NEGATIVE) {
                    c = negStrandColor;
                } else if (fragmentStrand == Strand.POSITIVE) {
                    c = posStrandColor;
                }
                break;
            case READ_ORDER:
                if (alignment.isPaired()) {
                    if (alignment.isFirstOfPair() && !alignment.isSecondOfPair()) {
                        c = firstOfPairColor;
                    } else if (!alignment.isFirstOfPair() && alignment.isSecondOfPair()) {
                        c = secondOfPairColor;
                    } else if (alignment.isFirstOfPair() && alignment.isSecondOfPair()) {
                        c = firstAndSecondofPairColor;
                    } else {
                        c = neitherForOrSecondOfPair;
                    }
                }
                break;
            case READ_GROUP:
                String rg = alignment.getReadGroup();
                if (rg != null) {
                    c = readGroupColors.get(rg);
                }
                break;
            case SAMPLE:
                String sample = alignment.getSample();
                if (sample != null) {
                    c = sampleColors.get(sample);
                }
                break;
            case LIBRARY:
                String library = alignment.getLibrary();
                if (library != null) {
                    c = sampleColors.get(library);
                }
                break;
            case MOVIE:
                readNameParts = alignment.getReadName().split("/");
                if (readNameParts.length >= 3) {
                    c = movieColors.get(readNameParts[0]);
                }
                break;
            case ZMW:
                readNameParts = alignment.getReadName().split("/");
                if (readNameParts.length >= 3) {
                    c = zmwColors.get(readNameParts[0] + "/" + readNameParts[1]);
                }
                break;
            case TAG:
                final String tag = renderOptions.getColorByTag();
                if (tag != null) {
                    Object tagValue = !colorByTagValueList.handlesTag(tag) ? alignment.getAttribute(tag) : colorByTagValueList.getValueForColorByTag(alignment, tag);
                    if (tagValue != null) {

                        ColorTable ctable;
                        String ctableKey;

                        String groupByTag = renderOptions.getGroupByTag();
                        if (groupByTag == null) {
                            ctable = defaultTagColors;

                        } else {
                            Object g = alignment.getAttribute(groupByTag);
                            String group = g == null ? "" : g.toString();
                            ctableKey = groupByTag + ":" + group;
                            ctable = tagValueColors.get(ctableKey);
                            if (ctable == null) {

                                if (groupByTag.equals("HP")) {
                                    ctable = getTenXColorTable(group);
                                } else {
                                    ctable = defaultTagColors;
                                }

                                tagValueColors.put(group, ctable);
                            }
                        }

                        c = ctable.get(tagValue.toString());

                    }
                }
                break;
            //     case LINK_STRAND:
            //         if (alignment instanceof LinkedAlignment && ((LinkedAlignment) alignment).getStrand() == Strand.NONE) {
            //             c = LL_COLOR;
            //         }
            //         break;


            default:
//                if (renderOptions.shadeCenters && center >= alignment.getStart() && center <= alignment.getEnd()) {
//                    if (locScale < 1) {
//                        c = grey2;
//                    }
//                }

        }
        if (c == null) c = defaultColor;

        c = shadeByMappingQuality(c, renderOptions, alignment.getMappingQuality());
        return c;

    }

    private Color shadeByMappingQuality(final Color initialColor, final AlignmentTrack.RenderOptions renderOptions, final int mappingQuality) {
        final float minAlpha = 0.15f;
        final float maxAlpha = 1;
        int maxMapQCutoff = renderOptions.getMappingQualityHigh();
        int minMapQCutoff = renderOptions.getMappingQualityLow();
        maxMapQCutoff = Ints.constrainToRange(maxMapQCutoff, 1, 255);
        minMapQCutoff = Ints.constrainToRange(minMapQCutoff, 0, maxMapQCutoff - 1);

        int clippedMQ = Ints.constrainToRange(mappingQuality, minMapQCutoff, maxMapQCutoff);
        float alphaRange = maxAlpha - minAlpha;
        float normalizedMQ = (float) (clippedMQ - minMapQCutoff) / (float) (maxMapQCutoff - minMapQCutoff);
        // Assuming white background TODO -- this should probably be passed in
        final Color backgroundColor = Color.white;

        // MQ of zero has special meaning, and pre-empts shading by mapping quality if "flag zero quality" is set
        if (mappingQuality == 0 && renderOptions.isFlagZeroQualityAlignments()) {
            return ColorUtilities.getCompositeColor(backgroundColor, initialColor, minAlpha);
        } else {
            switch (renderOptions.getShadeAlignmentsOption()) {
                case NONE:
                    return initialColor;
                case MAPPING_QUALITY_LOW:
                    normalizedMQ = 1 - normalizedMQ; //invert this and fall through.
                case MAPPING_QUALITY_HIGH:
                    float actualAlpha = minAlpha + alphaRange * normalizedMQ;
                    return ColorUtilities.getCompositeColor(backgroundColor, initialColor, actualAlpha);
                default:
                    return initialColor;
            }
        }
    }

    private ColorTable getTenXColorTable(String group) {
        ColorTable ctable;
        if (group.equals("1")) {
            ctable = tenXColorTable1;

        } else if (group.equals("2")) {
            ctable = tenXColorTable2;
        } else {
            ctable = tenXColorTable3;
        }
        return ctable;
    }

    public static PEStats getPEStats(Alignment alignment, AlignmentTrack.RenderOptions renderOptions) {
        String lb = alignment.getLibrary();
        if (lb == null) lb = "null";
        PEStats peStats = null;
        if (renderOptions.peStats != null) {
            peStats = renderOptions.peStats.get(lb);
        }
        return peStats;
    }

    /**
     * Determine if alignment insert size is outside max or min
     * range for outliers.
     *
     * @param alignment
     * @param renderOptions
     * @return -1 if unknown (stats not computed), 0 if not
     * an outlier, 1 if outlier
     */
    private int getOutlierStatus(Alignment alignment, AlignmentTrack.RenderOptions renderOptions) {
        PEStats peStats = getPEStats(alignment, renderOptions);
        if (renderOptions.isComputeIsizes() && peStats != null) {
            int minThreshold = peStats.getMinOutlierInsertSize();
            int maxThreshold = peStats.getMaxOutlierInsertSize();
            int dist = Math.abs(alignment.getInferredInsertSize());
            if (dist >= minThreshold || dist <= maxThreshold) {
                return 1;
            } else {
                return 0;
            }
        } else {
            return -1;
        }
    }

    /**
     * Returns -1 if alignment distance is less than minimum,
     * 0 if within bounds, and +1 if above maximum.
     *
     * @param alignment
     * @return
     */
    private int compareToBounds(Alignment alignment, AlignmentTrack.RenderOptions renderOptions) {
        int minThreshold = renderOptions.getMinInsertSize();
        int maxThreshold = renderOptions.getMaxInsertSize();
        PEStats peStats = getPEStats(alignment, renderOptions);
        if (renderOptions.isComputeIsizes() && peStats != null) {
            minThreshold = peStats.getMinThreshold();
            maxThreshold = peStats.getMaxThreshold();
        }

        int dist = Math.abs(alignment.getInferredInsertSize());
        if (dist < minThreshold) return -1;
        if (dist > maxThreshold) return +1;
        return 0;
    }

    /**
     * Returns a color to flag unexpected pair orientations.  Expected orientations (e.g. FR for Illumina) get the
     * neutral grey color
     *
     * @param alignment
     * @param peStats
     * @return
     */
    private Color getOrientationColor(Alignment alignment, PEStats peStats) {
        AlignmentTrack.OrientationType type = getOrientationType(alignment, peStats);
        Color c = typeToColorMap.get(type);
        return c == null ? track.getColor() : c;
    }

    static AlignmentTrack.OrientationType getOrientationType(Alignment alignment, PEStats peStats) {
        AlignmentTrack.OrientationType type = null;
        if (alignment.isPaired()) { // && !alignment.isProperPair()) {
            final String pairOrientation = alignment.getPairOrientation();
            if (peStats != null) {
                PEStats.Orientation libraryOrientation = peStats.getOrientation();
                switch (libraryOrientation) {
                    case FR:
                        //if (!alignment.isSmallInsert()) {
                        // if the isize < read length the reads might overlap, invalidating this test
                        type = frOrientationTypes.get(pairOrientation);
                        //}
                        break;
                    case RF:
                        type = rfOrientationTypes.get(pairOrientation);
                        break;
                    case F1F2:
                        type = f1f2OrientationTypes.get(pairOrientation);
                        break;
                    case F2F1:
                        type = f2f1OrientationTypes.get(pairOrientation);
                        break;
                }
            } else {
                // No peStats for this library, just guess
                if (alignment.getAttribute("CS") != null) {
                    type = f2f1OrientationTypes.get(pairOrientation);
                } else {
                    type = frOrientationTypes.get(pairOrientation);
                }
            }
        }

        return type;
    }

}