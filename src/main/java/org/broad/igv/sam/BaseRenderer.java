package org.broad.igv.sam;

import org.apache.batik.svggen.SVGGraphics2D;
import org.broad.igv.renderer.GraphicUtils;
import org.broad.igv.renderer.SequenceRenderer;
import org.broad.igv.sam.mods.BaseModficationFilter;
import org.broad.igv.sam.mods.BaseModificationRenderer;
import org.broad.igv.sam.mods.BaseModificationSet;
import org.broad.igv.track.RenderContext;
import org.broad.igv.ui.FontManager;
import org.broad.igv.ui.color.ColorUtilities;

import java.awt.*;
import java.awt.geom.Rectangle2D;
import java.awt.image.BufferedImage;
import java.awt.image.DataBufferInt;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

/**
 * Some static methods for rendering bases in the context of alignments.  *
 */

public class BaseRenderer {

    /**
     * Reusable, single-paint-context scratch for one color per pixel column.
     * The supplied graphics must be dedicated to the strip; its composite,
     * transform and clip are used directly, without pre-compositing the colors.
     */
    public static final class ColorStrip {

        private Graphics2D graphics;
        private boolean vector;
        private int startX;
        private int width;
        private int first;
        private int last;
        private BufferedImage image;
        private int[] pixels;
        private Color[] vectorColors;

        public void reset(Graphics2D graphics, int startX, int width) {
            this.graphics = graphics;
            this.startX = startX;
            this.width = width;
            first = width;
            last = -1;
            vector = graphics instanceof SVGGraphics2D;

            if (vector) {
                if (vectorColors == null || vectorColors.length < width) {
                    vectorColors = new Color[width];
                } else {
                    Arrays.fill(vectorColors, 0, width, null);
                }
            } else {
                if (image == null || image.getWidth() < width) {
                    image = new BufferedImage(width, 1, BufferedImage.TYPE_INT_ARGB);
                    pixels = ((DataBufferInt) image.getRaster().getDataBuffer()).getData();
                } else {
                    Arrays.fill(pixels, 0, width, 0);
                }
                graphics.setRenderingHint(RenderingHints.KEY_INTERPOLATION,
                        RenderingHints.VALUE_INTERPOLATION_NEAREST_NEIGHBOR);
            }
        }

        public void setColor(int pixelX, Color color) {
            int index = pixelX - startX;
            if (index < 0 || index >= width) return;

            if (vector) {
                vectorColors[index] = color;
            } else {
                pixels[index] = color == null ? 0 : color.getRGB();
            }
            if (color != null && color.getAlpha() != 0) {
                first = Math.min(first, index);
                last = Math.max(last, index);
            }
        }

        public void draw(int y, int height) {
            if (height <= 0 || last < first) return;

            if (vector) {
                int index = first;
                while (index <= last) {
                    Color color = vectorColors[index];
                    if (color == null || color.getAlpha() == 0) {
                        index++;
                        continue;
                    }
                    int rgba = color.getRGB();
                    int end = index + 1;
                    while (end <= last && vectorColors[end] != null &&
                            vectorColors[end].getRGB() == rgba) {
                        end++;
                    }
                    graphics.setColor(color);
                    graphics.fillRect(startX + index, y, end - index, height);
                    index = end;
                }
            } else {
                // Transparent endpoints can result from overwriting a populated column.
                while (first <= last && (pixels[first] >>> 24) == 0) first++;
                while (last >= first && (pixels[last] >>> 24) == 0) last--;
                if (last < first) return;
                graphics.drawImage(image, startX + first, y, startX + last + 1, y + height,
                        first, 0, last + 1, 1, null);
            }
        }
    }

    /**
     * Draw the base using either a letter or character, using the given color,
     * depending on size and bisulfite status
     *
     * @param g
     * @param color
     * @param c
     * @param pX
     * @param pY
     * @param dX
     * @param dY
     * @param bisulfiteMode
     * @param bisstatus
     */
    static void drawBase(Graphics2D g, Color color, char c, int pX, int pY, int dX, int dY, boolean bisulfiteMode,
                         BisulfiteBaseInfo.DisplayStatus bisstatus) {

        int fontSize = Math.min(Math.min(dX, dY), 12);
        if (fontSize >= 8 && (!bisulfiteMode || (bisulfiteMode && bisstatus.equals(BisulfiteBaseInfo.DisplayStatus.CHARACTER)))) {
            Font f = FontManager.getFont(Font.BOLD, fontSize);
            g.setFont(f);
            g.setColor(color);
            GraphicUtils.drawCenteredText(new char[]{c}, pX, pY, dX, dY, g);
        } else {

            // If bisulfite mode, we expand the rectangle to make it more visible
            if (bisulfiteMode && bisstatus.equals(BisulfiteBaseInfo.DisplayStatus.RECTANGLE)) {
                if (dX < 3) {
                    int expansion = dX;
                    pX -= expansion;
                    dX += (2 * expansion);
                }
            }

            if (color != null) {
                g.setColor(color);
                g.fillRect(pX, pY, dX, dY);
            }
        }
    }

    public static void drawExpandedInsertions(InsertionMarker insertionMarker,
                                              List<Alignment> alignments,
                                              RenderContext context,
                                              Rectangle rect,
                                              boolean leaveMargin,
                                              AlignmentTrack.RenderOptions renderOptions) {

        Graphics2D g = null;
        try {
            g = (Graphics2D) context.getGraphics().create();
            g.setRenderingHint(RenderingHints.KEY_ANTIALIASING, RenderingHints.VALUE_ANTIALIAS_ON);
            double dX = 1 / context.getScale();
            int fontSize = (int) Math.min(dX, 12);
            if (fontSize >= 8) {
                Font f = FontManager.getFont(Font.BOLD, fontSize);
                g.setFont(f);
            }
            AlignmentTrack.ColorOption colorOption = renderOptions.getColorOption();

            for (Alignment alignment : alignments) {

                if (alignment.getEnd() < insertionMarker.position) continue;
                if (alignment.getStart() > insertionMarker.position) break;

                AlignmentBlock insertion = alignment.getInsertionAt(insertionMarker.position);
                if (insertion != null && insertion.getBasesLength() > 0) {

                    double origin = context.getOrigin();
                    double locScale = context.getScale();

                    // Compute the start and end of the expanded insertion in pixels
                    int pixelStart = (int) ((insertion.getStart() - origin) / locScale);
                    int pixelEnd = (int) (pixelStart + insertion.getLength() / locScale);

                    // Skip if insertion is out of clipping rectangle -- this probably shouldn't happen
                    if (pixelEnd < rect.x || pixelStart > rect.getMaxX()) {
                        continue;
                    }

                    ByteSubarray bases = insertion.getBases();
                    int padding = insertion.getPadding();
                    final int size = insertion.getBasesLength() + padding;
                    for (int p = 0; p < size; p++) {

                        double pX = (pixelStart + (p / locScale));

                        // Don't draw out of clipping rect
                        if (pX > rect.getMaxX()) break;
                        else if (pX + dX < rect.getX()) continue;

                        // idx can be out of range due to (1) padding, (2) no read sequence (sequence == *)
                        int idx = p - padding;
                        char c = idx >= 0 && idx < bases.length ? (char) bases.getByte(p - padding) :
                                padding > 0 ? '-' : '*';

                        Color color = null;
                        // TODO -- support bisulfite mode?  Probably not possible as that depends on the reference
                        //if (bisulfiteMode) {
                        //     color = bisinfo.getDisplayColor(idx);
                        // } else
                        if (colorOption.isBaseMod() ||
                                colorOption.isSMRTKinetics()) {
                            color = Color.GRAY;
                        } else {
                            //color = nucleotideColors.get(c);
                            color = SequenceRenderer.nucleotideColors.get(c);
                        }
                        if (color == null) {
                            color = Color.black;
                        }

                        if (renderOptions.getShadeBasesOption() && p >= padding) {
                            byte qual = p < padding ? (byte) 126 : insertion.getQuality(p - padding);
                            color = BaseRenderer.getShadedColor(color, qual, renderOptions.getBaseQualityMin(), renderOptions.getBaseQualityMax());
                        }

                        if (dX < 8) {
                            g.setColor(color);
                            g.fill(new Rectangle2D.Double(pX, rect.y, dX, rect.height));
                        } else {
                            drawBase(g, color, c, (int) pX, rect.y, (int) dX, rect.height - (leaveMargin ? 2 : 0), false, null);
                        }
                    }
                    insertion.setPixelRange(context.translateX + pixelStart, context.translateX + pixelEnd);


                    if (colorOption.isBaseMod()) {
                        List<BaseModificationSet> baseModificationSets = alignment.getBaseModificationSets();
                        if (baseModificationSets != null) {
                            BaseModificationRenderer.drawBlock(origin, locScale, rect, g, renderOptions, baseModificationSets, insertion);
                        }
                    }
                }
            }
        } finally {
            if (g != null) g.dispose();
        }

    }

    public static Color getShadedColor(Color color, byte qual, int minQ, int maxQ) {

        if (qual >= maxQ) return color;

        float alpha;
        if (qual < minQ) {
            alpha = 0.2f;
        } else {
            alpha = Math.max(0.2f, Math.min(1.0f, 0.1f + 0.9f * (qual - minQ) / (maxQ - minQ)));
        }

        // Round alpha to nearest 0.1
        alpha = ((int) (alpha * 10 + 0.5f)) / 10.0f;
        String key = ColorUtilities.colorToString(color) + "_" + alpha;
        Color c = shadedColorCache.get(key);
        if (c == null) {
            c = ColorUtilities.modifyAlpha(color, (int) (alpha * 255));
            shadedColorCache.put(key, c);
        }
        return c;
    }

    static Map<String, Color> shadedColorCache = Collections.synchronizedMap(new HashMap<>());

}
