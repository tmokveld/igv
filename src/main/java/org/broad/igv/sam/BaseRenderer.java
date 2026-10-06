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
import java.util.IdentityHashMap;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

/**
 * Some static methods for rendering bases in the context of alignments.  *
 */

public class BaseRenderer {

    /**
     * Immutable palette/options snapshot. Resolve once per row (or insertion
     * paint), never per block/base. Source identity matters for unshaded and
     * high-quality returns; shaded colors deliberately share by RGB alone.
     */
    static final class PaletteColors {
        private final Map<Character, Color> palette;
        final boolean shade;
        final int minQ;
        final int maxQ;
        private final Color[] originals = new Color[256];
        private final Color[][] shaded;

        PaletteColors(Map<Character, Color> palette, boolean shade, int minQ, int maxQ) {
            this.palette = new HashMap<>(palette);
            this.shade = shade;
            this.minQ = minQ;
            this.maxQ = maxQ;
            shaded = shade ? new Color[256][] : null;
            int[] buckets = shade ? new int[256] : null;
            if (shade) {
                for (int i = 0; i < 256; i++) buckets[i] = qualityBucket((byte) i, minQ, maxQ);
            }
            Map<Color, Color[]> byIdentity = shade ? new IdentityHashMap<>() : null;
            for (int i = 0; i < 256; i++) {
                // Preserve the renderer's signed-byte-to-char conversion.
                Color color = this.palette.get((char) (byte) i);
                if (color == null) color = Color.BLACK;
                originals[i] = color;
                if (shade) {
                    Color[] qualities = byIdentity.get(color);
                    if (qualities == null) {
                        qualities = new Color[256];
                        Color[] alphas = new Color[11];
                        for (int q = 0; q < 256; q++) {
                            int bucket = buckets[q];
                            if (bucket < 0) {
                                qualities[q] = color;
                            } else {
                                if (alphas[bucket] == null) alphas[bucket] = shadedColor(color, bucket);
                                qualities[q] = alphas[bucket];
                            }
                        }
                        byIdentity.put(color, qualities);
                    }
                    shaded[i] = qualities;
                }
            }
        }

        boolean matches(Map<Character, Color> palette, boolean shade, int minQ, int maxQ) {
            if (this.shade != shade || this.minQ != minQ || this.maxQ != maxQ ||
                    this.palette.size() != palette.size()) return false;
            for (Map.Entry<Character, Color> entry : this.palette.entrySet()) {
                if (!palette.containsKey(entry.getKey()) || palette.get(entry.getKey()) != entry.getValue()) return false;
            }
            return true;
        }

        Color getColor(byte base, byte quality) {
            return shade ? shaded[base & 0xff][quality & 0xff] : originals[base & 0xff];
        }
    }

    private static volatile PaletteColors insertionPaletteColors;

    /**
     * Reusable scratch for ordered SRC_OVER overlays. Columns may be revisited;
     * each retains full-precision premultiplied channels until the strip is drawn.
     * The dedicated graphics' extra alpha applies to every layer, not the strip.
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
        private int column;
        private double[] channels;
        private double extraAlpha;
        private double alpha;
        private double red;
        private double green;
        private double blue;

        public void reset(Graphics2D graphics, int startX, int width) {
            this.graphics = graphics;
            this.startX = startX;
            this.width = width;
            first = width;
            last = -1;
            column = -1;
            alpha = red = green = blue = 0;
            extraAlpha = ((AlphaComposite) graphics.getComposite()).getAlpha();
            if (channels == null || channels.length < width * 4) {
                channels = new double[width * 4];
            } else {
                Arrays.fill(channels, 0, width * 4, 0);
            }
            vector = graphics instanceof SVGGraphics2D;

            if (vector) {
                if (pixels == null || pixels.length < width) {
                    pixels = new int[width];
                } else {
                    Arrays.fill(pixels, 0, width, 0);
                }
            } else {
                if (image == null || image.getWidth() < width) {
                    image = new BufferedImage(width, 1, BufferedImage.TYPE_INT_ARGB);
                }
                pixels = ((DataBufferInt) image.getRaster().getDataBuffer()).getData();
                Arrays.fill(pixels, 0, width, 0);
                graphics.setRenderingHint(RenderingHints.KEY_INTERPOLATION,
                        RenderingHints.VALUE_INTERPOLATION_NEAREST_NEIGHBOR);
            }
        }

        public void blendColor(int pixelX, Color color) {
            int index = pixelX - startX;
            if (index < 0 || index >= width || color == null) return;
            double sourceAlpha = color.getAlpha() / 255.0 * extraAlpha;
            if (sourceAlpha == 0) return;
            if (index != column) {
                finishColumn();
                column = index;
                int offset = index * 4;
                alpha = channels[offset];
                red = channels[offset + 1];
                green = channels[offset + 2];
                blue = channels[offset + 3];
            }
            double remaining = 1 - sourceAlpha;
            red = color.getRed() * sourceAlpha + red * remaining;
            green = color.getGreen() * sourceAlpha + green * remaining;
            blue = color.getBlue() * sourceAlpha + blue * remaining;
            alpha = sourceAlpha + alpha * remaining;
        }

        private void finishColumn() {
            if (column < 0) return;
            int offset = column * 4;
            channels[offset] = alpha;
            channels[offset + 1] = red;
            channels[offset + 2] = green;
            channels[offset + 3] = blue;
            // Keep premultiplied channels at full precision until the column is
            // complete. Rounding once can differ slightly from repeated Java2D fills.
            pixels[column] = ((int) Math.round(alpha * 255) << 24) |
                    ((int) Math.round(red / alpha) << 16) |
                    ((int) Math.round(green / alpha) << 8) |
                    (int) Math.round(blue / alpha);
            first = Math.min(first, column);
            last = Math.max(last, column);
        }

        public void draw(int y, int height) {
            finishColumn();
            if (height <= 0 || last < first) return;
            Composite composite = graphics.getComposite();
            graphics.setComposite(AlphaComposite.SrcOver);
            try {
                if (vector) {
                    int index = first;
                    while (index <= last) {
                        int rgba = pixels[index];
                        if ((rgba >>> 24) == 0) {
                            index++;
                            continue;
                        }
                        int end = index + 1;
                        while (end <= last && pixels[end] == rgba) end++;
                        graphics.setColor(new Color(rgba, true));
                        graphics.fillRect(startX + index, y, end - index, height);
                        index = end;
                    }
                } else {
                    graphics.drawImage(image, startX + first, y, startX + last + 1, y + height,
                            first, 0, last + 1, 1, null);
                }
            } finally {
                graphics.setComposite(composite);
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
            boolean shadeBases = renderOptions.getShadeBasesOption();
            int minQ = renderOptions.getBaseQualityMin();
            int maxQ = renderOptions.getBaseQualityMax();
            Map<Character, Color> palette = SequenceRenderer.getNucleotideColors();
            PaletteColors colors = insertionPaletteColors;
            if (colors == null || !colors.matches(palette, shadeBases, minQ, maxQ)) {
                colors = new PaletteColors(palette, shadeBases, minQ, maxQ);
                insertionPaletteColors = colors;
            }

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
                            byte quality = shadeBases && p >= padding ? insertion.getQuality(p - padding) : (byte) 126;
                            // Padding is never quality-shaded, even with unusual thresholds.
                            color = p < padding ? colors.originals[(byte) c & 0xff] : colors.getColor((byte) c, quality);
                        }
                        if (color == null) {
                            color = Color.black;
                        }

                        if (shadeBases && p >= padding && (colorOption.isBaseMod() || colorOption.isSMRTKinetics())) {
                            byte qual = insertion.getQuality(p - padding);
                            color = BaseRenderer.getShadedColor(color, qual, minQ, maxQ);
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
        int bucket = qualityBucket(qual, minQ, maxQ);
        return bucket < 0 ? color : shadedColor(color, bucket);
    }

    private static int qualityBucket(byte qual, int minQ, int maxQ) {
        if (qual >= maxQ) return -1;

        float alpha;
        if (qual < minQ) {
            alpha = 0.2f;
        } else {
            alpha = Math.max(0.2f, Math.min(1.0f, 0.1f + 0.9f * (qual - minQ) / (maxQ - minQ)));
        }
        return (int) (alpha * 10 + 0.5f);
    }

    // Bounded, direct-mapped RGB/tenth cache for arbitrary-color callers.
    // Palette lookups never acquire this lock. Collisions only evict a value;
    // neither source alpha nor thresholds belong in the shaded equivalence key.
    private static final Color[] shadedColors = new Color[256 * 11];

    private static synchronized Color shadedColor(Color color, int bucket) {
        int rgb = (color.getRed() << 16) | (color.getGreen() << 8) | color.getBlue();
        int slot = ((rgb ^ (rgb >>> 12)) & 0xff) * 11 + bucket;
        Color shaded = shadedColors[slot];
        if (shaded == null || (shaded.getRGB() & 0xffffff) != rgb) {
            float alpha = bucket / 10.0f;
            shaded = ColorUtilities.modifyAlpha(color, (int) (alpha * 255));
            shadedColors[slot] = shaded;
        }
        return shaded;
    }

}
