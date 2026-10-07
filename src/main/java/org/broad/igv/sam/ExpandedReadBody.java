package org.broad.igv.sam;

import java.awt.*;
import java.awt.geom.AffineTransform;
import java.awt.geom.Path2D;
import java.util.Arrays;

/**
 * Paint-owned scratch for one alignment's ordinary expanded rectangles. Each
 * coverage layer is a set of disjoint integer rectangles painted with the original
 * graphics, so repeated SRC_OVER paints retain Java2D's per-layer rounding.
 */
public final class ExpandedReadBody {
    private static final int DEVICE_LIMIT = 1 << 20;
    private Graphics2D graphics;
    private int[] coverage;
    private Path2D.Float path;
    private int capacity;
    private int startX;
    private int width;
    private int y;
    private int height;
    private int first;
    private int last;
    private double scale;
    private double translateX;

    /** Reset before each alignment; null disables scratch for an ineligible body. */
    public boolean reset(Graphics2D graphics, Rectangle row, int height) {
        this.graphics = null;
        if (graphics == null || !graphics.getClass().getName().equals("sun.java2d.SunGraphics2D") ||
                graphics.getPaint().getClass() != Color.class ||
                !(graphics.getComposite() instanceof AlphaComposite composite) ||
                composite.getRule() != AlphaComposite.SRC_OVER ||
                graphics.getRenderingHint(RenderingHints.KEY_ANTIALIASING) == RenderingHints.VALUE_ANTIALIAS_ON ||
                height <= 1 || row.width <= 0) {
            return false;
        }
        AffineTransform transform = graphics.getTransform();
        scale = transform.getScaleX();
        translateX = transform.getTranslateX();
        double translateY = transform.getTranslateY();
        long right = (long) row.x + row.width;
        long bottom = (long) row.y + height;
        if (transform.getShearX() != 0 || transform.getShearY() != 0 ||
                (scale != 1 && scale != 2) || transform.getScaleY() != scale ||
                translateX != Math.rint(translateX) || translateY != Math.rint(translateY) ||
                !deviceCoordinate(translateX) || !deviceCoordinate(translateY) ||
                !deviceCoordinate(row.y * scale + translateY) || !deviceCoordinate(bottom * scale + translateY)) {
            return false;
        }
        // Do not guess a raster's bounds with a null clip. Keep the original clip
        // (including curved damage); its bounds only limit horizontal scratch.
        Rectangle clip = graphics.getClipBounds();
        if (clip == null) return false;
        long left = Math.max(0L, clip.x);
        right = Math.min(right, (long) clip.x + clip.width);
        if (right <= left || right - left > DEVICE_LIMIT ||
                !deviceCoordinate(left * scale + translateX) || !deviceCoordinate(right * scale + translateX)) {
            return false;
        }
        startX = (int) left;
        width = (int) (right - left);
        if (capacity < width) {
            coverage = new int[width + 1];
            // At most ceil(width / 2) disjoint runs per layer, five segments each.
            path = new Path2D.Float(Path2D.WIND_NON_ZERO, ((width + 1) / 2) * 5);
            capacity = width;
        } else {
            Arrays.fill(coverage, 0, width + 1, 0);
        }
        this.graphics = graphics;
        y = row.y;
        this.height = height;
        first = width;
        last = 0;
        return true;
    }

    private static boolean deviceCoordinate(double value) {
        return Math.abs(value) <= DEVICE_LIMIT;
    }

    /** Add the original exclusive-right integer rectangle, without per-block work arrays. */
    public boolean fill(int start, int end) {
        if (graphics == null || end < start ||
                !deviceCoordinate(start * scale + translateX) || !deviceCoordinate(end * scale + translateX)) {
            return false;
        }
        int left = (int) Math.max(0L, (long) start - startX);
        int right = (int) Math.min(width, (long) end - startX);
        if (left >= right) return true; // Includes zero-width rectangles: no line caps.
        coverage[left]++;
        coverage[right]--;
        first = Math.min(first, left);
        last = Math.max(last, right);
        return true;
    }

    /** Flush before any other paint, then allow another batch in the same alignment. */
    public void draw() {
        if (graphics == null || first >= last) return;
        int count = 0;
        int layers = 0;
        for (int i = first; i < last; i++) {
            count += coverage[i];
            coverage[i] = count;
            layers = Math.max(layers, count);
        }
        for (int layer = 1; layer <= layers; layer++) {
            path.reset();
            int i = first;
            while (i < last) {
                if (coverage[i] < layer) {
                    i++;
                    continue;
                }
                int left = i++;
                while (i < last && coverage[i] >= layer) i++;
                path.moveTo(startX + left, y);
                path.lineTo(startX + i, y);
                path.lineTo(startX + i, y + height);
                path.lineTo(startX + left, y + height);
                path.closePath();
            }
            graphics.fill(path);
        }
        Arrays.fill(coverage, first, last + 1, 0);
        first = width;
        last = 0;
    }
}
