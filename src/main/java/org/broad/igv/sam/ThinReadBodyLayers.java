package org.broad.igv.sam;

import java.awt.*;
import java.awt.geom.Path2D;
import java.util.Arrays;

/**
 * Paint-owned layers of normalized thin-body device rectangles. Disjoint runs
 * are filled once per overdraw layer, retaining the original SRC_OVER rounding
 * and shape clip instead of flattening translucent contributions into an image.
 */
public final class ThinReadBodyLayers {
    private static final int DEVICE_LIMIT = 1 << 20;
    private Graphics2D graphics;
    private int[] coverage;
    private Path2D.Float path;
    private int capacity;
    private int startX;
    private int width;
    private int top;
    private int height;
    private int first;
    private int last;

    /** The caller has already restricted graphics to normalized device fills. */
    public boolean reset(Graphics2D graphics, Rectangle clip, int top, int height) {
        draw();
        this.graphics = null;
        // Use bounds only for scratch sizing. Graphics retains the complete
        // device shape clip, including curved damage and caps outside the row.
        if (clip == null || clip.width <= 0 || clip.width > DEVICE_LIMIT ||
                Math.abs((long) clip.x) > DEVICE_LIMIT ||
                Math.abs((long) clip.x + clip.width) > DEVICE_LIMIT ||
                Math.abs((long) top) > DEVICE_LIMIT ||
                Math.abs((long) top + height) > DEVICE_LIMIT) {
            return false;
        }
        startX = clip.x;
        width = clip.width;
        if (capacity < width) {
            coverage = new int[width + 1];
            // At most ceil(width / 2) disjoint runs, five path segments each.
            path = new Path2D.Float(Path2D.WIND_NON_ZERO, ((width + 1) / 2) * 5);
            capacity = width;
        }
        this.graphics = graphics;
        this.top = top;
        this.height = height;
        first = width;
        last = 0;
        return true;
    }

    /** Add the already normalized, exclusive-right device rectangle. */
    public void fill(int left, int right) {
        int start = (int) Math.max(0L, (long) left - startX);
        int end = (int) Math.min(width, (long) right - startX);
        if (start >= end) return;
        coverage[start]++;
        coverage[end]--;
        first = Math.min(first, start);
        last = Math.max(last, end);
    }

    /** Flush before changing color, painting another primitive or resetting. */
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
                path.moveTo(startX + left, top);
                path.lineTo(startX + i, top);
                path.lineTo(startX + i, top + height);
                path.lineTo(startX + left, top + height);
                path.closePath();
            }
            graphics.fill(path);
        }
        Arrays.fill(coverage, first, last + 1, 0);
        first = width;
        last = 0;
    }
}
