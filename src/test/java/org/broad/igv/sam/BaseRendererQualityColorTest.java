package org.broad.igv.sam;

import org.junit.Test;

import java.awt.Color;
import java.util.ArrayList;
import java.util.List;
import java.util.concurrent.CyclicBarrier;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.Future;
import java.util.concurrent.TimeUnit;

import static org.junit.Assert.*;

/** Public helper regressions against the pre-optimization arithmetic, not cache internals. */
public class BaseRendererQualityColorTest {
    private static final int[][] THRESHOLDS = {
            {5, 20}, {0, 60}, {-20, -5}, {-128, 127}, {10, 10}, {20, 5},
            {-300, 300}, {200, 300}, {-300, -200},
            {Integer.MIN_VALUE, Integer.MAX_VALUE}, {Integer.MIN_VALUE, 20},
            {-20, Integer.MAX_VALUE}, {Integer.MAX_VALUE, Integer.MIN_VALUE},
            {Integer.MIN_VALUE, Integer.MIN_VALUE}, {Integer.MAX_VALUE, Integer.MAX_VALUE}
    };

    @Test
    public void everySignedByteAndThresholdKeepsOldFloatRoundingAndTruncatedAlpha() {
        Color[] colors = {Color.BLACK, Color.WHITE, new Color(19, 173, 83),
                new Color(37, 71, 211, 0), new Color(37, 71, 211, 91),
                new Color(37, 71, 211, 255), new Color(254, 1, 127, 203)};
        for (int[] thresholds : THRESHOLDS) {
            for (Color color : colors) {
                for (int bits = 0; bits < 256; bits++) {
                    byte quality = (byte) bits;
                    Color actual = BaseRenderer.getShadedColor(color, quality, thresholds[0], thresholds[1]);
                    assertEquals("quality=" + quality + ", min=" + thresholds[0] + ", max=" + thresholds[1],
                            legacyColor(color, quality, thresholds[0], thresholds[1]), actual);
                    if (quality >= thresholds[1]) {
                        assertSame("Unshaded colors retain source identity and alpha", color, actual);
                    }
                }
            }
        }
    }

    @Test
    public void equalRgbSourcesIgnoreSourceAlphaOnlyWhenShaded() {
        Color transparent = new Color(37, 71, 211, 0);
        Color translucent = new Color(37, 71, 211, 91);
        for (int bits = 0; bits < 256; bits++) {
            byte quality = (byte) bits;
            Color first = BaseRenderer.getShadedColor(transparent, quality, 5, 20);
            Color second = BaseRenderer.getShadedColor(translucent, quality, 5, 20);
            if (quality < 20) {
                assertEquals(first, second);
                assertEquals(legacyColor(transparent, quality, 5, 20), first);
            } else {
                assertSame(transparent, first);
                assertSame(translucent, second);
                assertNotEquals(first, second);
            }
        }
        assertEquals(51, BaseRenderer.getShadedColor(translucent, (byte) 5, 5, 20).getAlpha());
        assertEquals(102, BaseRenderer.getShadedColor(translucent, (byte) 10, 5, 20).getAlpha());
        assertEquals(178, BaseRenderer.getShadedColor(translucent, (byte) 15, 5, 20).getAlpha());
        assertEquals(229, BaseRenderer.getShadedColor(translucent, (byte) 19, 5, 20).getAlpha());
    }

    @Test
    public void concurrentCallsWithManyCompetingRgbValuesNeverReturnAnotherColorsShade() throws Exception {
        int workers = 8;
        CyclicBarrier roundStart = new CyclicBarrier(workers);
        ExecutorService executor = Executors.newFixedThreadPool(workers);
        List<Future<?>> results = new ArrayList<>();
        try {
            for (int worker = 0; worker < workers; worker++) {
                final int lane = worker;
                results.add(executor.submit(() -> {
                    for (int round = 0; round < 256; round++) {
                        roundStart.await(30, TimeUnit.SECONDS);
                        // Distinct RGBs, shared low components, revisited in opposite orders.
                        // Thousands of competing keys exercise eviction/collisions without
                        // depending on cache size, indexing, or universal shaded identity.
                        int index = ((lane & 1) == 0 ? round : 255 - round) * workers + lane;
                        Color source = new Color((index >>> 3) & 255, index & 255, 83, (index * 31) & 255);
                        int[] thresholds = THRESHOLDS[round % THRESHOLDS.length];
                        for (int bits = 0; bits < 256; bits++) {
                            byte quality = (byte) bits;
                            Color actual = BaseRenderer.getShadedColor(source, quality, thresholds[0], thresholds[1]);
                            assertEquals(legacyColor(source, quality, thresholds[0], thresholds[1]), actual);
                            if (quality >= thresholds[1]) assertSame(source, actual);
                        }
                    }
                    return null;
                }));
            }
            for (Future<?> result : results) result.get(60, TimeUnit.SECONDS);
        } finally {
            executor.shutdownNow();
        }
    }

    // Copied unchanged from the old helper's branches, float expression, and rounding.
    // Only the removed String/Map lookup is replaced by direct Color construction.
    static Color legacyColor(Color color, byte qual, int minQ, int maxQ) {
        if (qual >= maxQ) return color;

        float alpha;
        if (qual < minQ) {
            alpha = 0.2f;
        } else {
            alpha = Math.max(0.2f, Math.min(1.0f, 0.1f + 0.9f * (qual - minQ) / (maxQ - minQ)));
        }

        // Round alpha to nearest 0.1
        alpha = ((int) (alpha * 10 + 0.5f)) / 10.0f;
        return new Color(color.getRed(), color.getGreen(), color.getBlue(), (int) (alpha * 255));
    }
}
