package org.broad.igv.sam;

import org.apache.batik.dom.GenericDOMImplementation;
import org.apache.batik.svggen.SVGGraphics2D;
import org.broad.igv.track.RenderContext;
import org.broad.igv.ui.panel.ReferenceFrame;
import org.junit.Test;

import java.awt.*;
import java.awt.geom.AffineTransform;
import java.awt.geom.Ellipse2D;
import java.awt.geom.Rectangle2D;
import java.awt.image.BufferedImage;
import java.awt.image.ColorModel;
import java.awt.image.DataBufferInt;
import java.util.List;
import java.util.function.Consumer;

import static org.junit.Assert.*;

/** Exact comparisons to separate, ordered Graphics2D.fillRect paints. */
public class ExpandedReadBodyTest {
    private static final Rectangle ROW = new Rectangle(0, 3, 100, 8);
    private static final int HEIGHT = 6;
    private static final Color INK = new Color(37, 149, 211, 113);
    private static final Color BACKGROUND = new Color(219, 183, 137, 197);
    private static final int[][] SPANS = {
            {0, 1}, {1, 4}, {4, 9}, {3, 7}, {7, 8}, {22, 23},
            {0, 1}, {4, 9}, {8, 8}, {99, 100}, {97, 100}, {22, 23}, {3, 7}
    };

    @Test
    public void adjacentOverlapRevisitedTinyZeroAndTerminalRectanglesKeepEveryAlphaLayer() {
        for (int type : new int[]{BufferedImage.TYPE_INT_ARGB, BufferedImage.TYPE_INT_ARGB_PRE}) {
            for (Color background : new Color[]{BACKGROUND, new Color(0, 0, 0, 0), new Color(241, 229, 201)}) {
                for (int alpha : new int[]{0, 1, 17, 51, 113, 173, 254, 255}) {
                    for (float extraAlpha : new float[]{0, 1, 0.75f, 0.37f}) {
                        for (int scale : new int[]{1, 2}) {
                            try (Pair pair = new Pair(type, background, g -> g.scale(scale, scale))) {
                                Color ink = new Color(37, 149, 211, alpha);
                                pair.color(ink, extraAlpha);
                                ExpandedReadBody body = new ExpandedReadBody();
                                assertTrue(body.reset(pair.actual, ROW, HEIGHT));
                                for (int[] span : SPANS) accept(body, pair.expected, ROW, HEIGHT, span[0], span[1]);
                                // Deep overlap exposes unioning or flattening of repeated alpha.
                                for (int layer = 0; layer < 37; layer++) accept(body, pair.expected, ROW, HEIGHT, 28, 34);
                                body.draw();
                                pair.assertExact();
                                body.draw();
                                pair.assertExact();
                            }
                        }
                    }
                }
            }
        }
    }

    @Test
    public void integerDeviceTranslationsAndDamageEdgesPreserveExclusiveRightCoverage() {
        for (int scale : new int[]{1, 2}) {
            for (int tx : new int[]{-3, 0, 5}) {
                for (int ty : new int[]{-2, 0, 7}) {
                    for (Shape damage : List.of(new Rectangle(0, 3, 19, 6),
                            new Rectangle(5, 4, 9, 3), new Rectangle(96, 2, 4, 8),
                            new Ellipse2D.Double(1.25, 2.25, 24.5, 9.5))) {
                        try (Pair pair = new Pair(BufferedImage.TYPE_INT_ARGB_PRE, BACKGROUND, g -> {
                            g.translate(tx, ty);
                            g.scale(scale, scale);
                            g.clip(damage);
                        })) {
                            pair.color(INK, 0.75f);
                            ExpandedReadBody body = new ExpandedReadBody();
                            assertTrue(body.reset(pair.actual, ROW, HEIGHT));
                            for (int[] span : new int[][]{{-6, 4}, {0, 19}, {4, 14}, {14, 15}, {96, 100}, {99, 100}, {0, 100}}) {
                                accept(body, pair.expected, ROW, HEIGHT, span[0], span[1]);
                            }
                            body.draw();
                            // Comparing the whole image also requires untouched pixels outside damage.
                            pair.assertExact();
                        }
                    }
                }
            }
        }
    }

    @Test
    public void flushThenResumePreservesDifferentGraphicsColorsAndPolygonOrder() {
        try (Pair pair = new Pair(BufferedImage.TYPE_INT_ARGB, BACKGROUND, g -> g.scale(2, 2))) {
            ExpandedReadBody body = new ExpandedReadBody();
            pair.color(INK, 0.37f);
            assertTrue(body.reset(pair.actual, ROW, HEIGHT));
            accept(body, pair.expected, ROW, HEIGHT, 2, 20);
            accept(body, pair.expected, ROW, HEIGHT, 6, 18);
            body.draw();
            Polygon decoration = new Polygon(new int[]{3, 15, 8}, new int[]{2, 9, 12}, 3);
            for (Graphics2D graphics : new Graphics2D[]{pair.actual, pair.expected}) {
                Graphics2D alternate = (Graphics2D) graphics.create();
                try {
                    alternate.setComposite(AlphaComposite.getInstance(AlphaComposite.SRC_OVER, 0.61f));
                    alternate.setColor(new Color(211, 43, 137, 179));
                    alternate.fill(decoration);
                } finally {
                    alternate.dispose();
                }
            }
            // draw() must clear coverage without disabling the current batch.
            accept(body, pair.expected, ROW, HEIGHT, 4, 19);
            accept(body, pair.expected, ROW, HEIGHT, 4, 5);
            body.draw();
            pair.assertExact();
            pair.color(new Color(211, 43, 137, 179), 0.61f);
            assertTrue(body.reset(pair.actual, ROW, HEIGHT));
            accept(body, pair.expected, ROW, HEIGHT, 1, 22);
            body.draw();
            pair.assertExact();
        }
    }

    @Test
    public void wideNarrowEmptyAndDiscardedResetsNeverLeakScratch() {
        try (Pair pair = new Pair(BufferedImage.TYPE_INT_ARGB_PRE, BACKGROUND, g -> {
            g.scale(2, 2);
            g.clip(new Rectangle(0, 0, 100, 40));
        })) {
            pair.color(INK, 0.75f);
            ExpandedReadBody body = new ExpandedReadBody();
            // The row's size cannot dictate the allocation: damage is only 100 columns.
            Rectangle wide = new Rectangle(0, 3, Integer.MAX_VALUE - 1, 8);
            assertTrue(body.reset(pair.actual, wide, HEIGHT));
            accept(body, pair.expected, wide, HEIGHT, 0, 100);
            accept(body, pair.expected, wide, HEIGHT, 90, 100);
            body.draw();
            Rectangle narrow = new Rectangle(0, 14, 7, 5);
            assertTrue(body.reset(pair.actual, narrow, 3));
            accept(body, pair.expected, narrow, 3, 0, 7);
            accept(body, pair.expected, narrow, 3, 6, 7);
            body.draw();
            // Reset discards unfinished work rather than replaying it in a later row.
            assertTrue(body.reset(pair.actual, ROW, HEIGHT));
            assertTrue(body.fill(0, 90));
            assertTrue(body.reset(pair.actual, new Rectangle(0, 22, 100, 8), HEIGHT));
            accept(body, pair.expected, new Rectangle(0, 22, 100, 8), HEIGHT, 12, 13);
            body.draw();
            pair.assertExact();
            pair.actual.setClip(new Rectangle(0, 0, 0, 0));
            pair.expected.setClip(new Rectangle(0, 0, 0, 0));
            // Empty damage may reject the batch or accept it as empty; neither paints.
            if (body.reset(pair.actual, ROW, HEIGHT)) assertTrue(body.fill(0, 100));
            body.draw();
            pair.assertExact();
        }
    }

    @Test
    public void hugeVerticalTranslationsRetainOriginalIntegerRectanglePaints() {
        for (int type : new int[]{BufferedImage.TYPE_INT_ARGB, BufferedImage.TYPE_INT_ARGB_PRE}) {
            for (int scale : new int[]{1, 2}) {
                for (int offset : new int[]{-25000000, 25000000}) {
                    try (Pair pair = new Pair(type, BACKGROUND, g -> {
                        g.scale(scale, scale);
                        g.translate(0, -offset);
                    })) {
                        pair.color(INK, 0.75f);
                        Rectangle row = new Rectangle(0, offset + 3, 100, 8);
                        ExpandedReadBody body = new ExpandedReadBody();
                        // Device y is small, but logical floats would lose integer edges.
                        assertFalse(body.reset(pair.actual, row, HEIGHT));
                        for (int[] span : new int[][]{{2, 20}, {7, 23}, {8, 8}, {30, 31}}) {
                            if (!body.fill(span[0], span[1])) {
                                body.draw();
                                pair.actual.fillRect(span[0], row.y, span[1] - span[0], HEIGHT);
                            }
                            pair.expected.fillRect(span[0], row.y, span[1] - span[0], HEIGHT);
                        }
                        body.draw();
                        pair.assertExact();
                    }
                }
            }
        }
    }

    @Test
    public void reversedEndpointsRejectWithoutConsumingEarlierOrLaterWork() {
        try (Pair pair = new Pair(BufferedImage.TYPE_INT_ARGB, BACKGROUND, g -> { })) {
            pair.color(INK, 0.75f);
            ExpandedReadBody body = new ExpandedReadBody();
            assertTrue(body.reset(pair.actual, ROW, HEIGHT));
            accept(body, pair.expected, ROW, HEIGHT, 2, 9);
            assertFalse(body.fill(9, 3));
            assertFalse(body.fill(Integer.MIN_VALUE, 3));
            assertFalse(body.fill(3, Integer.MAX_VALUE));
            body.draw();
            pair.actual.fillRect(9, ROW.y, 3 - 9, HEIGHT);
            pair.expected.fillRect(9, ROW.y, 3 - 9, HEIGHT);
            accept(body, pair.expected, ROW, HEIGHT, 4, 11);
            accept(body, pair.expected, ROW, HEIGHT, 11, 11);
            body.draw();
            pair.assertExact();
        }
    }

    @Test
    public void nullResetDisablesPendingScratchWithoutSuppressingOriginalPaints() {
        try (Pair pair = new Pair(BufferedImage.TYPE_INT_ARGB_PRE, BACKGROUND, g -> g.scale(2, 2))) {
            pair.color(INK, 0.75f);
            ExpandedReadBody body = new ExpandedReadBody();
            assertTrue(body.reset(pair.actual, ROW, HEIGHT));
            assertTrue(body.fill(2, 20));
            assertFalse(body.reset(null, ROW, HEIGHT));
            assertFalse(body.fill(2, 20));
            body.draw();
            pair.actual.fillRect(7, ROW.y, 13, HEIGHT);
            pair.expected.fillRect(7, ROW.y, 13, HEIGHT);
            pair.assertExact();
            assertTrue(body.reset(pair.actual, ROW, HEIGHT));
            accept(body, pair.expected, ROW, HEIGHT, 4, 11);
            body.draw();
            pair.assertExact();
        }
    }

    @Test
    public void copiedContextsHaveIndependentInterleavedPendingRectangles() {
        try (Pair pair = new Pair(BufferedImage.TYPE_INT_ARGB_PRE, BACKGROUND, g -> g.scale(2, 2))) {
            ReferenceFrame frame = new ReferenceFrame("expanded-scratch-copy");
            frame.setBounds(0, 100);
            RenderContext parent = new RenderContext(null, pair.actual, frame, new Rectangle(0, 0, 100, 40));
            RenderContext child = new RenderContext(parent);
            try {
                ExpandedReadBody first = parent.getExpandedReadBody();
                ExpandedReadBody second = child.getExpandedReadBody();
                Graphics2D firstGraphics = parent.getGraphics2D("pending-parent");
                Graphics2D secondGraphics = child.getGraphics2D("pending-child");
                firstGraphics.setColor(INK);
                secondGraphics.setColor(new Color(211, 43, 137, 179));
                firstGraphics.setComposite(AlphaComposite.getInstance(AlphaComposite.SRC_OVER, 0.37f));
                secondGraphics.setComposite(AlphaComposite.getInstance(AlphaComposite.SRC_OVER, 0.61f));
                assertTrue(first.reset(firstGraphics, ROW, HEIGHT));
                assertTrue(first.fill(2, 19));
                assertTrue(second.reset(secondGraphics, ROW, HEIGHT));
                assertTrue(second.fill(4, 16));
                assertTrue(first.fill(6, 13));
                assertTrue(second.fill(8, 21));
                first.draw();
                second.draw();
                pair.expected.setColor(INK);
                pair.expected.setComposite(firstGraphics.getComposite());
                pair.expected.fillRect(2, ROW.y, 17, HEIGHT);
                pair.expected.fillRect(6, ROW.y, 7, HEIGHT);
                pair.expected.setColor(secondGraphics.getColor());
                pair.expected.setComposite(secondGraphics.getComposite());
                pair.expected.fillRect(4, ROW.y, 12, HEIGHT);
                pair.expected.fillRect(8, ROW.y, 13, HEIGHT);
                pair.assertExact();
            } finally {
                child.dispose();
                child.getGraphics().dispose();
                parent.dispose();
            }
        }
    }

    @Test
    public void customColorContextsKeepOriginalPerRectanglePaints() {
        for (int type : new int[]{BufferedImage.TYPE_INT_ARGB, BufferedImage.TYPE_INT_ARGB_PRE}) {
            for (int scale : new int[]{1, 2}) {
                try (Pair pair = new Pair(type, BACKGROUND, g -> g.scale(scale, scale))) {
                    pair.color(INK, 0.75f);
                    pair.expected.setColor(new ContextColor());
                    pair.actual.setColor(new ContextColor());
                    ExpandedReadBody body = new ExpandedReadBody();
                    assertFalse(body.reset(pair.actual, ROW, HEIGHT));
                    for (int[] span : new int[][]{{2, 11}, {11, 20}}) {
                        if (!body.fill(span[0], span[1])) {
                            body.draw();
                            pair.actual.fillRect(span[0], ROW.y, span[1] - span[0], HEIGHT);
                        }
                        pair.expected.fillRect(span[0], ROW.y, span[1] - span[0], HEIGHT);
                    }
                    body.draw();
                    assertNotEquals(pair.expectedImage.getRGB(3 * scale, 4 * scale),
                            pair.expectedImage.getRGB(12 * scale, 4 * scale));
                    pair.assertExact();
                }
            }
        }
    }

    @Test
    public void unsupportedRasterStatesRejectAndClearPriorPendingWork() {
        List<Consumer<Graphics2D>> unsupported = List.of(
                g -> g.setRenderingHint(RenderingHints.KEY_ANTIALIASING, RenderingHints.VALUE_ANTIALIAS_ON),
                g -> g.setComposite(AlphaComposite.Src),
                g -> g.setComposite(AlphaComposite.Clear),
                g -> g.setComposite(new DelegatingComposite()),
                g -> g.setPaint(new GradientPaint(0, 0, Color.RED, 10, 0, Color.BLUE)),
                g -> g.scale(1.25, 1.25),
                g -> g.scale(1, 2),
                g -> g.translate(0.25, 0),
                g -> g.translate(0, 0.5),
                g -> g.shear(0.1, 0),
                g -> g.rotate(0.1),
                g -> g.scale(-1, 1),
                g -> g.setClip(null));
        for (Consumer<Graphics2D> configure : unsupported) {
            try (Pair pair = new Pair(BufferedImage.TYPE_INT_ARGB_PRE, BACKGROUND, g -> { })) {
                pair.color(INK, 0.75f);
                ExpandedReadBody body = new ExpandedReadBody();
                assertTrue(body.reset(pair.actual, ROW, HEIGHT));
                assertTrue(body.fill(2, 20));
                configure.accept(pair.actual);
                assertFalse(body.reset(pair.actual, ROW, HEIGHT));
                assertFalse(body.fill(2, 20));
                body.draw();
                pair.assertExact();
            }
        }
    }

    @Test
    public void vectorGraphicsRejectWithoutRasterizingOrPainting() {
        SVGGraphics2D graphics = new SVGGraphics2D(GenericDOMImplementation.getDOMImplementation()
                .createDocument("http://www.w3.org/2000/svg", "svg", null));
        try {
            graphics.setClip(new Rectangle(0, 0, 100, 40));
            graphics.setColor(INK);
            graphics.setComposite(AlphaComposite.getInstance(AlphaComposite.SRC_OVER, 0.75f));
            graphics.setRenderingHint(RenderingHints.KEY_ANTIALIASING, RenderingHints.VALUE_ANTIALIAS_OFF);
            ExpandedReadBody body = new ExpandedReadBody();
            assertFalse(body.reset(graphics, ROW, HEIGHT));
            assertFalse(body.fill(2, 20));
            body.draw();
            assertEquals(0, graphics.getRoot().getElementsByTagName("image").getLength());
            assertEquals(0, graphics.getRoot().getElementsByTagName("path").getLength());
        } finally {
            graphics.dispose();
        }
    }

    @Test
    public void emptyRowsThinHeightsAndUnsupportedDeviceEdgesDisablePendingWork() {
        for (Rectangle row : List.of(new Rectangle(0, 3, 0, 8), new Rectangle(0, 3, -1, 8),
                new Rectangle(0, 1 << 21, 100, 8))) {
            try (Pair pair = new Pair(BufferedImage.TYPE_INT_ARGB, BACKGROUND, g -> { })) {
                pair.color(INK, 0.75f);
                ExpandedReadBody body = new ExpandedReadBody();
                assertTrue(body.reset(pair.actual, ROW, HEIGHT));
                assertTrue(body.fill(2, 20));
                assertFalse(body.reset(pair.actual, row, HEIGHT));
                body.draw();
                pair.assertExact();
            }
        }
        for (int height : new int[]{-1, 0, 1}) {
            try (Pair pair = new Pair(BufferedImage.TYPE_INT_ARGB, BACKGROUND, g -> { })) {
                ExpandedReadBody body = new ExpandedReadBody();
                assertFalse(body.reset(pair.actual, ROW, height));
                assertFalse(body.fill(2, 20));
                body.draw();
                pair.assertExact();
            }
        }
    }

    private static void accept(ExpandedReadBody body, Graphics2D expected, Rectangle row, int height, int start, int end) {
        assertTrue("Supported rectangle " + start + ".." + end, body.fill(start, end));
        expected.fillRect(start, row.y, end - start, height);
    }

    static void assertExactPixels(BufferedImage expected, BufferedImage actual) {
        assertEquals(expected.getWidth(), actual.getWidth());
        assertEquals(expected.getHeight(), actual.getHeight());
        int[] expectedData = ((DataBufferInt) expected.getRaster().getDataBuffer()).getData();
        int[] actualData = ((DataBufferInt) actual.getRaster().getDataBuffer()).getData();
        for (int i = 0; i < expectedData.length; i++) {
            if (expectedData[i] != actualData[i]) {
                assertEquals("Raw pixel " + (i % expected.getWidth()) + "," + (i / expected.getWidth()), expectedData[i], actualData[i]);
            }
        }
    }

    private static final class Pair implements AutoCloseable {
        final BufferedImage expectedImage;
        final BufferedImage actualImage;
        final Graphics2D expected;
        final Graphics2D actual;

        Pair(int type, Color background, Consumer<Graphics2D> configure) {
            expectedImage = new BufferedImage(220, 90, type);
            actualImage = new BufferedImage(220, 90, type);
            expected = expectedImage.createGraphics();
            actual = actualImage.createGraphics();
            for (Graphics2D graphics : new Graphics2D[]{expected, actual}) {
                graphics.setComposite(AlphaComposite.Src);
                graphics.setColor(background);
                graphics.fillRect(0, 0, 220, 90);
                graphics.setComposite(AlphaComposite.SrcOver);
                graphics.setRenderingHint(RenderingHints.KEY_ANTIALIASING, RenderingHints.VALUE_ANTIALIAS_OFF);
                graphics.setClip(new Rectangle(0, 0, 220, 90));
                configure.accept(graphics);
            }
        }

        void color(Color color, float alpha) {
            for (Graphics2D graphics : new Graphics2D[]{expected, actual}) {
                graphics.setColor(color);
                graphics.setComposite(AlphaComposite.getInstance(AlphaComposite.SRC_OVER, alpha));
            }
        }

        void assertExact() {
            assertExactPixels(expectedImage, actualImage);
        }

        @Override
        public void close() {
            expected.dispose();
            actual.dispose();
        }
    }

    private static final class ContextColor extends Color {
        private int contexts;

        private ContextColor() {
            super(30, 80, 140, 113);
        }

        @Override
        public PaintContext createContext(ColorModel model, Rectangle bounds, Rectangle2D userBounds,
                                          AffineTransform transform, RenderingHints hints) {
            // A valid Color subclass can depend on separate paint-context creation.
            return new Color(30 + 20 * ++contexts, 80, 140, 113)
                    .createContext(model, bounds, userBounds, transform, hints);
        }
    }

    private static final class DelegatingComposite implements Composite {
        @Override
        public CompositeContext createContext(ColorModel source, ColorModel destination, RenderingHints hints) {
            return AlphaComposite.SrcOver.createContext(source, destination, hints);
        }
    }
}
