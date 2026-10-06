package org.broad.igv.feature.genome;

import org.broad.igv.feature.Chromosome;
import org.broad.igv.feature.genome.load.GenomeConfig;
import org.junit.After;
import org.junit.Before;
import org.junit.Test;

import java.lang.reflect.Field;
import java.nio.charset.StandardCharsets;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;
import java.util.concurrent.CyclicBarrier;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.Future;
import java.util.concurrent.TimeUnit;

import static org.junit.Assert.*;

/** Public sequence contracts: independent literals and a frozen pre-patch differential oracle. */
public class SequenceWrapperContractTest {

    private static final String DNA = "ABCDEFGHIJKLMNOPQRSTUVWXYZABCDEFGHIJKLMNOPQRSTUVWXYZ";
    private boolean savedCache;
    private int savedTileSize;
    private boolean savedOracleCache;
    private int savedOracleTileSize;
    private SequenceWrapper settings;
    private PrePatchSequenceWrapper oracleSettings;

    @Before
    public void saveGlobalSettings() throws Exception {
        // Reflection is restricted to saving preferences; no cache representation is inspected.
        savedCache = (boolean) setting(SequenceWrapper.class, "cacheSequences");
        savedTileSize = (int) setting(SequenceWrapper.class, "tileSize");
        savedOracleCache = (boolean) setting(PrePatchSequenceWrapper.class, "cacheSequences");
        savedOracleTileSize = (int) setting(PrePatchSequenceWrapper.class, "tileSize");
        settings = new SequenceWrapper(new MemorySequence());
        oracleSettings = new PrePatchSequenceWrapper(new MemorySequence());
        caching(true);
        tileSize(4);
    }

    @After
    public void restoreGlobalSettings() {
        if (settings != null) {
            settings.setTileSize(savedTileSize);
            SequenceWrapper.setCacheSequences(savedCache);
        }
        if (oracleSettings != null) {
            oracleSettings.setTileSize(savedOracleTileSize);
            PrePatchSequenceWrapper.setCacheSequences(savedOracleCache);
        }
    }

    @Test
    public void withinTileColdAndWarmReturnExactCopies() {
        Pair pair = pair("chr", DNA);
        pair.range("chr", 1, 3, bytes("BC"));
        pair.calls("chr:0-4");
        pair.range("chr", 2, 4, bytes("CD"));
        // End is inclusive for tile fetching, even though public coordinates are half-open.
        pair.calls("chr:0-4", "chr:4-8");
        pair.range("chr", 1, 3, bytes("BC"));
        pair.calls("chr:0-4", "chr:4-8");
    }

    @Test
    public void acrossTilesColdAndWarmCoalesceContiguousMisses() {
        Pair pair = pair("chr", DNA);
        pair.range("chr", 2, 11, bytes("CDEFGHIJK"));
        pair.calls("chr:0-12");
        pair.range("chr", 3, 9, bytes("DEFGHI"));
        pair.calls("chr:0-12");
    }

    @Test
    public void exactTileBoundariesPreserveExtraEndTileFetch() {
        Pair pair = pair("chr", DNA);
        pair.range("chr", 0, 4, bytes("ABCD"));
        pair.calls("chr:0-8");
        pair.range("chr", 4, 8, bytes("EFGH"));
        pair.calls("chr:0-8", "chr:8-12");
    }

    @Test
    public void mixedHitsSplitMissesIntoContiguousUnderlyingRanges() {
        Pair pair = pair("chr", DNA);
        pair.range("chr", 4, 5, bytes("E"));
        pair.range("chr", 16, 17, bytes("Q"));
        pair.range("chr", 1, 23, bytes("BCDEFGHIJKLMNOPQRSTUVW"));
        pair.calls("chr:4-8", "chr:16-20", "chr:0-4", "chr:8-16", "chr:20-24");
        pair.range("chr", 1, 23, bytes("BCDEFGHIJKLMNOPQRSTUVW"));
        pair.calls("chr:4-8", "chr:16-20", "chr:0-4", "chr:8-16", "chr:20-24");
    }

    @Test
    public void emptyRangesStillFetchTheirTile() {
        Pair pair = pair("chr", DNA);
        pair.range("chr", 2, 2, new byte[0]);
        pair.calls("chr:0-4");
        pair.range("chr", 4, 4, new byte[0]);
        pair.calls("chr:0-4", "chr:4-8");
    }

    @Test
    public void chromosomeLeftEdgePadsWithZeroBytes() {
        Pair pair = pair("chr", DNA);
        pair.range("chr", -2, 2, new byte[]{0, 0, 'A', 'B'});
        pair.calls("chr:0-4");
        pair.range("chr", -1, 0, new byte[]{0});
        pair.calls("chr:0-4");
    }

    @Test
    public void missingRangeLoadsZeroTilesInsteadOfReturningNull() {
        Pair pair = pair("chr", DNA);
        pair.range("missing", 1, 7, new byte[6]);
        pair.calls("missing:0-8");
        pair.range("missing", 2, 3, new byte[1]);
        pair.base("missing", 1, (byte) 0);
        pair.calls("missing:0-8");
    }

    @Test
    public void missingGetBaseCachesNullTileWhichMakesFirstRangeNull() {
        Pair pair = pair("chr", DNA);
        pair.base("missing", 1, (byte) 0);
        pair.range("missing", 1, 3, null);
        pair.range("missing", 2, 2, null);
        pair.calls("missing:0-4");
        pair.candidate.clearCache();
        pair.oracle.clearCache();
        pair.range("missing", 1, 3, new byte[2]);
        pair.calls("missing:0-4", "missing:0-4");
    }

    @Test
    public void nullLaterTileRetainsExistingErrorClass() {
        Pair pair = pair("chr", DNA);
        pair.base("missing", 5, (byte) 0);
        pair.error("missing", 1, 6, NullPointerException.class);
        pair.calls("missing:4-8", "missing:0-4");
    }

    @Test
    public void shortTerminalTilesReturnAvailableBytesAndTrailingZeros() {
        Pair pair = pair("chr", "ABCDEFGHIJ");
        pair.range("chr", 7, 10, bytes("HIJ"));
        pair.calls("chr:4-12");
        pair.range("chr", 8, 11, new byte[]{'I', 'J', 0});
        pair.base("chr", 9, (byte) 'J');
        pair.base("chr", 10, (byte) 0);
        pair.calls("chr:4-12");
    }

    @Test
    public void shortTerminalOffsetPastAvailableBytesRetainsErrorClass() {
        Pair pair = pair("chr", "ABCDEFGHIJ");
        pair.error("chr", 11, 11, ArrayIndexOutOfBoundsException.class);
        pair.calls("chr:8-12");
    }

    @Test
    public void tooManyShortTerminalTilesRetainNegativeArrayError() {
        Pair pair = pair("chr", "ABCDEFGHIJ");
        pair.error("chr", 1, 13, NegativeArraySizeException.class);
        pair.calls("chr:0-16");
    }

    @Test
    public void reversedRangesFailBeforeAnyUnderlyingCall() {
        Pair pair = pair("chr", DNA);
        pair.error("chr", 4, 3, NegativeArraySizeException.class);
        pair.calls();
    }

    @Test
    public void uncachedPathsDelegateExactRangesNullsAndBases() {
        caching(false);
        Pair pair = pair("chr", DNA);
        pair.range("chr", 1, 5, bytes("BCDE"));
        pair.range("chr", 1, 5, bytes("BCDE"));
        pair.range("chr", 2, 2, new byte[0]);
        pair.range("missing", 1, 3, null);
        pair.base("chr", 0, (byte) 'A');
        pair.base("chr", 4, (byte) 'E');
        pair.base("missing", 1, (byte) 0);
        pair.calls("chr:1-5", "chr:1-5", "chr:2-2", "missing:1-3");
        assertEquals(Arrays.asList("chr:0", "chr:4", "missing:1"), pair.backing.baseCalls);
        assertEquals(pair.backing.baseCalls, pair.oracleBacking.baseCalls);
    }

    @Test
    public void cachedGetBaseKeepsHistoricalTileStartZeroAndSharesRangeCache() {
        Pair pair = pair("chr", DNA);
        pair.base("chr", 0, (byte) 0);
        pair.base("chr", 1, (byte) 'B');
        pair.range("chr", 1, 3, bytes("BC"));
        pair.range("chr", 4, 6, bytes("EF"));
        pair.base("chr", 4, (byte) 0);
        pair.base("chr", 5, (byte) 'F');
        pair.calls("chr:0-4", "chr:4-8");
        assertTrue(pair.backing.baseCalls.isEmpty());
    }

    @Test
    public void returnedArrayMutationNeverPoisonsSingleOrMultipleCachedTiles() {
        Pair pair = pair("chr", DNA);
        byte[] one = pair.candidate.getSequence("chr", 1, 3);
        byte[] oldOne = pair.oracle.getSequence("chr", 1, 3);
        one[0] = oldOne[0] = '!';
        pair.range("chr", 1, 3, bytes("BC"));
        pair.base("chr", 1, (byte) 'B');
        byte[] many = pair.candidate.getSequence("chr", 2, 9);
        byte[] oldMany = pair.oracle.getSequence("chr", 2, 9);
        Arrays.fill(many, (byte) '!');
        Arrays.fill(oldMany, (byte) '!');
        pair.range("chr", 2, 9, bytes("CDEFGHI"));
        assertNotSame(one, pair.candidate.getSequence("chr", 1, 3));
        // Keep differential call accounting aligned after the identity assertion.
        pair.oracle.getSequence("chr", 1, 3);
        pair.calls("chr:0-4", "chr:4-12");
    }

    @Test
    public void distinctChromosomeNamesNeverShareTilesIncludingHashCollisions() {
        MemorySequence backing = new MemorySequence().add("Aa", "ABCD").add("BB", "WXYZ")
                .add("chr1", "EFGH").add("1", "IJKL").add("chr11", DNA);
        Pair pair = new Pair(backing);
        pair.range(new String("Aa"), 1, 3, bytes("BC"));
        pair.range(new String("BB"), 1, 3, bytes("XY"));
        pair.range(new String("Aa"), 1, 3, bytes("BC"));
        pair.range("chr1", 1, 3, bytes("FG"));
        pair.range("1", 1, 3, bytes("JK"));
        pair.range("chr11", 9, 11, bytes("JK"));
        pair.calls("Aa:0-4", "BB:0-4", "chr1:0-4", "1:0-4", "chr11:8-12");
    }

    @Test
    public void wrappersWithSameChromosomeRemainIsolated() {
        Pair first = pair("chr", "ABCDEFGH");
        Pair second = pair("chr", "WXYZQRST");
        first.range("chr", 1, 3, bytes("BC"));
        second.range("chr", 1, 3, bytes("XY"));
        first.candidate.clearCache();
        first.oracle.clearCache();
        second.range("chr", 1, 3, bytes("XY"));
        first.range("chr", 1, 3, bytes("BC"));
        first.calls("chr:0-4", "chr:0-4");
        second.calls("chr:0-4");
    }

    @Test
    public void fiftyEntryLimitSpansChromosomesAndEvictsInsertionOrderNotLru() {
        MemorySequence backing = new MemorySequence();
        for (int i = 0; i <= 50; i++) backing.add("c" + i, "ABCD");
        Pair pair = new Pair(backing);
        for (int i = 0; i < 50; i++) pair.base("c" + i, 1, (byte) 'B');
        pair.range("c0", 1, 3, bytes("BC")); // A hit must not refresh insertion order.
        pair.base("c50", 1, (byte) 'B');
        assertEquals(51, pair.backing.rangeCalls.size());
        pair.range("c1", 1, 3, bytes("BC"));
        assertEquals(51, pair.backing.rangeCalls.size());
        pair.range("c0", 1, 3, bytes("BC"));
        assertEquals(52, pair.backing.rangeCalls.size());
        assertEquals("c0:0-4", pair.backing.rangeCalls.get(51));
        pair.base("c1", 1, (byte) 'B');
        assertEquals(53, pair.backing.rangeCalls.size());
        assertEquals("c1:0-4", pair.backing.rangeCalls.get(52));
        pair.sameCalls();
    }

    @Test
    public void getBaseAndRangeInsertionsShareFiftyEntryEvictionOrder() {
        Pair pair = pair("chr", String.join("", Collections.nCopies(16, "ABCDEFGHIJKLMNOP")));
        for (int tile = 0; tile < 50; tile++) {
            pair.base("chr", tile * 4 + 1, (byte) "BFJN".charAt(tile % 4));
        }
        pair.base("chr", 1, (byte) 'B');
        pair.range("chr", 200, 201, bytes("I"));
        assertEquals(51, pair.backing.rangeCalls.size());
        pair.range("chr", 1, 3, bytes("BC"));
        assertEquals(52, pair.backing.rangeCalls.size());
        assertEquals("chr:0-4", pair.backing.rangeCalls.get(51));
        pair.sameCalls();
    }

    @Test
    public void explicitClearForcesFreshUnderlyingLoads() {
        Pair pair = pair("chr", DNA);
        pair.range("chr", 1, 3, bytes("BC"));
        pair.candidate.clearCache();
        pair.oracle.clearCache();
        pair.range("chr", 1, 3, bytes("BC"));
        pair.calls("chr:0-4", "chr:0-4");
    }

    @Test
    public void changingTileSizeInvalidatesTheChangedWrapper() {
        Pair pair = pair("chr", DNA);
        pair.range("chr", 1, 3, bytes("BC"));
        pair.candidate.setTileSize(8);
        pair.oracle.setTileSize(8);
        pair.range("chr", 5, 7, bytes("FG"));
        pair.calls("chr:0-4", "chr:0-8");
        pair.candidate.setTileSize(8);
        pair.oracle.setTileSize(8);
        pair.range("chr", 1, 3, bytes("BC"));
        pair.calls("chr:0-4", "chr:0-8");
    }

    @Test
    public void globalTileSizeChangesInvalidateOtherWrappersLazily() {
        // The frozen oracle has a known stale-cache bug here; assert the corrected contract literally.
        MemorySequence firstBacking = new MemorySequence().add("chr", DNA);
        MemorySequence otherBacking = new MemorySequence().add("chr", DNA);
        SequenceWrapper first = new SequenceWrapper(firstBacking);
        SequenceWrapper other = new SequenceWrapper(otherBacking);
        assertArrayEquals(bytes("FG"), other.getSequence("chr", 5, 7));
        first.setTileSize(8);
        assertArrayEquals(bytes("JK"), other.getSequence("chr", 9, 11));
        assertEquals(Arrays.asList("chr:4-8", "chr:8-16"), otherBacking.rangeCalls);
        assertArrayEquals(bytes("BC"), first.getSequence("chr", 1, 3));
        assertEquals(Collections.singletonList("chr:0-8"), firstBacking.rangeCalls);
    }

    @Test
    public void tileSizeChangeWhileCachingDisabledInvalidatesOnReenable() {
        // Do not perpetuate the pre-patch stale-cache behavior when global preferences change.
        MemorySequence backing = new MemorySequence().add("chr", DNA);
        SequenceWrapper wrapper = new SequenceWrapper(backing);
        assertArrayEquals(bytes("FG"), wrapper.getSequence("chr", 5, 7));
        SequenceWrapper.setCacheSequences(false);
        wrapper.setTileSize(8);
        assertArrayEquals(bytes("FG"), wrapper.getSequence("chr", 5, 7));
        SequenceWrapper.setCacheSequences(true);
        assertEquals((byte) 'J', wrapper.getBase("chr", 9));
        assertArrayEquals(bytes("JK"), wrapper.getSequence("chr", 9, 11));
        assertEquals(Arrays.asList("chr:4-8", "chr:5-7", "chr:8-16"), backing.rangeCalls);
    }

    @Test
    public void togglingCacheWithoutTileChangeRetainsExistingCachedTiles() {
        Pair pair = pair("chr", DNA);
        pair.range("chr", 1, 3, bytes("BC"));
        caching(false);
        pair.range("chr", 1, 3, bytes("BC"));
        caching(true);
        pair.range("chr", 1, 3, bytes("BC"));
        pair.calls("chr:0-4", "chr:1-3");
    }

    @Test
    public void genomeCanonicalAliasesClampEndsAndRejectEmptyOrUnknownRanges() throws Exception {
        MemorySequence backing = new MemorySequence().add("chr1", "ABCDEFGHIJ");
        MemorySequence oldBacking = backing.copy();
        Genome genome = genome("contract", backing);
        PrePatchSequenceWrapper oracle = new PrePatchSequenceWrapper(oldBacking);
        byte[] expected = bytes("BCDEFGHIJ");
        assertArrayEquals(expected, genome.getSequence("one", 1, 99));
        assertArrayEquals(oracle.getSequence("chr1", 1, 10), genome.getSequence("chr1", 1, 10));
        assertNull(genome.getSequence("unknown", 1, 3));
        assertNull(genome.getSequence("one", 2, 2));
        assertNull(genome.getSequence("one", 11, 99));
        assertEquals(Collections.singletonList("chr1:0-12"), backing.rangeCalls);
        assertEquals(oldBacking.rangeCalls, backing.rangeCalls);
        assertEquals((byte) 0, genome.getReference("chr1", 0));
        assertEquals((byte) 'B', genome.getReference("chr1", 1));
        assertEquals(oracle.getBase("chr1", 1), genome.getReference("chr1", 1));
    }

    @Test
    public void genomeWrappersStayIsolatedAndUseCacheFlagPreservesCurrentBehavior() throws Exception {
        MemorySequence firstBacking = new MemorySequence().add("chr1", "ABCDEFGH");
        MemorySequence secondBacking = new MemorySequence().add("chr1", "WXYZQRST");
        Genome first = genome("first", firstBacking);
        Genome second = genome("second", secondBacking);
        PrePatchSequenceWrapper firstOracle = new PrePatchSequenceWrapper(firstBacking.copy());
        PrePatchSequenceWrapper secondOracle = new PrePatchSequenceWrapper(secondBacking.copy());
        assertArrayEquals(bytes("BC"), first.getSequence("one", 1, 3));
        assertArrayEquals(bytes("XY"), second.getSequence("one", 1, 3));
        assertArrayEquals(firstOracle.getSequence("chr1", 1, 3), first.getSequence("chr1", 1, 3, false));
        assertArrayEquals(secondOracle.getSequence("chr1", 1, 3), second.getSequence("chr1", 1, 3, false));
        assertEquals(Collections.singletonList("chr1:0-4"), firstBacking.rangeCalls);
        assertEquals(Collections.singletonList("chr1:0-4"), secondBacking.rangeCalls);
    }


    @Test(timeout = 10000)
    public void concurrentPublicReadsAndClearHaveDeterministicResults() throws Exception {
        SequenceWrapper wrapper = new SequenceWrapper(new MemorySequence().add("chr", DNA));
        ExecutorService workers = Executors.newFixedThreadPool(3);
        CyclicBarrier start = new CyclicBarrier(3);
        try {
            Future<?> ranges = workers.submit(() -> {
                start.await(5, TimeUnit.SECONDS);
                for (int i = 0; i < 100; i++) {
                    assertArrayEquals(bytes("CDEFGHI"), wrapper.getSequence("chr", 2, 9));
                    byte[] copy = wrapper.getSequence("chr", 1, 3);
                    copy[0] = '!';
                }
                return null;
            });
            Future<?> bases = workers.submit(() -> {
                start.await(5, TimeUnit.SECONDS);
                for (int i = 0; i < 100; i++) {
                    assertEquals((byte) 'B', wrapper.getBase("chr", 1));
                    assertEquals((byte) 0, wrapper.getBase("chr", 4));
                }
                return null;
            });
            Future<?> clears = workers.submit(() -> {
                start.await(5, TimeUnit.SECONDS);
                for (int i = 0; i < 100; i++) wrapper.clearCache();
                return null;
            });
            ranges.get(5, TimeUnit.SECONDS);
            bases.get(5, TimeUnit.SECONDS);
            clears.get(5, TimeUnit.SECONDS);
        } finally {
            workers.shutdownNow();
            assertTrue("Workers must finish before restoring global preferences", workers.awaitTermination(5, TimeUnit.SECONDS));
        }
        assertArrayEquals(bytes("BC"), wrapper.getSequence("chr", 1, 3));
    }

    private static Object setting(Class<?> type, String name) throws Exception {
        Field field = type.getDeclaredField(name);
        field.setAccessible(true);
        return field.get(null);
    }

    private void caching(boolean enabled) {
        SequenceWrapper.setCacheSequences(enabled);
        PrePatchSequenceWrapper.setCacheSequences(enabled);
    }

    private void tileSize(int size) {
        settings.setTileSize(size);
        oracleSettings.setTileSize(size);
    }

    private static byte[] bytes(String sequence) {
        return sequence.getBytes(StandardCharsets.US_ASCII);
    }

    private static Pair pair(String chromosome, String sequence) {
        return new Pair(new MemorySequence().add(chromosome, sequence));
    }

    private static Genome genome(String id, MemorySequence sequence) throws Exception {
        GenomeConfig config = new GenomeConfig();
        config.setId(id);
        config.setName(id);
        config.setSequence(sequence);
        config.setChromAliases(Collections.singletonList(Arrays.asList("chr1", "one")));
        return new Genome(config);
    }

    private static final class Pair {
        final MemorySequence backing;
        final MemorySequence oracleBacking;
        final SequenceWrapper candidate;
        final PrePatchSequenceWrapper oracle;

        Pair(MemorySequence backing) {
            this.backing = backing;
            this.oracleBacking = backing.copy();
            this.candidate = new SequenceWrapper(backing);
            this.oracle = new PrePatchSequenceWrapper(oracleBacking);
        }

        void range(String chr, int start, int end, byte[] literal) {
            byte[] old = oracle.getSequence(chr, start, end);
            byte[] actual = candidate.getSequence(chr, start, end);
            assertArrayEquals("Independent oracle bytes", old, actual);
            assertArrayEquals("Literal range " + chr + ":" + start + "-" + end, literal, actual);
            sameCalls();
        }

        void base(String chr, int position, byte literal) {
            byte old = oracle.getBase(chr, position);
            byte actual = candidate.getBase(chr, position);
            assertEquals("Independent oracle base", old, actual);
            assertEquals("Literal base " + chr + ":" + position, literal, actual);
            sameCalls();
        }

        void error(String chr, int start, int end, Class<? extends RuntimeException> literalClass) {
            Class<?> old = errorClass(oracle, chr, start, end);
            Class<?> actual = errorClass(candidate, chr, start, end);
            assertEquals("Independent oracle error", old, actual);
            assertEquals("Historical error class", literalClass, actual);
            sameCalls();
        }

        void calls(String... expected) {
            assertEquals("Literal underlying ranges", Arrays.asList(expected), backing.rangeCalls);
            sameCalls();
        }

        void sameCalls() {
            assertEquals("Independent underlying ranges", oracleBacking.rangeCalls, backing.rangeCalls);
            assertEquals("Independent underlying base calls", oracleBacking.baseCalls, backing.baseCalls);
        }
    }

    private static Class<?> errorClass(Sequence sequence, String chr, int start, int end) {
        try {
            sequence.getSequence(chr, start, end);
            fail("Expected range failure for " + chr + ":" + start + "-" + end);
            return null;
        } catch (RuntimeException error) {
            return error.getClass();
        }
    }

    /** Private, deterministic backing with observable public calls and no filesystem dependencies. */
    private static final class MemorySequence implements Sequence {
        private final Map<String, byte[]> chromosomes = new LinkedHashMap<>();
        final List<String> rangeCalls = Collections.synchronizedList(new ArrayList<>());
        final List<String> baseCalls = Collections.synchronizedList(new ArrayList<>());

        MemorySequence add(String chromosome, String sequence) {
            chromosomes.put(chromosome, bytes(sequence));
            return this;
        }

        MemorySequence copy() {
            MemorySequence copy = new MemorySequence();
            chromosomes.forEach((name, data) -> copy.chromosomes.put(name, data.clone()));
            return copy;
        }

        @Override
        public byte[] getSequence(String chr, int start, int end) {
            rangeCalls.add(chr + ":" + start + "-" + end);
            byte[] data = chromosomes.get(chr);
            if (data == null) return null;
            int left = Math.max(0, start);
            int right = Math.min(data.length, end);
            byte[] result = new byte[right - left];
            System.arraycopy(data, left, result, 0, result.length);
            return result;
        }

        @Override
        public byte getBase(String chr, int position) {
            baseCalls.add(chr + ":" + position);
            byte[] data = chromosomes.get(chr);
            return data == null || position < 0 || position >= data.length ? 0 : data[position];
        }

        @Override
        public List<String> getChromosomeNames() {
            return new ArrayList<>(chromosomes.keySet());
        }

        @Override
        public int getChromosomeLength(String chr) {
            byte[] data = chromosomes.get(chr);
            return data == null ? -1 : data.length;
        }

        @Override
        public List<Chromosome> getChromosomes() {
            List<Chromosome> result = new ArrayList<>();
            chromosomes.forEach((name, data) -> result.add(new Chromosome(result.size(), name, data.length)));
            return result;
        }
    }
}
