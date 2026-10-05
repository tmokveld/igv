package org.broad.igv.sam.mods;

import htsjdk.samtools.util.SequenceUtil;

import java.util.HashMap;
import java.util.Map;

public class BaseModificationSet {

    private final char base;
    private final char strand;
    private final String modification;
    private final char canonicalBase;

    // Zero marks an empty slot; read positions are stored as position + 1.
    // Separate keys keep absent calls distinct from likelihoods 0 and unsigned 255.
    // The parser pre-counts calls, keeping occupancy at most 75% without resizing.
    private final int[] positions;
    private final byte[] likelihoods;

    BaseModificationSet(char base, char strand, String modification, int expectedCalls) {
        this.base = base;
        this.modification = modification;
        this.strand = strand;
        int capacity = 8;
        while (expectedCalls > capacity - capacity / 4) {
            capacity <<= 1;
        }
        this.positions = new int[capacity];
        this.likelihoods = new byte[capacity];
        this.canonicalBase = strand == '+' ? base : (char) SequenceUtil.complement((byte) base);
    }

    public char getBase() {
        return base;
    }

    public char getCanonicalBase() {
        return canonicalBase;
    }

    public String getModification() {
        return modification;
    }

    public char getStrand() {
        return strand;
    }

    /**
     * Returns the unsigned likelihood (0..255), or -1 if the read position has no call.
     */
    public int getLikelihood(int pos) {
        if (pos < 0) return -1;
        int key = pos + 1;
        int slot = slot(key, positions.length - 1);
        while (positions[slot] != 0) {
            if (positions[slot] == key) {
                return Byte.toUnsignedInt(likelihoods[slot]);
            }
            slot = (slot + 1) & (positions.length - 1);
        }
        return -1;
    }

    public boolean containsPosition(int pos) {
        return getLikelihood(pos) >= 0;
    }

    // Only the MM parser writes calls; consumers cannot mutate the likelihood storage.
    void putLikelihood(int pos, byte likelihood) {
        int key = pos + 1;
        int slot = slot(key, positions.length - 1);
        while (positions[slot] != 0 && positions[slot] != key) {
            slot = (slot + 1) & (positions.length - 1);
        }
        positions[slot] = key;
        likelihoods[slot] = likelihood;
    }

    private static int slot(int key, int mask) {
        int hash = key * 0x9E3779B9;
        return (hash ^ (hash >>> 16)) & mask;
    }

    /**
     * Return a descriptive string for the modification at the given position of the read sequence*
     * @param pos - position in the read sequence  (left to right, as recorded in BAM record, not 5'->3')
     * @return
     */
    public String valueString(int pos) {
        int likelihood = getLikelihood(pos);
        if (likelihood < 0) throw new NullPointerException("No modification call at read position " + pos);
        int l = (int) (100.0 * likelihood / 255);
        return "Base modification: " +
                ((codeValues.containsKey(modification)) ? codeValues.get(modification) : modification) +  " (" + l + "%)";
    }

    static Map<String, String> codeValues;

    static {
        codeValues = new HashMap<>();
        codeValues.put("m", "5mC");
        codeValues.put("h", "5hmC");
        codeValues.put("f", "5fC");
        codeValues.put("c", "5caC");
        codeValues.put("g", "5hmU");
        codeValues.put("e", "5fU");
        codeValues.put("b", "5caU");
        codeValues.put("a", "6mA");
        codeValues.put("o", "8xoG");
        codeValues.put("n", "Xao");
        codeValues.put("C", "Unknown C");
        codeValues.put("T", "Unknown T");
        codeValues.put("A", "Unknown A");
        codeValues.put("G", "Unknown G");
        codeValues.put("N", "Unknown");

    }
}
