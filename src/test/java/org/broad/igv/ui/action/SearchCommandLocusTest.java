package org.broad.igv.ui.action;

import org.broad.igv.feature.Chromosome;
import org.broad.igv.feature.genome.Genome;
import org.junit.Test;

import java.util.Arrays;
import java.util.List;

import static org.junit.Assert.*;

public class SearchCommandLocusTest {

    private final Genome genome = new Genome("underscore-locus-test", Arrays.asList(
            new Chromosome(0, "chr1", 100000),
            new Chromosome(1, "chrUn_KI270442v1", 100000),
            new Chromosome(2, "chr1_100_200", 50000)));

    @Test
    public void parsesUnderscoreCoordinates() {
        SearchCommand.SearchResult result = search("chr1_16682_16774");
        assertEquals(SearchCommand.ResultType.LOCUS, result.getType());
        assertEquals("chr1", result.getChr());
        assertEquals(16681, result.getStart());
        assertEquals(16774, result.getEnd());
    }

    @Test
    public void ignoresSequenceSuffix() {
        SearchCommand.SearchResult result = search("chr1_16682_16774_TGGTGGGGG");
        assertEquals(SearchCommand.ResultType.LOCUS, result.getType());
        assertEquals("chr1", result.getChr());
        assertEquals(16681, result.getStart());
        assertEquals(16774, result.getEnd());
    }

    @Test
    public void preservesUnderscoresInChromosomeNames() {
        SearchCommand.SearchResult result = search("chrUn_KI270442v1_100_200_acgtn");
        assertEquals("chrUn_KI270442v1", result.getChr());
        assertEquals(99, result.getStart());
        assertEquals(200, result.getEnd());
    }

    @Test
    public void exactChromosomeNameTakesPrecedence() {
        SearchCommand.SearchResult result = search("chr1_100_200");
        assertEquals("chr1_100_200", result.getChr());
        assertEquals(0, result.getStart());
        assertEquals(50000, result.getEnd());
    }

    @Test
    public void acceptsAliasesInMultipleLocusSearch() {
        genome.addChrAliases(Arrays.asList(Arrays.asList("chr1", "1")));
        SearchCommand command = new SearchCommand(null, "", genome);
        List<SearchCommand.SearchResult> results = command.runSearch("1_16682_16774 chr1:20000-20100");
        assertEquals(2, results.size());
        assertEquals("chr1", results.get(0).getChr());
        assertEquals(16681, results.get(0).getStart());
        assertEquals(16774, results.get(0).getEnd());
        assertEquals("chr1", results.get(1).getChr());
        assertEquals(19999, results.get(1).getStart());
        assertEquals(20100, results.get(1).getEnd());
    }

    @Test
    public void usesExistingShortIntervalExpansion() {
        SearchCommand.SearchResult canonical = search("chr1:1-2");
        SearchCommand.SearchResult result = search("chr1_1_2");
        assertEquals(canonical.getChr(), result.getChr());
        assertEquals(canonical.getStart(), result.getStart());
        assertEquals(canonical.getEnd(), result.getEnd());
    }

    private SearchCommand.SearchResult search(String input) {
        SearchCommand command = new SearchCommand(null, input, genome);
        List<SearchCommand.SearchResult> results = command.runSearch(input);
        assertEquals(input, 1, results.size());
        return results.get(0);
    }
}
