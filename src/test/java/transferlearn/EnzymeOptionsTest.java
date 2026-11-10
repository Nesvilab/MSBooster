/*
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or
 * implied. See the License for the specific language governing
 * permissions and limitations under the License.
 */

package transferlearn;

import static org.junit.jupiter.api.Assertions.assertArrayEquals;
import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertFalse;
import static org.junit.jupiter.api.Assertions.assertThrows;
import static org.junit.jupiter.api.Assertions.assertTrue;

import java.io.File;
import java.nio.file.Path;
import java.sql.Connection;
import java.sql.DriverManager;
import java.sql.Statement;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;

import org.junit.jupiter.api.Test;
import org.junit.jupiter.api.io.TempDir;
import testutils.ParquetFooters;

// FragDIA reads the enzyme out of a library's Parquet footer and takes no digestion flags of its
// own, so every library builder takes the same options under the same rules and writes the same
// keys. What a builder lets through here is what FragDIA digests with - and a builder that quietly
// fixed up a value would record an enzyme nobody asked for.
public class EnzymeOptionsTest {

    /** The options a command line of flag/value pairs gives, as a predictor's argument loop takes them. */
    private static EnzymeOptions given(String... flagValuePairs) {
        EnzymeOptions options = new EnzymeOptions();
        for (int i = 0; i < flagValuePairs.length; i += 2) {
            assertTrue(options.accept(flagValuePairs[i], flagValuePairs[i + 1]), flagValuePairs[i]);
        }
        return options;
    }

    private static final String[] TRYPSIN = {"--enzyme-name", "stricttrypsin", "--cut", "KR",
            "--no-cut", "@", "--sense", "C", "--num-enzyme-termini", "2", "--clip-n-term-m", "1"};

    private static String[] with(String[] base, String... more) {
        List<String> all = new ArrayList<>(Arrays.asList(base));
        all.addAll(Arrays.asList(more));
        return all.toArray(new String[0]);
    }

    private static String rejection(String... flagValuePairs) {
        return assertThrows(IllegalArgumentException.class, () -> given(flagValuePairs).validate(),
                "accepted " + Arrays.toString(flagValuePairs)).getMessage();
    }

    // --- what is taken, and what is left to the predictor ---------------------------------------

    @Test
    public void onlyTheEnzymeFlagsAreTaken() {
        EnzymeOptions options = new EnzymeOptions();
        for (String flag : EnzymeOptions.FLAGS) {
            assertTrue(options.accept(flag, "x"), flag);
        }
        for (String flag : new String[] {"--fasta", "--basename", "--model", "--enzyme", "cut"}) {
            assertFalse(new EnzymeOptions().accept(flag, "x"), flag);
        }
    }

    // --- taken off a predictor's command line ------------------------------------------------
    //
    // Both predictors take the enzyme options off their command line with extract before their own
    // loops run, and those loops parse what is left exactly as they did before there were any.

    // Each option is taken with its value, wherever it sits among the predictor's own arguments,
    // and those are handed back untouched and in their order.
    @Test
    public void extractTakesEachOptionWithItsValueAndLeavesTheRest() {
        EnzymeOptions enzyme = new EnzymeOptions();
        String[] rest = EnzymeOptions.extract(new String[] {"--paramsList", "msbooster_params.txt",
                "--enzyme-name", "stricttrypsin", "--model", "model.zip", "--cut", "KR", "--no-cut", "@",
                "--fasta", "C:\\fasta\\human.fasta", "--sense", "C", "--num-enzyme-termini", "2",
                "--output-format", "parquet", "--clip-n-term-m", "1"}, enzyme);
        enzyme.validate();

        assertEquals(Arrays.asList(TRYPSIN), enzyme.builderArgs());
        assertEquals(Arrays.asList("--paramsList", "msbooster_params.txt", "--model", "model.zip",
                "--fasta", "C:\\fasta\\human.fasta", "--output-format", "parquet"), Arrays.asList(rest));
    }

    // A repeated --cut, --no-cut, --sense and --enzyme-name pair up by the order each was given in,
    // however the flags interleave with one another and with the predictor's own.
    @Test
    public void extractPairsRepeatedOptionsByOrder() {
        EnzymeOptions enzyme = new EnzymeOptions();
        String[] rest = EnzymeOptions.extract(new String[] {"--cut", "KR", "--enzyme-name", "trypsin",
                "--basename", "lib", "--sense", "C", "--cut", "D", "--no-cut", "P", "--enzyme-name", "aspn",
                "--no-cut", "@", "--min-charge", "1", "--sense", "N", "--num-enzyme-termini", "1",
                "--clip-n-term-m", "0"}, enzyme);
        enzyme.validate();

        Map<String, String> expected = new LinkedHashMap<>();
        expected.put("search_enzyme_name_1", "trypsin");
        expected.put("search_enzyme_cut_1", "KR");
        expected.put("search_enzyme_nocut_1", "P");
        expected.put("search_enzyme_sense_1", "C");
        expected.put("search_enzyme_name_2", "aspn");
        expected.put("search_enzyme_cut_2", "D");
        expected.put("search_enzyme_nocut_2", "@");
        expected.put("search_enzyme_sense_2", "N");
        expected.put("num_enzyme_termini", "1");
        expected.put("clip_nTerm_M", "0");
        assertEquals(new ArrayList<>(expected.entrySet()), new ArrayList<>(enzyme.footer("").entrySet()));
        assertEquals(Arrays.asList("--basename", "lib", "--min-charge", "1"), Arrays.asList(rest));
    }

    // A value is a value whatever it looks like. Predictor's loop visits values as well as flags, so
    // a value read as a flag there would take the argument after it for its own value: here an
    // --enzyme-name of "--cut" is the name and adds no enzyme, one of "--paramsList" does not reach
    // the predictor as its flag, and a predictor option whose value looks like an enzyme flag keeps it.
    @Test
    public void extractNeverReadsAValueAsAFlag() {
        EnzymeOptions enzyme = new EnzymeOptions();
        String[] rest = EnzymeOptions.extract(new String[] {"--enzyme-name", "--cut", "--cut", "KR",
                "--basename", "--sense", "--no-cut", "@", "--sense", "C", "--model", "--enzyme-name",
                "--num-enzyme-termini", "2", "--clip-n-term-m", "1"}, enzyme);
        enzyme.validate();

        Map<String, String> footer = enzyme.footer("");
        assertEquals("--cut", footer.get("search_enzyme_name_1"));
        assertEquals("KR", footer.get("search_enzyme_cut_1"));
        assertEquals("C", footer.get("search_enzyme_sense_1"));
        assertFalse(footer.containsKey("search_enzyme_cut_2"), "the name was taken for a second --cut");
        assertEquals(Arrays.asList("--basename", "--sense", "--model", "--enzyme-name"), Arrays.asList(rest));

        EnzymeOptions named = new EnzymeOptions();
        assertArrayEquals(new String[] {"--fasta", "x.fasta"}, EnzymeOptions.extract(
                new String[] {"--enzyme-name", "--paramsList", "--fasta", "x.fasta"}, named));
        String message = assertThrows(IllegalArgumentException.class, named::validate).getMessage();
        assertTrue(message.contains("--enzyme-name"), "the value was not taken as the name: " + message);
    }

    // A command line with no enzyme options is handed back whole and in order - each predictor's own,
    // with every argument it takes - so its loop parses exactly what it parsed before.
    @Test
    public void extractLeavesACommandLineWithoutAnEnzymeWhole() {
        String[] server = {"--paramsList", "msbooster_params.txt", "--url", "http://localhost:8001",
                "--api-key", "key", "--model", "model.zip", "--peptide-list-to-predict", "peptides.csv",
                "--custom-mods", "K,114.0429", "--ms2", "true", "--rt", "true", "--im", "false",
                "--basename", "lib", "--min-charge", "1", "--max-charge", "4", "--fasta", "human.fasta",
                "--keep-decoys", "0", "--decoy-prefix", "rev_", "--output-format", "librarytsv",
                "--output-dir", "out"};
        String[] local = {"--paramsList", "msbooster_params.txt", "--model", "FragCast.zip",
                "--peptide-list-to-predict", "peptides.parquet", "--fasta", "", "--output-dir", "out",
                "--basename", "lib", "--output-format", "parquet", "--keep-decoys", "1",
                "--decoy-prefix", "DECOY_", "--min-charge", "2", "--max-charge", "3",
                "--custom-mods", "K,114.0429"};
        for (String[] args : new String[][] {server, local}) {
            EnzymeOptions enzyme = new EnzymeOptions();
            assertArrayEquals(args, EnzymeOptions.extract(args, enzyme));
            enzyme.validate();
            assertFalse(enzyme.isSet());
        }
    }

    // A flag left without a value - which both predictors refuse before this is reached - is not
    // taken as an option with an empty value: it is handed back for the predictor to deal with.
    @Test
    public void extractLeavesALastFlagWithoutAValue() {
        EnzymeOptions enzyme = new EnzymeOptions();
        assertArrayEquals(new String[] {"--paramsList", "p.txt", "--cut"},
                EnzymeOptions.extract(new String[] {"--paramsList", "p.txt", "--cut"}, enzyme));
        assertFalse(enzyme.isSet());
    }

    // --- no enzyme ---------------------------------------------------------------------------

    // A builder given no enzyme options writes no enzyme keys, and FragDIA then falls back to its
    // default rule and says so. It is not an error to give none.
    @Test
    public void noEnzymeIsValidAndRecordsNothing() {
        EnzymeOptions options = given();
        options.validate();
        assertFalse(options.isSet());
        assertTrue(options.builderArgs().isEmpty());
        assertTrue(options.footer("").isEmpty());
    }

    // The FASTA is independent of the enzyme: given alone, it is the only key.
    @Test
    public void aFastaAloneIsRecordedAsTheDatabaseNameOnly() {
        EnzymeOptions options = given();
        options.validate();
        Map<String, String> expected = new LinkedHashMap<>();
        expected.put("database_name", "C:\\fasta\\2024-01-01-decoys-human.fasta");
        assertEquals(expected, options.footer("C:\\fasta\\2024-01-01-decoys-human.fasta"));
    }

    // There are no defaults to fill in: an option that describes an enzyme no --cut names is a
    // command line that means something other than what it says.
    @Test
    public void everyOtherOptionNeedsACut() {
        for (String[] pair : new String[][] {{"--enzyme-name", "trypsin"}, {"--no-cut", "P"},
                {"--sense", "C"}, {"--num-enzyme-termini", "2"}, {"--clip-n-term-m", "1"}}) {
            String message = rejection(pair);
            assertTrue(message.contains(pair[0]) && message.contains("--cut"), message);
        }
    }

    // --- one enzyme, two enzymes -------------------------------------------------------------

    @Test
    public void oneEnzymeIsForwardedAndRecordedExactlyAsGiven() {
        EnzymeOptions options = given(TRYPSIN);
        options.validate();
        assertTrue(options.isSet());
        assertEquals(Arrays.asList(TRYPSIN), options.builderArgs());

        Map<String, String> expected = new LinkedHashMap<>();
        expected.put("search_enzyme_name_1", "stricttrypsin");
        expected.put("search_enzyme_cut_1", "KR");
        expected.put("search_enzyme_nocut_1", "@");
        expected.put("search_enzyme_sense_1", "C");
        expected.put("num_enzyme_termini", "2");
        expected.put("clip_nTerm_M", "1");
        expected.put("database_name", "/data/human.fasta");
        assertEquals(new ArrayList<>(expected.entrySet()),
                new ArrayList<>(options.footer("/data/human.fasta").entrySet()),
                "the keys, their values and nothing else");
    }

    // Each enzyme's options pair up by the order they were given in, whatever order the flags
    // themselves came in; the second enzyme's keys exist only because a second --cut was given.
    @Test
    public void twoEnzymesPairTheirOptionsByOrder() {
        EnzymeOptions options = given("--cut", "KR", "--cut", "D", "--sense", "C", "--sense", "N",
                "--no-cut", "P", "--no-cut", "@", "--clip-n-term-m", "0", "--num-enzyme-termini", "1");
        options.validate();

        assertEquals(Arrays.asList("--cut", "KR", "--no-cut", "P", "--sense", "C",
                "--cut", "D", "--no-cut", "@", "--sense", "N",
                "--num-enzyme-termini", "1", "--clip-n-term-m", "0"), options.builderArgs());

        Map<String, String> expected = new LinkedHashMap<>();
        expected.put("search_enzyme_cut_1", "KR");
        expected.put("search_enzyme_nocut_1", "P");
        expected.put("search_enzyme_sense_1", "C");
        expected.put("search_enzyme_cut_2", "D");
        expected.put("search_enzyme_nocut_2", "@");
        expected.put("search_enzyme_sense_2", "N");
        expected.put("num_enzyme_termini", "1");
        expected.put("clip_nTerm_M", "0");
        assertEquals(new ArrayList<>(expected.entrySet()), new ArrayList<>(options.footer("").entrySet()),
                "no names were given, so no name key; no FASTA, so no database_name");
    }

    // "@" is a value in its own right: a --cut of @ cuts nothing (nonspecific), and is recorded as
    // such rather than dropped or turned into a default.
    @Test
    public void aNonspecificEnzymeIsRecordedAsGiven() {
        EnzymeOptions options = given("--enzyme-name", "nonspecific", "--cut", "@", "--no-cut", "@",
                "--sense", "C", "--num-enzyme-termini", "0", "--clip-n-term-m", "0");
        options.validate();
        assertEquals("@", options.footer("").get("search_enzyme_cut_1"));
        assertEquals("0", options.footer("").get("num_enzyme_termini"));
    }

    // A name is free text, recorded exactly: FragDIA uses it in its log line and nowhere else.
    @Test
    public void aNameIsRecordedExactlyAsWritten() {
        EnzymeOptions options = given(with(new String[] {"--enzyme-name", "Trypsin/P", "--cut", "KR",
                "--no-cut", "@", "--sense", "C"}, "--num-enzyme-termini", "3", "--clip-n-term-m", "1"));
        options.validate();
        assertEquals("Trypsin/P", options.footer("").get("search_enzyme_name_1"));
        assertTrue(options.describe().contains("Trypsin/P"), options.describe());
    }

    // --- counts -------------------------------------------------------------------------------

    @Test
    public void atMostTwoEnzymes() {
        String message = rejection("--cut", "K", "--cut", "R", "--cut", "D",
                "--no-cut", "@", "--no-cut", "@", "--no-cut", "@", "--sense", "C", "--sense", "C",
                "--sense", "C", "--num-enzyme-termini", "2", "--clip-n-term-m", "1");
        assertTrue(message.contains("--cut") && message.contains("2"), message);
    }

    @Test
    public void namesAreGivenForEveryEnzymeOrForNone() {
        String message = rejection("--enzyme-name", "trypsin", "--cut", "KR", "--cut", "D",
                "--no-cut", "@", "--no-cut", "@", "--sense", "C", "--sense", "N",
                "--num-enzyme-termini", "2", "--clip-n-term-m", "1");
        assertTrue(message.contains("--enzyme-name"), message);
    }

    @Test
    public void everyCutNeedsItsNoCutAndSense() {
        assertTrue(rejection("--cut", "KR", "--sense", "C", "--num-enzyme-termini", "2",
                "--clip-n-term-m", "1").contains("--no-cut"));
        assertTrue(rejection("--cut", "KR", "--no-cut", "@", "--num-enzyme-termini", "2",
                "--clip-n-term-m", "1").contains("--sense"));
        assertTrue(rejection("--cut", "KR", "--cut", "D", "--no-cut", "@", "--sense", "C",
                "--sense", "N", "--num-enzyme-termini", "2", "--clip-n-term-m", "1").contains("--no-cut"),
                "two cuts, one no-cut");
    }

    // There is no default for either: FragDIA would otherwise be digesting with a number nobody chose.
    @Test
    public void terminiAndClipAreRequiredOnceWithACut() {
        assertTrue(rejection("--cut", "KR", "--no-cut", "@", "--sense", "C", "--clip-n-term-m", "1")
                .contains("--num-enzyme-termini"));
        assertTrue(rejection("--cut", "KR", "--no-cut", "@", "--sense", "C", "--num-enzyme-termini", "2")
                .contains("--clip-n-term-m"));
        assertTrue(rejection(with(TRYPSIN, "--num-enzyme-termini", "2")).contains("--num-enzyme-termini"),
                "given twice");
        assertTrue(rejection(with(TRYPSIN, "--clip-n-term-m", "1")).contains("--clip-n-term-m"),
                "given twice");
    }

    // --- values: validated, then used exactly as given ---------------------------------------

    @Test
    public void residuesAreUppercaseLettersOrAnAtSign() {
        for (String bad : new String[] {"kr", "Kr", "KR@", "@@", "K R", " KR", "KR ", "K,R", "1", ""}) {
            assertTrue(rejection("--cut", bad, "--no-cut", "@", "--sense", "C", "--num-enzyme-termini", "2",
                    "--clip-n-term-m", "1").contains("--cut"), "--cut " + bad);
            assertTrue(rejection("--cut", "KR", "--no-cut", bad, "--sense", "C", "--num-enzyme-termini", "2",
                    "--clip-n-term-m", "1").contains("--no-cut"), "--no-cut " + bad);
        }
        for (String good : new String[] {"KR", "@", "ABCDEFGHIJKLMNOPQRSTUVWXYZ", "P"}) {
            given("--cut", good, "--no-cut", good, "--sense", "N", "--num-enzyme-termini", "2",
                    "--clip-n-term-m", "1").validate();
        }
    }

    // Uppercase only, and only one letter: FragPipe upper-cases what it passes, so a lower-case
    // sense here is a caller that skipped that step, not one to be forgiven.
    @Test
    public void senseIsAnUppercaseCOrN() {
        for (String bad : new String[] {"c", "n", "CN", "after", "", " C"}) {
            assertTrue(rejection("--cut", "KR", "--no-cut", "@", "--sense", bad, "--num-enzyme-termini", "2",
                    "--clip-n-term-m", "1").contains("--sense"), "--sense " + bad);
        }
    }

    @Test
    public void terminiAreZeroToThree() {
        for (String bad : new String[] {"4", "-1", "02", "two", "", "2 "}) {
            assertTrue(rejection("--cut", "KR", "--no-cut", "@", "--sense", "C", "--num-enzyme-termini", bad,
                    "--clip-n-term-m", "1").contains("--num-enzyme-termini"), "--num-enzyme-termini " + bad);
        }
        for (String good : new String[] {"0", "1", "2", "3"}) {
            given("--cut", "KR", "--no-cut", "@", "--sense", "C", "--num-enzyme-termini", good,
                    "--clip-n-term-m", "1").validate();
        }
    }

    // Always a value, and only these two: there is no bare form to read as "on".
    @Test
    public void clipIsZeroOrOne() {
        for (String bad : new String[] {"true", "false", "yes", "2", "", "on"}) {
            assertTrue(rejection("--cut", "KR", "--no-cut", "@", "--sense", "C", "--num-enzyme-termini", "2",
                    "--clip-n-term-m", bad).contains("--clip-n-term-m"), "--clip-n-term-m " + bad);
        }
    }

    @Test
    public void anEmptyNameIsRejected() {
        assertTrue(rejection("--enzyme-name", "", "--cut", "KR", "--no-cut", "@", "--sense", "C",
                "--num-enzyme-termini", "2", "--clip-n-term-m", "1").contains("--enzyme-name"));
    }

    // --- the footer as DuckDB writes it --------------------------------------------------------

    // The proof that DuckDB, at the version this build pins, writes KV_METADATA into a Parquet footer
    // and that a reader finds exactly those keys there - nothing more. The values are the awkward
    // ones a path brings: a Windows backslash, a space, a quote, a character outside ASCII.
    @Test
    public void duckDbWritesTheFooterAndAReaderFindsExactlyIt(@TempDir Path dir) throws Exception {
        EnzymeOptions options = given(TRYPSIN);
        options.validate();
        String fasta = "C:\\Users\\o'brien\\FASTA files\\prot\u00e9ome.fasta";
        Map<String, String> footer = options.footer(fasta);

        File parquet = dir.resolve("library.parquet").toFile();
        try (Connection conn = DriverManager.getConnection("jdbc:duckdb:");
             Statement stmt = conn.createStatement()) {
            stmt.execute("COPY (SELECT CAST(1 AS FLOAT) AS PrecursorMz) TO '" +
                    parquet.getAbsolutePath().replace("\\", "/") + "' (FORMAT PARQUET" +
                    EnzymeOptions.kvMetadataOption(footer) + ")");
        }

        Map<String, String> read = ParquetFooters.read(parquet);
        assertEquals(footer, read);
        assertEquals(fasta, read.get("database_name"), "the path exactly as given");
    }

    // With nothing to record the COPY must be the one it always was, not one with an empty clause.
    @Test
    public void anEmptyFooterAddsNothingToTheCopy() {
        assertEquals("", EnzymeOptions.kvMetadataOption(new LinkedHashMap<>()));
    }
}
