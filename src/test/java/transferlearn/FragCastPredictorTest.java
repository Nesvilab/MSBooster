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

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertFalse;
import static org.junit.jupiter.api.Assertions.assertTrue;

import java.io.File;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.util.HashMap;
import java.nio.file.Path;
import java.sql.Connection;
import java.sql.DriverManager;
import java.sql.Statement;
import java.util.ArrayList;
import java.util.List;

import org.junit.jupiter.api.AfterEach;
import org.junit.jupiter.api.BeforeEach;
import org.junit.jupiter.api.Test;
import org.junit.jupiter.api.io.TempDir;

/**
 * The peptide list arrives with each precursor's charge glued onto the end of its peptide, and the
 * only way to get it back off is to strip trailing digits. That is the same thing the server path
 * does, and it is the one place in the local predictor where a wrong answer is silent: a peptide
 * split in the wrong place still looks like a peptide and still predicts, just for the wrong
 * precursor.
 */
public class FragCastPredictorTest {

    @TempDir
    Path tmp;

    // These drive the real parameter machinery, which writes into public statics across the whole
    // of Constants - not just the handful this class names. Snapshotting every one of them by
    // reflection is what keeps a test here from changing the outcome of a test that runs later:
    // createPredFileOnly and requirePinMzml in particular change what MainClass does, and several of
    // the paths point into a @TempDir that JUnit deletes.
    private final java.util.Map<java.lang.reflect.Field, Object> constantsBefore = new java.util.HashMap<>();

    @BeforeEach
    public void saveConstants() throws Exception {
        constantsBefore.clear();
        for (java.lang.reflect.Field f : allconstants.Constants.class.getFields()) {
            if (!java.lang.reflect.Modifier.isFinal(f.getModifiers())
                    && java.lang.reflect.Modifier.isStatic(f.getModifiers())) {
                constantsBefore.put(f, f.get(null));
            }
        }
    }

    @AfterEach
    public void restoreConstants() throws Exception {
        for (java.util.Map.Entry<java.lang.reflect.Field, Object> e : constantsBefore.entrySet()) {
            e.getKey().set(null, e.getValue());
        }
    }

    /** A peptide-list parquet with the one column the conversion reads. */
    private File peptideList(String... peptides) throws Exception {
        return peptideListLabelled("sp|P12345|PROT_HUMAN", peptides);
    }

    /** As above, with the protein labels the list carries named outright. */
    private File peptideListLabelled(String proteins, String... peptides) throws Exception {
        File parquet = tmp.resolve("peptides.parquet").toFile();
        StringBuilder values = new StringBuilder();
        for (int i = 0; i < peptides.length; i++) {
            if (i > 0) {
                values.append(", ");
            }
            values.append("('").append(peptides[i]).append("', '").append(proteins).append("', false)");
        }
        try (Connection conn = DriverManager.getConnection("jdbc:duckdb:");
             Statement stmt = conn.createStatement()) {
            stmt.execute("COPY (SELECT * FROM (VALUES " + values + ") AS t(peptide, proteins, is_decoy)) TO '" +
                    parquet.getAbsolutePath().replace("\\", "/") + "' (FORMAT PARQUET)");
        }
        return parquet;
    }

    private List<String[]> convert(File parquet, int minCharge, int maxCharge) throws Exception {
        File out = tmp.resolve("fragcast_in.tsv").toFile();
        FragCastPredictor.writeFragCastInput(parquet.getAbsolutePath(), out, minCharge, maxCharge);

        List<String> lines = Files.readAllLines(out.toPath(), StandardCharsets.UTF_8);
        assertEquals("peptide\tcharge\tproteins\tis_decoy", lines.get(0),
                "FragCast reads peptide, charge and the proteins it labels the library from");

        List<String[]> rows = new ArrayList<>();
        for (int i = 1; i < lines.size(); i++) {
            if (!lines.get(i).trim().isEmpty()) {
                rows.add(lines.get(i).split("\t", -1));
            }
        }
        return rows;
    }

    @Test
    public void takesTheChargeOffTheEndOfThePeptide() throws Exception {
        List<String[]> rows = convert(peptideList("PEPTIDEK2", "SEQUENCER3"), 2, 3);
        assertEquals(2, rows.size());
        assertEquals("PEPTIDEK", rows.get(0)[0]);
        assertEquals("2", rows.get(0)[1]);
        assertEquals("SEQUENCER", rows.get(1)[0]);
        assertEquals("3", rows.get(1)[1]);
    }

    @Test
    public void keepsAModificationWhoseMassEndsInADigit() throws Exception {
        // the reason trailing-digit stripping is safe: a delta mass closes with a bracket, so the
        // scan stops at it and cannot eat into the modification
        List<String[]> rows = convert(peptideList("PEPTIDEK[229.1629]2"), 2, 3);
        assertEquals(1, rows.size());
        assertEquals("2", rows.get(0)[1]);
        assertTrue(rows.get(0)[0].contains("229.16"), rows.get(0)[0]);
        assertFalse(rows.get(0)[0].endsWith("2"), "the charge must not be left on the peptide");
    }

    @Test
    public void expandsAPeptideThatNamesNoChargeAcrossTheRange() throws Exception {
        List<String[]> rows = convert(peptideList("PEPTIDEK"), 2, 4);
        assertEquals(3, rows.size());
        assertEquals("2", rows.get(0)[1]);
        assertEquals("3", rows.get(1)[1]);
        assertEquals("4", rows.get(2)[1]);
        for (String[] row : rows) {
            assertEquals("PEPTIDEK", row[0]);
        }
    }

    @Test
    public void dropsPrecursorsFragCastCannotPredict() throws Exception {
        // FragCastCharges caps what the model was trained for; a charge past it would otherwise be
        // sent to FragCast and come back as a prediction for a precursor it cannot represent
        List<String[]> rows = convert(peptideList("PEPTIDEK99"), 2, 3);
        assertTrue(rows.isEmpty(), "an unrepresentable charge should be skipped, got " + rows.size());
    }

    @Test
    public void skipsEmptyRowsWithoutFailing() throws Exception {
        List<String[]> rows = convert(peptideList("PEPTIDEK2", "", "SEQUENCER2"), 2, 3);
        assertEquals(2, rows.size());
    }



    // --- the command line must beat the parameter file -------------------------------------

    /** A parameter file that sets, and would otherwise silently win, everything the CLI names. */
    private File paramsFileNaming(String... keyValues) throws Exception {
        File f = tmp.resolve("params.txt").toFile();
        StringBuilder sb = new StringBuilder();
        for (String kv : keyValues) {
            sb.append(kv).append("\n");
        }
        Files.write(f.toPath(), sb.toString().getBytes(StandardCharsets.UTF_8));
        return f;
    }

    private static void applyAsPredictorWould(String[] args) throws Exception {
        java.util.HashMap<String, String> map = new java.util.HashMap<>();
        mainsteps.ParameterUtils.processCommandLineInputs(args, map);
        mainsteps.ParameterUtils.updateConstants(map);
    }

    @Test
    public void theNamedModelBundleSurvivesTheParameterFile() throws Exception {
        // the shipped template carries "FragCastModelZip =", and an empty value is still a value:
        // assigning Constants before reading the file let the file blank it out and the fine-tuned
        // model was silently ignored
        File params = paramsFileNaming("FragCastModelZip =", "decoyPrefix = rev_", "outputDirectory = " +
                tmp.toAbsolutePath().toString().replace("\\", "/"));
        String zip = tmp.resolve("finetuned.zip").toAbsolutePath().toString();

        applyAsPredictorWould(FragCastPredictor.peptideListArgs(params.getAbsolutePath(), "1", zip, "DECOY_"));

        assertEquals(zip, allconstants.Constants.FragCastModelZip);
        assertEquals("DECOY_", allconstants.Constants.decoyPrefix);
        assertEquals(Integer.valueOf(1), allconstants.Constants.keepDecoys);
    }

    @Test
    public void thePinPathAlsoOutranksTheParameterFile() throws Exception {
        File params = paramsFileNaming("FragCastModelZip =", "spectraModel = DIA-NN", "keepDecoys = 1",
                "outputDirectory = " + tmp.toAbsolutePath().toString().replace("\\", "/"));
        String zip = tmp.resolve("finetuned.zip").toAbsolutePath().toString();

        applyAsPredictorWould(FragCastPredictor.pinPathArgs(params.getAbsolutePath(), "0", zip, "DECOY_", false));

        assertEquals(zip, allconstants.Constants.FragCastModelZip);
        assertEquals(allconstants.FragCastModels.CONFORMER, allconstants.Constants.spectraModel);
        assertEquals(Integer.valueOf(0), allconstants.Constants.keepDecoys);
    }

    // --- the parameter file's fast choice must survive the pin path's model overrides --------

    @Test
    public void theParameterFileCanAskForTheFastSpecModel() throws Exception {
        File params = paramsFileNaming("spectraModel = " + allconstants.FragCastModels.FAST);
        assertTrue(FragCastPredictor.requestsFastSpecModel(params.getAbsolutePath()));
    }

    @Test
    public void theFastRequestReadsLikeARescoringRunWould() throws Exception {
        // MainClass reads the same line through LowercaseModelMapper, so a hand-written
        // "fragcast-fast" runs fast when rescoring; this entry point must agree with it
        File params = paramsFileNaming("spectraModel = fragcast-fast");
        assertTrue(FragCastPredictor.requestsFastSpecModel(params.getAbsolutePath()));
    }

    @Test
    public void anyOtherSpecModelIsNotAFastRequest() throws Exception {
        assertFalse(FragCastPredictor.requestsFastSpecModel(
                paramsFileNaming("spectraModel = " + allconstants.FragCastModels.CONFORMER).getAbsolutePath()));
        assertFalse(FragCastPredictor.requestsFastSpecModel(
                paramsFileNaming("spectraModel = DIA-NN").getAbsolutePath()));
        assertFalse(FragCastPredictor.requestsFastSpecModel(
                paramsFileNaming("decoyPrefix = rev_").getAbsolutePath()),
                "a file naming no Spec model asks for nothing");
    }

    @Test
    public void thePinPathForwardsTheFastChoiceForSpectraOnly() throws Exception {
        // Without this the pin path flattened the file's choice to the Conformer, and the run
        // silently predicted with weights the user never picked. RT and IM stay under the plain
        // name: only the MS2 weights come in a fast variant.
        File params = paramsFileNaming("spectraModel = " + allconstants.FragCastModels.FAST,
                "outputDirectory = " + tmp.toAbsolutePath().toString().replace("\\", "/"));

        applyAsPredictorWould(FragCastPredictor.pinPathArgs(params.getAbsolutePath(), "0", "", "rev_", true));

        assertEquals(allconstants.FragCastModels.FAST, allconstants.Constants.spectraModel);
        assertEquals(allconstants.FragCastModels.CONFORMER, allconstants.Constants.rtModel);
        assertEquals(allconstants.FragCastModels.CONFORMER, allconstants.Constants.imModel);
        assertTrue(allconstants.FragCastModels.usesFastSpecModel(allconstants.Constants.spectraModel),
                "the prediction derives FragCast's --fast flag from what Constants now says");
    }

    @Test
    public void noBundleMeansNoOverride() throws Exception {
        // an empty --model must not blank out a path the parameter file legitimately set
        File params = paramsFileNaming("FragCastModelZip = /from/the/file.zip",
                "outputDirectory = " + tmp.toAbsolutePath().toString().replace("\\", "/"));

        applyAsPredictorWould(FragCastPredictor.peptideListArgs(params.getAbsolutePath(), "1", "", "rev_"));

        assertEquals("/from/the/file.zip", allconstants.Constants.FragCastModelZip);
    }
    @Test
    public void handsTheProteinsStraightToFragCast() throws Exception {
        // FragCast fills ProteinId, GeneName, AllMappedProteins and Proteotypic from this column, so
        // passing it is the whole of the protein annotation. It also has to survive the peptide being
        // one whose sequence FragCast rewrites: it drops the non-standard residues B/J/O/U/X/Z from
        // PeptideSequence, so anything that labelled the library afterwards by matching on that
        // sequence would miss exactly the selenocysteine peptides this covers.
        List<String[]> rows = convert(peptideList("PEPUTIDEK2"), 2, 3);
        assertEquals(1, rows.size());
        assertEquals(4, rows.get(0).length, "peptide, charge, proteins and is_decoy");
        assertEquals("sp|P12345|PROT_HUMAN", rows.get(0)[2]);
    }

    @Test
    public void everyExpandedChargeCarriesTheSameProteins() throws Exception {
        List<String[]> rows = convert(peptideList("PEPTIDEK"), 2, 3);
        assertEquals(2, rows.size());
        for (String[] row : rows) {
            assertEquals("sp|P12345|PROT_HUMAN", row[2]);
        }
    }

    // A peptide list is as free to terminate its labels with a ';' as MSFragger's pin is, and
    // FragCast reads one either way as "this peptide is shared" - so a list written that way would
    // predict a library with Proteotypic = 0 on every row.
    @Test
    public void aTerminatedProteinLabelDoesNotMakeThePeptideLookShared() throws Exception {
        List<String[]> rows = convert(peptideListLabelled("sp|P12345|PROT_HUMAN;", "PEPTIDEK2"), 2, 3);
        assertEquals(1, rows.size());
        assertEquals("sp|P12345|PROT_HUMAN", rows.get(0)[2]);
    }

    // The separator between two real labels is what says the peptide IS shared, and it survives.
    @Test
    public void aGenuinelySharedPeptideKeepsItsSeparator() throws Exception {
        List<String[]> rows = convert(
                peptideListLabelled("sp|P1|A_HUMAN;sp|P2|B_HUMAN;", "PEPTIDEK2"), 2, 3);
        assertEquals(1, rows.size());
        assertEquals("sp|P1|A_HUMAN;sp|P2|B_HUMAN", rows.get(0)[2]);
    }
    // Without --peptide-list the peptides come from this job's pin files, and the only file
    // carrying their proteins is the TSV written for FragCast. Reading it as Parquet - as the
    // server path can, because its own generated input is Parquet - fails after the whole
    // prediction has already been paid for.
    @Test
    public void theProteinsOfAPinDerivedRunAreReadFromTheFragCastInput() throws Exception {
        File list = tmp.resolve("spectraRT.tsv").toFile();
        Files.write(list.toPath(), String.join(System.lineSeparator(),
                "peptide	charge	proteins",
                "KWFSQHNHLK	2	sp|O75747|P3C2G_HUMAN;",
                "PEPTIDEK	2	sp|P12345|TEST_HUMAN").getBytes(StandardCharsets.UTF_8));

        HashMap<String, String> protToGene = new HashMap<>();
        protToGene.put("O75747", "PIK3C2G");
        protToGene.put("P12345", "TEST");

        HashMap<String, String> mapped = Helpers.mapProteinsListToGenes(list.getAbsolutePath(), protToGene);
        assertEquals("PIK3C2G", mapped.get("O75747"), "a TSV peptide list produced no gene label");
        assertEquals("TEST", mapped.get("P12345"));
    }

    // --- the library the converter reads, whichever predictor wrote it ----------------------

    /**
     * The schema everything downstream reads: the AlphaPeptDeep server's own columns, taken from a
     * real one of its outputs, plus the two only FragCast can fill. Those two sit where the
     * experimental library FragSpecLib writes puts them - after FragmentLossType, before Proteotypic -
     * so a predicted library and an experimental one line up column for column.
     */
    private static final String[][] LIBRARY_SCHEMA = {
            {"PrecursorMz", "FLOAT"}, {"ProductMz", "FLOAT"}, {"Annotation", "VARCHAR"},
            {"ProteinId", "VARCHAR"}, {"PeptideSequence", "VARCHAR"},
            {"ModifiedPeptideSequence", "VARCHAR"}, {"PrecursorCharge", "SMALLINT"},
            {"LibraryIntensity", "FLOAT"}, {"NormalizedRetentionTime", "FLOAT"},
            {"PrecursorIonMobility", "VARCHAR"}, {"FragmentType", "VARCHAR"},
            {"FragmentCharge", "SMALLINT"}, {"FragmentSeriesNumber", "SMALLINT"},
            {"FragmentLossType", "VARCHAR"},
            {"AverageExperimentalRetentionTime", "VARCHAR"}, {"AllMappedProteins", "VARCHAR"},
            {"Proteotypic", "SMALLINT"}};

    /** A copy of the committed real FragCast library, which carries FragCast's own 19 columns. */
    private File fragCastLibrary() throws Exception {
        File copy = tmp.resolve("lib.parquet").toFile();
        Files.copy(new File("src/test/resources/fragcast_lib.parquet").toPath(), copy.toPath());
        return copy;
    }

    @Test
    public void theLibraryEndsUpInTheServersSchema() throws Exception {
        File lib = fragCastLibrary();
        FragCastPredictor.toLibrarySchema(lib);

        try (Connection conn = DriverManager.getConnection("jdbc:duckdb:");
             Statement stmt = conn.createStatement();
             java.sql.ResultSet rs = stmt.executeQuery(
                     "SELECT * FROM read_parquet('" + lib.getAbsolutePath().replace("\\", "/") + "') LIMIT 0")) {
            java.sql.ResultSetMetaData m = rs.getMetaData();
            assertEquals(LIBRARY_SCHEMA.length, m.getColumnCount(),
                    "the gene columns go, the group columns stay");
            for (int i = 0; i < LIBRARY_SCHEMA.length; i++) {
                assertEquals(LIBRARY_SCHEMA[i][0], m.getColumnName(i + 1), "column " + (i + 1));
                assertEquals(LIBRARY_SCHEMA[i][1], m.getColumnTypeName(i + 1),
                        "type of " + m.getColumnName(i + 1));
            }
        }
    }

    @Test
    public void theConversionKeepsEveryRow() throws Exception {
        File lib = fragCastLibrary();
        long before = rowCount(lib);
        FragCastPredictor.toLibrarySchema(lib);
        assertEquals(before, rowCount(lib));
        assertTrue(before > 0);
    }

    @Test
    public void everyPropertyIsKept() throws Exception {
        File lib = fragCastLibrary();
        FragCastPredictor.toLibrarySchema(lib);
        try (Connection conn = DriverManager.getConnection("jdbc:duckdb:");
             Statement stmt = conn.createStatement();
             java.sql.ResultSet rs = stmt.executeQuery(
                     "SELECT COUNT(NormalizedRetentionTime) AS rt, " +
                     "COUNT(NULLIF(PrecursorIonMobility, '')) AS im FROM read_parquet('" +
                     lib.getAbsolutePath().replace("\\", "/") + "')")) {
            rs.next();
            // one build-library run predicts MS2, RT and IM together and nothing drops any of them
            assertTrue(rs.getLong("rt") > 0, "retention times must survive");
            assertTrue(rs.getLong("im") > 0, "ion mobilities must survive");
        }
    }

    private long rowCount(File parquet) throws Exception {
        try (Connection conn = DriverManager.getConnection("jdbc:duckdb:");
             Statement stmt = conn.createStatement();
             java.sql.ResultSet rs = stmt.executeQuery("SELECT COUNT(*) FROM read_parquet('" +
                     parquet.getAbsolutePath().replace("\\", "/") + "')")) {
            rs.next();
            return rs.getLong(1);
        }
    }
    @Test
    public void marksDecoysSoFragCastCanTellThemApart() throws Exception {
        // the peptide list is roughly half decoys; is_decoy is what lets FragCast keep them under
        // --keep-decoys and leave them out without it, so the marking must survive the conversion
        // truthfully either way
        File parquet = tmp.resolve("mixed.parquet").toFile();
        try (Connection conn = DriverManager.getConnection("jdbc:duckdb:");
             Statement stmt = conn.createStatement()) {
            stmt.execute("COPY (SELECT * FROM (VALUES " +
                    "('PEPTIDEK2', 'sp|P12345|PROT_HUMAN', false), " +
                    "('KEDITPEP2', 'rev_sp|P12345|PROT_HUMAN', true)) " +
                    "AS t(peptide, proteins, is_decoy)) TO '" +
                    parquet.getAbsolutePath().replace("\\", "/") + "' (FORMAT PARQUET)");
        }
        List<String[]> rows = convert(parquet, 2, 3);
        assertEquals(2, rows.size());
        assertEquals("false", rows.get(0)[3]);
        assertEquals("true", rows.get(1)[3]);
    }
}
