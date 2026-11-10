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
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashMap;
import java.nio.file.Path;
import java.sql.Connection;
import java.sql.DriverManager;
import java.sql.Statement;
import java.util.ArrayList;
import java.util.List;
import java.util.Map;

import blue.strategic.parquet.ParquetWriter;
import org.apache.parquet.schema.LogicalTypeAnnotation;
import org.apache.parquet.schema.MessageType;
import org.apache.parquet.schema.PrimitiveType.PrimitiveTypeName;
import org.apache.parquet.schema.Type;
import org.apache.parquet.schema.Types;
import org.junit.jupiter.api.AfterEach;
import org.junit.jupiter.api.BeforeEach;
import org.junit.jupiter.api.Test;
import org.junit.jupiter.api.io.TempDir;
import testutils.ParquetFooters;

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
     * The schema everything downstream reads: the AlphaPeptDeep server's own columns, in the types
     * its writer gives them - the ion mobility a float, which is what the server writes whenever it
     * predicts one, and the only type FragDIA takes one in - plus the two only FragCast can fill. Those two sit where the
     * experimental library FragSpecLib writes puts them - after FragmentLossType, before Proteotypic -
     * so a predicted library and an experimental one line up column for column. The experimental
     * retention time is a float, as FragSpecLib writes its measured one; a predicted library's is null.
     */
    private static final String[][] LIBRARY_SCHEMA = {
            {"PrecursorMz", "FLOAT"}, {"ProductMz", "FLOAT"}, {"Annotation", "VARCHAR"},
            {"ProteinId", "VARCHAR"}, {"PeptideSequence", "VARCHAR"},
            {"ModifiedPeptideSequence", "VARCHAR"}, {"PrecursorCharge", "SMALLINT"},
            {"LibraryIntensity", "FLOAT"}, {"NormalizedRetentionTime", "FLOAT"},
            {"PrecursorIonMobility", "FLOAT"}, {"FragmentType", "VARCHAR"},
            {"FragmentCharge", "SMALLINT"}, {"FragmentSeriesNumber", "SMALLINT"},
            {"FragmentLossType", "VARCHAR"},
            {"AverageExperimentalRetentionTime", "FLOAT"}, {"AllMappedProteins", "VARCHAR"},
            {"Proteotypic", "SMALLINT"}};

    private static final String AVERAGE_EXPERIMENTAL_RT = "AverageExperimentalRetentionTime";

    /** What a run that names neither an enzyme nor a FASTA records in a library's footer. */
    private static final Map<String, String> NO_FOOTER = Collections.emptyMap();

    /** The options FragPipe passes for strict trypsin. */
    private static EnzymeOptions strictTrypsin() {
        EnzymeOptions enzyme = new EnzymeOptions();
        String[] args = {"--enzyme-name", "stricttrypsin", "--cut", "KR", "--no-cut", "@", "--sense", "C",
                "--num-enzyme-termini", "2", "--clip-n-term-m", "1"};
        for (int i = 0; i < args.length; i += 2) {
            assertTrue(enzyme.accept(args[i], args[i + 1]));
        }
        enzyme.validate();
        return enzyme;
    }

    /**
     * A copy of fragcast_lib_narrow.parquet: a real build-library output, in the 17 columns of the
     * schema, from the FragCast before AverageExperimentalRetentionTime became a float - it holds the
     * text "" there.
     */
    private File olderFragCastLibrary() throws Exception {
        File copy = tmp.resolve("narrow.parquet").toFile();
        Files.copy(new File("src/test/resources/fragcast_lib_narrow.parquet").toPath(), copy.toPath());
        return copy;
    }

    /**
     * A library in the schema as FragCast writes it now, by a writer other than DuckDB: the rows of
     * the real build-library output above, written by parquet-java with every column required - as
     * parquet-rs writes FragCast's - except AverageExperimentalRetentionTime, an optional float that
     * holds nothing. A stand-in until a build-library output from that FragCast is committed.
     */
    private File inTheSchemaFromAnotherWriter() throws Exception {
        List<Type> fields = new ArrayList<>();
        for (String[] column : LIBRARY_SCHEMA) {
            Type.Repetition repetition = column[0].equals(AVERAGE_EXPERIMENTAL_RT)
                    ? Type.Repetition.OPTIONAL : Type.Repetition.REQUIRED;
            switch (column[1]) {
                case "FLOAT":
                    fields.add(Types.primitive(PrimitiveTypeName.FLOAT, repetition).named(column[0]));
                    break;
                case "SMALLINT":
                    fields.add(Types.primitive(PrimitiveTypeName.INT32, repetition)
                            .as(LogicalTypeAnnotation.intType(16, true)).named(column[0]));
                    break;
                default:
                    fields.add(Types.primitive(PrimitiveTypeName.BINARY, repetition)
                            .as(LogicalTypeAnnotation.stringType()).named(column[0]));
                    break;
            }
        }
        MessageType schema = new MessageType("library", fields);

        File out = tmp.resolve("in_the_schema.parquet").toFile();
        try (Connection conn = DriverManager.getConnection("jdbc:duckdb:");
             Statement stmt = conn.createStatement();
             java.sql.ResultSet rs = stmt.executeQuery("SELECT * REPLACE (CAST(NULL AS FLOAT) AS " +
                     AVERAGE_EXPERIMENTAL_RT + ") FROM read_parquet('" +
                     olderFragCastLibrary().getAbsolutePath().replace("\\", "/") + "')");
             ParquetWriter<Object[]> writer = ParquetWriter.writeFile(schema, out, (row, values) -> {
                 for (int i = 0; i < LIBRARY_SCHEMA.length; i++) {
                     if (row[i] != null) {
                         values.write(LIBRARY_SCHEMA[i][0], row[i]);
                     }
                 }
             })) {
            while (rs.next()) {
                Object[] row = new Object[LIBRARY_SCHEMA.length];
                for (int i = 0; i < row.length; i++) {
                    Object value = rs.getObject(i + 1);
                    row[i] = value instanceof Short ? Integer.valueOf((Short) value) : value;
                }
                writer.write(row);
            }
        }
        return out;
    }

    /** A copy of the committed real FragCast library, which carries FragCast's own 19 columns. */
    private File fragCastLibrary() throws Exception {
        File copy = tmp.resolve("lib.parquet").toFile();
        Files.copy(new File("src/test/resources/fragcast_lib.parquet").toPath(), copy.toPath());
        return copy;
    }

    private static void assertLibrarySchema(File lib) throws Exception {
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
    public void theLibraryEndsUpInTheServersSchema() throws Exception {
        File lib = fragCastLibrary();
        FragCastPredictor.toLibrarySchema(lib, NO_FOOTER);
        assertLibrarySchema(lib);
    }

    // FragCast now writes this schema itself, and then the conversion has nothing to do - which on a
    // FASTA-scale library is a whole pass over the file, minutes of it, not spent. The library has to
    // be left exactly as it was found: not rewritten into the same thing, not touched at all.
    @Test
    public void aLibraryAlreadyInTheSchemaIsNotRewritten() throws Exception {
        File lib = fragCastLibrary();
        FragCastPredictor.toLibrarySchema(lib, NO_FOOTER); //an older FragCast's 19 columns, brought into the schema
        assertLibrarySchema(lib); //or the second call below would be skipping for the wrong reason
        byte[] inTheSchema = Files.readAllBytes(lib.toPath());
        long longAgo = 1_000_000_000_000L; //2001: any rewrite, even to identical bytes, moves this
        assertTrue(lib.setLastModified(longAgo));

        FragCastPredictor.toLibrarySchema(lib, NO_FOOTER);

        assertEquals(longAgo, lib.lastModified(), "a library already in the schema was written again");
        assertArrayEquals(inTheSchema, Files.readAllBytes(lib.toPath()));
        assertFalse(new File(lib.getAbsolutePath() + ".apd.parquet").exists());
    }

    // The same, for the file this is really about: a library FragCast itself wrote in the schema.
    // The one above was written by DuckDB; FragCast's has non-nullable int16/float32 columns, one
    // nullable float, and a writer's own encodings, and has to be recognised just the same - or the
    // six minutes are spent after all. Whatever FragCast recorded in its footer stays with it.
    @Test
    public void aLibraryFragCastWroteInTheSchemaIsLeftAsItIs() throws Exception {
        File lib = inTheSchemaFromAnotherWriter();
        byte[] asWritten = Files.readAllBytes(lib.toPath());
        Map<String, String> ownFooter = ParquetFooters.read(lib);
        long longAgo = 1_000_000_000_000L;
        assertTrue(lib.setLastModified(longAgo));

        assertLibrarySchema(lib); //FragCast's own output is the schema, before anything touches it
        FragCastPredictor.toLibrarySchema(lib, strictTrypsin().footer("/data/human.fasta"));

        assertEquals(longAgo, lib.lastModified(), "FragCast's own library was written again");
        assertArrayEquals(asWritten, Files.readAllBytes(lib.toPath()));
        assertEquals(ownFooter, ParquetFooters.read(lib));
        assertFalse(new File(lib.getAbsolutePath() + ".apd.parquet").exists());
    }

    // The FragCast before this one wrote AverageExperimentalRetentionTime as the text "", so its
    // library - a real build-library output - is not the schema any more and is converted. DuckDB
    // refuses to cast '' to a float, so the empty text has to become null rather than fail the run;
    // nothing else about the library changes.
    @Test
    public void anOlderFragCastsEmptyTextRetentionTimeBecomesNull() throws Exception {
        File original = olderFragCastLibrary();
        File lib = tmp.resolve("converted.parquet").toFile();
        Files.copy(original.toPath(), lib.toPath());

        FragCastPredictor.toLibrarySchema(lib, NO_FOOTER);

        assertLibrarySchema(lib);
        assertEquals(rowCount(original), rowCount(lib));
        String others = "SELECT * EXCLUDE (" + AVERAGE_EXPERIMENTAL_RT + ") FROM read_parquet('%s')";
        try (Connection conn = DriverManager.getConnection("jdbc:duckdb:");
             Statement stmt = conn.createStatement();
             java.sql.ResultSet rs = stmt.executeQuery("SELECT " +
                     "(SELECT COUNT(" + AVERAGE_EXPERIMENTAL_RT + ") FROM read_parquet('" + path(lib) + "')) AS kept, " +
                     "(SELECT COUNT(*) FROM (" + String.format(others, path(original)) + " EXCEPT ALL " +
                     String.format(others, path(lib)) + ")) AS changed")) {
            rs.next();
            assertEquals(0, rs.getLong("kept"), "a predicted library has no experimental retention time");
            assertEquals(0, rs.getLong("changed"), "another column changed in the conversion");
        }
    }

    // Only the EMPTY text is read as null. Text that is a number is that number, and text that is
    // neither still fails the conversion, rather than quietly becoming null.
    @Test
    public void aTextRetentionTimeMustStillBeANumberOrEmpty() throws Exception {
        String narrow = path(olderFragCastLibrary());
        File numbers = tmp.resolve("numbers.parquet").toFile();
        File words = tmp.resolve("words.parquet").toFile();
        try (Connection conn = DriverManager.getConnection("jdbc:duckdb:");
             Statement stmt = conn.createStatement()) {
            stmt.execute("COPY (SELECT * REPLACE (CASE WHEN FragmentSeriesNumber = 2 THEN '12.5' ELSE '' END " +
                    "AS " + AVERAGE_EXPERIMENTAL_RT + ") FROM read_parquet('" + narrow + "')) TO '" +
                    path(numbers) + "' (FORMAT PARQUET)");
            stmt.execute("COPY (SELECT * REPLACE ('n/a' AS " + AVERAGE_EXPERIMENTAL_RT + ") FROM read_parquet('" +
                    narrow + "')) TO '" + path(words) + "' (FORMAT PARQUET)");
        }

        FragCastPredictor.toLibrarySchema(numbers, NO_FOOTER);
        assertLibrarySchema(numbers);
        try (Connection conn = DriverManager.getConnection("jdbc:duckdb:");
             Statement stmt = conn.createStatement();
             java.sql.ResultSet rs = stmt.executeQuery("SELECT COUNT(*) FILTER (WHERE FragmentSeriesNumber = 2) " +
                     "AS twos, COUNT(" + AVERAGE_EXPERIMENTAL_RT + ") AS kept, MIN(" + AVERAGE_EXPERIMENTAL_RT +
                     ") AS lo, MAX(" + AVERAGE_EXPERIMENTAL_RT + ") AS hi FROM read_parquet('" + path(numbers) + "')")) {
            rs.next();
            assertTrue(rs.getLong("twos") > 0);
            assertEquals(rs.getLong("twos"), rs.getLong("kept"));
            assertEquals(12.5f, rs.getFloat("lo"));
            assertEquals(12.5f, rs.getFloat("hi"));
        }

        assertThrows(java.sql.SQLException.class, () -> FragCastPredictor.toLibrarySchema(words, NO_FOOTER));
        assertFalse(new File(words.getAbsolutePath() + ".apd.parquet").exists());
    }

    // --- the library's footer --------------------------------------------------------------

    // A COPY carries no key-value metadata across, so a library the conversion rewrites would lose
    // the footer FragCast wrote into it - the enzyme FragDIA digests with among it. It is given the
    // same keys instead, from the same options, and nothing else: FragCast's own arrow schema goes.
    @Test
    public void aRewrittenLibraryCarriesTheEnzymeAndFasta() throws Exception {
        File lib = fragCastLibrary();
        assertTrue(ParquetFooters.read(lib).containsKey("ARROW:schema"), "the fixture changed");
        Map<String, String> footer = strictTrypsin().footer("C:\\fasta\\human.fasta");

        FragCastPredictor.toLibrarySchema(lib, footer);

        assertLibrarySchema(lib);
        assertEquals(footer, ParquetFooters.read(lib));
    }

    // With no enzyme and no FASTA there is nothing to record, and nothing is: the conversion writes
    // the library it always wrote.
    @Test
    public void aRewriteWithNothingToRecordRecordsNothing() throws Exception {
        File lib = fragCastLibrary();
        FragCastPredictor.toLibrarySchema(lib, NO_FOOTER);
        assertTrue(ParquetFooters.read(lib).isEmpty(), ParquetFooters.read(lib).toString());
    }

    // FragCast is told what to record exactly as this command line gave it, the FASTA after the
    // enzyme; told nothing when neither was given, so its command line is the one it always had.
    @Test
    public void fragCastIsHandedTheEnzymeAndTheFastaAsGiven() {
        EnzymeOptions trypsin = strictTrypsin();
        List<String> expected = new ArrayList<>(trypsin.builderArgs());
        expected.addAll(Arrays.asList("--fasta", "C:\\fasta\\human.fasta"));
        assertEquals(expected, FragCastPredictor.libraryArgs(trypsin, "C:\\fasta\\human.fasta"));

        assertEquals(trypsin.builderArgs(), FragCastPredictor.libraryArgs(trypsin, ""));
        assertEquals(Arrays.asList("--fasta", "/data/human.fasta"),
                FragCastPredictor.libraryArgs(new EnzymeOptions(), "/data/human.fasta"),
                "a FASTA is recorded with or without an enzyme");
        assertTrue(FragCastPredictor.libraryArgs(new EnzymeOptions(), "").isEmpty());
    }

    private static String path(File file) {
        return file.getAbsolutePath().replace("\\", "/");
    }

    // Nor is it enough for the schema to be in there somewhere. A library with a column more than the
    // schema, or one short of it, is not the file downstream expects, so it goes through the
    // conversion - which drops the extra one, and says so plainly about a missing one.
    @Test
    public void aColumnTooManyOrTooFewIsNotTheSchema() throws Exception {
        String narrow = new File("src/test/resources/fragcast_lib_narrow.parquet").getAbsolutePath().replace("\\", "/");
        File extra = tmp.resolve("extra.parquet").toFile();
        File shortOne = tmp.resolve("short.parquet").toFile();
        try (Connection conn = DriverManager.getConnection("jdbc:duckdb:");
             Statement stmt = conn.createStatement()) {
            stmt.execute("COPY (SELECT *, 'x' AS GeneName FROM read_parquet('" + narrow + "')) TO '" +
                    extra.getAbsolutePath().replace("\\", "/") + "' (FORMAT PARQUET)");
            stmt.execute("COPY (SELECT * EXCLUDE (Proteotypic) FROM read_parquet('" + narrow + "')) TO '" +
                    shortOne.getAbsolutePath().replace("\\", "/") + "' (FORMAT PARQUET)");
        }

        FragCastPredictor.toLibrarySchema(extra, NO_FOOTER);
        assertLibrarySchema(extra);

        java.sql.SQLException missing = org.junit.jupiter.api.Assertions.assertThrows(
                java.sql.SQLException.class, () -> FragCastPredictor.toLibrarySchema(shortOne, NO_FOOTER));
        assertTrue(missing.getMessage().contains("Proteotypic"), missing.getMessage());
        assertFalse(missing.getMessage().toLowerCase().contains("self"),
                "a missing column was reported as a self-reference: " + missing.getMessage());
        assertFalse(new File(shortOne.getAbsolutePath() + ".apd.parquet").exists());
    }

    // The right column names are not enough to skip it: what the conversion does to these columns
    // is cast them, so it is their types that say whether anything is left to do.
    @Test
    public void theRightColumnsInWiderTypesAreStillConverted() throws Exception {
        File wide = tmp.resolve("wide.parquet").toFile();
        List<String> names = new ArrayList<>();
        for (String[] column : LIBRARY_SCHEMA) {
            names.add(column[0]);
        }
        try (Connection conn = DriverManager.getConnection("jdbc:duckdb:");
             Statement stmt = conn.createStatement()) {
            //the schema's 17 columns by name, still as the doubles and ints an older FragCast wrote
            stmt.execute("COPY (SELECT " + String.join(", ", names) + " FROM read_parquet('" +
                    fragCastLibrary().getAbsolutePath().replace("\\", "/") + "')) TO '" +
                    wide.getAbsolutePath().replace("\\", "/") + "' (FORMAT PARQUET)");
        }
        long before = rowCount(wide);

        FragCastPredictor.toLibrarySchema(wide, NO_FOOTER);

        assertLibrarySchema(wide);
        assertEquals(before, rowCount(wide));
    }

    @Test
    public void theConversionKeepsEveryRow() throws Exception {
        File lib = fragCastLibrary();
        long before = rowCount(lib);
        FragCastPredictor.toLibrarySchema(lib, NO_FOOTER);
        assertEquals(before, rowCount(lib));
        assertTrue(before > 0);
    }

    @Test
    public void everyPropertyIsKept() throws Exception {
        File lib = fragCastLibrary();
        FragCastPredictor.toLibrarySchema(lib, NO_FOOTER);
        try (Connection conn = DriverManager.getConnection("jdbc:duckdb:");
             Statement stmt = conn.createStatement();
             java.sql.ResultSet rs = stmt.executeQuery(
                     "SELECT COUNT(NormalizedRetentionTime) AS rt, " +
                     "COUNT(PrecursorIonMobility) AS im FROM read_parquet('" +
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
