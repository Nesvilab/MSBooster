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

import java.io.File;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.sql.Connection;
import java.sql.DriverManager;
import java.sql.Statement;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;

import org.junit.jupiter.api.Test;
import org.junit.jupiter.api.io.TempDir;

// The gene map answers three different lookups, and only one of them asks about a whole protein
// group: the librarytsv join matches ProteinId, which FragCast writes as the FIRST protein of the
// group, and ParquetToSpecLib splits ProteinId and looks its parts up one at a time. A map holding
// group keys alone misses both, and a gene name degrades to the accession it was meant to replace.
public class HelpersTest {
    /** A peptide list in the form both predictors' inputs take: a proteins column, one row each. */
    private static File peptideList(Path dir, String... proteins) throws Exception {
        StringBuilder sb = new StringBuilder("peptide,proteins,is_decoy\n");
        for (int i = 0; i < proteins.length; i++) {
            sb.append("PEPTIDE").append((char) ('A' + i)).append("K,").append(proteins[i]).append(",false\n");
        }
        File file = dir.resolve("peptides.csv").toFile();
        Files.write(file.toPath(), sb.toString().getBytes(StandardCharsets.UTF_8));
        return file;
    }

    private static HashMap<String, String> fasta(String... protGenePairs) {
        HashMap<String, String> protToGene = new HashMap<>();
        for (int i = 0; i < protGenePairs.length; i += 2) {
            protToGene.put(protGenePairs[i], protGenePairs[i + 1]);
        }
        return protToGene;
    }

    // The histone case from a real run: every peptide of the group is shared, so no row ever names
    // P1 on its own and the only key was P1;P2. FragCast's ProteinId is P1, which matched nothing.
    @Test
    public void everyProteinOfASharedGroupIsAKeyOfItsOwn(@TempDir Path dir) throws Exception {
        HashMap<String, String> map = Helpers.mapProteinsListToGenes(
                peptideList(dir, "sp|P1|A_HUMAN;sp|P2|B_HUMAN").getAbsolutePath(),
                fasta("P1", "GENEA", "P2", "GENEB"));

        assertEquals("GENEA", map.get("P1"), "FragCast names ProteinId after the first protein only");
        assertEquals("GENEB", map.get("P2"), "the speclib writer looks the group up one protein at a time");
        assertEquals("GENEA;GENEB", map.get("P1;P2"), "the server writes the whole group as ProteinId");
    }

    // A peptide in one protein is the plain case and must keep answering as it always did.
    @Test
    public void aProteotypicPeptideStillMapsItsOneProtein(@TempDir Path dir) throws Exception {
        HashMap<String, String> map = Helpers.mapProteinsListToGenes(
                peptideList(dir, "sp|P1|A_HUMAN").getAbsolutePath(), fasta("P1", "GENEA"));

        assertEquals("GENEA", map.get("P1"));
        assertEquals(1, map.size(), "a single protein is not also a group: " + map);
    }

    // What a peptide states directly outranks what a group implies, whichever order they are read in.
    @Test
    public void aStandaloneProteinOutranksTheGroupsThatContainIt(@TempDir Path dir) throws Exception {
        for (String[] order : Arrays.asList(
                new String[] {"sp|P1|A_HUMAN", "sp|P1|A_HUMAN;sp|P2|B_HUMAN"},
                new String[] {"sp|P1|A_HUMAN;sp|P2|B_HUMAN", "sp|P1|A_HUMAN"})) {
            HashMap<String, String> map = Helpers.mapProteinsListToGenes(
                    peptideList(dir, order).getAbsolutePath(), fasta("P1", "GENEA", "P2", "GENEB"));
            assertEquals("GENEA", map.get("P1"), "read in the order " + Arrays.toString(order));
        }
    }

    // A protein the FASTA does not name keeps answering with its accession, as it did before - the
    // point is that the lookup now finds an entry at all, not that it invents a gene.
    @Test
    public void aProteinMissingFromTheFastaAnswersWithItsAccession(@TempDir Path dir) throws Exception {
        HashMap<String, String> map = Helpers.mapProteinsListToGenes(
                peptideList(dir, "sp|P1|A_HUMAN;sp|P9|UNKNOWN_HUMAN").getAbsolutePath(),
                fasta("P1", "GENEA"));

        assertEquals("GENEA", map.get("P1"));
        assertEquals("P9", map.get("P9"));
    }

    // Both predictors label a decoy by prefixing the accession they read out of its header, so the
    // map has to be keyed that way or every decoy's gene degrades to the accession. A peptide sitting
    // in ONE decoy is the case a group's own handling never reaches.
    @Test
    public void aDecoyOnItsOwnIsKeyedAndNamedLikeOneInAGroup(@TempDir Path dir) throws Exception {
        HashMap<String, String> alone = Helpers.mapProteinsListToGenes(
                peptideList(dir, "rev_sp|P1|A_HUMAN").getAbsolutePath(), fasta("P1", "GENEA"));
        assertEquals("rev_GENEA", alone.get("rev_P1"));

        HashMap<String, String> grouped = Helpers.mapProteinsListToGenes(
                peptideList(dir, "rev_sp|P1|A_HUMAN;sp|P2|B_HUMAN").getAbsolutePath(),
                fasta("P1", "GENEA", "P2", "GENEB"));
        assertEquals(alone.get("rev_P1"), grouped.get("rev_P1"),
                "one protein and a group of one protein must answer the same way");
    }

    // A label with no database token is already its own accession, prefix and all. Putting the prefix
    // back on it gives rev_rev_P1 - a key nothing matches, which is the very failure the prefix
    // handling exists to prevent. FASTAs without pipes in their headers are the ones at risk.
    @Test
    public void aDecoyLabelWithNoDatabaseTokenIsNotPrefixedTwice(@TempDir Path dir) throws Exception {
        HashMap<String, String> map = Helpers.mapProteinsListToGenes(
                peptideList(dir, "rev_P1").getAbsolutePath(), fasta("P1", "GENEA"));

        assertEquals(new HashSet<>(Arrays.asList("rev_P1")), map.keySet());
        assertEquals("rev_P1", map.get("rev_P1"));
    }

    // A decoy the FASTA does not name answers with its prefixed accession, never with the target's.
    @Test
    public void aDecoyMissingFromTheFastaKeepsItsPrefix(@TempDir Path dir) throws Exception {
        HashMap<String, String> map = Helpers.mapProteinsListToGenes(
                peptideList(dir, "rev_sp|P9|UNKNOWN_HUMAN").getAbsolutePath(), fasta("P1", "GENEA"));
        assertEquals("rev_P9", map.get("rev_P9"));
    }

    // The list FragCast is given is a TSV whose protein column is full of pipes, and a short one
    // gives a delimiter sniffer more pipes than tabs to go on. Guessing '|' makes the whole run die
    // on "Referenced column proteins not found" - and only for small jobs, which is the worst
    // possible size to fail at. The extension already says what the delimiter is.
    @Test
    public void aShortPeptideListIsReadAtTheDelimiterItsNameImplies(@TempDir Path dir) throws Exception {
        File tsv = dir.resolve("peptides.tsv").toFile();
        Files.write(tsv.toPath(), ("peptide\tcharge\tproteins\n"
                + "PEPTIDEK\t2\trev_sp|P1|A_HUMAN\n"
                + "SAMPLERK\t2\tsp|P1|A_HUMAN\n"
                + "ANOTHERR\t2\trev_sp|P2|B_HUMAN;sp|P3|C_HUMAN\n").getBytes(StandardCharsets.UTF_8));

        HashMap<String, String> map = Helpers.mapProteinsListToGenes(
                tsv.getAbsolutePath(), fasta("P1", "GENEA"));

        assertEquals("GENEA", map.get("P1"));
    }

    // --- the written library.tsv --------------------------------------------------------------

    /** The 19 columns FragSpecLib writes, in its order - what a predicted library must line up with. */
    private static final String FRAGSPECLIB_HEADER = String.join("\t",
            "PrecursorMz", "ProductMz", "Annotation", "ProteinId", "GeneName", "PeptideSequence",
            "ModifiedPeptideSequence", "PrecursorCharge", "LibraryIntensity",
            "NormalizedRetentionTime", "PrecursorIonMobility", "FragmentType", "FragmentCharge",
            "FragmentSeriesNumber", "FragmentLossType", "AverageExperimentalRetentionTime",
            "AllMappedProteins", "AllMappedGenes", "Proteotypic");

    /**
     * One library row in the schema {@code FragCastPredictor.toLibrarySchema} leaves behind, or the
     * server's own when {@code group} is null - the difference being the two columns only FragCast
     * fills.
     */
    private static File library(Path dir, String proteinId, String group) throws Exception {
        File parquet = dir.resolve(group == null ? "server.parquet" : "fragcast.parquet").toFile();
        String extra = group == null ? ""
                : "'' AS AverageExperimentalRetentionTime, '" + group + "' AS AllMappedProteins, ";
        try (Connection conn = DriverManager.getConnection("jdbc:duckdb:");
             Statement stmt = conn.createStatement()) {
            stmt.execute("COPY (SELECT " +
                    "CAST(464.73 AS FLOAT) AS PrecursorMz, CAST(702.36 AS FLOAT) AS ProductMz, " +
                    "'y6^1' AS Annotation, '" + proteinId + "' AS ProteinId, " +
                    "'PEPTIDEK' AS PeptideSequence, 'PEPTIDEK' AS ModifiedPeptideSequence, " +
                    "CAST(2 AS SMALLINT) AS PrecursorCharge, CAST(10000 AS FLOAT) AS LibraryIntensity, " +
                    "CAST(-17.75 AS FLOAT) AS NormalizedRetentionTime, " +
                    "CAST('' AS VARCHAR) AS PrecursorIonMobility, 'y' AS FragmentType, " +
                    "CAST(1 AS SMALLINT) AS FragmentCharge, CAST(6 AS SMALLINT) AS FragmentSeriesNumber, " +
                    "'' AS FragmentLossType, " + extra +
                    "CAST(0 AS SMALLINT) AS Proteotypic) TO '" +
                    parquet.getAbsolutePath().replace("\\", "/") + "' (FORMAT PARQUET)");
        }
        return parquet;
    }

    private static String[] writeAndRead(Path dir, File parquet, HashMap<String, String> map,
                                         List<String> header) throws Exception {
        String tsv = dir.resolve(parquet.getName() + ".tsv").toString();
        try (Connection conn = DriverManager.getConnection("jdbc:duckdb:")) {
            Helpers.convertParquetToLibraryTsv(parquet.getAbsolutePath(), tsv, map, conn);
        }
        List<String> lines = Files.readAllLines(Paths.get(tsv), StandardCharsets.UTF_8);
        header.add(lines.get(0));
        return lines.get(1).split("\t", -1);
    }

    // A FragCast library must line up with the experimental one column for column, and both of its
    // gene columns must hold GENES: FragCast reads its own off the protein header, so it would write
    // the entry mnemonics A;B where the genes are GENEA;GENEB.
    @Test
    public void aFragCastLibraryKeepsItsGroupAndGetsRealGeneNames(@TempDir Path dir) throws Exception {
        String group = "sp|P1|A_HUMAN;sp|P2|B_HUMAN";
        HashMap<String, String> map = Helpers.mapProteinsListToGenes(
                peptideList(dir, group).getAbsolutePath(), fasta("P1", "GENEA", "P2", "GENEB"));

        List<String> header = new ArrayList<>();
        String[] row = writeAndRead(dir, library(dir, "P1", group), map, header);

        assertEquals(FRAGSPECLIB_HEADER, header.get(0));
        assertEquals(19, row.length);
        assertEquals("P1", row[3], "ProteinId names the leading protein, as FragSpecLib's does");
        assertEquals("GENEA", row[4]);
        assertEquals(group, row[16], "the group is the only place a shared peptide's other proteins live");
        assertEquals("GENEA;GENEB", row[17]);
        assertEquals("0", row[18]);
    }

    // A protein group the map has no entry for still answers with accessions rather than mnemonics,
    // and never with the pipe-delimited headers.
    @Test
    public void anUnmappedGroupFallsBackToItsAccessions(@TempDir Path dir) throws Exception {
        List<String> header = new ArrayList<>();
        String[] row = writeAndRead(dir, library(dir, "P1", "sp|P1|A_HUMAN;sp|P2|B_HUMAN"),
                new HashMap<>(), header);

        assertEquals("P1;P2", row[17]);
        assertEquals("P1", row[4], "GeneName falls back the same way it always did");
    }

    // The server returns no group column, and none may be invented for it - its library keeps the 16
    // columns it has always had.
    @Test
    public void aServerLibraryGainsNoGroupColumns(@TempDir Path dir) throws Exception {
        List<String> header = new ArrayList<>();
        String[] row = writeAndRead(dir, library(dir, "P1", null), fasta("P1", "GENEA"), header);

        assertEquals(16, row.length);
        assertFalse(header.get(0).contains("AllMapped"), header.get(0));
        assertEquals("GENEA", row[4]);
    }
}
