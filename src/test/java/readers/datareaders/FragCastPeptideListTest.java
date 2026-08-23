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

package readers.datareaders;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertFalse;
import static org.junit.jupiter.api.Assertions.assertTrue;

import java.io.File;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashSet;
import java.util.List;
import java.util.Set;
import java.util.concurrent.ConcurrentHashMap;

import allconstants.FragCastCharges;
import org.junit.jupiter.api.AfterEach;
import org.junit.jupiter.api.BeforeEach;
import org.junit.jupiter.api.Test;
import org.junit.jupiter.api.io.TempDir;

// The peptide list decides how large the predicted library is. FragCast expands only the rows that
// omit a charge, so naming the charge per row is what keeps the library to the precursors the search
// actually found - and a named charge outside 1-6 is a hard error that fails the whole prediction,
// not a dropped row.
public class FragCastPeptideListTest {
    private static final String TAB = "\t";

    @BeforeEach
    public void setUp() {
        FragCastCharges.reset();
    }

    @AfterEach
    public void tearDown() {
        FragCastCharges.reset();
    }

    /** A minimal pin holding one row per (peptide, charge, scan) triple. */
    private static File pin(Path dir, String... specs) throws Exception {
        List<String> lines = new ArrayList<>();
        lines.add(String.join(TAB, "SpecId", "Label", "ScanNr", "rank", "Peptide",
                "retentiontime", "log10_evalue"));
        int scan = 1000;
        for (String spec : specs) {
            String[] parts = spec.split("@");
            //SpecId's last period-delimited field is the charge, which is where getPep reads it from
            lines.add(String.join(TAB, "run.1." + scan + "." + parts[1], "1", String.valueOf(scan),
                    "1", "K." + parts[0] + ".R", "30.0", "-8.0"));
            scan++;
        }
        File file = dir.resolve("test.pin").toFile();
        Files.write(file.toPath(), (String.join(System.lineSeparator(), lines)
                + System.lineSeparator()).getBytes(StandardCharsets.UTF_8));
        return file;
    }

    private static Set<String> collect(File pinFile) throws Exception {
        Set<String> hSetHits = ConcurrentHashMap.newKeySet();
        PinReader pin = new PinReader(pinFile.getCanonicalPath());
        try {
            pin.createFragCastList(hSetHits);
        } finally {
            pin.close();
        }
        return hSetHits;
    }

    // Every row names its charge, so FragCast predicts these three precursors and no others. Written
    // without the charge, the same peptides would be expanded across a charge range instead.
    @Test
    public void everyPrecursorCarriesItsOwnCharge(@TempDir Path dir) throws Exception {
        Set<String> hits = collect(pin(dir, "PEPTIDEK@2", "PEPTIDEK@3", "AFTFVAKGGR@4"));
        assertEquals(new HashSet<>(Arrays.asList(
                "PEPTIDEK" + TAB + "2", "PEPTIDEK" + TAB + "3", "AFTFVAKGGR" + TAB + "4")), hits);
    }

    // Nine PSMs of three precursors are three lines: the list is deduplicated on the precursor, so
    // nothing repeats and nothing is predicted twice.
    @Test
    public void repeatedPsmsOfOnePrecursorWriteOneLine(@TempDir Path dir) throws Exception {
        Set<String> hits = collect(pin(dir,
                "PEPTIDEK@2", "PEPTIDEK@3", "PEPTIDEK@4",
                "PEPTIDEK@2", "PEPTIDEK@3", "PEPTIDEK@4",
                "PEPTIDEK@2", "PEPTIDEK@3", "PEPTIDEK@4"));
        assertEquals(3, hits.size(), "the peptide list repeats itself: " + hits);
    }

    // Naming a charge FragCast cannot represent is a hard error that fails the whole prediction, so
    // those precursors have to be left out of the list entirely.
    @Test
    public void chargesFragCastCannotRepresentAreLeftOut(@TempDir Path dir) throws Exception {
        Set<String> hits = collect(pin(dir, "PEPTIDEK@2", "BIGCHARGEPK@8", "ANOTHERPEPR@7"));
        assertEquals(1, hits.size(), "an unpredictable charge reached the list: " + hits);
        assertTrue(hits.contains("PEPTIDEK" + TAB + "2"));
        assertEquals(2, FragCastCharges.unrepresentableCount(),
                "the precursors left out were not counted, so the run cannot say why");
    }

    // Modified and unmodified forms are different peptides to FragCast, so they must not collapse -
    // that would silently drop one of them from the predictions entirely.
    @Test
    public void modifiedFormsAreNotCollapsedIntoTheUnmodifiedOne(@TempDir Path dir) throws Exception {
        Set<String> hits = collect(pin(dir,
                "PEPTC[57.0215]IDEK@2", "PEPTC[57.0215]IDEK@2", "PEPTCIDEK@2"));
        assertEquals(2, hits.size(), "a modified form was merged with the unmodified one: " + hits);
        assertTrue(hits.contains("PEPTCIDEK" + TAB + "2"));
    }

    // Two pin files of one job share the collector, so a precursor found in both is still one line.
    @Test
    public void aprecursorSeenInTwoFilesIsStillOneLine(@TempDir Path dir) throws Exception {
        Set<String> hSetHits = ConcurrentHashMap.newKeySet();
        for (String name : new String[] {"a", "b"}) {
            Path sub = dir.resolve(name);
            sub.toFile().mkdirs();
            PinReader pin = new PinReader(pin(sub, "PEPTIDEK@2", "SHAREDPEPK@3").getCanonicalPath());
            try {
                pin.createFragCastList(hSetHits);
            } finally {
                pin.close();
            }
        }
        assertEquals(2, hSetHits.size(), "the same precursors from two files were written twice: " + hSetHits);
    }

    // FragCast refuses a peptide field holding a comma, space, semicolon, pipe or quote - it means
    // the wrong delimiter or a protein accession. MSBooster's base format can hold none of them.
    @Test
    public void thepeptidesWrittenHoldNothingFragCastRefuses(@TempDir Path dir) throws Exception {
        Set<String> hits = collect(pin(dir, "PEPTC[57.0215]IDEK@2", "n[42.0106]PEPTIDEK@3"));
        assertFalse(hits.isEmpty());
        for (String hit : hits) {
            String peptide = hit.split(TAB, -1)[0];
            for (char c : new char[] {',', ' ', ';', '|', '"'}) {
                assertFalse(peptide.indexOf(c) >= 0,
                        "peptide " + peptide + " holds " + c + ", which FragCast rejects");
            }
        }
    }
    // --- the proteins column, which is what lets FragCast label a library from pin files ---------

    /**
     * A pin that carries protein labels, written exactly as MSFragger writes them: every entry is
     * TERMINATED by a ';', so even a peptide in a single protein reads {@code sp|P12345|AAA_HUMAN;}.
     * That terminator is the whole point of the fixture - both predictors call a peptide
     * non-proteotypic when this string holds a ';', so passing the field through unchanged returned
     * a library with {@code Proteotypic = 0} on every row.
     */
    private static File pinWithProteins(Path dir, String proteins) throws Exception {
        List<String> lines = new ArrayList<>();
        lines.add(String.join(TAB, "SpecId", "Label", "ScanNr", "rank", "Peptide",
                "retentiontime", "log10_evalue", "Proteins"));
        lines.add(String.join(TAB, "run.1.1000.2", "1", "1000", "1", "K.PEPTIDEK.R", "30.0", "-8.0",
                proteins));
        File file = dir.resolve("prot.pin").toFile();
        Files.write(file.toPath(), (String.join(System.lineSeparator(), lines)
                + System.lineSeparator()).getBytes(StandardCharsets.UTF_8));
        return file;
    }

    private static File pinWithProteins(Path dir) throws Exception {
        return pinWithProteins(dir, "sp|P12345|AAA_HUMAN;");
    }

    @Test
    public void rescoringInputKeepsItsTwoColumns(@TempDir Path dir) throws Exception {
        // rescoring never reads a protein, so its input is left exactly as it was
        allconstants.Constants.createPredFileOnly = false;
        try {
            assertEquals(new HashSet<>(Arrays.asList("PEPTIDEK" + TAB + "2")),
                    collect(pinWithProteins(dir)));
        } finally {
            allconstants.Constants.createPredFileOnly = false;
        }
    }

    @Test
    public void libraryInputCarriesTheProteins(@TempDir Path dir) throws Exception {
        // FragCast fills ProteinId, GeneName and Proteotypic from this column, which is the only way
        // a library predicted from pin files can be annotated. MSFragger's terminating ';' is not
        // part of any protein name, and leaving it on made every peptide read as shared.
        allconstants.Constants.createPredFileOnly = true;
        try {
            assertEquals(new HashSet<>(Arrays.asList(
                            "PEPTIDEK" + TAB + "2" + TAB + "sp|P12345|AAA_HUMAN")),
                    collect(pinWithProteins(dir)));
        } finally {
            allconstants.Constants.createPredFileOnly = false;
        }
    }

    // A peptide in more than one protein still reads as shared: the separators BETWEEN entries are
    // what says so, and only the empty entries around them are dropped.
    @Test
    public void aSharedPeptideKeepsEverySeparatorBetweenItsProteins(@TempDir Path dir) throws Exception {
        allconstants.Constants.createPredFileOnly = true;
        try {
            assertEquals(new HashSet<>(Arrays.asList(
                            "PEPTIDEK" + TAB + "2" + TAB + "sp|P1|A_HUMAN;sp|P2|B_HUMAN")),
                    collect(pinWithProteins(dir, "sp|P1|A_HUMAN;sp|P2|B_HUMAN;")));
        } finally {
            allconstants.Constants.createPredFileOnly = false;
        }
    }

    // A pin may spread its protein IDs across extra tab-separated columns the header does not name.
    // Reading only the named one would drop every protein after the first, which calls a shared
    // peptide proteotypic - the same defect as the terminator, in the other direction.
    @Test
    public void proteinsSpreadAcrossExtraColumnsAreAllKept(@TempDir Path dir) throws Exception {
        allconstants.Constants.createPredFileOnly = true;
        try {
            assertEquals(new HashSet<>(Arrays.asList(
                            "PEPTIDEK" + TAB + "2" + TAB + "sp|P1|A_HUMAN;sp|P2|B_HUMAN")),
                    collect(pinWithProteins(dir, "sp|P1|A_HUMAN" + TAB + "sp|P2|B_HUMAN")));
        } finally {
            allconstants.Constants.createPredFileOnly = false;
        }
    }

    @Test
    public void aPinWithNoProteinColumnStillPredicts(@TempDir Path dir) throws Exception {
        // an unlabelled library beats a crash: getColumn would index at -1 and throw
        allconstants.Constants.createPredFileOnly = true;
        try {
            assertEquals(new HashSet<>(Arrays.asList("PEPTIDEK" + TAB + "2" + TAB + "")),
                    collect(pin(dir, "PEPTIDEK@2")));
        } finally {
            allconstants.Constants.createPredFileOnly = false;
        }
    }
}
