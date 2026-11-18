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

package koinaclasses;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertNull;
import static org.junit.jupiter.api.Assertions.assertThrows;
import static org.junit.jupiter.api.Assertions.assertTrue;

import java.io.ByteArrayOutputStream;
import java.io.PrintStream;
import java.nio.charset.StandardCharsets;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.LinkedHashSet;
import java.util.List;
import java.util.Set;
import java.util.TreeSet;
import koinaclasses.KoinaModelSpec.Input;
import koinaclasses.KoinaModelSpec.Kind;
import koinaclasses.KoinaModelSpec.Limits;
import org.junit.jupiter.api.Test;
import org.junit.jupiter.api.io.TempDir;
import peptideptmformatting.PeptideFormatter;
import peptideptmformatting.PeptideSkipper;
import readers.datareaders.PinReader;
import testutils.KoinaGolden;
import testutils.KoinaGolden.Case;
import testutils.KoinaGolden.Psm;

/**
 * Before a request is written, each PSM is checked against the model's spec: a PSM the model cannot take is not
 * sent, a mod it cannot encode is left out of the peptide sent, and both are counted for one summary per model. A
 * model that can take none of the PSMs ends the run (D28).
 */
class KoinaPreCheckTest {
    private static final String MS2_MODEL = "Prosit_2020_intensity_HCD";
    private static final String RT_MODEL = "Prosit_2019_irt";
    private static final String LONG = "ADGLRTEFLKADGLRTEFLKADGLRTEFLKAK"; //32 residues

    private static final List<Psm> PSMS = List.of(
            new Psm("K.LGEPTIDEK.A", 2),
            new Psm("K.LGEPTIDEK.A", 2), //the same precursor in another scan: one more PSM, no more precursor
            new Psm("R.GS[79.9663]PEDLHVTAR.L", 2),
            new Psm("R.GS[79.9663]PEDLHVTAR.L", 3),
            new Psm("R.S[79.9663]PEDLS[79.9663]HVTAR.L", 2), //two of the same mod: counted once for the PSM
            new Psm("K.n[42.0106]ADVTGLPHK.E", 2),
            new Psm("R.AVLDQGSTWPEM[15.9949]AGK.S", 2),
            new Psm("K.GDUVLPAEK.A", 2),
            new Psm("K." + LONG + ".T", 3),
            new Psm("K.LGEPTIDEK.A", 7),
            new Psm("K.TEC[57.0215]LDHFAVNPR.G", 3));

    @TempDir
    Path tmp;

    @Test
    void psmsAreCountedByWhyTheyAreNotSentAndByTheModsLeftOut() throws Exception {
        Checked checked = check(KoinaModelRegistry.require(MS2_MODEL), PSMS);
        KoinaPreCheck preCheck = checked.preCheck;

        assertEquals(11, preCheck.psms());
        assertEquals(3, preCheck.skippedPsms());
        //most first, then by name
        assertEquals("{charge above 6=1, longer than 30 residues=1, residue U=1}",
                preCheck.skippedByReason().toString());
        assertEquals("{S[79.9663]=3, n[42.0106]=1}", preCheck.droppedModPsms().toString());
        assertEquals(List.of("LGEPTIDEK,2", "GSPEDLHVTAR,2", "GSPEDLHVTAR,3", "SPEDLSHVTAR,2", "ADVTGLPHK,2",
                "AVLDQGSTWPEM[UNIMOD:35]AGK,2", "TEC[UNIMOD:4]LDHFAVNPR,3"), precursors(checked.hits));

        assertEquals(List.of(
                "[INFO] - " + MS2_MODEL + " is sent 8 of 11 PSMs, as 7 precursors",
                "[INFO] -   3 PSMs not sent, as the model cannot take them: charge above 6 (1), "
                        + "longer than 30 residues (1), residue U (1)",
                "[INFO] -   mods the model cannot encode, left out of the peptides sent (PSMs): S[79.9663] (3), "
                        + "n[42.0106] (1); the fragments that hold one are shifted by its mass"),
                report(preCheck, checked.hits.size()));
        preCheck.requirePredictable(checked.hits.size());
    }

    @Test
    void anRtModelIsPredictedAsIfTheModsLeftOutWereNotThere() throws Exception {
        Checked checked = check(KoinaModelRegistry.require(RT_MODEL), PSMS);
        List<String> report = report(checked.preCheck, checked.hits.size());
        assertEquals(3, report.size(), String.join("\n", report));
        assertTrue(report.get(2).endsWith(": S[79.9663] (3), n[42.0106] (1); predicted as if unmodified there"),
                report.get(2));
    }

    @Test
    void aModelThatIsSentEveryPsmWhollyReportsOneLine() throws Exception {
        Checked checked = check(KoinaModelRegistry.require(MS2_MODEL), PSMS.subList(0, 2));
        assertEquals(List.of("[INFO] - " + MS2_MODEL + " is sent 2 of 2 PSMs, as 1 precursor"),
                report(checked.preCheck, checked.hits.size()));
    }

    //D28
    @Test
    void aModelThatCanTakeNoneOfThePsmsEndsTheRunWithWhy() throws Exception {
        Checked checked = check(KoinaModelRegistry.require(MS2_MODEL),
                List.of(new Psm("K.GDUVLPAEK.A", 2), new Psm("K.LGEPTIDEK.A", 7), new Psm("K.GDUVLPAEK.A", 3)));
        assertTrue(checked.hits.isEmpty());
        KoinaException e = assertThrows(KoinaException.class,
                () -> checked.preCheck.requirePredictable(checked.hits.size()));
        assertEquals("Koina model " + MS2_MODEL + " cannot take any of the 3 PSMs: residue U (2), charge above 6 (1), "
                + "so it has nothing to predict. Choose another model", e.getMessage());

        KoinaPreCheck none = new KoinaPreCheck(KoinaModelRegistry.require(MS2_MODEL));
        e = assertThrows(KoinaException.class, () -> none.requirePredictable(0));
        assertEquals("Koina model " + MS2_MODEL + " has nothing to predict: there are no PSMs", e.getMessage());
    }

    //the Prosit-PTM models are sent every mod Koina's hand-written table or zero-shot encoding takes
    //(KoinaPrositPtm), the peptides of their golden run here; the others are left out and counted, C-term mods too
    @Test
    void thePrositPtmModelsLeaveOutTheModsKoinaCannotEncode() throws Exception {
        for (String model : List.of("Prosit_2025_intensity_22PTM", "Prosit_2025_irt_22PTM")) {
            List<Psm> psms = Case.load(model).psms;
            Checked checked = check(KoinaModelRegistry.require(model), psms);
            KoinaPreCheck preCheck = checked.preCheck;

            assertEquals(28, preCheck.psms(), model);
            assertEquals("{charge above 6=1, longer than 30 residues=1, residue U=1}",
                    preCheck.skippedByReason().toString(), model);
            assertEquals("{N[203.0794]=1, S[-18.0106]=1, T[146.0579]=1, c[-0.9840]=1, n[-17.0265]=1}",
                    preCheck.droppedModPsms().toString(), model);
            List<String> sent = precursors(checked.hits);
            assertTrue(sent.contains("[UNIMOD:737]-GS[UNIMOD:737]PEDLHVTAK[UNIMOD:737],2"), model + " " + sent);
            assertTrue(sent.contains("SHCVAELLGPEK,2"), model + ": a free cysteine is sent as a bare C");
            assertTrue(sent.contains("GSPEDLHVTAR,2"), model + ": the dehydration is left out");

            List<String> report = report(preCheck, checked.hits.size());
            assertEquals("[INFO] - " + model + " is sent 25 of 28 PSMs, as 24 precursors", report.get(0));
            assertEquals("[INFO] -   mods the model cannot encode, left out of the peptides sent (PSMs): "
                    + "N[203.0794] (1), S[-18.0106] (1), T[146.0579] (1), c[-0.9840] (1), n[-17.0265] (1); "
                    + (model.contains("irt") ? "predicted as if unmodified there"
                    : "the fragments that hold one are shifted by its mass"), report.get(2));
        }
    }

    @Test
    void limitsCanSkipPeptidesWithAFreeCysteine() {
        KoinaModelSpec spec = KoinaModelSpec.builder("NoFreeCys", Kind.MS2)
                .limits(Limits.NONE.withSkipFreeCys(true)).build();
        assertNull(PeptideSkipper.skipReason(new PeptideFormatter("TEC[57.0215]LDHFAVNPR", "3", "base"), spec));
        assertNull(PeptideSkipper.skipReason(new PeptideFormatter("LGEPTIDEK", "2", "base"), spec));
        assertEquals(Limits.FREE_CYS,
                PeptideSkipper.skipReason(new PeptideFormatter("TECLDHFAVNPR", "3", "base"), spec));
        assertEquals(Limits.FREE_CYS,
                PeptideSkipper.skipReason(new PeptideFormatter("C[57.0215]PEPTIDEC", "2", "base"), spec));
        assertThrows(IllegalStateException.class, () -> spec.limits().skipReason("TECLDHFAVNPR", 3, null),
                "limits that skip free cysteines need the peptide with its mods");

        //the models that read every C as carbamidomethylated, or take no other
        Set<String> skipping = new TreeSet<>();
        for (KoinaModelSpec model : KoinaModelRegistry.all()) {
            if (model.limits().skipFreeCys()) {
                skipping.add(model.name());
            }
        }
        assertEquals(Set.of("Altimeter_2024_intensities", "Prosit_2025_intensity_lac", "Prosit_2025_irt_lac"),
                skipping);
    }

    //the models of D12 part B, on the PSMs of their golden runs: what each is not sent, the mods left out, and the
    //precursors sent (PSMs of the same peptide and charge, or of the same peptide for a retention time model, are one)
    @Test
    void theLactylationChronologerAndAltimeterModelsSayWhatTheyLeaveOut() throws Exception {
        String[][] expected = {
                {"Prosit_2025_intensity_lac", "20", "{charge above 6=1, free cysteine=1, longer than 30 residues=1, "
                        + "residue U=1}", "{C[39.9950]=1, K[42.0106]=1, K[56.0262]=1, Q[-17.0265]=1, S[79.9663]=1, "
                        + "c[72.0211]=1, n[42.0106]=1}", "14"},
                {"Prosit_2025_irt_lac", "20", "{charge above 6=1, free cysteine=1, longer than 30 residues=1, "
                        + "residue U=1}", "{C[39.9950]=1, K[42.0106]=1, K[56.0262]=1, Q[-17.0265]=1, S[79.9663]=1, "
                        + "c[72.0211]=1, n[42.0106]=1}", "14"},
                {"Chronologer_RT", "30", "{longer than 50 residues=1, residue U=1}",
                        "{K[100.0160]=1, N[0.9840]=1, c[-0.9840]=1, c[42.0106]=1, n[-17.0265]=1}", "26"},
                {"Altimeter_2024_intensities", "18", "{charge above 7=1, free cysteine=1, longer than 40 residues=1, "
                        + "residue U=1, shorter than 6 residues=1}",
                        "{C[39.9950]=1, S[79.9663]=1, c[-0.9840]=1, n[42.0106]=1}", "12"}};
        for (String[] e : expected) {
            String model = e[0];
            List<Psm> psms = Case.load(model).psms;
            Checked checked = check(KoinaModelRegistry.require(model), psms);
            assertEquals(Long.parseLong(e[1]), checked.preCheck.psms(), model);
            assertEquals(e[2], checked.preCheck.skippedByReason().toString(), model);
            assertEquals(e[3], checked.preCheck.droppedModPsms().toString(), model);
            assertEquals(Integer.parseInt(e[4]), checked.hits.size(), model);
        }
    }

    @Test
    void theFragmentationIsSentAsTheSpecWritesIt() throws Exception {
        KoinaModelSpec spec = KoinaModelSpec.builder("ExactFragmentation", Kind.MS2)
                .inputs(Input.PRECURSOR_CHARGES, Input.FRAGMENTATION_TYPES)
                .fragmentation("HCD", "HCD_EXACT").fragmentations("CID").build();
        assertEquals("HCD_EXACT", spec.fragmentationSentAs("HCD"));
        assertEquals("CID", spec.fragmentationSentAs("CID"));
        assertEquals("ETD", spec.fragmentationSentAs("ETD"), "a type the spec does not name is sent as it is");

        Checked checked = check(spec, PSMS.subList(0, 1));
        assertEquals(List.of("LGEPTIDEK,2,27,LUMOS,HCD_EXACT,9"), new ArrayList<>(checked.hits));
        //the models of today send the run's fragmentation as it is
        for (KoinaModelSpec model : KoinaModelRegistry.all()) {
            assertEquals("HCD", model.fragmentationSentAs("HCD"), model.name());
        }
    }

    //what PinReader.createJSON made of PSMs for a model
    private static final class Checked {
        KoinaPreCheck preCheck;
        Set<String> hits = new LinkedHashSet<>();
    }

    private Checked check(KoinaModelSpec spec, List<Psm> psms) throws Exception {
        Path pin = KoinaGolden.writePin(psms, tmp.resolve(spec.name() + ".pin"));
        Checked checked = new Checked();
        checked.preCheck = new KoinaPreCheck(spec);
        KoinaGolden.Settings settings = KoinaGolden.Settings.apply(Case.of(spec.name(), "LUMOS", "HCD", "27", psms));
        try {
            PinReader reader = new PinReader(pin.toString());
            try {
                reader.createJSON(checked.hits, checked.preCheck);
            } finally {
                reader.close();
            }
        } finally {
            settings.restore();
        }
        return checked;
    }

    //the peptide and charge of each request entry
    private static List<String> precursors(Set<String> hits) {
        List<String> precursors = new ArrayList<>();
        for (String hit : hits) {
            String[] fields = hit.split(",");
            precursors.add(fields[0] + "," + fields[1]);
        }
        return precursors;
    }

    //the lines KoinaPreCheck.report prints, without their time
    private static List<String> report(KoinaPreCheck preCheck, int precursors) {
        PrintStream out = System.out;
        ByteArrayOutputStream printed = new ByteArrayOutputStream();
        try {
            System.setOut(new PrintStream(printed, true, StandardCharsets.UTF_8));
            preCheck.report(precursors);
        } finally {
            System.setOut(out);
        }
        List<String> lines = new ArrayList<>();
        for (String line : printed.toString(StandardCharsets.UTF_8).split("\\R")) {
            int at = line.indexOf("[INFO]");
            if (at >= 0) {
                lines.add(line.substring(at));
            }
        }
        System.out.print(printed.toString(StandardCharsets.UTF_8));
        return lines;
    }
}
