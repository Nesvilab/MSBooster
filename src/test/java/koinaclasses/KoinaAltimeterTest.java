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

import static org.junit.jupiter.api.Assertions.assertArrayEquals;
import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertFalse;
import static org.junit.jupiter.api.Assertions.assertNull;
import static org.junit.jupiter.api.Assertions.assertSame;
import static org.junit.jupiter.api.Assertions.assertThrows;
import static org.junit.jupiter.api.Assertions.assertTrue;

import allconstants.Constants;
import com.google.gson.JsonArray;
import com.google.gson.JsonElement;
import com.google.gson.JsonObject;
import com.google.gson.JsonParser;
import java.util.ArrayList;
import java.util.LinkedHashSet;
import java.util.List;
import java.util.Set;
import java.util.TreeSet;
import koinaclasses.KoinaModelSpec.AnnotationParser;
import koinaclasses.KoinaModelSpec.Input;
import koinaclasses.KoinaModelSpec.Limits;
import org.junit.jupiter.api.Test;
import peptideptmformatting.PeptideFormatter;
import peptideptmformatting.PeptideSkipper;
import predictions.PredictionEntry;
import predictions.PredictionEntryHashMap;
import testutils.KoinaGolden;
import testutils.KoinaGolden.Case;

/**
 * Altimeter_2024_intensities: oxidized M and carbamidomethylated C only, every C carbamidomethylated, 6 to 40
 * residues and charges 1 to 7 (models/Altimeter/Altimeter_2024_preprocess_sequence in Koina, commit d176169). It
 * answers 200 fragments per precursor with their m/z and mzPAF labels: b and y ions, the precursor and immonium ions.
 */
class KoinaAltimeterTest {
    private static final String MODEL = "Altimeter_2024_intensities";
    private static final KoinaModelSpec ALTIMETER = KoinaModelRegistry.require(MODEL);

    @Test
    void theSpecKeepsKoinasMzAndReadsItsLabels() {
        assertEquals(KoinaModelSpec.Kind.MS2, ALTIMETER.kind());
        assertEquals(List.of(Input.PRECURSOR_CHARGES, Input.COLLISION_ENERGIES), new ArrayList<>(ALTIMETER.inputs()));
        assertTrue(ALTIMETER.calibrateCe());
        assertEquals(20, ALTIMETER.minCe());
        assertEquals(40, ALTIMETER.maxCe());
        assertEquals(List.of("intensities", "mz", "annotations"), ALTIMETER.outputs());
        assertEquals(AnnotationParser.ALTIMETER, ALTIMETER.annotationParser());
        assertTrue(ALTIMETER.trustKoinaMz());
        assertEquals(100, ALTIMETER.batchSize());
        assertTrue(ALTIMETER.allowedFragmentationTypes().isEmpty(), "it is sent no fragmentation");
        assertTrue(ALTIMETER.peptideFormat().camOnBareCys());
        assertEquals(Set.of("C4", "M35"), ALTIMETER.modRule().allowedMods());
        assertEquals(MODEL, ALTIMETER.citationKey());
    }

    @Test
    void theModelTakes6To40ResiduesChargesUpTo7AndNoFreeCysteine() {
        assertEquals("shorter than 6 residues", skip("ADGLR", 1));
        assertNull(skip("ADGLRK", 1));
        String forty = "ELVSAPLGTNEIYDHFAQWTGLPVEDAKLGRSTEPVNGHA";
        assertNull(skip(forty, 4));
        assertEquals("longer than 40 residues", skip(forty + "K", 4));
        assertNull(skip("TEALDHFAVNPRGSLEWQKPEDVTAR", 7));
        assertEquals("charge above 7", skip("TEALDHFAVNPRGSLEWQKPEDVTAR", 8));
        assertEquals("charge below 1", skip("TEALDHFAVNPRGSLEWQKPEDVTAR", 0));
        assertEquals("residue U", skip("GDUVLPAEK", 2));
        assertEquals(Limits.FREE_CYS, skip("SHCVAELLGPEK", 2));
        //a C-term mod is written after the last residue, and does not make a last C a modified one
        assertEquals(Limits.FREE_CYS, skip("LGEPTIDEC[-0.9840]cterm", 2));
        assertNull(skip("LGEPTIDEC[57.0215][-0.9840]cterm", 2));
        assertNull(skip("TEC[57.0215]LDHFAVNPR", 3));
        assertNull(skip("C[39.9950]ADLSGFTPEVIAR", 2));
    }

    //pin peptide, the peptide the model is sent, the mods left out of it, and the peptide read back
    private static final String[][] ROUND_TRIPS = {
            {"K.LGEPTIDEK.A", "LGEPTIDEK", "", "LGEPTIDEK"},
            {"R.AVLDQGSTWPEM[15.9949]AGK.S", "AVLDQGSTWPEM[UNIMOD:35]AGK", "", "AVLDQGSTWPEM[15.9949]AGK"},
            {"K.TEC[57.0215]LDHFAVNPR.G", "TEC[UNIMOD:4]LDHFAVNPR", "", "TEC[57.0215]LDHFAVNPR"},
            //sent as carbamidomethylated: the -17.0265 of the cyclization is left out, and added back to the m/z of
            //the fragments that hold the cysteine
            {"K.C[39.9950]ADLSGFTPEVIAR.Q", "C[UNIMOD:4]ADLSGFTPEVIAR", "C[39.9950]", "C[57.0215]ADLSGFTPEVIAR"},
            {"K.n[42.0106]ADVTGLPHK.E", "ADVTGLPHK", "n[42.0106]", "ADVTGLPHK"},
            {"R.GS[79.9663]PEDLHVTAR.L", "GSPEDLHVTAR", "S[79.9663]", "GSPEDLHVTAR"},
            {"K.PEPTM[31.9898]IDEK.A", "PEPTMIDEK", "M[31.9898]", "PEPTMIDEK"},
            {"K.LGEPTIDEKc[-0.9840].A", "LGEPTIDEK", "c[-0.9840]", "LGEPTIDEK"}};

    @Test
    void peptidesGoToTheModelAndBack() {
        for (String[] trip : ROUND_TRIPS) {
            String what = trip[0];
            PeptideFormatter pf = new PeptideFormatter(trip[0], 2, "pin");
            assertNull(PeptideSkipper.skipReason(pf, ALTIMETER), what);
            List<String> dropped = new ArrayList<>();
            String sent = pf.getModel(ALTIMETER, dropped);
            assertEquals(trip[1], sent, what);
            assertEquals(trip[2], String.join(",", dropped), what);
            PeptideFormatter back = new PeptideFormatter(sent, 2, ALTIMETER);
            assertEquals(trip[3], back.getBase(), what);
            assertEquals(sent, back.getModel(ALTIMETER), what);
        }
    }

    @Test
    void theParserReadsEachKindOfLabel() {
        assertParsed("y5", 1, "y", 5);
        assertParsed("b12", 1, "b", 12);
        assertParsed("y12^3", 3, "y", 12);
        assertParsed("b3^7", 7, "b", 3);
        //a loss on a b or y ion makes a b-NL or y-NL, as in UniSpec's labels
        assertParsed("y5-H2O", 1, "y-NL", 5);
        assertParsed("b4-NH3^2", 2, "b-NL", 4);
        //the precursor is p, with or without a loss
        assertParsed("p", 1, "p", 0);
        assertParsed("p^4", 4, "p", 0);
        assertParsed("p-H2O^2", 2, "p", 0);
        //an immonium ion is imm, with its mod, losses and gains
        assertParsed("IY", 1, "imm", 0);
        assertParsed("IC(Carbamidomethyl)", 1, "imm", 0);
        assertParsed("IM(Oxidation)", 1, "imm", 0);
        assertParsed("IY-NH3", 1, "imm", 0);
        assertParsed("IH-NH3+CO", 1, "imm", 0);
        assertParsed("IH+2CO", 1, "imm", 0);
        assertParsed("IR-N3H3+O2", 1, "imm", 0);
        assertParsed("IC(Carbamidomethyl)-NH3", 1, "imm", 0);
    }

    //a label of another kind, or one written otherwise, fails the response rather than being read as a wrong ion
    @Test
    void anUnknownLabelFailsLoudly() {
        for (String label : new String[] {"a3", "c3^2", "x2", "z5", "y3+1", "y5+i", "y5^2+i", "b", "y", "p3",
                "p^", "y5^", "Int/AE/3", "IMA", "I", "iy", "IY+", "IC(Carbamidomethyl", "y-H2O", "Y5", ""}) {
            IllegalArgumentException e = assertThrows(IllegalArgumentException.class,
                    () -> AnnotationParser.ALTIMETER.parse(label, new ArrayList<>(), new ArrayList<>(),
                            new ArrayList<>()), label);
            assertTrue(e.getMessage().contains("\"" + label + "\""), e.getMessage());
            assertTrue(e.getMessage().contains("Altimeter"), e.getMessage());
        }
    }

    /**
     * Every label of the recorded golden response (src/test/resources/koina/Altimeter_2024_intensities) is read: its
     * b and y ions up to charge 7, its precursors and its immonium ions with their mods, losses and gains. The
     * recording has no b or y ion with a loss, which the labels above cover. A label the parser does not know still
     * fails beside them.
     */
    @Test
    void theParserReadsEveryLabelKoinaAnswered() throws Exception {
        JsonObject response = JsonParser.parseString(Case.load(MODEL).read(KoinaGolden.RESPONSE)).getAsJsonObject();
        JsonArray annotations = null;
        for (JsonElement output : response.getAsJsonArray("outputs")) {
            if (output.getAsJsonObject().get("name").getAsString().equals("annotations")) {
                annotations = output.getAsJsonObject().getAsJsonArray("data");
            }
        }
        assertTrue(annotations != null && annotations.size() == 12 * 200, "12 precursors of 200 fragments");

        Set<String> labels = new LinkedHashSet<>();
        for (JsonElement annotation : annotations) {
            if (!annotation.getAsString().isEmpty()) {
                labels.add(annotation.getAsString());
            }
        }
        Set<String> ionTypes = new TreeSet<>();
        Set<Integer> charges = new TreeSet<>();
        for (String label : labels) {
            List<Integer> labelCharges = new ArrayList<>();
            List<String> labelIonTypes = new ArrayList<>();
            List<Integer> fragNums = new ArrayList<>();
            AnnotationParser.ALTIMETER.parse(label, labelCharges, labelIonTypes, fragNums);
            assertEquals(1, labelIonTypes.size(), label);
            String ionType = labelIonTypes.get(0);
            ionTypes.add(ionType);
            charges.add(labelCharges.get(0));
            assertEquals(label.startsWith("I") ? "imm" : label.substring(0, 1), ionType, label);
            assertEquals(ionType.equals("b") || ionType.equals("y"), fragNums.get(0) > 0, label);
        }
        assertEquals(Set.of("b", "imm", "p", "y"), ionTypes);
        assertEquals(Set.of(1, 2, 3, 4, 5, 6, 7), charges);
        for (String label : new String[] {"IC(Carbamidomethyl)", "IM(Oxidation)", "IR-N3H3+O2", "IK+CO2+H2N2", "p",
                "p^7", "b12^5", "y1"}) {
            assertTrue(labels.contains(label), label + " is in the recorded response");
        }

        IllegalArgumentException e = assertThrows(IllegalArgumentException.class,
                () -> AnnotationParser.ALTIMETER.parse("a2", new ArrayList<>(), new ArrayList<>(), new ArrayList<>()));
        assertTrue(e.getMessage().contains("\"a2\""), e.getMessage());
    }

    private static void assertParsed(String label, int charge, String ionType, int fragNum) {
        List<Integer> charges = new ArrayList<>();
        List<String> ionTypes = new ArrayList<>();
        List<Integer> fragNums = new ArrayList<>();
        AnnotationParser.ALTIMETER.parse(label, charges, ionTypes, fragNums);
        assertEquals(List.of(charge), charges, label);
        assertEquals(List.of(ionType), ionTypes, label);
        assertEquals(List.of(fragNum), fragNums, label);
    }

    private static String skip(String base, int charge) {
        return PeptideSkipper.skipReason(new PeptideFormatter(base, charge, "base"), ALTIMETER);
    }

    /**
     * Altimeter as the auxiliary spectra model, filling in the spectrum of a PSM the main model has none for: like
     * UniSpec's, its m/z are exact, so the spectrum is matched in ppm, and only monoisotopic peaks are kept. Its
     * entries, read back from the MGF file, carry no isotope numbers, as it predicts no isotope peak, so all their
     * peaks are kept. A UniSpec entry with isotope peaks keeps its monoisotopic ones; a model whose m/z are not
     * kept fills it in to be matched in Daltons.
     */
    @Test
    void asAuxiliaryModelItFillsInSpectraMatchedInPpm() {
        String saved = Constants.auxSpectraModel;
        try {
            Constants.auxSpectraModel = MODEL;
            PredictionEntry altimeter = merged(aux(new int[0]));
            assertEquals(3, altimeter.numFragments());
            assertArrayEquals(new float[] {300f, 400f, 500f}, altimeter.getMzs());
            assertFalse(altimeter.daltonMatching);

            Constants.auxSpectraModel = "UniSpec";
            PredictionEntry unispec = merged(aux(new int[] {0, 1, 0}));
            assertArrayEquals(new float[] {300f, 500f}, unispec.getMzs());
            assertArrayEquals(new int[] {0, 0}, unispec.isotopes);
            assertFalse(unispec.daltonMatching);
            PredictionEntry unispecMonoisotopic = merged(aux(new int[0]));
            assertEquals(3, unispecMonoisotopic.numFragments(), "an entry without isotope peaks keeps its peaks");
            assertFalse(unispecMonoisotopic.daltonMatching);

            Constants.auxSpectraModel = "Prosit_2020_intensity_HCD";
            PredictionEntry prosit = merged(aux(new int[0]));
            assertEquals(3, prosit.numFragments());
            assertTrue(prosit.daltonMatching);
        } finally {
            Constants.auxSpectraModel = saved;
        }
    }

    //an auxiliary prediction of three fragments, with these isotope numbers (none, as read from an MGF file without)
    private static PredictionEntry aux(int[] isotopes) {
        PredictionEntry pe = new PredictionEntry(new float[] {300f, 400f, 500f}, new float[] {1f, 0.5f, 0.25f},
                new int[0], new int[0], new String[] {"y", "imm", "b"});
        pe.isotopes = isotopes;
        return pe;
    }

    //the entry of a PSM with no main spectrum once the auxiliary prediction is merged into it
    private static PredictionEntry merged(PredictionEntry auxPrediction) {
        PredictionEntryHashMap library = new PredictionEntryHashMap();
        PredictionEntry main = new PredictionEntry();
        library.put("LGEPTIDEK|2", main);
        PredictionEntryHashMap auxLibrary = new PredictionEntryHashMap();
        auxLibrary.put("LGEPTIDEK|2", auxPrediction);
        auxLibrary.mergeIntoLibrary(library, "auxSpectra");
        PredictionEntry entry = library.get("LGEPTIDEK|2");
        assertSame(auxPrediction, entry.auxSpectra);
        return entry;
    }
}
