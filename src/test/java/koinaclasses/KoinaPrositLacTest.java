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
import static org.junit.jupiter.api.Assertions.assertFalse;
import static org.junit.jupiter.api.Assertions.assertNull;
import static org.junit.jupiter.api.Assertions.assertTrue;

import allconstants.Constants;
import java.util.ArrayList;
import java.util.List;
import java.util.Map;
import java.util.Set;
import koinaclasses.KoinaModelSpec.EncodableMods;
import koinaclasses.KoinaModelSpec.Input;
import koinaclasses.KoinaModelSpec.InstrumentMap;
import koinaclasses.KoinaModelSpec.Limits;
import koinaclasses.KoinaModelSpec.ModMassMap;
import koinaclasses.KoinaModelSpec.ModRule;
import koinaclasses.KoinaModelSpec.PeptideFormat;
import org.junit.jupiter.api.Test;
import peptideptmformatting.PTMhandler;
import peptideptmformatting.PeptideFormatter;
import peptideptmformatting.PeptideSkipper;
import utils.InstrumentUtils;

/**
 * The Prosit lactylation models (Prosit_2025_intensity_lac, Prosit_2025_irt_lac): Koina reads K[UNIMOD:2114],
 * [UNIMOD:2114]-, M[UNIMOD:35] and C[UNIMOD:4], and every C as carbamidomethylated
 * (models/Prosit/Prosit_Preprocess_peptide_lac/1/sequence_conversion.py, commit d176169).
 */
class KoinaPrositLacTest {
    private static final List<String> LAC_MODELS = List.of("Prosit_2025_intensity_lac", "Prosit_2025_irt_lac");

    @Test
    void theTwoModelsShareTheirRules() {
        for (String model : LAC_MODELS) {
            KoinaModelSpec spec = KoinaModelRegistry.require(model);
            assertEquals(PeptideFormat.PROSIT, spec.peptideFormat(), model);
            assertTrue(spec.peptideFormat().camOnBareCys(), model + ": a C is sent as carbamidomethylated");
            assertEquals(ModMassMap.LAC, spec.modMassMap(), model);
            assertEquals(ModRule.LAC, spec.modRule(), model);
            assertEquals(EncodableMods.NO_CTERM, spec.encodableMods(), model);
            assertTrue(spec.limits().skipFreeCys(), model + ": a free cysteine would be read as carbamidomethylated");
            assertEquals(Limits.NONSTANDARD_RESIDUES, spec.limits().forbiddenResidues(), model);
            assertEquals(30, spec.limits().maxLength(), model);
            assertEquals(1, spec.limits().minCharge(), model);
            assertEquals(6, spec.limits().maxCharge(), model);
            assertEquals(1000, spec.batchSize(), model);
            assertEquals(model, spec.citationKey(), model);
            assertFalse(spec.singleThreadFormatting(), model);
        }
        KoinaModelSpec intensity = KoinaModelRegistry.require("Prosit_2025_intensity_lac");
        assertEquals(List.of(Input.PRECURSOR_CHARGES, Input.COLLISION_ENERGIES, Input.INSTRUMENT_TYPES,
                Input.FRAGMENTATION_TYPES), new ArrayList<>(intensity.inputs()));
        assertTrue(intensity.calibrateCe());
        assertEquals(20, intensity.minCe());
        assertEquals(40, intensity.maxCe());
        assertEquals(Set.of("HCD", "CID"), intensity.allowedFragmentationTypes());
        assertEquals("HCD", intensity.fragmentationSentAs("HCD"));
        assertEquals("CID", intensity.fragmentationSentAs("CID"));
        assertEquals(InstrumentMap.LAC, intensity.instrumentMap());
        assertEquals(List.of("intensities", "mz", "annotation"), intensity.outputs());
        assertEquals(List.of("irt"), KoinaModelRegistry.require("Prosit_2025_irt_lac").outputs());
    }

    //D21: Koina reads LUMOS, ASTRAL and ECLIPSE; MSBooster's Astral is THERMOTOF, and an Eclipse is detected as LUMOS
    @Test
    void anAstralIsSentAsAstralAndAnythingElseAsLumos() {
        KoinaModelSpec spec = KoinaModelRegistry.require("Prosit_2025_intensity_lac");
        String saved = Constants.instrument;
        try {
            for (String instrument : new String[] {"THERMOTOF", "thermotof", "ThermoTOF"}) {
                for (String charge : new String[] {"1", "2", "6"}) {
                    assertEquals("ASTRAL", spec.instrumentMap().map(instrument, charge), instrument);
                    Constants.instrument = instrument;
                    assertEquals("ASTRAL", InstrumentUtils.mapInstrumentToModelSpecific(spec, charge), instrument);
                }
            }
            for (String instrument : new String[] {"", "LUMOS", "lumos", "QE", "QEHFX", "TIMSTOF", "SCIEXTOF",
                    "ELITE", "VELOS", "NONE", "ASTRAL", "Orbitrap"}) {
                for (String charge : new String[] {"1", "2", "6"}) {
                    assertEquals("LUMOS", spec.instrumentMap().map(instrument, charge), instrument);
                    Constants.instrument = instrument;
                    assertEquals("LUMOS", InstrumentUtils.mapInstrumentToModelSpecific(spec, charge), instrument);
                }
            }
        } finally {
            Constants.instrument = saved;
        }
    }

    @Test
    void theMassTableIsLactylOxidationAndCarbamidomethyl() {
        assertEquals(Map.of("K2114", 72.021129, "[2114", 72.021129, "M35", PTMhandler.oxidationMass,
                "C4", PTMhandler.carbamidomethylationMass), PTMhandler.AAunimodToModMassLac);
        assertEquals(Map.of(2114, 72.021129, 35, PTMhandler.oxidationMass, 4, PTMhandler.carbamidomethylationMass),
                PTMhandler.unimodToModMassLac);
        //the ALL table gives a K of the lactyl mass another code, which the models do not read
        assertFalse(PTMhandler.AAunimodToModMassAll.containsKey("K2114"));
        assertTrue(PTMhandler.AAunimodToModMassAll.containsKey("K378"));
    }

    //pin peptide, the peptide the model is sent (null if not sent), the mods left out of it, and the peptide read back
    private static final String[][] ROUND_TRIPS = {
            {"K.LGEPTIDEK.A", "LGEPTIDEK", "", "LGEPTIDEK"},
            {"K.ADK[72.0211]LGEPTIDER.S", "ADK[UNIMOD:2114]LGEPTIDER", "", "ADK[72.0211]LGEPTIDER"},
            {"K.n[72.0211]ADVTGLPHK.E", "[UNIMOD:2114]-ADVTGLPHK", "", "[72.0211]ADVTGLPHK"},
            {"K.n[72.0211]LGEPTK[72.0211]IDEK.A", "[UNIMOD:2114]-LGEPTK[UNIMOD:2114]IDEK", "",
                    "[72.0211]LGEPTK[72.0211]IDEK"},
            {"K.M[15.9949]DPSLVTAEGK[72.0211].A", "M[UNIMOD:35]DPSLVTAEGK[UNIMOD:2114]", "",
                    "M[15.9949]DPSLVTAEGK[72.0211]"},
            {"K.TEC[57.0215]LDHFAVNPR.G", "TEC[UNIMOD:4]LDHFAVNPR", "", "TEC[57.0215]LDHFAVNPR"},
            //the pyro-carbamidomethyl cysteine is sent as carbamidomethylated, its -17.0265 left out
            {"K.C[39.9950]ADLSGFTPEVIAR.Q", "C[UNIMOD:4]ADLSGFTPEVIAR", "C[39.9950]", "C[57.0215]ADLSGFTPEVIAR"},
            //propionyl (UNIMOD:58) is not lactyl, and the models have no code for it
            {"K.ADK[56.0262]LGEPTIDER.S", "ADKLGEPTIDER", "K[56.0262]", "ADKLGEPTIDER"},
            {"K.ADK[42.0106]LGEPTIDER.S", "ADKLGEPTIDER", "K[42.0106]", "ADKLGEPTIDER"},
            {"K.n[42.0106]ADVTGLPHK.E", "ADVTGLPHK", "n[42.0106]", "ADVTGLPHK"},
            {"R.GS[79.9663]PEDLHVTAR.L", "GSPEDLHVTAR", "S[79.9663]", "GSPEDLHVTAR"},
            {"R.Q[-17.0265]YLDTEGAIPAK.D", "QYLDTEGAIPAK", "Q[-17.0265]", "QYLDTEGAIPAK"},
            //a C-term mod of the lactyl mass is not written after the last residue, where it would read as K-lactyl
            {"K.LGEPTIDEKc[72.0211].A", "LGEPTIDEK", "c[72.0211]", "LGEPTIDEK"},
            {"K.LGEPTIDEK[72.0211]c[-0.9840].A", "LGEPTIDEK[UNIMOD:2114]", "c[-0.9840]", "LGEPTIDEK[72.0211]"},
            {"K.SHCVAELLGPEK.T", null, "", null}};

    //base -> model -> base: what is sent reads back as the PSM's peptide without the mods left out
    @Test
    void peptidesGoToTheModelAndBack() {
        for (String model : LAC_MODELS) {
            KoinaModelSpec spec = KoinaModelRegistry.require(model);
            for (String[] trip : ROUND_TRIPS) {
                String what = model + " " + trip[0];
                PeptideFormatter pf = new PeptideFormatter(trip[0], 2, "pin");
                if (trip[1] == null) {
                    assertEquals(Limits.FREE_CYS, PeptideSkipper.skipReason(pf, spec), what);
                    continue;
                }
                assertNull(PeptideSkipper.skipReason(pf, spec), what);
                List<String> dropped = new ArrayList<>();
                String sent = pf.getModel(spec, dropped);
                assertEquals(trip[1], sent, what);
                assertEquals(trip[2], String.join(",", dropped), what);
                PeptideFormatter back = new PeptideFormatter(sent, 2, spec);
                assertEquals(trip[3], back.getBase(), what);
                //read back, it is sent as it was
                assertEquals(sent, back.getModel(spec), what);
            }
        }
    }

    //propionyl (K[UNIMOD:58]) is a mod the lactylation models leave out, not lactyl; the lactyl mass is lactyl for
    //them only, and the models that read any UniMod code keep their own codes for both
    @Test
    void lactylIsWrittenForTheLactylationModelsOnly() {
        KoinaModelSpec lac = KoinaModelRegistry.require("Prosit_2025_intensity_lac");
        KoinaModelSpec alphapept = KoinaModelRegistry.require("AlphaPept_ms2_generic");
        List<String> dropped = new ArrayList<>();
        assertEquals("ADKLGEPTIDER",
                new PeptideFormatter("ADK[56.0262]LGEPTIDER", 2, "base").getModel(lac, dropped));
        assertEquals(List.of("K[56.0262]"), dropped);
        assertEquals("ADK[UNIMOD:2114]LGEPTIDER",
                new PeptideFormatter("ADK[72.0211]LGEPTIDER", 2, "base").getModel(lac));
        assertEquals("ADKLGEPTIDER", new PeptideFormatter("ADK[72.0211]LGEPTIDER", 2, "base")
                .getModel(KoinaModelRegistry.require("Prosit_2020_intensity_HCD")));

        Set<String> found = Set.copyOf(PeptideFormatter.foundUnimods);
        PeptideFormatter.foundUnimods.clear();
        try {
            for (String base : new String[] {"ADK[56.0262]LGEPTIDER", "ADK[72.0211]LGEPTIDER"}) {
                String sent = new PeptideFormatter(base, 2, "base").getModel(alphapept);
                assertTrue(sent.startsWith("ADK[UNIMOD:") && !sent.contains("2114"), sent);
                assertEquals(base, new PeptideFormatter(sent, 2, alphapept).getBase());
            }
        } finally {
            PeptideFormatter.foundUnimods.clear();
            PeptideFormatter.foundUnimods.addAll(found);
        }
    }

    @Test
    void theModelsTakePeptidesUpTo30ResiduesAndCharges1To6WithoutAFreeCysteine() {
        for (String model : LAC_MODELS) {
            KoinaModelSpec spec = KoinaModelRegistry.require(model);
            assertNull(skip(spec, "ELVSAPLGTNEIYDHFAQWTGLPVEDAKLG", 2), model);
            assertEquals("longer than 30 residues", skip(spec, "ELVSAPLGTNEIYDHFAQWTGLPVEDAKLGR", 2), model);
            assertNull(skip(spec, "LGEPTIDEK", 1), model);
            assertNull(skip(spec, "LGEPTIDEK", 6), model);
            assertEquals("charge above 6", skip(spec, "LGEPTIDEK", 7), model);
            assertEquals("residue U", skip(spec, "GDUVLPAEK", 2), model);
            assertEquals(Limits.FREE_CYS, skip(spec, "SHCVAELLGPEK", 2), model);
            assertEquals(Limits.FREE_CYS, skip(spec, "C[57.0215]PEPTIDEC", 2), model);
            //a C-term mod is written after the last residue, and does not make a last C a modified one
            assertEquals(Limits.FREE_CYS, skip(spec, "LGEPTIDEC[-0.9840]cterm", 2), model);
            assertNull(skip(spec, "LGEPTIDEC[57.0215][-0.9840]cterm", 2), model);
            assertNull(skip(spec, "C[39.9950]PEPTIDEC[57.0215]", 2), model);
        }
    }

    private static String skip(KoinaModelSpec spec, String base, int charge) {
        return PeptideSkipper.skipReason(new PeptideFormatter(base, charge, "base"), spec);
    }
}
