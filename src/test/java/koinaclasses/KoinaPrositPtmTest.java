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
import static org.junit.jupiter.api.Assertions.assertNotNull;
import static org.junit.jupiter.api.Assertions.assertNull;
import static org.junit.jupiter.api.Assertions.assertThrows;
import static org.junit.jupiter.api.Assertions.assertTrue;

import java.util.ArrayList;
import java.util.List;
import java.util.Map;
import koinaclasses.KoinaModelSpec.EncodableMods;
import koinaclasses.KoinaModelSpec.ModMassMap;
import koinaclasses.KoinaModelSpec.ModRule;
import koinaclasses.KoinaModelSpec.PeptideFormat;
import org.junit.jupiter.api.Test;
import peptideptmformatting.PTMhandler;
import peptideptmformatting.PeptideFormatter;

/**
 * The mods Koina's Prosit-PTM models (Prosit_2025_intensity_22PTM, Prosit_2025_intensity_40PTM,
 * Prosit_2025_irt_22PTM) can be sent, and how their peptides are written and read back. The pairs below were sent to
 * the live models one at a time, as synthetic peptides, when the rule was written: Koina predicted every pair the
 * rule takes, and answered HTTP 400 for every pair it refuses.
 */
class KoinaPrositPtmTest {
    private static final List<String> PTM_MODELS = List.of(
            "Prosit_2025_intensity_22PTM", "Prosit_2025_intensity_40PTM", "Prosit_2025_irt_22PTM");

    @Test
    void handTablePairsAreTaken() {
        String[][] pairs = {{"S", "21"}, {"T", "21"}, {"Y", "21"}, {"P", "21"}, {"n", "1"}, {"K", "1"},
                {"n", "737"}, {"K", "737"}, {"n", "2016"}, {"C", "4"}, {"K", "4"}, {"K", "1263"}, {"M", "35"},
                {"Q", "28"}, {"E", "27"}, {"E", "34"}, {"N", "7"}, {"R", "7"}, {"K", "121"}, {"K", "36"},
                {"S", "43"}};
        for (String[] pair : pairs) {
            char site = pair[0].charAt(0);
            assertTrue(KoinaPrositPtm.inHandTable(site, pair[1]), pair[0] + pair[1]);
            assertNull(KoinaPrositPtm.refusal(site, pair[1]), pair[0] + pair[1]);
        }
    }

    //a pair outside the hand-written table: the site's attachment group plus the delta composition
    @Test
    void theZeroShotEncodingTakesPairsWithAnAttachmentGroupAndNoNegativeCount() {
        String[][] taken = {
                {"S", "737"}, //TMT on S: attachment group H O
                {"K", "747"}, //malonyl
                {"M", "425"}, //dioxidation
                {"C", "26"}, //pyro-carbamidomethyl
                {"Y", "354"}, //nitro, H(-1) N O(2): attachment group H O plus it leaves no H, which is not negative
                {"n", "122"}, //formyl on the N-terminus: attachment group H(2) N
                {"K", "259"}}; //13C(6) 15N(2): isotopes count as their element, so it adds no atom
        for (String[] pair : taken) {
            char site = pair[0].charAt(0);
            assertFalse(KoinaPrositPtm.inHandTable(site, pair[1]), pair[0] + pair[1]);
            assertNull(KoinaPrositPtm.refusal(site, pair[1]), pair[0] + pair[1]);
        }
    }

    @Test
    void pairsKoinaCannotEncodeAreRefusedWithItsReason() {
        //dehydration on S: H O plus H(-2) O(-1) leaves -1 H
        assertRefused('S', "23", "residue S plus UNIMOD:23 (H(-2) O(-1)) has a negative H count");
        //dHex is a sugar Koina lists, but its token pattern cannot read a symbol that starts in lower case
        assertRefused('T', "295", "cannot read composition token 'dHex' of 'dHex'");
        assertRefused('T', "142", "cannot read composition token 'dHex' of 'dHex HexNAc'");
        assertRefused('N', "43", "no attachment group for residue N");
        assertRefused('E', "1", "no attachment group for residue E");
        assertRefused('A', "34", "no attachment group for residue A");
        //pyro-Glu is in the hand-written table on Q, not on the N-terminus
        assertRefused('n', "28", "the N-terminus plus UNIMOD:28 (H(-3) N(-1)) has a negative H count");
        assertRefused('D', "30", "element 'Na' cannot be encoded (only H C N O P S)");
        assertRefused('K', "2114", "UNIMOD:2114 is not in Koina's unimod.obo");
        //in the hand-written table, but Koina's lookup of the code fails before it reaches the table
        assertRefused('K', "12118", "UNIMOD:12118 is not in Koina's unimod.obo");
        assertRefused('K', "Mymod", "UNIMOD:Mymod is not in Koina's unimod.obo");
    }

    private static void assertRefused(char site, String code, String reason) {
        assertEquals(reason, KoinaPrositPtm.refusal(site, code), site + code);
        assertFalse(KoinaPrositPtm.encodes((site == KoinaPrositPtm.NTERM ? "[" : String.valueOf(site)) + code),
                site + code);
    }

    @Test
    void everyHandTablePairIsEncodableIfKoinaHasItsCode() {
        List<String> notInKoinasObo = new ArrayList<>();
        for (String mod : KoinaPrositPtm.handTableMods()) {
            String code = mod.substring(1);
            if (KoinaUnimod.composition(code) == null) {
                notInKoinasObo.add(mod);
                assertFalse(KoinaPrositPtm.encodes(mod), mod);
            } else {
                assertTrue(KoinaPrositPtm.encodes(mod), mod);
            }
        }
        assertEquals(List.of("K5634", "K12118", "K19903", "K129317"), notInKoinasObo);
        assertEquals(64, KoinaPrositPtm.handTableMods().size());
        assertTrue(KoinaPrositPtm.handTableMods().containsAll(List.of("[1", "[737", "[2016", "C4", "C1263", "P21")));
        assertFalse(KoinaPrositPtm.preferredMods().contains("C1263"), "carbamidomethyl, not Gly, is written");
        assertFalse(KoinaPrositPtm.preferredMods().contains("K1263"), "carbamidomethyl, not Gly, is written");
        assertEquals(KoinaPrositPtm.handTableMods().size() - 2, KoinaPrositPtm.preferredMods().size());
    }

    @Test
    void compositionsAreReadAsTheZeroShotEncodingReadsThem() {
        //H C N O P S
        assertArrayEquals(new int[] {20, 12, 2, 2, 0, 0}, KoinaPrositPtm.atoms("H(20) C(8) 13C(4) N 15N O(2)"));
        assertArrayEquals(new int[] {33, 20, 1, 15, 0, 0}, KoinaPrositPtm.atoms("Hex(2) HexNAc"));
        assertArrayEquals(new int[] {17, 11, 1, 9, 0, 0}, KoinaPrositPtm.atoms("NeuGc"));
        assertArrayEquals(new int[] {1, 0, 0, 3, 1, 0}, KoinaPrositPtm.atoms("H O(3) P"));
        assertThrows(IllegalArgumentException.class, () -> KoinaPrositPtm.atoms("dHex HexNAc"));
        assertThrows(IllegalArgumentException.class, () -> KoinaPrositPtm.atoms("H(-1) Na"));
        assertThrows(IllegalArgumentException.class, () -> KoinaPrositPtm.atoms("Hep(2)"),
                "a sugar Koina does not list");
    }

    @Test
    void aCtermModIsNeverEncodable() {
        assertTrue(EncodableMods.PROSIT_PTM.encodes("[1"));
        assertTrue(EncodableMods.PROSIT_PTM.encodes("[1", false));
        assertFalse(EncodableMods.PROSIT_PTM.encodes("[1", true),
                "Koina would read a C-terminal acetyl as one on the N-terminus");
        assertFalse(EncodableMods.PROSIT_PTM.encodes("[2", true));
        assertTrue(EncodableMods.PROSIT_PTM.encodes("S737"));
        assertFalse(EncodableMods.PROSIT_PTM.encodes("S23"));
        //the other rules do not look at the terminus
        assertTrue(EncodableMods.ANY.encodes("[2", true));
        assertTrue(EncodableMods.KOINA_ELEMENTS.encodes("[2", true));
    }

    @Test
    void theThreeModelsShareTheirRules() {
        for (String model : PTM_MODELS) {
            KoinaModelSpec spec = KoinaModelRegistry.require(model);
            assertEquals(PeptideFormat.PROSIT_PTM, spec.peptideFormat(), model);
            assertEquals(ModMassMap.PROSIT_PTM, spec.modMassMap(), model);
            assertEquals(ModRule.PROSIT_PTM, spec.modRule(), model);
            assertEquals(EncodableMods.PROSIT_PTM, spec.encodableMods(), model);
            assertFalse(spec.peptideFormat().camOnBareCys(), model + ": a bare C is a free cysteine");
            assertFalse(spec.limits().skipFreeCys(), model);
            assertFalse(spec.singleThreadFormatting(), model);
            assertEquals(30, spec.limits().maxLength(), model);
            assertEquals(1, spec.limits().minCharge(), model);
            assertEquals(6, spec.limits().maxCharge(), model);
            assertEquals(1000, spec.batchSize(), model);
            assertEquals(model, spec.citationKey(), model);
        }
    }

    @Test
    void theMassTableAddsPhosphoOnPAndGlyOnC() {
        Map<String, Double> table = PTMhandler.getAAunimodToModMassPrositPtm();
        assertFalse(PTMhandler.AAunimodToModMassAll.containsKey("P21"));
        assertFalse(PTMhandler.AAunimodToModMassAll.containsKey("C1263"));
        assertEquals(PTMhandler.unimodToModMassAll.get(21), table.get("P21"));
        assertEquals(PTMhandler.unimodToModMassAll.get(1263), table.get("C1263"));
        assertEquals(PTMhandler.AAunimodToModMassAll.size() + 2, table.size());
        assertTrue(table.entrySet().containsAll(PTMhandler.AAunimodToModMassAll.entrySet()));
        assertEquals(table, PTMhandler.getAAunimodToModMassPrositPtm(), "made once for one ALL table");
    }

    //pin peptide, the peptide the model is sent, the mods left out of it, and the peptide read back from it
    private static final String[][] ROUND_TRIPS = {
            {"K.LGEPTIDEK.A", "LGEPTIDEK", "", "LGEPTIDEK"},
            {"K.n[42.0106]ADVTGLPHK.E", "[UNIMOD:1]-ADVTGLPHK", "", "[42.0106]ADVTGLPHK"},
            //propionyl is K[UNIMOD:58] of the hand-written table, not its isobaric K[UNIMOD:206] or K[UNIMOD:5634]
            {"K.ADK[56.0262]LGEPTIDER.S", "ADK[UNIMOD:58]LGEPTIDER", "", "ADK[56.0262]LGEPTIDER"},
            {"K.n[56.0262]ADVTGLPHK.E", "[UNIMOD:58]-ADVTGLPHK", "", "[56.0262]ADVTGLPHK"},
            {"K.n[42.0106]M[15.9949]DPSLVTAEGK.A", "[UNIMOD:1]-M[UNIMOD:35]DPSLVTAEGK", "",
                    "[42.0106]M[15.9949]DPSLVTAEGK"},
            {"K.n[229.1629]LGEPTIDEK[229.1629].A", "[UNIMOD:737]-LGEPTIDEK[UNIMOD:737]", "",
                    "[229.1629]LGEPTIDEK[229.1629]"},
            {"K.n[229.1629]GS[229.1629]PEDLHVTAK[229.1629].A", "[UNIMOD:737]-GS[UNIMOD:737]PEDLHVTAK[UNIMOD:737]",
                    "", "[229.1629]GS[229.1629]PEDLHVTAK[229.1629]"},
            {"K.n[27.9949]LGEPTIDEK.A", "[UNIMOD:122]-LGEPTIDEK", "", "[27.9949]LGEPTIDEK"},
            {"K.SHCVAELLGPEK.T", "SHCVAELLGPEK", "", "SHCVAELLGPEK"},
            {"K.TEC[57.0215]LDHFAVNPR.G", "TEC[UNIMOD:4]LDHFAVNPR", "", "TEC[57.0215]LDHFAVNPR"},
            {"K.ADK[57.0215]LGEPTIDER.S", "ADK[UNIMOD:4]LGEPTIDER", "", "ADK[57.0215]LGEPTIDER"},
            {"K.C[39.9950]ADLSGFTPEVIAR.Q", "C[UNIMOD:26]ADLSGFTPEVIAR", "", "C[39.9949]ADLSGFTPEVIAR"},
            {"R.GS[79.9663]PEDLHVTAR.L", "GS[UNIMOD:21]PEDLHVTAR", "", "GS[79.9663]PEDLHVTAR"},
            {"R.GSP[79.9663]EDLHVTAR.L", "GSP[UNIMOD:21]EDLHVTAR", "", "GSP[79.9663]EDLHVTAR"},
            {"K.LGVEK[114.0429]AYNDGR.F", "LGVEK[UNIMOD:121]AYNDGR", "", "LGVEK[114.0429]AYNDGR"},
            {"K.ADK[86.0004]LGEPTIDER.S", "ADK[UNIMOD:747]LGEPTIDER", "", "ADK[86.0004]LGEPTIDER"},
            {"K.PEPTM[31.9898]IDEK.A", "PEPTM[UNIMOD:425]IDEK", "", "PEPTM[31.9898]IDEK"},
            {"R.Q[-17.0265]YLDTEGAIPAK.D", "Q[UNIMOD:28]YLDTEGAIPAK", "", "Q[-17.0265]YLDTEGAIPAK"},
            {"K.VFLEN[0.9840]GTLR.A", "VFLEN[UNIMOD:7]GTLR", "", "VFLEN[0.9840]GTLR"},
            {"R.GS[-18.0106]PEDLHVTAR.L", "GSPEDLHVTAR", "S[-18.0106]", "GSPEDLHVTAR"},
            {"R.GT[146.0579]PEDLHVTAR.L", "GTPEDLHVTAR", "T[146.0579]", "GTPEDLHVTAR"},
            {"K.VFLEN[203.0794]GTLR.A", "VFLENGTLR", "N[203.0794]", "VFLENGTLR"},
            {"K.n[-17.0265]QYLDTEGAIPAK.D", "QYLDTEGAIPAK", "n[-17.0265]", "QYLDTEGAIPAK"},
            {"K.LGEPTIDEKc[-0.9840].A", "LGEPTIDEK", "c[-0.9840]", "LGEPTIDEK"},
            //a methyl ester: Koina would take methyl on the N-terminus, and read this one there
            {"K.ADK[42.0106]LGEPTIDEKc[14.0157].A", "ADK[UNIMOD:1]LGEPTIDEK", "c[14.0157]", "ADK[42.0106]LGEPTIDEK"}};

    //base -> model -> base: what is sent reads back as the PSM's peptide without the mods left out
    @Test
    void peptidesGoToTheModelAndBack() {
        for (String model : PTM_MODELS) {
            KoinaModelSpec spec = KoinaModelRegistry.require(model);
            for (String[] trip : ROUND_TRIPS) {
                String what = model + " " + trip[0];
                PeptideFormatter pf = new PeptideFormatter(trip[0], 2, "pin");
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

    //ms2pip_timsTOF2023 is the ms2pip_timsTOF2024 spec under its own name
    @Test
    void msPipTimsTof2023IsTheTimsTof2024Spec() {
        KoinaModelSpec tof2023 = KoinaModelRegistry.require("ms2pip_timsTOF2023");
        KoinaModelSpec tof2024 = KoinaModelRegistry.require("ms2pip_timsTOF2024");
        KoinaModelSpec renamed = new KoinaModelSpec(tof2023.name(), tof2024.kind(), tof2024.inputs(),
                tof2024.outputs(), tof2024.imConversion(), tof2024.peptideFormat(), tof2024.annotationParser(),
                tof2024.trustKoinaMz(), tof2024.modMassMap(), tof2024.modRule(), tof2024.encodableMods(),
                tof2024.limits(), tof2024.batchSize(), tof2024.minCe(), tof2024.maxCe(), tof2024.calibrateCe(),
                tof2024.instrumentMap(), tof2024.fragmentations(), tof2024.isobaricLabel(),
                tof2024.takesIsobaricLabels(), tof2024.sortBatchByLength(), tof2024.singleThreadFormatting(),
                tof2023.name());
        assertEquals(renamed, tof2023);
        assertEquals("ms2pip_timsTOF2023", tof2023.citationKey());

        String[][] trips = {
                {"R.AVLDQGSTWPEM[15.9949]AGK.S", "AVLDQGSTWPEM[UNIMOD:35]AGK", "", "AVLDQGSTWPEM[15.9949]AGK"},
                {"K.TEC[57.0215]LDHFAVNPR.G", "TEC[UNIMOD:4]LDHFAVNPR", "", "TEC[57.0215]LDHFAVNPR"},
                {"K.SHCVAELLGPEK.T", "SHCVAELLGPEK", "", "SHCVAELLGPEK"},
                {"R.GS[79.9663]PEDLHVTAR.L", "GSPEDLHVTAR", "S[79.9663]", "GSPEDLHVTAR"},
                {"K.n[42.0106]ADVTGLPHK.E", "ADVTGLPHK", "n[42.0106]", "ADVTGLPHK"}};
        for (String[] trip : trips) {
            PeptideFormatter pf = new PeptideFormatter(trip[0], 2, "pin");
            List<String> dropped = new ArrayList<>();
            String sent = pf.getModel(tof2023, dropped);
            assertEquals(trip[1], sent, trip[0]);
            assertEquals(trip[2], String.join(",", dropped), trip[0]);
            assertEquals(trip[3], new PeptideFormatter(sent, 2, tof2023).getBase(), trip[0]);
            assertEquals(sent, new PeptideFormatter(trip[0], 2, "pin").getModel(tof2024), trip[0]);
        }
        assertNotNull(KoinaModelRegistry.find("MS2PIP_TIMSTOF2023"));
    }
}
