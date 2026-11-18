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

import java.util.ArrayList;
import java.util.List;
import java.util.Set;
import java.util.TreeSet;
import koinaclasses.KoinaModelSpec.EncodableMods;
import koinaclasses.KoinaModelSpec.Limits;
import koinaclasses.KoinaModelSpec.ModMassMap;
import koinaclasses.KoinaModelSpec.ModRule;
import koinaclasses.KoinaModelSpec.PeptideFormat;
import org.junit.jupiter.api.Test;
import peptideptmformatting.PTMhandler;
import peptideptmformatting.PeptideFormatter;
import peptideptmformatting.PeptideSkipper;

/**
 * Chronologer_RT: Koina reads the residue and N-term pairs of its table
 * (models/Chronologer/Chronologer_Preprocess/1/sequence_conversion.py, commit d176169), a bare C as a free cysteine,
 * and the first 50 residues of a peptide. It reads an unknown pair, or the letters of a C-term mod, as something
 * else without saying so, so MSBooster sends neither.
 */
class KoinaChronologerTest {
    private static final KoinaModelSpec CHRONOLOGER = KoinaModelRegistry.require("Chronologer_RT");

    @Test
    void theSpecIsAnRtModelWithItsOwnTable() {
        assertEquals(KoinaModelSpec.Kind.RT, CHRONOLOGER.kind());
        assertEquals(List.of("rt"), CHRONOLOGER.outputs());
        assertTrue(CHRONOLOGER.inputs().isEmpty());
        assertEquals(PeptideFormat.CHRONOLOGER, CHRONOLOGER.peptideFormat());
        assertFalse(CHRONOLOGER.peptideFormat().camOnBareCys(), "a bare C is a free cysteine");
        assertEquals("-", CHRONOLOGER.peptideFormat().ntermSeparator());
        assertEquals(ModMassMap.CHRONOLOGER, CHRONOLOGER.modMassMap());
        assertEquals(ModRule.CHRONOLOGER, CHRONOLOGER.modRule());
        assertEquals(EncodableMods.NO_CTERM, CHRONOLOGER.encodableMods());
        assertEquals(Limits.NONSTANDARD_RESIDUES, CHRONOLOGER.limits().forbiddenResidues());
        assertEquals(50, CHRONOLOGER.limits().maxLength());
        assertFalse(CHRONOLOGER.limits().skipFreeCys());
        assertEquals(1000, CHRONOLOGER.batchSize());
        assertEquals("Chronologer_RT", CHRONOLOGER.citationKey());
    }

    //the pairs of Koina's amino_acid_unimod_map and n_terminal_mod_map, but succinyl (K-125) and the N-term pyro
    //forms (27, 28 and 385 before residue 1)
    @Test
    void theTableIsKoinasTableWithoutSuccinylAndNtermPyroForms() {
        assertEquals(new TreeSet<>(Set.of("C4", "M35", "S21", "T21", "Y21", "K1", "K121", "K34", "K36", "K37",
                        "R34", "R36", "K739", "K737", "C385", "E27", "Q28", "[1", "[737", "[739")),
                new TreeSet<>(PTMhandler.AAunimodToModMassChronologer.keySet()));
        assertEquals(PTMhandler.pyroCarbamidomethylMass, PTMhandler.AAunimodToModMassChronologer.get("C385"));
        assertEquals(PTMhandler.pyroCarbamidomethylMass, PTMhandler.unimodToModMassChronologer.get(385));
        for (String aaCode : PTMhandler.AAunimodToModMassChronologer.keySet()) {
            int code = Integer.parseInt(aaCode.substring(1));
            assertEquals(PTMhandler.AAunimodToModMassChronologer.get(aaCode),
                    PTMhandler.unimodToModMassChronologer.get(code), aaCode + " reads back with the mass it is sent for");
        }
    }

    //pin peptide, the peptide the model is sent (null if not sent), the mods left out of it, and the peptide read back
    private static final String[][] ROUND_TRIPS = {
            {"K.LGEPTIDEK.A", "LGEPTIDEK", "", "LGEPTIDEK"},
            {"R.AVLDQGSTWPEM[15.9949]AGK.S", "AVLDQGSTWPEM[UNIMOD:35]AGK", "", "AVLDQGSTWPEM[15.9949]AGK"},
            {"K.TEC[57.0215]LDHFAVNPR.G", "TEC[UNIMOD:4]LDHFAVNPR", "", "TEC[57.0215]LDHFAVNPR"},
            {"K.SHCVAELLGPEK.T", "SHCVAELLGPEK", "", "SHCVAELLGPEK"},
            {"K.C[39.9950]ADLSGFTPEVIAR.Q", "C[UNIMOD:385]ADLSGFTPEVIAR", "", "C[39.9950]ADLSGFTPEVIAR"},
            {"K.n[42.0106]ADVTGLPHK.E", "[UNIMOD:1]-ADVTGLPHK", "", "[42.0106]ADVTGLPHK"},
            {"K.n[42.0106]M[15.9949]DPSLVTAEGK.A", "[UNIMOD:1]-M[UNIMOD:35]DPSLVTAEGK", "",
                    "[42.0106]M[15.9949]DPSLVTAEGK"},
            {"R.GS[79.9663]PEDLHVTAR.L", "GS[UNIMOD:21]PEDLHVTAR", "", "GS[79.9663]PEDLHVTAR"},
            {"R.NDT[79.9663]PLM[15.9949]EAGR.V", "NDT[UNIMOD:21]PLM[UNIMOD:35]EAGR", "",
                    "NDT[79.9663]PLM[15.9949]EAGR"},
            {"R.GEDY[79.9663]LPTAK.S", "GEDY[UNIMOD:21]LPTAK", "", "GEDY[79.9663]LPTAK"},
            {"K.ADK[42.0106]LGEPTIDER.S", "ADK[UNIMOD:1]LGEPTIDER", "", "ADK[42.0106]LGEPTIDER"},
            {"K.LGVEK[114.0429]AYNDGR.F", "LGVEK[UNIMOD:121]AYNDGR", "", "LGVEK[114.0429]AYNDGR"},
            {"K.ADK[14.0157]LGEPTIDER.S", "ADK[UNIMOD:34]LGEPTIDER", "", "ADK[14.0157]LGEPTIDER"},
            {"K.ADK[28.0313]LGEPTIDER.S", "ADK[UNIMOD:36]LGEPTIDER", "", "ADK[28.0313]LGEPTIDER"},
            {"K.ADK[42.0470]LGEPTIDER.S", "ADK[UNIMOD:37]LGEPTIDER", "", "ADK[42.0470]LGEPTIDER"},
            {"R.GR[14.0157]PEDLHVTAK.L", "GR[UNIMOD:34]PEDLHVTAK", "", "GR[14.0157]PEDLHVTAK"},
            {"R.GR[28.0313]PEDLHVTAK.L", "GR[UNIMOD:36]PEDLHVTAK", "", "GR[28.0313]PEDLHVTAK"},
            {"K.n[229.1629]LGEPTIDEK[229.1629].A", "[UNIMOD:737]-LGEPTIDEK[UNIMOD:737]", "",
                    "[229.1629]LGEPTIDEK[229.1629]"},
            {"K.n[224.1525]LGEPTIDEK[224.1525].A", "[UNIMOD:739]-LGEPTIDEK[UNIMOD:739]", "",
                    "[224.1525]LGEPTIDEK[224.1525]"},
            {"R.Q[-17.0265]YLDTEGAIPAK.D", "Q[UNIMOD:28]YLDTEGAIPAK", "", "Q[-17.0265]YLDTEGAIPAK"},
            {"K.E[-18.0106]WLSPQDGAVR.T", "E[UNIMOD:27]WLSPQDGAVR", "", "E[-18.0106]WLSPQDGAVR"},
            //succinyl: Koina's table has it as K-125, not UniMod's 64, and MSBooster leaves it out
            {"K.ADK[100.0160]LGEPTIDER.S", "ADKLGEPTIDER", "K[100.0160]", "ADKLGEPTIDER"},
            {"K.VFLEN[0.9840]GTLR.A", "VFLENGTLR", "N[0.9840]", "VFLENGTLR"},
            {"K.ADS[42.0106]LGEPTIDER.S", "ADSLGEPTIDER", "S[42.0106]", "ADSLGEPTIDER"},
            {"K.n[-17.0265]QYLDTEGAIPAK.D", "QYLDTEGAIPAK", "n[-17.0265]", "QYLDTEGAIPAK"},
            {"K.LGEPTIDEKc[-0.9840].A", "LGEPTIDEK", "c[-0.9840]", "LGEPTIDEK"},
            //a C-term mod of the N-term acetyl mass: Koina would read the letters of [UNIMOD:1] after K as residues
            {"K.LGEPTIDEKc[42.0106].A", "LGEPTIDEK", "c[42.0106]", "LGEPTIDEK"},
            {"K.LGEPTIDEK[42.0106]c[-0.9840].A", "LGEPTIDEK[UNIMOD:1]", "c[-0.9840]", "LGEPTIDEK[42.0106]"}};

    //base -> model -> base: what is sent reads back as the PSM's peptide without the mods left out
    @Test
    void peptidesGoToTheModelAndBack() {
        for (String[] trip : ROUND_TRIPS) {
            String what = trip[0];
            PeptideFormatter pf = new PeptideFormatter(trip[0], 2, "pin");
            assertNull(PeptideSkipper.skipReason(pf, CHRONOLOGER), what);
            List<String> dropped = new ArrayList<>();
            String sent = pf.getModel(CHRONOLOGER, dropped);
            assertEquals(trip[1], sent, what);
            assertEquals(trip[2], String.join(",", dropped), what);
            PeptideFormatter back = new PeptideFormatter(sent, 2, CHRONOLOGER);
            assertEquals(trip[3], back.getBase(), what);
            //read back, it is sent as it was
            assertEquals(sent, back.getModel(CHRONOLOGER), what);
        }
    }

    //Koina reads a longer peptide's first 50 residues and drops the rest, so MSBooster does not send it
    @Test
    void peptidesOfMoreThan50ResiduesAreNotSent() {
        String fifty = "ELVSAPLGTNEIYDHFAQWTGLPVEDAKLGRSTEPVNGHAKYLDTEGAIP";
        assertEquals(50, fifty.length());
        assertNull(skip(fifty, 3));
        assertEquals("longer than 50 residues", skip(fifty + "K", 3));
        assertEquals("residue U", skip("GDUVLPAEK", 2));
        assertNull(skip("TEALDHFAVNPRGSLEWQKPEDVTAR", 7), "an RT model takes any charge");
        assertNull(skip("SHCVAELLGPEK", 2), "a free cysteine is sent as a bare C");
    }

    //TMT has the run's TMT mass, as in the LIMITED table: a TMTpro run is sent its label as TMT
    @Test
    void tmtFollowsTheIsobaricLabelOfTheRun() {
        try {
            PTMhandler.setIsobaricLabel(PTMhandler.tmtproMass);
            assertEquals(PTMhandler.tmtproMass, PTMhandler.AAunimodToModMassChronologer.get("K737"));
            assertEquals(PTMhandler.tmtproMass, PTMhandler.AAunimodToModMassChronologer.get("[737"));
            assertEquals(PTMhandler.tmtproMass, PTMhandler.unimodToModMassChronologer.get(737));
            String sent = new PeptideFormatter("[304.2071]LGEPTIDEK[304.2071]", 2, "base").getModel(CHRONOLOGER);
            assertEquals("[UNIMOD:737]-LGEPTIDEK[UNIMOD:737]", sent);
            assertEquals("[304.2072]LGEPTIDEK[304.2072]", new PeptideFormatter(sent, 2, CHRONOLOGER).getBase());
        } finally {
            PTMhandler.resetIsobaricLabel();
        }
        assertEquals(PTMhandler.tmt6_10_11Mass, PTMhandler.AAunimodToModMassChronologer.get("K737"));
        assertEquals(PTMhandler.tmtZeroMass, PTMhandler.AAunimodToModMassChronologer.get("K739"));
    }

    private static String skip(String base, int charge) {
        return PeptideSkipper.skipReason(new PeptideFormatter(base, charge, "base"), CHRONOLOGER);
    }
}
