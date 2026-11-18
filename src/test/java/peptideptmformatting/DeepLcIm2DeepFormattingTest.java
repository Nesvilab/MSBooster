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

package peptideptmformatting;

import static org.junit.jupiter.api.Assertions.assertEquals;

import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;
import java.util.Set;
import koinaclasses.KoinaModelRegistry;
import koinaclasses.KoinaModelSpec;
import koinaclasses.KoinaPreCheck;
import org.junit.jupiter.api.AfterEach;
import org.junit.jupiter.api.BeforeEach;
import org.junit.jupiter.api.Test;

/**
 * Koina's DeepLC and IM2Deep preprocessing read the atoms of one mod per residue, from the composition of its UniMod
 * code, and nothing written before the first residue. So the peptides sent to Deeplc_hela_hf and IM2Deep carry an
 * N-term mod on the first residue and a C-term mod on the last; a terminal residue with a mod of its own keeps only
 * that one (D24); and a mod whose composition the model cannot read is left out: one with an element other than C,
 * H, N, O, S and P for both, and also one with isotope atoms for IM2Deep. Read back, the peptides name the same
 * precursors.
 */
class DeepLcIm2DeepFormattingTest {
    private static final KoinaModelSpec DEEPLC = KoinaModelRegistry.require("Deeplc_hela_hf");
    private static final KoinaModelSpec IM2DEEP = KoinaModelRegistry.require("IM2Deep");
    private static final KoinaModelSpec ALPHAPEPT = KoinaModelRegistry.require("AlphaPept_rt_generic");

    private final Set<String> savedFoundUnimods = new HashSet<>();

    @BeforeEach
    void clearFoundUnimods() {
        savedFoundUnimods.addAll(PeptideFormatter.foundUnimods);
        PeptideFormatter.foundUnimods.clear();
    }

    @AfterEach
    void restoreFoundUnimods() {
        PeptideFormatter.foundUnimods.clear();
        PeptideFormatter.foundUnimods.addAll(savedFoundUnimods);
    }

    @Test
    void anNtermAcetylGoesOnResidue1AndReadsBackAsTheNtermMod() {
        for (KoinaModelSpec spec : List.of(DEEPLC, IM2DEEP)) {
            Sent sent = send("K.n[42.0106]ADVTGM[15.9949]LPHK.E", spec);
            assertEquals("A[UNIMOD:1]DVTGM[UNIMOD:35]LPHK", sent.peptide, spec.name());
            assertEquals(List.of(), sent.dropped, spec.name());
            assertEquals("[42.0106]ADVTGM[15.9949]LPHK", readBack(sent.peptide, spec), spec.name());
        }
    }

    @Test
    void anNtermTmtGoesOnResidue1() {
        Sent sent = send("K.n[229.1629]LGEPTIDEK[229.1629].A", DEEPLC);
        assertEquals("L[UNIMOD:737]GEPTIDEK[UNIMOD:737]", sent.peptide);
        assertEquals("[229.1629]LGEPTIDEK[229.1629]", readBack(sent.peptide, DEEPLC));
    }

    //a code the mass table has as a mod of residue 1 reads back as that, as Koina reads both alike
    @Test
    void aResidue1ModReadsBackAsItself() {
        Sent sent = send("K.K[42.0106]LGEPTIDEK.A", DEEPLC);
        assertEquals("K[UNIMOD:1]LGEPTIDEK", sent.peptide);
        assertEquals("K[42.0106]LGEPTIDEK", readBack(sent.peptide, DEEPLC));
    }

    //D24: the residue's own mod is kept and the terminal mod is left out, and counted
    @Test
    void anNtermModOnAModifiedResidue1IsLeftOutAndCounted() {
        for (KoinaModelSpec spec : List.of(DEEPLC, IM2DEEP)) {
            Sent sent = send("K.n[42.0106]M[15.9949]DPSLVTAEGK.A", spec);
            assertEquals("M[UNIMOD:35]DPSLVTAEGK", sent.peptide, spec.name());
            assertEquals(List.of("n[42.0106]"), sent.dropped, spec.name());
            assertEquals("M[15.9949]DPSLVTAEGK", readBack(sent.peptide, spec), spec.name());

            KoinaPreCheck preCheck = new KoinaPreCheck(spec);
            assertEquals("M[UNIMOD:35]DPSLVTAEGK",
                    preCheck.modelPeptide(new PeptideFormatter("K.n[42.0106]M[15.9949]DPSLVTAEGK.A", 2, "pin")));
            assertEquals("{n[42.0106]=1}", preCheck.droppedModPsms().toString(), spec.name());
        }
    }

    @Test
    void aCtermModGoesOnTheLastResidueUnlessItHasAModOfItsOwn() {
        Sent amidated = send("K.LGEPTIDEKc[-0.9840].A", DEEPLC);
        assertEquals("LGEPTIDEK[UNIMOD:2]", amidated.peptide);
        assertEquals(List.of(), amidated.dropped);

        Sent both = send("K.LGEPTIDEK[229.1629]c[-0.9840].A", DEEPLC);
        assertEquals("LGEPTIDEK[UNIMOD:737]", both.peptide);
        assertEquals(List.of("c[-0.9840]"), both.dropped);
    }

    //SILAC Lys (UNIMOD:259) is 13C(6) 15N(2), which IM2Deep's preprocessing cannot read
    @Test
    void im2deepLeavesOutAnIsotopeLabel() {
        Sent im2deep = send("K.LGEPTIDEK[8.0142].A", IM2DEEP);
        assertEquals("LGEPTIDEK", im2deep.peptide);
        assertEquals(List.of("K[8.0142]"), im2deep.dropped);
        assertEquals("LGEPTIDEK[UNIMOD:259]", send("K.LGEPTIDEK[8.0142].A", DEEPLC).peptide,
                "DeepLC reads an isotope as its element");
    }

    //a residue 1 whose own mod is left out takes the N-term mod
    @Test
    void anNtermModTakesResidue1WhenItsOwnModIsLeftOut() {
        Sent sent = send("R.n[42.0106]K[8.0142]LGEPTIDEK.A", IM2DEEP);
        assertEquals("K[UNIMOD:1]LGEPTIDEK", sent.peptide);
        assertEquals(List.of("K[8.0142]"), sent.dropped);
    }

    //Na (Cation:Na, UNIMOD:30) and a sugar (HexNAc, UNIMOD:43) are not among the elements either model reads
    @Test
    void bothLeaveOutModsWithOtherElements() {
        for (KoinaModelSpec spec : List.of(DEEPLC, IM2DEEP)) {
            Sent sent = send("K.LGE[21.9819]PTIDES[203.0794]K.A", spec);
            assertEquals("LGEPTIDESK", sent.peptide, spec.name());
            assertEquals(List.of("S[203.0794]", "E[21.9819]"), sent.dropped, spec.name());
        }
        assertEquals("LGE[UNIMOD:30]PTIDES[UNIMOD:43]K", send("K.LGE[21.9819]PTIDES[203.0794]K.A", ALPHAPEPT).peptide,
                "AlphaPeptDeep is sent every mod it knows");
    }

    //the codes a model leaves out are not remembered as found, so the next peptide is written the same way
    @Test
    void aCodeLeftOutIsNotRememberedAsFound() {
        send("K.LGEPTIDEK[8.0142].A", IM2DEEP);
        assertEquals(Set.of(), new HashSet<>(PeptideFormatter.foundUnimods));
        send("K.LGEPTIDEK[8.0142].A", DEEPLC);
        assertEquals(Set.of("K259"), new HashSet<>(PeptideFormatter.foundUnimods));
        assertEquals("LGEPTIDEK", send("K.LGEPTIDEK[8.0142].A", IM2DEEP).peptide,
                "a code another model found is still one IM2Deep cannot read");
    }

    private record Sent(String peptide, List<String> dropped) {}

    private static Sent send(String pinPeptide, KoinaModelSpec spec) {
        List<String> dropped = new ArrayList<>();
        String peptide = new PeptideFormatter(pinPeptide, 2, "pin").getModel(spec, dropped);
        return new Sent(peptide, dropped);
    }

    private static String readBack(String modelPeptide, KoinaModelSpec spec) {
        return new PeptideFormatter(modelPeptide, 2, spec).getBase();
    }
}
