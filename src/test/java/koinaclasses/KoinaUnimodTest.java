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
import static org.junit.jupiter.api.Assertions.assertThrows;
import static org.junit.jupiter.api.Assertions.assertTrue;

import java.util.List;
import koinaclasses.KoinaModelSpec.EncodableMods;
import koinaclasses.KoinaUnimod.Token;
import org.junit.jupiter.api.Test;

/**
 * The UniMod compositions Koina reads, and the mods DeepLC and IM2Deep can be sent by them.
 */
class KoinaUnimodTest {
    @Test
    void theCompositionsAreKoinas() {
        assertEquals("H(2) C(2) O", KoinaUnimod.composition(1));
        assertEquals("H(20) C(8) 13C(4) N 15N O(2)", KoinaUnimod.composition(737));
        assertEquals("HexNAc", KoinaUnimod.composition(43));
        assertEquals("H(2) C(2) O", KoinaUnimod.composition("1"));
        assertNull(KoinaUnimod.composition(0), "the root term has no composition");
        assertNull(KoinaUnimod.composition(2114), "lactylation is newer than Koina's unimod.obo");
        assertNull(KoinaUnimod.composition("Mymod"), "a custom mod without a UniMod code");
        assertNull(KoinaUnimod.composition(""));
    }

    @Test
    void tokensAreReadAsKoinasZeroShotEncodingReadsThem() {
        assertEquals(List.of(new Token("", "C", -6), new Token("13", "C", 6), new Token("", "N", -2),
                new Token("15", "N", 2)), KoinaUnimod.tokens("C(-6) 13C(6) N(-2) 15N(2)"));
        assertEquals(List.of(new Token("", "Hex", 5), new Token("", "HexNAc", 2)),
                KoinaUnimod.tokens("Hex(5) HexNAc(2)"));
        assertEquals(List.of(new Token("", "H", 1), new Token("", "Na", 1)), KoinaUnimod.tokens("H Na"));
        assertTrue(KoinaUnimod.tokens("2H(4) C").get(0).isIsotope());
        assertThrows(IllegalArgumentException.class, () -> KoinaUnimod.tokens("H(2) c(2)"));
    }

    //every composition of Koina's file is read but those with dHex, a symbol the token pattern of Koina's zero-shot
    //encoding does not read either, as it starts with a lower-case letter
    @Test
    void everyCompositionOfTheFileIsReadableButThoseWithDHex() {
        int read = 0;
        for (int id = 0; id < 3000; id++) {
            String composition = KoinaUnimod.composition(id);
            if (composition == null) {
                continue;
            }
            if (composition.matches("(.* )?dHex(\\(\\d+\\))?( .*)?")) {
                assertThrows(IllegalArgumentException.class, () -> KoinaUnimod.tokens(composition), "UNIMOD:" + id);
            } else {
                assertFalse(KoinaUnimod.tokens(composition).isEmpty(), "UNIMOD:" + id);
                read++;
            }
        }
        assertTrue(read > 1000, read + " compositions read");
        assertFalse(EncodableMods.KOINA_ELEMENTS.encodes("N142"), "HexNAc(1)dHex(1)");
    }

    @Test
    void deepLcTakesTheSixElementsAndTheirIsotopesAndIm2DeepOnlyTheElements() {
        //residue (or "[") plus code, as in the mod sets
        String[] both = {"[1", "M35", "C4", "S21", "K121", "Q28", "E27", "C26", "N7", "[2"};
        String[] deepLcOnly = {"K737", "[737", "K259", "R267", "K2016", "[214"};
        String[] neither = {"E30", "S43", "N43", "K999999", "KMymod", "M0"};
        for (String code : both) {
            assertTrue(EncodableMods.KOINA_ELEMENTS.encodes(code), code);
            assertTrue(EncodableMods.KOINA_ELEMENTS_NO_ISOTOPES.encodes(code), code);
        }
        for (String code : deepLcOnly) {
            assertTrue(EncodableMods.KOINA_ELEMENTS.encodes(code), code);
            assertFalse(EncodableMods.KOINA_ELEMENTS_NO_ISOTOPES.encodes(code), code);
        }
        for (String code : neither) {
            assertFalse(EncodableMods.KOINA_ELEMENTS.encodes(code), code);
            assertFalse(EncodableMods.KOINA_ELEMENTS_NO_ISOTOPES.encodes(code), code);
            assertTrue(EncodableMods.ANY.encodes(code), code);
        }
    }

    @Test
    void theModelsHaveTheirRules() {
        assertEquals(EncodableMods.KOINA_ELEMENTS, KoinaModelRegistry.require("Deeplc_hela_hf").encodableMods());
        assertEquals(EncodableMods.KOINA_ELEMENTS_NO_ISOTOPES, KoinaModelRegistry.require("IM2Deep").encodableMods());
        for (KoinaModelSpec spec : KoinaModelRegistry.all()) {
            if (spec.modRule() == KoinaModelSpec.ModRule.PROSIT_PTM) { //see KoinaPrositPtmTest
                assertEquals(EncodableMods.PROSIT_PTM, spec.encodableMods(), spec.name());
            } else if (spec.modRule() == KoinaModelSpec.ModRule.LAC
                    || spec.modRule() == KoinaModelSpec.ModRule.CHRONOLOGER) { //KoinaPrositLacTest, KoinaChronologerTest
                assertEquals(EncodableMods.NO_CTERM, spec.encodableMods(), spec.name());
            } else if (!spec.name().equals("Deeplc_hela_hf") && !spec.name().equals("IM2Deep")) {
                assertEquals(EncodableMods.ANY, spec.encodableMods(), spec.name());
            }
            assertEquals(spec.name().equals("Deeplc_hela_hf") || spec.name().equals("IM2Deep"),
                    spec.peptideFormat().termModsOnResidues(), spec.name());
        }
    }
}
