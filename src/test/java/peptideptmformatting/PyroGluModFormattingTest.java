/*
 * This file is part of MSBooster.
 *
 * MSBooster is free software: you can redistribute it and/or modify it under
 * the terms of the GNU Lesser General Public License as published by the Free
 * Software Foundation, either version 3 of the License, or (at your
 * option) any later version.
 *
 * MSBooster is distributed in the hope that it will be useful, but WITHOUT
 * ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
 * FITNESS FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public License
 * for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with MSBooster. If not, see <https://www.gnu.org/licenses/>.
 */

package peptideptmformatting;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertFalse;
import static org.junit.jupiter.api.Assertions.assertTrue;

import java.util.HashSet;
import java.util.Map;

import org.junit.jupiter.api.BeforeEach;
import org.junit.jupiter.api.Test;

/**
 * Regression tests for N-terminal pyro-Glu ProForma serialization.
 *
 * <p>Historically (v1.4.14) the alphapeptdeep-name parser built its
 * amino-acid/unimod lookup keys with a string replacement on the localization
 * suffix: {@code "Q^Any_N-term".replace("Any_N-term", "[")} yielded {@code "Q^["},
 * so the map key for Gln->pyro-Glu (UNIMOD:28) became {@code "Q^[28"} instead of
 * {@code "[28"}. When that key was chosen for a peptide, the serializer emitted the
 * malformed token {@code Q[UNIMOD:^[28]...}, an unclosed ProForma group that Koina
 * (and downstream parsing) could not read. The closely related Glu->pyro-Glu
 * (UNIMOD:27) only escaped by luck of hash ordering. Fixed in v1.4.16 by collapsing
 * any "term" localization to a bare {@code "["}.
 *
 * <p>These tests pin the corrected behaviour and guard the underlying map invariant
 * so the malformed key can never reappear.
 */
class PyroGluModFormattingTest {

    @BeforeEach
    void resetFoundUnimods() {
        // foundUnimods is static/shared; clear it so each case exercises the
        // map-driven lookup path deterministically, independent of test order.
        PeptideFormatter.foundUnimods.clear();
    }

    @Test
    void pyroGluQ_serializesToUnimod28() {
        PeptideFormatter pf = new PeptideFormatter("Q[-17.0265]IDNARLAADDF", 2, "base");
        assertEquals("Q[UNIMOD:28]IDNARLAADDF", pf.getAlphapept());
    }

    @Test
    void pyroGluE_serializesToUnimod27() {
        PeptideFormatter pf = new PeptideFormatter("E[-18.0106]IDNARLAADDF", 2, "base");
        assertEquals("E[UNIMOD:27]IDNARLAADDF", pf.getAlphapept());
    }

    @Test
    void pyroGluQ_hasNoStrayCaretOrBracket() {
        // The exact shape reported in the issue: Q[UNIMOD:^[28]GVNDNEEGFFSAR
        PeptideFormatter pf = new PeptideFormatter("Q[-17.0265]GVNDNEEGFFSAR", 2, "base");
        String alphapept = pf.getAlphapept();

        assertEquals("Q[UNIMOD:28]GVNDNEEGFFSAR", alphapept);
        assertFalse(alphapept.contains("^"), "stray localization separator leaked: " + alphapept);
        assertFalse(alphapept.contains("[UNIMOD:^"), "malformed UNIMOD token: " + alphapept);
        assertBalancedBrackets(alphapept);
    }

    @Test
    void formatBaseToSpecific_directCall_isWellFormed() {
        // Exercise the serializer directly with a fresh, empty foundUnimods set so the
        // result depends only on the AAunimodToModMassAll map, not on prior state.
        String base = "Q[-17.0265]IDNARLAADDF";
        String[] result = PTMhandler.formatPeptideBaseToSpecific(
                base, base.indexOf('['), base.indexOf(']'), "alphapept", new HashSet<>(), false);

        assertEquals("Q[UNIMOD:28]IDNARLAADDF", result[0]);
        assertEquals("Q28", result[1]);
    }

    @Test
    void aaUnimodMap_keysNeverContainLocalizationSeparator() {
        // Root-cause invariant: an alphapeptdeep localization separator ('^') must
        // never survive into a lookup key. If it does, serialization produces the
        // malformed [UNIMOD:^[...] token again.
        for (String key : PTMhandler.AAunimodToModMassAll.keySet()) {
            assertFalse(key.contains("^"),
                    "AAunimodToModMassAll contains a malformed key: [" + key + "]");
        }
    }

    @Test
    void pyroGluQ_roundTripsBackToBaseWithoutCrashing() {
        // The reported client-side IndexOutOfBoundsException happened while parsing the
        // malformed peptide string back to base form. A well-formed token must survive
        // the reverse conversion and recover the modification mass.
        Map<Integer, Double> modMap = PTMhandler.unimodToModMassAll;
        String base = PTMhandler.formatPeptideSpecificToBase("Q[UNIMOD:28]IDNARLAADDF", modMap, "[]");

        assertTrue(base.startsWith("Q[-17.026"), "unexpected round-trip base: " + base);
        assertTrue(base.endsWith("]IDNARLAADDF"), "unexpected round-trip base: " + base);
        assertFalse(base.contains("^"), "stray localization separator after round-trip: " + base);
        assertBalancedBrackets(base);
    }

    private static void assertBalancedBrackets(String peptide) {
        int depth = 0;
        for (int i = 0; i < peptide.length(); i++) {
            char c = peptide.charAt(i);
            if (c == '[') {
                depth++;
            } else if (c == ']') {
                depth--;
            }
            assertTrue(depth >= 0 && depth <= 1,
                    "unbalanced/nested brackets at index " + i + " in: " + peptide);
        }
        assertEquals(0, depth, "unclosed bracket group in: " + peptide);
    }
}
