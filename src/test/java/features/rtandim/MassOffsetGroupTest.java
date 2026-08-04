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

package features.rtandim;

import static org.junit.jupiter.api.Assertions.assertArrayEquals;
import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertFalse;
import static org.junit.jupiter.api.Assertions.assertTrue;

import org.junit.jupiter.api.Test;

public class MassOffsetGroupTest {

    @Test
    public void parsesDeltaMassesFromBaseFormatName() {
        assertArrayEquals(new double[]{57.0215, -14.0157},
                MassOffsetGroup.deltaMasses("AC[57.0215]DEK[-14.0157]|2"));
    }

    @Test
    public void unmodifiedPeptideHasNoDeltaMasses() {
        assertEquals(0, MassOffsetGroup.deltaMasses("PEPTIDEK|2").length);
        assertEquals(0, MassOffsetGroup.deltaMasses(null).length);
    }

    @Test
    public void skipsBracketedTextThatIsNotAMass() {
        assertArrayEquals(new double[]{15.9949},
                MassOffsetGroup.deltaMasses("PEM[Oxidation]TID[15.9949]EK|2"));
    }

    @Test
    public void toleratesUnclosedBracket() {
        assertArrayEquals(new double[]{57.0215},
                MassOffsetGroup.deltaMasses("AC[57.0215]DEK[42.010|2"));
    }

    @Test
    public void matchesTheOffsetItWasBuiltWith() {
        MassOffsetGroup group = MassOffsetGroup.of("422.1624");
        assertTrue(group.matches("GSLVVVSS[422.1624]LLGR|2"));
        assertFalse(group.matches("GSLVVVSSLLGR|2"));
    }

    /**
     * The reason this class exists: MSFragger rounds the fragger.params offset and the localized
     * delta mass independently, so a pin can hold [79.9663] and [79.9664] for the same offset.
     */
    @Test
    public void matchesAcrossLastDecimalRoundingDisagreement() {
        MassOffsetGroup group = MassOffsetGroup.of("79.9663");
        assertTrue(group.matches("PEPS[79.9663]TIDEK|2"));
        assertTrue(group.matches("PEPS[79.9664]TIDEK|2"));
    }

    @Test
    public void doesNotMatchADistinctNearbyOffset() {
        //79.9568 and 79.9663 are separate entries of a real detailed offset list
        MassOffsetGroup group = MassOffsetGroup.of("79.9663");
        assertFalse(group.matches("PEPS[79.9568]TIDEK|2"));
    }

    /** A substring match would have let the positive offset claim the negative one's PSMs. */
    @Test
    public void signIsRespected() {
        assertFalse(MassOffsetGroup.of("14.0157").matches("PEPC[-14.0157]TIDEK|2"));
        assertTrue(MassOffsetGroup.of("-14.0157").matches("PEPC[-14.0157]TIDEK|2"));
    }

    /** A substring match would have let a short mass match inside a longer one. */
    @Test
    public void shortMassDoesNotMatchInsideLongerMass() {
        assertFalse(MassOffsetGroup.of("2.0157").matches("PEPK[42.0157]TIDEK|2"));
    }

    @Test
    public void ampersandGroupMatchesAnyOfItsMasses() {
        MassOffsetGroup group = MassOffsetGroup.of("541.0611&484.0396&-116.0586");
        assertTrue(group.matches("PEPD[541.0611]TIDEK|2"));
        assertTrue(group.matches("PEPC[-116.0586]TIDEK|2"));
        assertFalse(group.matches("PEPM[15.9949]TIDEK|2"));
    }

    @Test
    public void emptyKeyMatchesEveryPeptide() {
        MassOffsetGroup group = MassOffsetGroup.of("");
        assertTrue(group.matches("PEPTIDEK|2"));
        assertTrue(group.matches("PEPS[422.1624]TIDEK|2"));
    }

    @Test
    public void othersKeyMatchesNothingDirectly() {
        MassOffsetGroup group = MassOffsetGroup.of(MassOffsetGroup.OTHERS);
        assertFalse(group.matches("PEPTIDEK|2"));
        assertFalse(group.matches("PEPS[422.1624]TIDEK|2"));
    }

    /** A hand-written non-numeric massesForLoessCalibration keeps its historical substring match. */
    @Test
    public void nonNumericKeyFallsBackToSubstringMatch() {
        MassOffsetGroup group = MassOffsetGroup.of("Oxidation");
        assertTrue(group.matches("PEM[Oxidation]TIDEK|2"));
        assertFalse(group.matches("PEPTIDEK|2"));
    }

    /**
     * The tolerance has to sit between the two scales it separates: the 1e-4 disagreement between
     * MSFragger's two roundings of one offset, and the ~8e-3 smallest gap between distinct offsets
     * in a realistic list. Exactly 1e-3 is deliberately not asserted — it is not representable in
     * binary floating point, so it is not a case that can arise from real 4-decimal pin values.
     */
    @Test
    public void toleranceSeparatesRoundingJitterFromDistinctOffsets() {
        MassOffsetGroup group = MassOffsetGroup.of("100.0000");
        assertTrue(group.matches("PEPK[100.0001]TIDEK|2"));
        assertTrue(group.matches("PEPK[99.9999]TIDEK|2"));
        assertFalse(group.matches("PEPK[100.0080]TIDEK|2"));
        assertFalse(group.matches("PEPK[99.9920]TIDEK|2"));
    }
}
