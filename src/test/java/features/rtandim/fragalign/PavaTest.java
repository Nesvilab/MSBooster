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

package features.rtandim.fragalign;

import static org.junit.jupiter.api.Assertions.assertArrayEquals;
import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertTrue;

import org.junit.jupiter.api.Test;

/**
 * Transcribed from FragAlign's {@code src/pava.rs} unit tests so the Java port
 * is held to the same behaviour as the Rust reference.
 */
public class PavaTest {

    @Test
    public void monotoneInputPassesThroughUnchanged() {
        double[] y = {1.0, 2.0, 3.0, 4.0, 5.0};
        assertArrayEquals(new double[]{1.0, 2.0, 3.0, 4.0, 5.0},
                Pava.poolAdjacentViolators(y, null), 0.0);
    }

    @Test
    public void singleViolationPoolsAdjacentPairToTheirMean() {
        // 3.0 followed by 1.0 violates monotonicity, so the pair pools to 2.0.
        double[] out = Pava.poolAdjacentViolators(new double[]{1.0, 3.0, 1.0, 4.0}, null);
        assertArrayEquals(new double[]{1.0, 2.0, 2.0, 4.0}, out, 0.0);
        for (int i = 1; i < out.length; i++) {
            assertTrue(out[i] >= out[i - 1], "output must be non-decreasing");
        }
    }

    @Test
    public void descendingRunCollapsesToASingleBlock() {
        double[] out = Pava.poolAdjacentViolators(new double[]{4.0, 3.0, 2.0, 1.0}, null);
        assertArrayEquals(new double[]{2.5, 2.5, 2.5, 2.5}, out, 0.0);
    }

    @Test
    public void smallWeightBarelyMovesThePooledMean() {
        double[] y = {10.0, 0.0};
        double[] out = Pava.poolAdjacentViolators(y, new double[]{1.0, 1e-6});

        double expected = (10.0 * 1.0 + 0.0 * 1e-6) / (1.0 + 1e-6);
        assertEquals(expected, out[0], 1e-12);
        assertEquals(expected, out[1], 1e-12);
        assertTrue(out[0] > 9.99, "heavy anchor must dominate the merge");

        // The unweighted case pools to the plain midpoint instead.
        assertArrayEquals(new double[]{5.0, 5.0},
                Pava.poolAdjacentViolators(y, null), 0.0);
    }

    @Test
    public void allZeroWeightsFallBackToUnweightedPooling() {
        double[] out = Pava.poolAdjacentViolators(new double[]{3.0, 1.0}, new double[]{0.0, 0.0});
        assertArrayEquals(new double[]{2.0, 2.0}, out, 0.0);
    }

    @Test
    public void emptyInputYieldsEmptyOutput() {
        assertEquals(0, Pava.poolAdjacentViolators(new double[0], null).length);
    }
}
