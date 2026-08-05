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

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertSame;
import static org.junit.jupiter.api.Assertions.assertTrue;

import org.junit.jupiter.api.Test;

/**
 * Covers the Fritsch-Carlson monotone-preserving interpolant and the
 * monotonicity guard from FragAlign's {@code regression.rs}
 * ({@code monotone_pchip} and {@code enforce_monotonic}).
 */
public class MonotonePchipTest {

    private static void assertNonDecreasing(HermiteSpline s, double x0, double x1) {
        double previous = Double.NEGATIVE_INFINITY;
        for (int k = 0; k < 400; k++) {
            double x = x0 + (x1 - x0) * k / 399.0;
            double v = s.evaluate(x);
            assertTrue(v >= previous - 1e-9, "descends at x = " + x);
            previous = v;
        }
    }

    @Test
    public void interpolantPassesThroughItsKnots() {
        double[] kx = {0.0, 1.0, 2.0, 3.0, 4.0};
        double[] ky = {0.0, 3.0, 3.0, 8.0, 20.0};
        HermiteSpline s = MonotonePchip.through(kx, ky);
        for (int i = 0; i < kx.length; i++) {
            assertEquals(ky[i], s.evaluate(kx[i]), 1e-9);
        }
    }

    @Test
    public void interpolantDoesNotOvershootBetweenKnots() {
        // A flat segment between two steep ones is where a plain cubic Hermite
        // would overshoot; the Fritsch-Carlson limiter must prevent it.
        double[] kx = {0.0, 1.0, 2.0, 3.0};
        double[] ky = {0.0, 10.0, 10.0, 20.0};
        assertNonDecreasing(MonotonePchip.through(kx, ky), 0.0, 3.0);
    }

    @Test
    public void alreadyMonotoneSplineIsReturnedUnchanged() {
        HermiteSpline line = new HermiteSpline(new double[]{0.0, 10.0}, new double[]{1.0, 2.0, 21.0, 2.0});
        assertSame(line, MonotonePchip.enforce(line, 0.0, 10.0));
    }

    @Test
    public void nonMonotoneSplineIsRebuiltAsItsRunningMax() {
        // Rises to 10 at x = 5, then falls back to 5 at x = 10 — the spurious
        // tail fall-back the guard exists to remove.
        HermiteSpline falling = new HermiteSpline(
                new double[]{0.0, 5.0, 10.0},
                new double[]{0.0, 2.0, 10.0, 0.0, 5.0, -2.0});

        HermiteSpline fixed = MonotonePchip.enforce(falling, 0.0, 10.0);

        assertNonDecreasing(fixed, 0.0, 10.0);
        // The running max is taken over 512 discrete samples, so the captured
        // peak sits a sampling interval below the true continuous one.
        assertEquals(10.0, fixed.evaluate(10.0), 1e-4, "tail must hold at the running max");
        assertEquals(10.0, fixed.evaluate(5.0), 1e-4, "the rise up to the peak is preserved");
    }

    @Test
    public void degenerateRangeLeavesTheSplineAlone() {
        HermiteSpline line = new HermiteSpline(new double[]{0.0, 10.0}, new double[]{1.0, 2.0, 21.0, 2.0});
        assertSame(line, MonotonePchip.enforce(line, 5.0, 5.0));
    }
}
