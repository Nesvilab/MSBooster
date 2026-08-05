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

import org.junit.jupiter.api.Test;

/**
 * Transcribed from FragAlign's {@code src/spline.rs} unit tests, plus the batch
 * coverage from {@code tests/rt_regression.rs}.
 *
 * <p>The fixture is the line y = 2x + 1 over knots at x = 0 and x = 10: value 1
 * at x = 0, value 21 at x = 10, slope 2 at both.
 */
public class HermiteSplineTest {

    private static HermiteSpline line() {
        return new HermiteSpline(new double[]{0.0, 10.0}, new double[]{1.0, 2.0, 21.0, 2.0});
    }

    @Test
    public void evaluatesLinearlyInsideTheFittedRange() {
        HermiteSpline s = line();
        assertEquals(11.0, s.evaluate(5.0), 1e-9);
        assertEquals(1.0, s.evaluate(0.0), 1e-9);
        assertEquals(21.0, s.evaluate(10.0), 1e-9);
    }

    @Test
    public void extrapolatesLinearlyOutsideTheFittedRange() {
        HermiteSpline s = line();
        assertEquals(-9.0, s.evaluate(-5.0), 1e-9);
        assertEquals(31.0, s.evaluate(15.0), 1e-9);
    }

    @Test
    public void invertRecoversTheAbscissaForAKnownValue() {
        assertEquals(5.0, line().invert(11.0), 1e-3);
    }

    @Test
    public void invertClampsOutsideTheFittedRange() {
        HermiteSpline s = line();
        assertEquals(0.0, s.invert(-100.0), 0.0);
        assertEquals(10.0, s.invert(1000.0), 0.0);
    }

    @Test
    public void exposesTheFittedRangeBounds() {
        HermiteSpline s = line();
        assertEquals(0.0, s.xMin(), 0.0);
        assertEquals(10.0, s.xMax(), 0.0);
    }

    @Test
    public void batchEvaluateMatchesTheScalarLoop() {
        HermiteSpline s = line();
        double[] xs = {-3.0, 0.0, 2.5, 7.5, 10.0, 12.0};
        double[] expected = new double[xs.length];
        for (int i = 0; i < xs.length; i++) {
            expected[i] = s.evaluate(xs[i]);
        }
        assertArrayEquals(expected, s.evaluateAll(xs), 0.0);
    }

    @Test
    public void batchInvertMatchesTheScalarLoop() {
        HermiteSpline s = line();
        double[] values = {1.0, 6.0, 11.0, 21.0};
        double[] expected = new double[values.length];
        for (int i = 0; i < values.length; i++) {
            expected[i] = s.invert(values[i]);
        }
        assertArrayEquals(expected, s.invertAll(values), 0.0);
    }
}
