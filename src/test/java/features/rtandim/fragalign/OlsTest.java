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
import static org.junit.jupiter.api.Assertions.assertNull;

import org.junit.jupiter.api.Test;

/**
 * Transcribed from FragAlign's {@code src/ols.rs} unit tests.
 *
 * <p>A {@code null} return is the Java equivalent of the reference's
 * {@code None}: it tells the caller to fall back to the LOESS path.
 */
public class OlsTest {

    @Test
    public void recoversAnExactLinearRelation() {
        // Fit y = 2x + 3 with design rows [1, x] — the intercept is supplied
        // explicitly as a constant column, since the solve adds none.
        double[] xs = {0.0, 1.0, 2.0, 3.0, 4.0, 5.0};
        double[][] a = new double[xs.length][];
        double[] b = new double[xs.length];
        for (int i = 0; i < xs.length; i++) {
            a[i] = new double[]{1.0, xs[i]};
            b[i] = 2.0 * xs[i] + 3.0;
        }

        double[] coeffs = Ols.solveNoIntercept(a, b);
        assertEquals(2, coeffs.length);
        assertEquals(3.0, coeffs[0], 1e-9);
        assertEquals(2.0, coeffs[1], 1e-9);
    }

    @Test
    public void fitsASingleColumnThroughTheOrigin() {
        double[] xs = {1.0, 2.0, 3.0, 4.0};
        double[][] a = new double[xs.length][];
        double[] b = new double[xs.length];
        for (int i = 0; i < xs.length; i++) {
            a[i] = new double[]{xs[i]};
            b[i] = 4.0 * xs[i];
        }

        double[] coeffs = Ols.solveNoIntercept(a, b);
        assertEquals(1, coeffs.length);
        assertEquals(4.0, coeffs[0], 1e-9);
    }

    @Test
    public void rejectsAnEmptyDesignMatrix() {
        assertNull(Ols.solveNoIntercept(new double[0][], new double[0]));
    }

    @Test
    public void rejectsZeroWidthRows() {
        assertNull(Ols.solveNoIntercept(new double[][]{{}, {}}, new double[]{1.0, 2.0}));
    }

    @Test
    public void rejectsRaggedRows() {
        double[][] a = {{1.0, 2.0}, {1.0}};
        assertNull(Ols.solveNoIntercept(a, new double[]{1.0, 2.0}));
    }

    @Test
    public void rejectsATargetLengthMismatch() {
        double[][] a = {{1.0, 0.0}, {1.0, 1.0}};
        assertNull(Ols.solveNoIntercept(a, new double[]{1.0, 2.0, 3.0}));
    }

    @Test
    public void rejectsANonFiniteTarget() {
        // A NaN target propagates into every coefficient; the reference rejects
        // the solve rather than returning a spline full of NaN.
        double[][] a = {{1.0, 0.0}, {1.0, 1.0}, {1.0, 2.0}};
        assertNull(Ols.solveNoIntercept(a, new double[]{1.0, Double.NaN, 3.0}));
    }
}
