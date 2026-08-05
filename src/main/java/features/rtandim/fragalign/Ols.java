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

import org.apache.commons.math3.linear.Array2DRowRealMatrix;
import org.apache.commons.math3.linear.RealMatrix;
import org.apache.commons.math3.linear.SingularValueDecomposition;

/**
 * Ordinary least squares with no intercept term.
 *
 * <p>Port of FragAlign's {@code src/ols.rs}. The design matrix is used exactly
 * as supplied; no constant column is appended.
 */
public final class Ols {

    /**
     * Singular values at or below this are treated as zero and their direction
     * dropped from the solution. Matches the cutoff the reference passes to
     * nalgebra's SVD solve.
     */
    private static final double SINGULAR_VALUE_CUTOFF = 1e-12;

    private Ols() {}

    /**
     * Solves {@code min_x ||A x - b||} by truncated singular value decomposition.
     *
     * @param a design matrix as {@code n} rows of length {@code p}
     * @param b target vector of length {@code n}
     * @return the {@code p} coefficients, or {@code null} when the solve cannot
     *         be trusted — an empty or ragged matrix, zero-width rows, a target
     *         length mismatch, a decomposition failure, or any non-finite input
     *         or coefficient. The caller falls back to the LOESS path.
     */
    public static double[] solveNoIntercept(double[][] a, double[] b) {
        int n = a.length;
        if (n == 0) {
            return null;
        }
        int p = a[0].length;
        if (p == 0) {
            return null;
        }
        for (double[] row : a) {
            if (row.length != p) {
                return null;
            }
        }
        if (n != b.length) {
            return null;
        }
        if (!allFinite(b) || !allFinite(a)) {
            return null;
        }

        double[] x = new double[p];
        try {
            SingularValueDecomposition svd =
                    new SingularValueDecomposition(new Array2DRowRealMatrix(a, false));
            RealMatrix u = svd.getU();
            RealMatrix v = svd.getV();
            double[] singularValues = svd.getSingularValues();

            // Truncated pseudo-inverse: x = sum_i (u_i . b / s_i) v_i, skipping
            // directions whose singular value falls at or below the cutoff.
            for (int i = 0; i < singularValues.length; i++) {
                double s = singularValues[i];
                if (s <= SINGULAR_VALUE_CUTOFF) {
                    continue;
                }
                double dot = 0.0;
                for (int r = 0; r < n; r++) {
                    dot += u.getEntry(r, i) * b[r];
                }
                double scale = dot / s;
                for (int j = 0; j < p; j++) {
                    x[j] += scale * v.getEntry(j, i);
                }
            }
        } catch (RuntimeException e) {
            return null;
        }

        return allFinite(x) ? x : null;
    }

    private static boolean allFinite(double[] values) {
        for (double v : values) {
            if (!Double.isFinite(v)) {
                return false;
            }
        }
        return true;
    }

    private static boolean allFinite(double[][] rows) {
        for (double[] row : rows) {
            if (!allFinite(row)) {
                return false;
            }
        }
        return true;
    }
}
