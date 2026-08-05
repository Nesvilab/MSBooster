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

import java.util.Arrays;
import java.util.Comparator;
import java.util.stream.IntStream;

/**
 * Calibration-curve regression: fits a piecewise cubic Hermite spline to
 * {@code (x, y)} anchor pairs.
 *
 * <p>Port of FragAlign's {@code src/regression.rs}. Two paths are tried in
 * order:
 *
 * <ol>
 *   <li><b>PAVA + Hermite-basis least squares.</b> Weighted PAVA monotonises the
 *       values, knots are placed at equally-spaced data positions, and a
 *       Hermite-basis design matrix is solved by {@link Ols}. Edges extrapolate
 *       linearly rather than going flat.</li>
 *   <li><b>LOESS + PAV + PCHIP.</b> Used when the least-squares solve fails; see
 *       {@link LoessFallback}.</li>
 * </ol>
 *
 * <p>Two robustness steps wrap the fit: {@link OutlierMask} drops gradient-end
 * misassignments before fitting, and {@link MonotonePchip#enforce} guarantees
 * the returned curve is non-decreasing.
 *
 * <p>FragAlign names the axes iRT and RT. MSBooster fits in both directions, so
 * the axes are neutral here — {@code x} is the independent variable and
 * {@code y} the dependent one, whichever way round the caller supplies them.
 */
public final class FragAlignRegression {

    /** Below this many anchors no fit is attempted. */
    private static final int MIN_ANCHORS = 5;

    /** Upper bound on the number of spline segments. */
    private static final int MAX_SEGMENTS = 20;

    /**
     * Quantile of the anchor distribution the outermost knots are pinned to.
     *
     * <p>Knots otherwise sit at interval midpoints, which leaves the fitted
     * domain strictly inside the data at both ends and truncates an inverse
     * lookup, since {@link HermiteSpline#invert} clamps to that domain. Pinning
     * to the literal extremes fixes the truncation but hands the ceiling to
     * whatever single anchor is furthest out. Pinning just inside them keeps the
     * domain on the body of the distribution: on the CCRCC run the top anchor
     * (predRT 169.8) sat 15 units clear of its nearest neighbour (154.5), and
     * this quantile excludes it while keeping the genuine late-eluting cluster.
     */
    private static final double KNOT_EDGE_QUANTILE = 0.999;

    private FragAlignRegression() {}

    /**
     * Sorts the anchors by ascending {@code x}, drops backward misassignments,
     * de-ties the abscissae and fits a spline.
     *
     * @param x       independent variable, any order
     * @param y       dependent variable, parallel to {@code x}
     * @param weights per-anchor non-negative weights, or {@code null}
     * @return the fitted spline, or {@code null} when fewer than
     *         {@value #MIN_ANCHORS} anchors are supplied or both fitting paths
     *         fail
     * @throws IllegalArgumentException if the array lengths disagree
     */
    public static HermiteSpline fit(double[] x, double[] y, double[] weights) {
        if (x.length != y.length) {
            throw new IllegalArgumentException(
                    "x.length (" + x.length + ") != y.length (" + y.length + ")");
        }
        if (weights != null && weights.length != x.length) {
            throw new IllegalArgumentException(
                    "weights.length (" + weights.length + ") != x.length (" + x.length + ")");
        }

        int n = x.length;
        if (n < MIN_ANCHORS) {
            return null;
        }

        // Stable sort by x ascending. Double.compare gives the same total order
        // the reference gets from total_cmp, so NaN abscissae land at the end
        // deterministically.
        int[] order = IntStream.range(0, n).boxed()
                .sorted(Comparator.comparingDouble(i -> x[i]))
                .mapToInt(Integer::intValue).toArray();

        double[] xs = new double[n];
        double[] ys = new double[n];
        double[] ws = weights == null ? null : new double[n];
        for (int i = 0; i < n; i++) {
            xs[i] = x[order[i]];
            ys[i] = y[order[i]];
            if (ws != null) {
                ws[i] = weights[order[i]];
            }
        }

        // Drop anchors whose value sits far below the local high level. At the
        // end of the gradient the relation is multivalued, and those points
        // would otherwise average the fit down through the cloud.
        boolean[] keep = OutlierMask.keep(ys);
        int kept = 0;
        for (boolean k : keep) {
            if (k) {
                kept++;
            }
        }
        if (kept >= MIN_ANCHORS && kept < n) {
            double[] fx = new double[kept];
            double[] fy = new double[kept];
            double[] fw = ws == null ? null : new double[kept];
            int j = 0;
            for (int i = 0; i < n; i++) {
                if (!keep[i]) {
                    continue;
                }
                fx[j] = xs[i];
                fy[j] = ys[i];
                if (fw != null) {
                    fw[j] = ws[i];
                }
                j++;
            }
            xs = fx;
            ys = fy;
            ws = fw;
            n = kept;
        }

        // Enforce strictly increasing abscissae.
        for (int i = 1; i < n; i++) {
            if (xs[i] <= xs[i - 1]) {
                xs[i] = xs[i - 1] + 1e-8;
            }
        }

        return fitSorted(xs, ys, ws);
    }

    /**
     * Fits pre-sorted, strictly-increasing anchors. Prefer {@link #fit} unless
     * the caller has already sorted and de-tied the input.
     *
     * @return the fitted spline, or {@code null} if both paths fail
     */
    public static HermiteSpline fitSorted(double[] x, double[] y, double[] weights) {
        HermiteSpline spline = fitPavaLeastSquares(x, y, weights);
        if (spline == null) {
            spline = LoessFallback.fit(x, y, weights);
        }
        if (spline == null) {
            return null;
        }
        return MonotonePchip.enforce(spline, x[0], x[x.length - 1]);
    }

    /**
     * PAVA monotonisation followed by a Hermite-basis least-squares solve.
     *
     * <p>Edges use a linear basis (value and slope) rather than a cubic one, so
     * the fitted curve extrapolates instead of kinking flat at the boundary.
     *
     * @return the fitted spline, or {@code null} if the solve fails
     */
    private static HermiteSpline fitPavaLeastSquares(double[] x, double[] y, double[] weights) {
        int n = x.length;
        if (n < MIN_ANCHORS) {
            return null;
        }

        double[] w = normalizeWeights(weights, n);
        double[] pava = Pava.poolAdjacentViolators(y, weights == null ? null : w);

        double[] points = knotPositions(x, n);
        int m = points.length;
        int nCoeffs = 2 * m;

        double[][] a = new double[n][nCoeffs];
        for (int i = 0; i < n; i++) {
            double xi = x[i];
            int k = 0;
            while (k < m && points[k] < xi) {
                k++;
            }

            if (k == 0 && xi <= points[0]) {
                // Left edge: linear basis.
                a[i][0] = 1.0;
                a[i][1] = xi - points[0];
            } else if (k == m) {
                // Right edge: linear basis.
                a[i][(m - 1) * 2] = 1.0;
                a[i][(m - 1) * 2 + 1] = xi - points[m - 1];
            } else {
                double segW = points[k] - points[k - 1];
                if (segW < 1e-12) {
                    a[i][(k - 1) * 2] = 1.0;
                    continue;
                }
                double u = (xi - points[k - 1]) / segW;
                double v = u - 1.0;
                a[i][(k - 1) * 2] = (1.0 + 2.0 * u) * v * v;
                a[i][(k - 1) * 2 + 1] = segW * u * v * v;
                a[i][k * 2] = u * u * (1.0 - 2.0 * v);
                a[i][k * 2 + 1] = segW * u * u * v;
            }
        }

        // Weighted-least-squares reduction: scale each row and target by sqrt(w).
        boolean nonUniform = false;
        for (double wi : w) {
            if (Math.abs(wi - 1.0) > 1e-12) {
                nonUniform = true;
                break;
            }
        }
        double[] target = pava;
        if (nonUniform) {
            target = new double[n];
            double[][] scaled = new double[n][nCoeffs];
            for (int i = 0; i < n; i++) {
                double s = Math.sqrt(w[i]);
                target[i] = pava[i] * s;
                for (int j = 0; j < nCoeffs; j++) {
                    scaled[i][j] = a[i][j] * s;
                }
            }
            a = scaled;
        }

        double[] coeffs = Ols.solveNoIntercept(a, target);
        return coeffs == null ? null : new HermiteSpline(points, coeffs);
    }

    /**
     * Places knots at the midpoints of equally-spaced <em>data positions</em>,
     * so knot density follows the density of the anchors rather than the span of
     * x. Duplicate knots are collapsed.
     */
    private static double[] knotPositions(double[] x, int n) {
        int segments = Math.min(MAX_SEGMENTS, Math.max(1, (int) (2.0 * Math.sqrt(n / 20.0))));
        int m = segments;

        double split = (n - 1) / (double) m;
        double[] seg = new double[m + 1];
        for (int i = 0; i < m; i++) {
            seg[i] = x[Math.min((int) (split * i), n - 1)];
        }
        seg[m] = x[n - 1];

        double[] points = new double[m];
        for (int i = 0; i < m; i++) {
            points[i] = 0.5 * (seg[i] + seg[i + 1]);
        }

        // Pin the outermost knots near the data extremes. Interval midpoints
        // leave the fitted domain strictly inside the data at both ends — on a
        // 300-point range spanning 0..299 the first knot lands at 21 — which
        // truncates the curve. evaluate() hides this by extrapolating linearly,
        // but invert() clamps to the fitted range, so an inverse lookup
        // saturates before it reaches the last anchor. On the CCRCC run that
        // cost the top 24 iRT units of the late-eluting cluster. The pins sit
        // just inside the extremes: see KNOT_EDGE_QUANTILE.
        if (m >= 2) {
            double span = n - 1;
            int lo = (int) ((1.0 - KNOT_EDGE_QUANTILE) * span);
            int hi = (int) (KNOT_EDGE_QUANTILE * span);
            points[0] = x[Math.min(lo, n - 1)];
            points[m - 1] = x[Math.min(hi, n - 1)];
        }

        // Collapse duplicates, re-testing the shifted-in value at the same index.
        int i = 1;
        while (i < m) {
            if (points[i] - points[i - 1] < 1e-10) {
                System.arraycopy(points, i + 1, points, i, m - i - 1);
                m--;
            } else {
                i++;
            }
        }
        return Arrays.copyOf(points, m);
    }

    /**
     * Scales weights to a maximum of one. A {@code null} array, or one whose
     * maximum is non-positive, yields all ones.
     */
    static double[] normalizeWeights(double[] weights, int n) {
        double[] w = new double[n];
        if (weights == null) {
            Arrays.fill(w, 1.0);
            return w;
        }

        double wMax = 0.0;
        for (double v : weights) {
            if (v > wMax) {
                wMax = v;
            }
        }
        if (wMax <= 0.0) {
            Arrays.fill(w, 1.0);
            return w;
        }

        for (int i = 0; i < n; i++) {
            w[i] = Math.max(0.0, weights[i]) / wMax;
        }
        return w;
    }
}
