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

import org.apache.commons.math3.analysis.interpolation.LoessInterpolator;
import org.apache.commons.math3.analysis.polynomials.PolynomialSplineFunction;

/**
 * LOESS + PAV + monotone PCHIP fitting path, used when the Hermite-basis
 * least-squares solve in {@link FragAlignRegression} fails.
 *
 * <p>Ports {@code fit_spline_from_sorted_arrays}, {@code select_loess_bandwidth},
 * {@code loess_cv_mse} and {@code fix_boundary_plateaus} from FragAlign's
 * {@code src/regression.rs}.
 *
 * <p>FragAlign carries its own LOESS smoother in {@code src/loess.rs} because
 * Rust has no commons-math3; that file is explicitly a port of commons-math3's
 * {@code LoessInterpolator} and {@code SplineInterpolator}, down to the
 * truncating cast on the bandwidth and the 1e-12 accuracy. Java calls the
 * original directly, which is the same numerics from the same source.
 */
public final class LoessFallback {

    /** Bandwidth candidates for cross-validated selection, as a fraction of the data. */
    private static final double[] BANDWIDTH_CANDIDATES = {0.05, 0.10, 0.15, 0.20, 0.30};

    /** Robustness re-weighting passes in each LOESS fit. */
    private static final int ROBUSTNESS_ITERS = 2;

    /** Folds used to cross-validate the bandwidth. */
    private static final int CV_FOLDS = 3;

    /** Uniformly-spaced control points sampled from the smoothed curve. */
    private static final int CONTROL_POINTS = 25;

    /** Below this central slope the data is too pathological to synthesise a straight edge. */
    private static final double MIN_CENTRAL_SLOPE = 1e-3;

    /** A control point within this much of the edge value still counts as plateau. */
    private static final double PLATEAU_THRESHOLD = 0.3;

    private LoessFallback() {}

    /**
     * Fits pre-sorted anchors by cross-validated LOESS, PAV monotonisation and a
     * monotone-preserving interpolant.
     *
     * @return the fitted spline, or {@code null} if the smooth fails
     */
    public static HermiteSpline fit(double[] x, double[] y, double[] weights) {
        int n = x.length;
        if (n < 5) {
            return null;
        }

        double[] w = FragAlignRegression.normalizeWeights(weights, n);

        double[] smoothed;
        try {
            smoothed = new LoessInterpolator(selectBandwidth(x, y), ROBUSTNESS_ITERS).smooth(x, y, w);
        } catch (RuntimeException e) {
            return null;
        }

        double[] monotone = Pava.poolAdjacentViolators(smoothed, null);

        double[] kx = new double[CONTROL_POINTS];
        double[] ky = new double[CONTROL_POINTS];
        double xMin = x[0];
        double xMax = x[n - 1];
        for (int k = 0; k < CONTROL_POINTS; k++) {
            kx[k] = xMin + (double) k / (CONTROL_POINTS - 1) * (xMax - xMin);
            ky[k] = interpolateLinear(x, monotone, kx[k]);
        }

        return fixLeadingPlateau(MonotonePchip.through(kx, ky), x, y);
    }

    /**
     * Picks the bandwidth candidate minimising cross-validated mean squared
     * error, floored so each local window keeps at least three points.
     */
    private static double selectBandwidth(double[] x, double[] y) {
        double minBandwidth = Math.max(3.0 / x.length, 0.01);
        double best = 0.15;
        double bestMse = Double.MAX_VALUE;
        for (double candidate : BANDWIDTH_CANDIDATES) {
            double bandwidth = Math.max(candidate, minBandwidth);
            double mse = crossValidatedMse(x, y, bandwidth);
            if (mse < bestMse) {
                bestMse = mse;
                best = bandwidth;
            }
        }
        return best;
    }

    /**
     * Mean squared error of a LOESS fit at the given bandwidth, estimated by
     * {@value #CV_FOLDS}-fold cross-validation over interleaved anchors.
     *
     * @return the error, or {@link Double#MAX_VALUE} if any fold fails
     */
    private static double crossValidatedMse(double[] x, double[] y, double bandwidth) {
        int n = x.length;
        double totalSquaredError = 0.0;
        int tested = 0;

        for (int fold = 0; fold < CV_FOLDS; fold++) {
            int nTrain = 0;
            for (int i = 0; i < n; i++) {
                if (i % CV_FOLDS != fold) {
                    nTrain++;
                }
            }
            if (nTrain < 4) {
                continue;
            }
            int nTest = n - nTrain;

            double[] xTrain = new double[nTrain];
            double[] yTrain = new double[nTrain];
            double[] xTest = new double[nTest];
            double[] yTest = new double[nTest];
            int trainIdx = 0;
            int testIdx = 0;
            for (int i = 0; i < n; i++) {
                if (i % CV_FOLDS == fold) {
                    xTest[testIdx] = x[i];
                    yTest[testIdx] = y[i];
                    testIdx++;
                } else {
                    xTrain[trainIdx] = x[i];
                    yTrain[trainIdx] = y[i];
                    trainIdx++;
                }
            }
            for (int i = 1; i < nTrain; i++) {
                if (xTrain[i] <= xTrain[i - 1]) {
                    xTrain[i] = xTrain[i - 1] + 1e-8;
                }
            }

            try {
                PolynomialSplineFunction f =
                        new LoessInterpolator(Math.max(bandwidth, 3.0 / nTrain), ROBUSTNESS_ITERS)
                                .interpolate(xTrain, yTrain);
                double trainMin = xTrain[0];
                double trainMax = xTrain[nTrain - 1];
                for (int i = 0; i < nTest; i++) {
                    double xi = Math.min(Math.max(xTest[i], trainMin), trainMax);
                    double diff = f.value(xi) - yTest[i];
                    totalSquaredError += diff * diff;
                    tested++;
                }
            } catch (RuntimeException e) {
                return Double.MAX_VALUE;
            }
        }

        return tested > 0 ? totalSquaredError / tested : Double.MAX_VALUE;
    }

    /**
     * Replaces a PAV-induced plateau at the low-x edge with a straight line at
     * the central-region slope, joined to the interior fit.
     *
     * <p>When anchors at the low end are noisy, PAV collapses them into a flat
     * constant, which would assign an extreme-x peptide a common mid-range
     * value. Trailing plateaus are deliberately left alone — in FragAlign's
     * iRT-to-RT orientation they reflect the physical RT cap at the end of the
     * run.
     */
    private static HermiteSpline fixLeadingPlateau(HermiteSpline spline, double[] x, double[] y) {
        double[] points = spline.points();
        double[] coeff = spline.coefficients();
        int nCtrl = points.length;
        if (nCtrl < 3 || x.length < 8) {
            return spline;
        }

        int n = x.length;
        int from = n / 4;
        int to = 3 * n / 4;
        double centralSlope = regressionSlope(
                Arrays.copyOfRange(x, from, to), Arrays.copyOfRange(y, from, to));
        if (!Double.isFinite(centralSlope) || centralSlope < MIN_CENTRAL_SLOPE) {
            return spline;
        }

        int firstRisingKnot = 0;
        for (int k = 1; k < nCtrl; k++) {
            if (coeff[2 * k] - coeff[0] > PLATEAU_THRESHOLD) {
                firstRisingKnot = k;
                break;
            }
        }
        if (firstRisingKnot == 0) {
            return spline;
        }

        double anchorX = points[firstRisingKnot];
        double anchorY = coeff[2 * firstRisingKnot];
        for (int k = 0; k < firstRisingKnot; k++) {
            coeff[2 * k] = anchorY - centralSlope * (anchorX - points[k]);
            coeff[2 * k + 1] = centralSlope;
        }
        coeff[2 * firstRisingKnot + 1] = centralSlope;
        return spline.withCoefficients(coeff);
    }

    /** Linear interpolation of {@code (x, y)} at {@code xi}, clamped to the end values. */
    private static double interpolateLinear(double[] x, double[] y, double xi) {
        int n = x.length;
        if (xi <= x[0]) {
            return y[0];
        }
        if (xi >= x[n - 1]) {
            return y[n - 1];
        }
        int idx = Arrays.binarySearch(x, xi);
        if (idx >= 0) {
            return y[idx];
        }
        int insertion = -idx - 1;
        double t = (xi - x[insertion - 1]) / (x[insertion] - x[insertion - 1]);
        return y[insertion - 1] + t * (y[insertion] - y[insertion - 1]);
    }

    /** Least-squares slope of {@code ys} on {@code xs}, or NaN if undefined. */
    private static double regressionSlope(double[] xs, double[] ys) {
        int count = xs.length;
        if (count < 2) {
            return Double.NaN;
        }
        double meanX = 0.0;
        double meanY = 0.0;
        for (int i = 0; i < count; i++) {
            meanX += xs[i];
            meanY += ys[i];
        }
        meanX /= count;
        meanY /= count;

        double sxx = 0.0;
        double sxy = 0.0;
        for (int i = 0; i < count; i++) {
            double dx = xs[i] - meanX;
            sxx += dx * dx;
            sxy += dx * (ys[i] - meanY);
        }
        return sxx == 0.0 ? Double.NaN : sxy / sxx;
    }
}
