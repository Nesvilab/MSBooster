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

/**
 * Fritsch-Carlson monotone-preserving cubic Hermite interpolation, and the
 * monotonicity guard applied to every fitted spline.
 *
 * <p>Ports {@code monotone_pchip} and {@code enforce_monotonic} from FragAlign's
 * {@code src/regression.rs}.
 */
public final class MonotonePchip {

    /** Samples used to test a fitted spline for monotonicity. */
    private static final int MONOTONICITY_SAMPLES = 512;

    private MonotonePchip() {}

    /**
     * Builds a monotone-preserving cubic Hermite spline through the given knots.
     *
     * @param kx knot positions, sorted ascending
     * @param ky knot values, non-decreasing
     */
    public static HermiteSpline through(double[] kx, double[] ky) {
        int nCtrl = kx.length;
        int m = nCtrl - 1;

        double[] delta = new double[m];
        for (int k = 0; k < m; k++) {
            double dx = kx[k + 1] - kx[k];
            delta[k] = dx > 1e-12 ? (ky[k + 1] - ky[k]) / dx : 0.0;
        }

        double[] tangent = new double[nCtrl];
        tangent[0] = delta[0];
        tangent[m] = delta[m - 1];
        for (int k = 1; k < m; k++) {
            tangent[k] = 0.5 * (delta[k - 1] + delta[k]);
        }

        // Fritsch-Carlson limiter: clamp tangents into the circle of radius 3
        // so the interpolant cannot overshoot and break monotonicity.
        for (int k = 0; k < m; k++) {
            if (Math.abs(delta[k]) < 1e-12) {
                tangent[k] = 0.0;
                tangent[k + 1] = 0.0;
                continue;
            }
            double alpha = tangent[k] / delta[k];
            double beta = tangent[k + 1] / delta[k];
            double s = alpha * alpha + beta * beta;
            if (s > 9.0) {
                double tau = 3.0 / Math.sqrt(s);
                tangent[k] = tau * alpha * delta[k];
                tangent[k + 1] = tau * beta * delta[k];
            }
        }

        double[] coeff = new double[2 * nCtrl];
        for (int k = 0; k < nCtrl; k++) {
            coeff[2 * k] = ky[k];
            coeff[2 * k + 1] = tangent[k];
        }
        return new HermiteSpline(kx, coeff);
    }

    /**
     * Returns {@code spline} unchanged when it is already non-decreasing across
     * {@code [x0, x1]}; otherwise rebuilds it as a monotone interpolant through
     * the running maximum of its densely-sampled values.
     *
     * <p>A calibration curve must be non-decreasing, but the Hermite-basis solve
     * can overshoot and then dip across a sparse, steeply-rising tail. The
     * running maximum removes that spurious fall-back while preserving
     * everything up to it. Well-behaved fits are returned untouched.
     */
    public static HermiteSpline enforce(HermiteSpline spline, double x0, double x1) {
        if (x1 <= x0) {
            return spline;
        }

        int n = MONOTONICITY_SAMPLES;
        double[] xs = new double[n];
        double[] ys = new double[n];
        for (int k = 0; k < n; k++) {
            xs[k] = x0 + (x1 - x0) * k / (double) (n - 1);
            ys[k] = spline.evaluate(xs[k]);
        }

        boolean descends = false;
        for (int k = 1; k < n; k++) {
            if (ys[k] < ys[k - 1] - 1e-6) {
                descends = true;
                break;
            }
        }
        if (!descends) {
            return spline;
        }

        for (int k = 1; k < n; k++) {
            if (ys[k] < ys[k - 1]) {
                ys[k] = ys[k - 1];
            }
        }
        return through(xs, ys);
    }
}
