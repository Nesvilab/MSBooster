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

import java.util.function.DoubleUnaryOperator;

/**
 * Piecewise cubic Hermite spline, the fitted calibration curve.
 *
 * <p>Port of FragAlign's {@code src/spline.rs}. FragAlign names the axes iRT and
 * RT; MSBooster fits in both directions (experimental RT to predicted RT in
 * {@code MzmlReader.setLOESS}, the reverse in the best-model search), so the
 * axes are named neutrally here: {@link #evaluate} maps x to y and
 * {@link #invert} maps y back to x.
 *
 * <p>{@code points} holds the strictly-increasing control-point positions
 * (length m + 1). {@code coeff} holds, for each control point, the alternating
 * value and derivative (length 2(m + 1)): {@code coeff[2k]} is the value at knot
 * k, {@code coeff[2k + 1]} the slope there. Outside the fitted range the spline
 * extrapolates linearly from the edge knot's value and slope rather than going
 * flat.
 *
 * <p>Instances are immutable and safe to share across threads.
 */
public final class HermiteSpline implements DoubleUnaryOperator {

    private final double[] points;
    private final double[] coeff;

    /**
     * @param points strictly-increasing knot positions, length m + 1
     * @param coeff  alternating value/derivative coefficients, length 2(m + 1)
     */
    public HermiteSpline(double[] points, double[] coeff) {
        if (points.length == 0) {
            throw new IllegalArgumentException("points must not be empty");
        }
        if (coeff.length != 2 * points.length) {
            throw new IllegalArgumentException(
                    "coeff.length (" + coeff.length + ") != 2 * points.length (" + 2 * points.length + ")");
        }
        this.points = points.clone();
        this.coeff = coeff.clone();
    }

    /** Knot positions. Returns a copy; the spline itself is immutable. */
    public double[] points() {
        return points.clone();
    }

    /** Value/derivative coefficients. Returns a copy; the spline is immutable. */
    public double[] coefficients() {
        return coeff.clone();
    }

    /**
     * Returns a spline with the same knots and the supplied coefficients. Used
     * by the boundary-plateau fix, which the reference performs as an in-place
     * mutation of {@code coeff}.
     */
    public HermiteSpline withCoefficients(double[] newCoeff) {
        return new HermiteSpline(points, newCoeff);
    }

    /**
     * Evaluates the spline at {@code x}, extrapolating linearly outside the
     * fitted range.
     */
    public double evaluate(double x) {
        int m = points.length - 1;

        if (x <= points[0]) {
            return coeff[0] + coeff[1] * (x - points[0]);
        }
        if (x >= points[m]) {
            return coeff[2 * m] + coeff[2 * m + 1] * (x - points[m]);
        }

        int k = m - 1;
        for (int s = 0; s < m; s++) {
            if (x <= points[s + 1]) {
                k = s;
                break;
            }
        }

        double w = points[k + 1] - points[k];
        if (w < 1e-12) {
            return coeff[2 * k];
        }

        double u = (x - points[k]) / w;
        double v = u - 1.0;

        return (1.0 + 2.0 * u) * v * v * coeff[2 * k]
                + w * u * v * v * coeff[2 * k + 1]
                + u * u * (1.0 - 2.0 * v) * coeff[2 * (k + 1)]
                + w * u * u * v * coeff[2 * (k + 1) + 1];
    }

    /** Batch {@link #evaluate}. {@code xs} need not be sorted. */
    public double[] evaluateAll(double[] xs) {
        double[] out = new double[xs.length];
        for (int i = 0; i < xs.length; i++) {
            out[i] = evaluate(xs[i]);
        }
        return out;
    }

    /**
     * Inverse of {@link #evaluate}: the x whose fitted value is {@code y}, found
     * by bisection on the monotone-increasing spline.
     *
     * <p>The search domain is clamped to the fitted range, so a {@code y} below
     * the left edge returns {@link #xMin} and one above the right edge returns
     * {@link #xMax}.
     */
    public double invert(double y) {
        double lo = points[0];
        double hi = points[points.length - 1];
        if (y <= evaluate(lo)) {
            return lo;
        }
        if (y >= evaluate(hi)) {
            return hi;
        }
        for (int i = 0; i < 60; i++) {
            double mid = 0.5 * (lo + hi);
            double yMid = evaluate(mid);
            if (Math.abs(yMid - y) < 1e-4) {
                return mid;
            }
            if (yMid < y) {
                lo = mid;
            } else {
                hi = mid;
            }
        }
        return 0.5 * (lo + hi);
    }

    /** Batch {@link #invert}. */
    public double[] invertAll(double[] ys) {
        double[] out = new double[ys.length];
        for (int i = 0; i < ys.length; i++) {
            out[i] = invert(ys[i]);
        }
        return out;
    }

    /** Lower bound of the fitted range. */
    public double xMin() {
        return points[0];
    }

    /** Upper bound of the fitted range. */
    public double xMax() {
        return points[points.length - 1];
    }

    /** Lets a fitted spline be used directly wherever a calibration curve is expected. */
    @Override
    public double applyAsDouble(double x) {
        return evaluate(x);
    }
}
