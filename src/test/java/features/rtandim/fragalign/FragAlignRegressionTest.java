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
import static org.junit.jupiter.api.Assertions.assertNotNull;
import static org.junit.jupiter.api.Assertions.assertNull;
import static org.junit.jupiter.api.Assertions.assertThrows;
import static org.junit.jupiter.api.Assertions.assertTrue;

import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.io.BufferedReader;
import java.nio.charset.StandardCharsets;
import java.util.ArrayList;
import java.util.List;

import org.junit.jupiter.api.Test;
import org.junit.jupiter.params.ParameterizedTest;
import org.junit.jupiter.params.provider.ValueSource;

/**
 * Holds the Java port to FragAlign's committed reference output.
 *
 * <p>The four {@code *.tsv} anchor sets and their {@code *_fit.tsv} fitted
 * curves are copied verbatim from FragAlign's {@code tests/data}. The fitted
 * files were produced by the Rust {@code dump_fits} binary: 200 samples of
 * {@code spline.evaluate} spanning the anchor range. Reproducing them is a
 * direct fidelity check of the port across linear, nonlinear, outlier-laden and
 * sparse-edge inputs, including the weighted path.
 *
 * <p>The remaining tests are transcribed from {@code tests/rt_regression.rs}.
 */
public class FragAlignRegressionTest {

    /** Columns of one anchor TSV. {@code weights} is null when the file has no weight column. */
    private static final class Anchors {
        final double[] x;
        final double[] y;
        final double[] weights;

        Anchors(double[] x, double[] y, double[] weights) {
            this.x = x;
            this.y = y;
            this.weights = weights;
        }
    }

    private static List<String[]> readRows(String resource) throws IOException {
        List<String[]> rows = new ArrayList<>();
        try (InputStream in = FragAlignRegressionTest.class.getClassLoader().getResourceAsStream(resource)) {
            assertNotNull(in, "missing test resource: " + resource);
            BufferedReader reader = new BufferedReader(new InputStreamReader(in, StandardCharsets.UTF_8));
            String line;
            boolean headerSeen = false;
            while ((line = reader.readLine()) != null) {
                line = line.trim();
                if (line.isEmpty() || line.startsWith("#")) {
                    continue;
                }
                if (!headerSeen) {
                    headerSeen = true;
                    continue;
                }
                rows.add(line.split("\t"));
            }
        }
        return rows;
    }

    private static Anchors readAnchors(String scenario) throws IOException {
        List<String[]> rows = readRows("fragalign/" + scenario + ".tsv");
        int n = rows.size();
        double[] x = new double[n];
        double[] y = new double[n];
        boolean weighted = rows.get(0).length >= 3;
        double[] weights = weighted ? new double[n] : null;
        for (int i = 0; i < n; i++) {
            x[i] = Double.parseDouble(rows.get(i)[0]);
            y[i] = Double.parseDouble(rows.get(i)[1]);
            if (weighted) {
                weights[i] = Double.parseDouble(rows.get(i)[2]);
            }
        }
        return new Anchors(x, y, weights);
    }

    @ParameterizedTest
    @ValueSource(strings = {"linear_gradient", "nonlinear_gradient", "with_outliers", "sparse_edges"})
    public void reproducesTheRustReferenceCurve(String scenario) throws IOException {
        Anchors anchors = readAnchors(scenario);
        HermiteSpline spline = FragAlignRegression.fit(anchors.x, anchors.y, anchors.weights);
        assertNotNull(spline, "fit returned null for " + scenario);

        List<String[]> expected = readRows("fragalign/" + scenario + "_fit.tsv");
        assertEquals(200, expected.size(), "unexpected sample count for " + scenario);

        for (String[] row : expected) {
            double x = Double.parseDouble(row[0]);
            double want = Double.parseDouble(row[1]);
            double got = spline.evaluate(x);
            // The reference files carry six decimals, so the floor on agreement
            // is file rounding rather than the solve.
            double tolerance = Math.max(1e-5, 1e-6 * Math.abs(want));
            assertEquals(want, got, tolerance, scenario + " diverges at x = " + x);
        }
    }

    @Test
    public void reconstructsALinearGradient() {
        int n = 200;
        double[] x = new double[n];
        double[] y = new double[n];
        for (int i = 0; i < n; i++) {
            x[i] = i;
            y[i] = 0.05 * x[i] + 5.0;
        }

        HermiteSpline spline = FragAlignRegression.fit(x, y, null);
        assertNotNull(spline);

        double sumSq = 0.0;
        for (int i = 0; i < n; i++) {
            double residual = y[i] - spline.evaluate(x[i]);
            sumSq += residual * residual;
        }
        assertTrue(Math.sqrt(sumSq / n) < 0.5, "RMSE too high on a linear gradient");
    }

    @Test
    public void reconstructsANonlinearGradient() {
        int n = 500;
        double[] x = new double[n];
        double[] y = new double[n];
        for (int i = 0; i < n; i++) {
            double t = i / (double) (n - 1);
            x[i] = -50.0 + 150.0 * t;
            y[i] = 5.0 + 60.0 * (t - 0.5 * t * t + 0.2 * t * t * t);
        }

        HermiteSpline spline = FragAlignRegression.fit(x, y, null);
        assertNotNull(spline);

        double sumSq = 0.0;
        for (int i = 0; i < n; i++) {
            double residual = y[i] - spline.evaluate(x[i]);
            sumSq += residual * residual;
        }
        assertTrue(Math.sqrt(sumSq / n) < 0.3, "RMSE too high on a nonlinear gradient");
    }

    @Test
    public void sortsUnsortedInputInternally() {
        double[] x = {50.0, 10.0, 30.0, 70.0, 90.0, 20.0, 40.0, 60.0, 80.0, 0.0};
        double[] y = new double[x.length];
        for (int i = 0; i < x.length; i++) {
            y[i] = 0.1 * x[i] + 2.0;
        }

        HermiteSpline spline = FragAlignRegression.fit(x, y, null);
        assertNotNull(spline);
        assertEquals(7.0, spline.evaluate(50.0), 0.5);
    }

    @Test
    public void fittedDomainReachesTheEndsOfTheAnchorRange() {
        // Interval midpoints alone leave the outermost knots well inside the
        // data — on this 0..299 range the first knot landed at 21 — which
        // truncates the curve, because invert() clamps to the fitted domain
        // even though evaluate() extrapolates. On the real CCRCC run that put
        // the ceiling at predRT 146.2 while the top anchor was 169.8.
        //
        // The pins sit at the 99.9th percentile rather than the literal
        // extremes, so the domain reaches to within one anchor of each end
        // without handing the boundary to a lone outlier.
        int n = 300;
        double[] x = new double[n];
        double[] y = new double[n];
        for (int i = 0; i < n; i++) {
            x[i] = i;
            y[i] = 0.05 * i + 5.0;
        }

        HermiteSpline spline = FragAlignRegression.fit(x, y, null);
        assertNotNull(spline);
        assertTrue(spline.xMin() <= x[1],
                "fit starts at " + spline.xMin() + ", short of the first anchors");
        assertTrue(spline.xMax() >= x[n - 2],
                "fit ends at " + spline.xMax() + ", short of the last anchors");
    }

    @Test
    public void outermostKnotsIgnoreALoneExtremeAnchor() {
        // One anchor far beyond the rest must not become the clamp ceiling.
        int n = 300;
        double[] x = new double[n];
        double[] y = new double[n];
        for (int i = 0; i < n - 1; i++) {
            x[i] = i;
            y[i] = 0.05 * i + 5.0;
        }
        x[n - 1] = 5000.0;                 // isolated outlier, far past the body
        y[n - 1] = 0.05 * 299.0 + 5.0;

        HermiteSpline spline = FragAlignRegression.fit(x, y, null);
        assertNotNull(spline);
        assertTrue(spline.xMax() < 1000.0,
                "ceiling " + spline.xMax() + " was handed to the lone outlier");
        assertTrue(spline.xMax() >= 290.0,
                "ceiling " + spline.xMax() + " fell short of the anchor body");
    }

    @Test
    public void returnsNullBelowFiveAnchors() {
        assertNull(FragAlignRegression.fit(
                new double[]{1.0, 2.0, 3.0, 4.0},
                new double[]{5.0, 6.0, 7.0, 8.0},
                null));
    }

    @Test
    public void rejectsMismatchedInputLengths() {
        assertThrows(IllegalArgumentException.class, () -> FragAlignRegression.fit(
                new double[]{1.0, 2.0, 3.0, 4.0, 5.0},
                new double[]{1.0, 2.0},
                null));
    }

    @Test
    public void rejectsMismatchedWeightLength() {
        assertThrows(IllegalArgumentException.class, () -> FragAlignRegression.fit(
                new double[]{1.0, 2.0, 3.0, 4.0, 5.0},
                new double[]{1.0, 2.0, 3.0, 4.0, 5.0},
                new double[]{1.0, 1.0}));
    }

    @Test
    public void fittedCurveIsAlwaysNonDecreasing() {
        // A rising ramp whose tail collapses into a multivalued cloud — the
        // shape that flattens the top right of MSBooster's calibration curve.
        int n = 300;
        double[] x = new double[n];
        double[] y = new double[n];
        for (int i = 0; i < n; i++) {
            x[i] = i;
            y[i] = i < 260 ? i : (i % 3 == 0 ? 40.0 : 265.0);
        }

        HermiteSpline spline = FragAlignRegression.fit(x, y, null);
        assertNotNull(spline);

        double previous = Double.NEGATIVE_INFINITY;
        for (int k = 0; k <= 600; k++) {
            double xi = x[0] + (x[n - 1] - x[0]) * k / 600.0;
            double v = spline.evaluate(xi);
            assertTrue(v >= previous - 1e-9, "curve descends at x = " + xi);
            previous = v;
        }
    }

    @Test
    public void theTopOfTheCurveDoesNotGoFlat() {
        // The regression this whole change exists for: with a noisy multivalued
        // tail the fitted curve must keep rising, not collapse to a constant.
        int n = 400;
        double[] x = new double[n];
        double[] y = new double[n];
        for (int i = 0; i < n; i++) {
            x[i] = i;
            y[i] = i;
        }
        // Gradient-end cloud: every third anchor past 360 is a gross misassignment.
        for (int i = 360; i < n; i += 3) {
            y[i] = 50.0;
        }

        HermiteSpline spline = FragAlignRegression.fit(x, y, null);
        assertNotNull(spline);

        double terminalSlope = (spline.evaluate(399.0) - spline.evaluate(359.0)) / 40.0;
        double middleSlope = (spline.evaluate(300.0) - spline.evaluate(100.0)) / 200.0;
        assertTrue(terminalSlope > 0.25 * middleSlope,
                "terminal slope " + terminalSlope + " collapsed against middle slope " + middleSlope);
    }
}
