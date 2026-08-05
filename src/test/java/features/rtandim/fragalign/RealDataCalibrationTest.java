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

import static org.junit.jupiter.api.Assertions.assertNotNull;
import static org.junit.jupiter.api.Assertions.assertTrue;

import java.io.BufferedReader;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.nio.charset.StandardCharsets;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

import org.junit.jupiter.api.BeforeAll;
import org.junit.jupiter.api.Test;

/**
 * Regression test against the dataset that motivated replacing the old
 * LOESS + PAV fitter: 5000 RT calibration anchors from a CPTAC CCRCC LUMOS DDA
 * run, dumped from {@code MzmlReader.expAndPredRTs}.
 *
 * <p>The old fitter collapsed the top right of this curve into a flat line: the
 * final Pool Adjacent Violators block absorbed the whole sparse tail and every
 * point in it received the block mean, so the curve stopped at predRT 134.2 even
 * though anchors around expRT 98-100 have a median predRT of 140.9 rising to
 * 144.9. The tests below pin the two properties that failure violated.
 *
 * <p>The fit is performed the way {@code MzmlReader.setLOESS} performs it:
 * predicted RT as the independent variable, experimental RT as the dependent
 * one, with {@link HermiteSpline#invert} supplying the experimental-to-predicted
 * direction the RT features consume.
 */
public class RealDataCalibrationTest {

    private static final String RESOURCE = "fragalign/ccrcc_lumos_f01_rt_anchors.tsv";

    /** Highest predicted RT the superseded fitter reached on this dataset. */
    private static final double OLD_FITTER_CEILING = 134.21;

    private static double[] expRT;
    private static double[] predRT;
    private static HermiteSpline spline;

    @BeforeAll
    public static void fitOnce() throws IOException {
        List<double[]> rows = new ArrayList<>();
        try (InputStream in = RealDataCalibrationTest.class.getClassLoader().getResourceAsStream(RESOURCE)) {
            assertNotNull(in, "missing test resource: " + RESOURCE);
            BufferedReader reader = new BufferedReader(new InputStreamReader(in, StandardCharsets.UTF_8));
            String line = reader.readLine(); // header
            while ((line = reader.readLine()) != null) {
                if (line.trim().isEmpty()) {
                    continue;
                }
                String[] parts = line.split("\t");
                rows.add(new double[]{Double.parseDouble(parts[0]), Double.parseDouble(parts[1])});
            }
        }

        expRT = new double[rows.size()];
        predRT = new double[rows.size()];
        for (int i = 0; i < rows.size(); i++) {
            expRT[i] = rows.get(i)[0];
            predRT[i] = rows.get(i)[1];
        }

        spline = FragAlignRegression.fit(predRT, expRT, null);
        assertNotNull(spline, "fit returned null on real calibration anchors");
    }

    @Test
    public void curveReachesTheTopOfTheDenseAnchorRidge() {
        double maxFitted = Double.NEGATIVE_INFINITY;
        double lo = Arrays.stream(expRT).min().orElseThrow(IllegalStateException::new);
        double hi = Arrays.stream(expRT).max().orElseThrow(IllegalStateException::new);
        for (int k = 0; k <= 1000; k++) {
            maxFitted = Math.max(maxFitted, spline.invert(lo + (hi - lo) * k / 1000.0));
        }

        // Anchors at expRT 98-100 have a median predRT of 140.9. A fitter that
        // truncates the ridge, as the old one did at 134.2, fails here.
        assertTrue(maxFitted > OLD_FITTER_CEILING + 5.0,
                "curve tops out at " + maxFitted + ", short of the anchor ridge");
        assertTrue(maxFitted >= 140.0,
                "curve tops out at " + maxFitted + ", below the local anchor median of 140.9");
    }

    @Test
    public void curveIsNonDecreasingAcrossTheRun() {
        double lo = Arrays.stream(expRT).min().orElseThrow(IllegalStateException::new);
        double hi = Arrays.stream(expRT).max().orElseThrow(IllegalStateException::new);
        double previous = Double.NEGATIVE_INFINITY;
        for (int k = 0; k <= 2000; k++) {
            double x = lo + (hi - lo) * k / 2000.0;
            double v = spline.invert(x);
            assertTrue(v >= previous - 1e-6, "calibration curve descends at expRT " + x);
            previous = v;
        }
    }

    @Test
    public void curveTracksTheAnchorsThroughTheBodyOfTheGradient() {
        // Median absolute deviation between the fit and its anchors, over the
        // well-populated middle of the run.
        List<Double> deviations = new ArrayList<>();
        for (int i = 0; i < expRT.length; i++) {
            if (expRT[i] >= 20.0 && expRT[i] <= 95.0) {
                deviations.add(Math.abs(spline.invert(expRT[i]) - predRT[i]));
            }
        }
        assertTrue(deviations.size() > 3000, "unexpectedly few anchors in the gradient body");
        double[] sorted = deviations.stream().mapToDouble(Double::doubleValue).sorted().toArray();
        double median = sorted[sorted.length / 2];
        assertTrue(median < 5.0, "median deviation from anchors is " + median);
    }
}
