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

import static org.junit.jupiter.api.Assertions.assertFalse;
import static org.junit.jupiter.api.Assertions.assertTrue;

import org.junit.jupiter.api.Test;

/**
 * Covers FragAlign's {@code backward_outlier_mask} contract: drop anchors whose
 * value sits far below the local high level, because at the gradient end the
 * relation becomes multivalued and those points drag a monotone fit down through
 * the cloud.
 */
public class OutlierMaskTest {

    /** A clean ramp of {@code n} points, ordered by ascending x. */
    private static double[] ramp(int n) {
        double[] y = new double[n];
        for (int i = 0; i < n; i++) {
            y[i] = i;
        }
        return y;
    }

    private static void assertAllKept(boolean[] keep) {
        for (int i = 0; i < keep.length; i++) {
            assertTrue(keep[i], "index " + i + " should have been kept");
        }
    }

    @Test
    public void keepsEverythingOnCleanSingleValuedData() {
        assertAllKept(OutlierMask.keep(ramp(200)));
    }

    @Test
    public void isANoOpBelowTheMinimumAnchorCount() {
        // 49 anchors is under the threshold, so even a gross outlier survives —
        // a local quantile is meaningless on that few points.
        double[] y = ramp(49);
        y[40] = -1000.0;
        assertAllKept(OutlierMask.keep(y));
    }

    @Test
    public void isANoOpWhenEveryValueIsIdentical() {
        double[] y = new double[100];
        java.util.Arrays.fill(y, 7.0);
        assertAllKept(OutlierMask.keep(y));
    }

    @Test
    public void dropsTheGradientEndCloudButKeepsTheRamp() {
        // Reproduces the shape of the MSBooster failure: a clean rising ramp
        // with a cluster of far-below anchors appended at the high-x end.
        int rampSize = 200;
        int cloudSize = 10;
        double[] y = new double[rampSize + cloudSize];
        System.arraycopy(ramp(rampSize), 0, y, 0, rampSize);
        for (int i = 0; i < cloudSize; i++) {
            y[rampSize + i] = 20.0;
        }

        boolean[] keep = OutlierMask.keep(y);

        for (int i = 0; i < rampSize; i++) {
            assertTrue(keep[i], "ramp anchor " + i + " should have been kept");
        }
        for (int i = rampSize; i < y.length; i++) {
            assertFalse(keep[i], "cloud anchor " + i + " should have been dropped");
        }
    }

    @Test
    public void keepsAModerateDipButDropsADeepOne() {
        // The anchor is compared against the 75th percentile of its window, not
        // against its own trend value. In a 200-point ramp that quantile sits
        // roughly 10 above the anchor, so the effective tolerance is a little
        // under 0.18 of the range rather than exactly it.
        double range = 199.0;

        double[] shallow = ramp(200);
        shallow[150] -= 0.10 * range;
        assertTrue(OutlierMask.keep(shallow)[150], "a 10% dip is noise and must survive");

        double[] deep = ramp(200);
        deep[150] -= 0.30 * range;
        assertFalse(OutlierMask.keep(deep)[150], "a 30% dip is a misassignment and must be dropped");
    }
}
