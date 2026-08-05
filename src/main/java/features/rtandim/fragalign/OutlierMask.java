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

/**
 * Flags "backward" misassignment anchors before fitting.
 *
 * <p>Port of FragAlign's {@code backward_outlier_mask}. An anchor is dropped
 * when its value sits more than {@link #MARGIN_FRACTION} of the value range
 * below the local high quantile of a centred rank window. At the end of the
 * gradient the relation becomes multivalued — the same experimental RT carries
 * predicted values spanning most of the scale — and those low anchors drag a
 * monotone fit down through the cloud. On clean, single-valued data nothing is
 * flagged.
 *
 * <p>The windows are rank-based, so the input must already be ordered by
 * ascending x. The reference takes x as a parameter but reads only its length.
 */
public final class OutlierMask {

    /** Below this many anchors a local quantile is meaningless, so nothing is dropped. */
    private static final int MIN_ANCHORS = 50;

    /** Floor for the centred rank-window half-width. */
    private static final int WINDOW_HALF = 20;

    /** The window half-width is the larger of {@link #WINDOW_HALF} and n / this. */
    private static final int WINDOW_DIVISOR = 80;

    /** Quantile defining the "local high level" an anchor is compared against. */
    private static final double HIGH_QUANTILE = 0.75;

    /** An anchor this far below the local high level, as a fraction of the value range, is dropped. */
    private static final double MARGIN_FRACTION = 0.18;

    private OutlierMask() {}

    /**
     * Returns a per-anchor keep mask.
     *
     * @param y anchor values, ordered by ascending x
     * @return {@code keep[i]} false where anchor i is a backward misassignment;
     *         all true for fewer than {@link #MIN_ANCHORS} anchors or a
     *         degenerate value range
     */
    public static boolean[] keep(double[] y) {
        int n = y.length;
        boolean[] keep = new boolean[n];
        Arrays.fill(keep, true);
        if (n < MIN_ANCHORS) {
            return keep;
        }

        double yMin = Double.POSITIVE_INFINITY;
        double yMax = Double.NEGATIVE_INFINITY;
        for (double v : y) {
            yMin = Math.min(yMin, v);
            yMax = Math.max(yMax, v);
        }
        double margin = MARGIN_FRACTION * (yMax - yMin);
        if (!Double.isFinite(margin) || margin <= 0.0) {
            return keep;
        }

        int half = Math.max(WINDOW_HALF, n / WINDOW_DIVISOR);
        double[] window = new double[2 * half + 1];
        for (int i = 0; i < n; i++) {
            int lo = Math.max(0, i - half);
            int hi = Math.min(n, i + half + 1);
            int len = hi - lo;
            System.arraycopy(y, lo, window, 0, len);
            Arrays.sort(window, 0, len);
            int quantileIdx = (int) ((len - 1) * HIGH_QUANTILE);
            if (y[i] < window[quantileIdx] - margin) {
                keep[i] = false;
            }
        }
        return keep;
    }
}
