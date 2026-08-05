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
 * Pool Adjacent Violators Algorithm (PAVA).
 *
 * <p>Port of FragAlign's {@code src/pava.rs}. The unweighted and weighted
 * variants collapse into a single method: passing {@code null} weights
 * reproduces the unweighted behaviour exactly.
 */
public final class Pava {

    private Pava() {}

    /**
     * Returns a non-decreasing version of {@code y} minimising the weighted sum
     * of squared deviations from the input.
     *
     * <p>Merging two adjacent blocks uses the weighted mean of the block values,
     * each block's accumulated weight being the sum of its members'. This stops
     * low-quality anchors from dragging the monotone-enforced curve down through
     * unweighted averaging with the high-quality signal.
     *
     * <p>{@code y} is not modified; a fresh array is returned.
     *
     * @param y       values to monotonise
     * @param weights per-value weights, or {@code null} for unweighted pooling
     * @throws IllegalArgumentException if {@code weights} is non-null and its
     *                                  length differs from {@code y}'s
     */
    public static double[] poolAdjacentViolators(double[] y, double[] weights) {
        int n = y.length;
        if (weights != null && weights.length != n) {
            throw new IllegalArgumentException(
                    "weights.length (" + weights.length + ") != y.length (" + n + ")");
        }

        double[] result = new double[n];
        if (n == 0) {
            return result;
        }

        double[] w0 = initialWeights(weights, n);

        // One block per input element, pooled in place on a stack. value/weight
        // hold the running weighted mean and accumulated weight of each block;
        // runLength counts the elements pooled into it.
        double[] value = new double[n];
        double[] weight = new double[n];
        int[] runLength = new int[n];
        int top = 0;

        for (int i = 0; i < n; i++) {
            value[top] = y[i];
            weight[top] = w0[i];
            runLength[top] = 1;
            top++;

            // Merge while the newest block violates monotonicity against its
            // predecessor. Merging can re-violate further left, so this
            // cascades — the backtracking of the reference implementation.
            while (top > 1 && value[top - 1] < value[top - 2]) {
                double v1 = value[top - 2];
                double w1 = weight[top - 2];
                double v2 = value[top - 1];
                double w2 = weight[top - 1];
                top--;
                value[top - 1] = (v1 * w1 + v2 * w2) / (w1 + w2);
                weight[top - 1] = w1 + w2;
                runLength[top - 1] += runLength[top];
            }
        }

        int pos = 0;
        for (int k = 0; k < top; k++) {
            for (int j = 0; j < runLength[k]; j++) {
                result[pos] = value[k];
                pos++;
            }
        }
        return result;
    }

    /**
     * Initial per-element weights.
     *
     * <p>{@code null} weights, or weights whose non-negative part sums to zero,
     * fall back to all ones; otherwise each entry is floored at {@code 1e-12} so
     * a merge can never divide by zero.
     */
    private static double[] initialWeights(double[] weights, int n) {
        double[] w0 = new double[n];
        if (weights == null) {
            java.util.Arrays.fill(w0, 1.0);
            return w0;
        }

        double total = 0.0;
        for (double v : weights) {
            total += Math.max(0.0, v);
        }
        if (total <= 0.0) {
            java.util.Arrays.fill(w0, 1.0);
            return w0;
        }

        for (int i = 0; i < n; i++) {
            w0[i] = Math.max(1e-12, weights[i]);
        }
        return w0;
    }
}
