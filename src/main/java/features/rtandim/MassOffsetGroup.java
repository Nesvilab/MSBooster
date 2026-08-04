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

package features.rtandim;

import java.util.ArrayList;
import java.util.Arrays;

/**
 * One RT/IM calibration group: the set of mass shifts whose PSMs share a calibration curve.
 *
 * A group key comes from {@code Constants.massesForLoessCalibration}. Several masses may share one
 * curve, written joined by "&amp;" (this is how all offsets of a single MSFragger mass_offsets line
 * end up on one curve). Two keys are special:
 * <ul>
 *   <li>{@code ""} — no grouping was requested, so every peptide belongs to the one curve</li>
 *   <li>{@code "others"} — the catch-all for peptides in no other group. It never matches directly;
 *       callers fall back to it when nothing else matched</li>
 * </ul>
 *
 * Membership is decided on the delta masses themselves rather than on the text of the peptide name.
 * MSFragger writes the mass from fragger.params and the localized delta mass independently, each
 * rounded to 4 decimals, and the two can disagree in the last digit — a pin holds both
 * {@code [79.9663]} and {@code [79.9664]} for the same offset — so a string match drops PSMs.
 * Comparing numerically also stops a mass from matching inside a longer one, e.g. {@code 14.0157}
 * inside {@code [-14.0157]}.
 *
 * Instances are immutable and safe to share across threads.
 */
public final class MassOffsetGroup {
    /**
     * How far a peptide's delta mass may sit from the group's mass and still be considered the same
     * offset. Comfortably above the 1e-4 rounding disagreement described above, and well below the
     * smallest gap between distinct offsets in a realistic list (~8e-3 Da in a 85-offset search).
     */
    public static final double TOLERANCE = 1e-3;

    public static final String OTHERS = "others";

    private static final double[] NO_MASSES = new double[0];

    private final String key;
    private final String[] tokens;
    private final double[] masses; //NaN where the matching token is not a number

    private MassOffsetGroup(String key, String[] tokens, double[] masses) {
        this.key = key;
        this.tokens = tokens;
        this.masses = masses;
    }

    /** Builds the group for one key from {@code Constants.massesForLoessCalibration}. */
    public static MassOffsetGroup of(String key) {
        String[] tokens = key.split("&");
        double[] masses = new double[tokens.length];
        for (int i = 0; i < tokens.length; i++) {
            masses[i] = parseOrNaN(tokens[i]);
        }
        return new MassOffsetGroup(key, tokens, masses);
    }

    private static double parseOrNaN(String s) {
        try {
            return Double.parseDouble(s.trim());
        } catch (NumberFormatException e) {
            return Double.NaN;
        }
    }

    public String getKey() {
        return key;
    }

    /**
     * The delta masses a peptide carries, read from the "[]" annotations of a base-format peptide
     * name such as {@code "AC[57.0215]DEK[-14.0157]|2"}. Bracketed text that is not a number is
     * skipped, since it cannot be a mass offset.
     */
    public static double[] deltaMasses(String peptideName) {
        if (peptideName == null || peptideName.indexOf('[') < 0) {
            return NO_MASSES;
        }
        ArrayList<Double> found = new ArrayList<>(4);
        int from = 0;
        while (true) {
            int open = peptideName.indexOf('[', from);
            if (open < 0) {
                break;
            }
            int close = peptideName.indexOf(']', open + 1);
            if (close < 0) {
                break;
            }
            double mass = parseOrNaN(peptideName.substring(open + 1, close));
            if (!Double.isNaN(mass)) {
                found.add(mass);
            }
            from = close + 1;
        }
        double[] result = new double[found.size()];
        for (int i = 0; i < result.length; i++) {
            result[i] = found.get(i);
        }
        return result;
    }

    /**
     * Whether a peptide belongs to this group. Pass the peptide's delta masses from
     * {@link #deltaMasses}, parsed once and reused across all groups.
     *
     * A token that is not a number falls back to the historical substring match on the peptide name,
     * so a hand-written non-numeric {@code massesForLoessCalibration} keeps working. That fallback is
     * also what gives the two special keys their behaviour: {@code ""} is contained in every name, and
     * {@code "others"} is contained in none.
     */
    public boolean matches(String peptideName, double[] peptideDeltaMasses) {
        for (int i = 0; i < tokens.length; i++) {
            if (Double.isNaN(masses[i])) {
                if (peptideName.contains(tokens[i])) {
                    return true;
                }
            } else {
                for (double delta : peptideDeltaMasses) {
                    if (Math.abs(delta - masses[i]) <= TOLERANCE) {
                        return true;
                    }
                }
            }
        }
        return false;
    }

    /** Convenience for callers that do not already hold the parsed masses. */
    public boolean matches(String peptideName) {
        return matches(peptideName, deltaMasses(peptideName));
    }

    @Override
    public String toString() {
        return "MassOffsetGroup(" + key + " -> " + Arrays.toString(masses) + ")";
    }
}
