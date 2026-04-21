/*
 * This file is part of MSBooster.
 *
 * MSBooster is free software: you can redistribute it and/or modify it under
 * the terms of the GNU Lesser General Public License as published by the Free
 * Software Foundation, either version 3 of the License, or (at your
 * option) any later version.
 *
 * MSBooster is distributed in the hope that it will be useful, but WITHOUT
 * ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
 * FITNESS FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public License
 * for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with MSBooster. If not, see <https://www.gnu.org/licenses/>.
 */

package predictions;

import allconstants.Constants;
import allconstants.FragmentIonConstants;

import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;

public class PredictionEntry {
    // All 5 per-fragment fields are packed into a single int[] to eliminate 4 array headers
    // and reduce per-element cost from 20 bytes to 12 bytes per fragment (~13 GB savings at 22.6M precursors).
    //
    // Layout — STRIDE = 3 ints per fragment at index i:
    //   fragments[STRIDE*i]     = Float.floatToRawIntBits(mz)
    //   fragments[STRIDE*i + 1] = Float.floatToRawIntBits(intensity)
    //   fragments[STRIDE*i + 2] = (fragNum & 0xFF) << 16
    //                           | (ionTypeByte & 0xFF) << 8
    //                           | (charge & 0xFF)
    // ionTypeByte is an index into FragmentIonConstants.ION_NAMES.
    public static final int STRIDE = 3;

    public int[] fragments = new int[0];
    public String[] fullAnnotations = new String[0];
    public int[] isotopes = new int[0];
    public float RT;
    public float IM;
    public PredictionEntry auxSpectra = null;
    public HashMap<String, Float[]> scores = null;
    public double precursorMz = 0d;
    private static final float maxIntensity = 1f;
    public boolean daltonMatching = false;

    // ---- Scalar accessors ----
    public int numFragments() { return fragments.length / STRIDE; }

    public float getMz(int i)           { return Float.intBitsToFloat(fragments[STRIDE * i]); }
    public float getIntensity(int i)    { return Float.intBitsToFloat(fragments[STRIDE * i + 1]); }
    public int   getFragNum(int i)      { return (fragments[STRIDE * i + 2] >> 16) & 0xFF; }
    public int   getIonTypeByte(int i)  { return (fragments[STRIDE * i + 2] >>  8) & 0xFF; }
    public String getIonTypeString(int i) { return FragmentIonConstants.ION_NAMES[getIonTypeByte(i)]; }
    public int   getCharge(int i)       { return  fragments[STRIDE * i + 2]         & 0xFF; }

    // ---- Bulk decoders — allocate on demand ----
    public float[] getMzs() {
        int n = numFragments();
        float[] out = new float[n];
        for (int i = 0; i < n; i++) out[i] = getMz(i);
        return out;
    }
    public float[] getIntensities() {
        int n = numFragments();
        float[] out = new float[n];
        for (int i = 0; i < n; i++) out[i] = getIntensity(i);
        return out;
    }
    public int[] getFragNums() {
        int n = numFragments();
        int[] out = new int[n];
        for (int i = 0; i < n; i++) out[i] = getFragNum(i);
        return out;
    }
    public int[] getCharges() {
        int n = numFragments();
        int[] out = new int[n];
        for (int i = 0; i < n; i++) out[i] = getCharge(i);
        return out;
    }
    public String[] getFragmentIonTypes() {
        int n = numFragments();
        String[] out = new String[n];
        for (int i = 0; i < n; i++) out[i] = getIonTypeString(i);
        return out;
    }

    // ---- Constructors ----
    public PredictionEntry() {}

    public PredictionEntry(float[] mzs, float[] intensities, int[] fragNums, int[] charges,
                           String[] fragmentIonTypes) {
        this.fullAnnotations = new String[0];
        this.isotopes = new int[0];
        int[] sortedIndices = sortedIndicesByMz(mzs);
        this.fragments = buildFragments(mzs, intensities, fragNums, charges, fragmentIonTypes, sortedIndices);
    }

    public PredictionEntry(float[] mzs, float[] intensities, int[] fragNums, int[] charges,
                           String[] fragmentIonTypes, String[] fullAnnotations) {
        this.isotopes = new int[0];
        int[] sortedIndices = sortedIndicesByMz(mzs);
        this.fragments = buildFragments(mzs, intensities, fragNums, charges, fragmentIonTypes, sortedIndices);
        setFullAnnotations(fullAnnotations, sortedIndices);
    }

    // ---- Internal pack helpers ----
    private static int[] buildFragments(float[] mzs, float[] intensities, int[] fragNums,
                                        int[] charges, String[] fragmentIonTypes,
                                        int[] sortedIndices) {
        int n = mzs.length;
        int[] packed = new int[STRIDE * n];
        for (int i = 0; i < n; i++) {
            int src = sortedIndices[i];
            int ionByte = (fragmentIonTypes != null && fragmentIonTypes.length != 0)
                    ? FragmentIonConstants.ION_INDEX.getOrDefault(fragmentIonTypes[src],
                                                                   FragmentIonConstants.UNKNOWN_ION_BYTE)
                    : FragmentIonConstants.UNKNOWN_ION_BYTE;
            packSlot(packed, i,
                    mzs[src],
                    intensities[src],
                    (fragNums != null && fragNums.length != 0)   ? fragNums[src]  : 0,
                    ionByte,
                    (charges  != null && charges.length  != 0)   ? charges[src]   : 0);
        }
        return packed;
    }

    static void packSlot(int[] arr, int i, float mz, float intensity,
                         int fragNum, int ionByte, int charge) {
        arr[STRIDE * i]     = Float.floatToRawIntBits(mz);
        arr[STRIDE * i + 1] = Float.floatToRawIntBits(intensity);
        arr[STRIDE * i + 2] = (fragNum & 0xFF) << 16 | (ionByte & 0xFF) << 8 | (charge & 0xFF);
    }

    // Sort indices by ascending m/z without boxing (valid for positive floats).
    private static int[] sortedIndicesByMz(float[] mzs) {
        long[] keys = new long[mzs.length];
        for (int i = 0; i < mzs.length; i++) {
            keys[i] = ((long) Float.floatToRawIntBits(mzs[i]) << 32) | i;
        }
        Arrays.sort(keys);
        int[] indices = new int[keys.length];
        for (int i = 0; i < keys.length; i++) {
            indices[i] = (int) keys[i];
        }
        return indices;
    }

    public void preprocessFragments(HashSet<String> fragmentIonTypesSet, int numTopFragments) {
        int n = numFragments();
        if (n == 0) return;

        mergeCloseMzs();
        n = numFragments();

        // Decode intensities once for all filtering steps
        float[] intensities = new float[n];
        for (int i = 0; i < n; i++) intensities[i] = getIntensity(i);

        // Normalize to max = 1
        float maxInt = 0f;
        for (float v : intensities) if (v > maxInt) maxInt = v;
        if (maxInt != maxIntensity) for (int i = 0; i < n; i++) intensities[i] /= maxInt;

        int potentialFragments = n;
        float[] tmpInts = new float[n];
        System.arraycopy(intensities, 0, tmpInts, 0, n);

        // Filter by allowed fragment ion types
        if (!fragmentIonTypesSet.isEmpty()) {
            for (int i = 0; i < n; i++) {
                if (!fragmentIonTypesSet.contains(getIonTypeString(i))) {
                    tmpInts[i] = 0f;
                    potentialFragments--;
                }
            }
        }

        // Filter by percent base peak
        if (Constants.useBasePeak && Constants.percentBasePeak < 100) {
            maxInt = 0f;
            for (float v : tmpInts) if (v > maxInt) maxInt = v;
            float cutoff = Constants.percentBasePeak / 100f * maxInt;
            for (int i = 0; i < n; i++) {
                if (tmpInts[i] < cutoff && tmpInts[i] != 0f) {
                    tmpInts[i] = 0f;
                    potentialFragments--;
                }
            }
        }

        // Select top N fragments
        if (potentialFragments > numTopFragments && Constants.useTopFragments) {
            potentialFragments = numTopFragments;
            for (int k = 0; k < numTopFragments; k++) {
                maxInt = tmpInts[0];
                int index = 0;
                for (int j = 1; j < n; j++) {
                    if (tmpInts[j] > maxInt) { maxInt = tmpInts[j]; index = j; }
                }
                tmpInts[index] = -1f;
            }
        } else {
            for (int i = 0; i < n; i++) {
                if (tmpInts[i] != 0f) tmpInts[i] = -1f;
            }
        }

        // Rebuild packed array with normalized intensities for kept fragments
        int[] newFragments = new int[STRIDE * potentialFragments];
        String[] newFullAnnotations = fullAnnotations.length > 0 ? new String[potentialFragments] : fullAnnotations;
        int[] newIsotopes = isotopes.length > 0 ? new int[potentialFragments] : isotopes;
        int addIdx = 0;
        for (int i = 0; i < n; i++) {
            boolean keep = Constants.useTopFragments ? (tmpInts[i] == -1f) : (tmpInts[i] != 0f);
            if (keep) {
                newFragments[STRIDE * addIdx]     = fragments[STRIDE * i];                     // mz bits
                newFragments[STRIDE * addIdx + 1] = Float.floatToRawIntBits(intensities[i]);   // normalized intensity
                newFragments[STRIDE * addIdx + 2] = fragments[STRIDE * i + 2];                 // meta
                if (fullAnnotations.length > 0) newFullAnnotations[addIdx] = fullAnnotations[i];
                if (isotopes.length > 0)        newIsotopes[addIdx] = isotopes[i];
                addIdx++;
            }
        }
        fragments = newFragments;
        fullAnnotations = newFullAnnotations;
        isotopes = newIsotopes;
    }

    private void mergeCloseMzs() {
        int n = numFragments();
        HashSet<Integer> excludedIdx = new HashSet<>();
        for (int i = 0; i < n - 1; i++) {
            if (getMz(i + 1) - getMz(i) < 0.00001f) {
                excludedIdx.add(i);

                float mergedMz = (getMz(i + 1) + getMz(i)) / 2;
                float mergedIntensity = getIntensity(i + 1) + getIntensity(i);

                // Pick winner by fragment ion hierarchy
                int metaI    = fragments[STRIDE * i + 2];
                int metaNext = fragments[STRIDE * (i + 1) + 2];
                int winnerMeta = metaNext;
                for (String fit : FragmentIonConstants.fragmentIonHierarchy) {
                    if (fit.equals(getIonTypeString(i))) {
                        winnerMeta = metaI;
                        break;
                    } else if (fit.equals(getIonTypeString(i + 1))) {
                        break;
                    }
                }

                fragments[STRIDE * (i + 1)]     = Float.floatToRawIntBits(mergedMz);
                fragments[STRIDE * (i + 1) + 1] = Float.floatToRawIntBits(mergedIntensity);
                fragments[STRIDE * (i + 1) + 2] = winnerMeta;

                if (fullAnnotations.length > 0 && winnerMeta == metaI)
                    fullAnnotations[i + 1] = fullAnnotations[i];
                if (isotopes.length > 0 && winnerMeta == metaI)
                    isotopes[i + 1] = isotopes[i];
            }
        }
        if (!excludedIdx.isEmpty()) {
            int newN = n - excludedIdx.size();
            int[] newFragments = new int[STRIDE * newN];
            String[] newFullAnnotations = fullAnnotations.length > 0 ? new String[newN] : fullAnnotations;
            int[] newIsotopes = isotopes.length > 0 ? new int[newN] : isotopes;
            int idx = 0;
            for (int j = 0; j < n; j++) {
                if (!excludedIdx.contains(j)) {
                    newFragments[STRIDE * idx]     = fragments[STRIDE * j];
                    newFragments[STRIDE * idx + 1] = fragments[STRIDE * j + 1];
                    newFragments[STRIDE * idx + 2] = fragments[STRIDE * j + 2];
                    if (fullAnnotations.length > 0) newFullAnnotations[idx] = fullAnnotations[j];
                    if (isotopes.length > 0)        newIsotopes[idx] = isotopes[j];
                    idx++;
                }
            }
            fragments = newFragments;
            fullAnnotations = newFullAnnotations;
            isotopes = newIsotopes;
        }
    }

    public void setRT(float RT) { this.RT = RT; }
    public float getRT() { return RT; }
    public void setIM(float IM) { this.IM = IM; }
    public float getIM() { return IM; }
    public int[] getIsotopes() { return isotopes; }
    public String[] getFullAnnotations() { return fullAnnotations; }

    public void setFullAnnotations(String[] fa, int[] sortedIndices) {
        this.fullAnnotations = new String[fa.length];
        this.isotopes = new int[fa.length];

        for (int i = 0; i < sortedIndices.length; i++) {
            this.fullAnnotations[i] = fa[sortedIndices[i]];
        }

        int i = 0;
        for (String annotation : this.fullAnnotations) {
            if (annotation.startsWith("Int")) {
                annotation = annotation.split("/")[1];
                this.fullAnnotations[i] = annotation;
            }
            if (annotation.endsWith("i")) {
                char num = annotation.charAt(annotation.length() - 2);
                isotopes[i] = (num == '+') ? 1 : Integer.parseInt(String.valueOf(num));
            } else {
                isotopes[i] = 0;
            }
            i++;
        }
    }
}
