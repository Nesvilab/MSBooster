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

import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.stream.IntStream;

public class PredictionEntry {
    public float[] mzs = new float[0];
    public float[] intensities = new float[0];
    public int[] fragNums = new int[0];
    public int[] charges = new int[0];
    public String[] fragmentIonTypes = new String[0];
    public String[] fullAnnotations = new String[0];
    public int[] isotopes = new int[0];
    public float RT;
    public float IM;
    public PredictionEntry auxSpectra = null;
    public HashMap<String, Float[]> scores = new HashMap<>();
    public ArrayList<Integer> times = new ArrayList<>();
    public double precursorMz = 0d;
    //public boolean filtered = false; //work with this if filtering step is unnecessarily run multiple times
    private static final float maxIntensity = 1f;
    public boolean daltonMatching = false;

    public PredictionEntry() {}

    public PredictionEntry(float[] mzs, float[] intensities, int[] fragNums, int[] charges,
                           String[] fragmentIonTypes) {

        this.mzs = new float[mzs.length];
        this.intensities = new float[intensities.length];
        this.fragNums = new int[fragNums.length];
        this.charges = new int[charges.length];
        this.fragmentIonTypes = new String[fragmentIonTypes.length];
        this.isotopes = new int[mzs.length];
        this.fullAnnotations = new String[mzs.length];

        int[] sortedIndices = IntStream.range(0, mzs.length)
                .boxed().sorted((k, j) -> Float.compare(mzs[k], mzs[j]))
                .mapToInt(ele -> ele).toArray();

        for (int i = 0; i < sortedIndices.length; i++) {
            this.mzs[i] = mzs[sortedIndices[i]];
            this.intensities[i] = intensities[sortedIndices[i]];
            if (fragNums.length != 0) {
                this.fragNums[i] = fragNums[sortedIndices[i]];
            }
            if (charges.length != 0) {
                this.charges[i] = charges[sortedIndices[i]];
            }
            if (fragmentIonTypes.length != 0) {
                this.fragmentIonTypes[i] = fragmentIonTypes[sortedIndices[i]];
            }
        }
    }

    public PredictionEntry(float[] mzs, float[] intensities, int[] fragNums, int[] charges,
                           String[] fragmentIonTypes, String[] fullAnnotations) {

        this.mzs = new float[mzs.length];
        this.intensities = new float[intensities.length];
        this.fragNums = new int[fragNums.length];
        this.charges = new int[charges.length];
        this.fragmentIonTypes = new String[fragmentIonTypes.length];

        int[] sortedIndices = IntStream.range(0, mzs.length)
                .boxed().sorted((k, j) -> Float.compare(mzs[k], mzs[j]))
                .mapToInt(ele -> ele).toArray();

        for (int i = 0; i < sortedIndices.length; i++) {
            this.mzs[i] = mzs[sortedIndices[i]];
            this.intensities[i] = intensities[sortedIndices[i]];
            if (fragNums.length != 0) {
                this.fragNums[i] = fragNums[sortedIndices[i]];
            }
            if (charges.length != 0) {
                this.charges[i] = charges[sortedIndices[i]];
            }
            if (fragmentIonTypes.length != 0) {
                this.fragmentIonTypes[i] = fragmentIonTypes[sortedIndices[i]];
            }
        }

        setFullAnnotations(fullAnnotations, sortedIndices);
    }

    public void preprocessFragments(HashSet<String> fragmentIonTypesSet, int numTopFragments) {
        if (intensities.length != 0) {
            mergeCloseMzs();

            //set max intensity as 1
            //do it before fragment ion filtering so it considers truly highest intensity peak as intensity 1
            float maxInt = 0f;
            for (int i = 0; i < intensities.length; i++) {
                if (intensities[i] > maxInt) {
                    maxInt = intensities[i];
                }
            }
            if (maxInt != maxIntensity) {
                for (int i = 0; i < intensities.length; i++) {
                    intensities[i] /= maxInt;
                }
            }

            int potentialFragments = intensities.length; //number of fragments to consider
            float[] tmpInts = new float[potentialFragments]; //indicates if fragment should be used (0 means no)
            System.arraycopy(intensities, 0, tmpInts, 0, potentialFragments);

            //filter for allowed fragment ion types
            if (fragmentIonTypes.length > 0 && !fragmentIonTypesSet.isEmpty()) {
                for (int i = 0; i < fragmentIonTypes.length; i++) {
                    if (!fragmentIonTypesSet.contains(fragmentIonTypes[i])) {
                        tmpInts[i] = 0f;
                        potentialFragments--;
                    }
                }
            }

            //above intensity threshold
            if (Constants.useBasePeak && Constants.percentBasePeak < 100) {
                //get max intensity
                maxInt = 0f;
                for (float f : tmpInts) {
                    if (f > maxInt) {
                        maxInt = f;
                    }
                }

                //make cutoff by percentage
                float cutoff = Constants.percentBasePeak / 100f * maxInt;
                for (int i = 0; i < tmpInts.length; i++) {
                    if (tmpInts[i] < cutoff && tmpInts[i] != 0f) {
                        tmpInts[i] = 0f;
                        potentialFragments--;
                    }
                }
            }

            //only use N top fragments
            if ((potentialFragments > numTopFragments) && Constants.useTopFragments) {
                potentialFragments = numTopFragments;

                //setting highest intensities to -1
                for (int i = 0; i < numTopFragments; i++) {
                    maxInt = tmpInts[0];
                    int index = 0;
                    for (int j = 1; j < tmpInts.length; j++) {
                        float thisInt = tmpInts[j];
                        if (thisInt > maxInt) {
                            maxInt = thisInt;
                            index = j;
                        }
                    }
                    tmpInts[index] = -1f; //no longer highest intensity peak
                }
            } else {
                for (int i = 0; i < tmpInts.length; i++) {
                    if (tmpInts[i] != 0f) {
                        tmpInts[i] = -1f;
                    }
                }
            }

            //assigning final values
            float[] predIntensities = new float[potentialFragments];
            float[] predMZs = new float[potentialFragments];
            int[] pfragNums = new int[potentialFragments];
            int[] pcharges = new int[potentialFragments];
            String[] pfragmentIonTypes = new String[potentialFragments];
            String[] pfullAnnotations = new String[potentialFragments];
            int[] pisotopes = new int[potentialFragments];
            int addIdx = 0;
            for (int i = 0; i < tmpInts.length; i++) {
                if (Constants.useTopFragments) {
                    if (tmpInts[i] == -1f) { //TODO: could make the steps below a method
                        predIntensities[addIdx] = intensities[i];
                        predMZs[addIdx] = mzs[i];
                        if (fragNums.length > 0) {
                            pfragNums[addIdx] = fragNums[i];
                        }
                        if (charges.length > 0) {
                            pcharges[addIdx] = charges[i];
                        }
                        if (fragmentIonTypes.length > 0) {
                            pfragmentIonTypes[addIdx] = fragmentIonTypes[i];
                        }
                        if (fullAnnotations.length > 0) {
                            pfullAnnotations[addIdx] = fullAnnotations[i];
                        }
                        if (isotopes.length > 0) {
                            pisotopes[addIdx] = isotopes[i];
                        }
                        addIdx++;
                    }
                } else {
                    if (tmpInts[i] != 0f) {
                        predIntensities[addIdx] = intensities[i];
                        predMZs[addIdx] = mzs[i];
                        if (fragNums.length > 0) {
                            pfragNums[addIdx] = fragNums[i];
                        }
                        if (charges.length > 0) {
                            pcharges[addIdx] = charges[i];
                        }
                        if (fragmentIonTypes.length > 0) {
                            pfragmentIonTypes[addIdx] = fragmentIonTypes[i];
                        }
                        if (fullAnnotations.length > 0) {
                            pfullAnnotations[addIdx] = fullAnnotations[i];
                        }
                        if (isotopes.length > 0) {
                            pisotopes[addIdx] = isotopes[i];
                        }
                        addIdx++;
                    }
                }
            }

            mzs = predMZs;
            intensities = predIntensities;
            if (fragNums.length > 0) {
                fragNums = pfragNums;
            }
            if (charges.length > 0) {
                charges = pcharges;
            }
            if (fragmentIonTypes.length > 0) {
                fragmentIonTypes = pfragmentIonTypes;
            }
            if (fullAnnotations.length > 0) {
                fullAnnotations = pfullAnnotations;
            }
            if (isotopes.length > 0) {
                isotopes = pisotopes;
            }
        }
    }

    private void mergeCloseMzs() {
        //first, need to merge mzs that are too close to distinguish
        HashSet<Integer> excludedIdx = new HashSet<>();
        for (int i = 0; i < mzs.length - 1; i++) {
            if (mzs[i + 1] - mzs[i] < 0.00001f) {
                excludedIdx.add(i);

                mzs[i + 1] = (mzs[i + 1] + mzs[i]) / 2;
                intensities[i + 1] += intensities[i];

                //now need to decide how to annotate it
                //currently only choosing whichever annotation is first in fragment ion hierarchy
                for (String fit : FragmentIonConstants.fragmentIonHierarchy) {
                    if (fit.equals(fragmentIonTypes[i])) {
                        fragmentIonTypes[i + 1] = fragmentIonTypes[i];
                        if (fragNums.length > 0) {
                            fragNums[i + 1] = fragNums[i];
                        }
                        if (charges.length > 0) {
                            charges[i + 1] = charges[i];
                        }
                        if (fullAnnotations.length > 0) {
                            fullAnnotations[i + 1] = fullAnnotations[i];
                        }
                        if (isotopes.length > 0) {
                            isotopes[i + 1] = isotopes[i];
                        }
                        break;
                    } else if (fit.equals(fragmentIonTypes[i + 1])) {
                        break;
                    }
                }
            }
        }
        if (!excludedIdx.isEmpty()) {
            int newLength = mzs.length - excludedIdx.size();
            float[] predMZs = new float[newLength];
            float[] predIntensities = new float[newLength];
            int[] pfragNums = new int[newLength];
            int[] pcharges = new int[newLength];
            String[] pfragmentIonTypes = new String[newLength];
            String[] pfullAnnotations = new String[newLength];
            int[] pisotopes = new int[newLength];

            int i = 0;
            for (int j = 0; j < mzs.length; j++) {
                if (!excludedIdx.contains(j)) {
                    predMZs[i] = mzs[j];
                    predIntensities[i] = intensities[j];
                    if (fragNums.length > 0) {
                        pfragNums[i] = fragNums[j];
                    }
                    if (charges.length > 0) {
                        pcharges[i] = charges[j];
                    }
                    if (fragmentIonTypes.length > 0) {
                        pfragmentIonTypes[i] = fragmentIonTypes[j];
                    }
                    if (isotopes.length > 0) {
                        pisotopes[i] = isotopes[j];
                    }
                    if (fullAnnotations.length > 0) {
                        pfullAnnotations[i] = fullAnnotations[j];
                    }
                    i++;
                }
            }

            mzs = predMZs;
            intensities = predIntensities;
            if (fragNums.length > 0) {
                fragNums = pfragNums;
            }
            if (charges.length > 0) {
                charges = pcharges;
            }
            if (fragmentIonTypes.length > 0) {
                fragmentIonTypes = pfragmentIonTypes;
            }
            if (fullAnnotations.length > 0) {
                fullAnnotations = pfullAnnotations;
            }
            if (isotopes.length > 0) {
                isotopes = pisotopes;
            }
        }
    }

    public float[] getMzs() {return mzs;}
    public float[] getIntensities() {return intensities;}
    public int[] getFragNums() {return fragNums;}
    public int[] getCharges() {return charges;}
    public void setRT(float RT) {this.RT = RT;}
    public float getRT() {return RT;}
    public void setIM(float IM) {this.IM = IM;}
    public float getIM() {return IM;}
    public String[] getFragmentIonTypes() {return fragmentIonTypes;}
    public int[] getIsotopes() {return isotopes;}
    public void setFullAnnotations(String[] fa, int[] sortedIndices) {
        this.fullAnnotations = new String[fa.length];
        this.isotopes = new int[fa.length];

        for (int i = 0; i < sortedIndices.length; i++) {
            this.fullAnnotations[i] = fa[sortedIndices[i]];
        }

        //set isotopes
        int i = 0;
        for (String annotation : this.fullAnnotations) {
            if (annotation.startsWith("Int")) {
                annotation = annotation.split("/")[1];
                this.fullAnnotations[i] = annotation;
            }

            if (annotation.endsWith("i")) {
                char num = annotation.charAt(annotation.length() - 2);
                if (num == '+') {
                    isotopes[i] = 1;
                } else {
                    isotopes[i] = Integer.parseInt(String.valueOf(num));
                }
            } else {
                isotopes[i] = 0;
            }
            i++;
        }
    }
    public String[] getFullAnnotations() {return fullAnnotations;}
}
