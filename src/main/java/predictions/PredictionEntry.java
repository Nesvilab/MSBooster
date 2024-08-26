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
import features.spectra.MassCalculator;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.stream.IntStream;

import static utils.Print.printError;

public class PredictionEntry {
    public float[] mzs = new float[0];
    public float[] intensities = new float[0];
    public int[] fragNums = new int[0];
    public int[] flags = new int[0];
    public int[] charges = new int[0];
    public String[] fragmentIonTypes = new String[0];
    public float RT;
    public float IM;
    float detectability;
    int counter;
    public HashMap<String, Float[]> scores = new HashMap<>();
    public ArrayList<Integer> times = new ArrayList<>();
    public double precursorMz = 0d;
    public boolean filtered = false;

    public PredictionEntry() {}

    public PredictionEntry(float[] mzs, float[] intensities, int[] fragNums, int[] charges,
                           String[] fragmentIonTypes, int[] flags) {

        this.mzs = new float[mzs.length];
        this.intensities = new float[intensities.length];
        this.fragNums = new int[fragNums.length];
        this.charges = new int[charges.length];
        this.fragmentIonTypes = new String[fragmentIonTypes.length];
        this.flags = new int[flags.length];

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
            if (flags.length != 0) {
                this.flags[i] = flags[sortedIndices[i]];
            }
        }
        if (fragmentIonTypes.length != 0 && flags.length == 0) {
            setFlags();
        } else if (flags.length != 0 && fragmentIonTypes.length == 0) {
            setFragmentIonTypes();
        }
    }

    public void filterFragments() {
        if (!filtered) {
            filtered = true;
            if (intensities.length != 0) {
                if (Constants.useBasePeak || Constants.useTopFragments) {
                    int potentialFragments = intensities.length;
                    float[] tmpInts = new float[potentialFragments];
                    for (int i = 0; i < potentialFragments; i++) {
                        tmpInts[i] = intensities[i];
                    }
                    if (Constants.useBasePeak) {
                        if (Constants.percentBasePeak < 100) {
                            //get max intensity
                            float maxIntensity = 0f;
                            for (float f : intensities) {
                                if (f > maxIntensity) {
                                    maxIntensity = f;
                                }
                            }

                            //make cutoff by percentage
                            float cutoff = Constants.percentBasePeak / 100f * maxIntensity;
                            for (int i = 0; i < tmpInts.length; i++) {
                                if (tmpInts[i] < cutoff) {
                                    tmpInts[i] = 0f;
                                    potentialFragments--;
                                }
                            }
                        }
                    }

                    if ((potentialFragments > Constants.topFragments) &&
                            (Constants.useTopFragments)) {
                        potentialFragments = Constants.topFragments;

                        //setting highest intensities to -1
                        for (int i = 0; i < Constants.topFragments; i++) {
                            float maxInt = tmpInts[0];
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
                    int[] pflags = new int[potentialFragments];
                    int[] pcharges = new int[potentialFragments];
                    String[] pfragmentIonTypes = new String[potentialFragments];
                    int addIdx = 0;
                    for (int i = 0; i < tmpInts.length; i++) {
                        if (Constants.useTopFragments) {
                            if (tmpInts[i] == -1f) {
                                predIntensities[addIdx] = intensities[i];
                                predMZs[addIdx] = mzs[i];
                                if (fragNums.length > 0) {
                                    pfragNums[addIdx] = fragNums[i];
                                }
                                if (flags.length > 0) {
                                    pflags[addIdx] = flags[i];
                                }
                                if (charges.length > 0) {
                                    pcharges[addIdx] = charges[i];
                                }
                                if (fragmentIonTypes.length > 0) {
                                    pfragmentIonTypes[addIdx] = fragmentIonTypes[i];
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
                                if (flags.length > 0) {
                                    pflags[addIdx] = flags[i];
                                }
                                if (charges.length > 0) {
                                    pcharges[addIdx] = charges[i];
                                }
                                if (fragmentIonTypes.length > 0) {
                                    pfragmentIonTypes[addIdx] = fragmentIonTypes[i];
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
                    if (flags.length > 0) {
                        flags = pflags;
                    }
                    if (charges.length > 0) {
                        charges = pcharges;
                    }
                    if (fragmentIonTypes.length > 0) {
                        fragmentIonTypes = pfragmentIonTypes;
                    }
                }
            }
        }
    }

    public float[] getMzs() {return mzs;}
    public float[] getIntensities() {return intensities;}
    public int[] getFragNums() {return fragNums;}
    private void setFlags() {
        this.flags = new int[fragmentIonTypes.length];
        for (int i = 0; i < fragmentIonTypes.length; i++) {
            this.flags[i] = MassCalculator.ionTOflag.get(fragmentIonTypes[i]);
        }
    }
    public int[] getFlags() {return flags;}
    public int[] getCharges() {return charges;}
    public void setRT(float RT) {this.RT = RT;}
    public float getRT() {return RT;}
    public void setIM(float IM) {this.IM = IM;}
    public float getIM() {return IM;}
    public void setDetectability(float detectability) {this.detectability = detectability;}
    private void setFragmentIonTypes() {
        this.fragmentIonTypes = new String[flags.length];
        for (int i = 0; i < flags.length; i++) {
            this.fragmentIonTypes[i] = MassCalculator.flagTOion.get(flags[i]);
        }
    }
    public String[] getFragmentIonTypes() {return fragmentIonTypes;}
}
