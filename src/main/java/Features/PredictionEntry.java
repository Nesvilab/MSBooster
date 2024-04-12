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

package Features;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.stream.IntStream;

public class PredictionEntry {
    float[] mzs = new float[0];
    float[] intensities = new float[0];
    int[] fragNums = new int[0];
    int[] flags = new int[0];
    int[] charges = new int[0];
    String[] fragmentIonTypes = new String[0];
    float RT;
    float IM;
    float detectability;
    int counter;
    HashMap<String, Float[]> scores = new HashMap<>();
    ArrayList<Integer> times = new ArrayList<>();
    double precursorMz = 0d;
    boolean filtered = false;

    public PredictionEntry() {}

    public PredictionEntry(float[] mzs, float[] intensities, int[] fragNums, int[] charges,
                           String[] fragmentIonTypes) {

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
        if (fragmentIonTypes.length != 0) {
            setFlags();
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

    public void setMzs(float[] mzs) {this.mzs = mzs;}
    public float[] getMzs() {return mzs;}

    public void setIntensities(float[] intensities) {this.intensities = intensities;}
    public float[] getIntensities() {return intensities;}

    public void setFragNums(int[] fragNums) {this.fragNums = fragNums;}
    public int[] getFragNums() {return fragNums;}

    public void setFlags(int[] flags) {
        this.flags = flags;
        this.fragmentIonTypes = new String[flags.length];
        for (int i = 0; i < flags.length; i++) {
            this.fragmentIonTypes[i] = Constants.flagTOion.get(flags[i]);
        }
    }
    public void setFlags() {
        this.flags = new int[fragmentIonTypes.length];
        for (int i = 0; i < fragmentIonTypes.length; i++) {
            this.flags[i] = Constants.ionTOflag.get(fragmentIonTypes[i]);
        }
    }
    public int[] getFlags() {return flags;}

    public void setCharges(int[] charges) {this.charges = charges;}
    public int[] getCharges() {return charges;}

    public void setRT(float RT) {this.RT = RT;}
    public float getRT() {return RT;}

    public void setIM(float IM) {this.IM = IM;}
    public float getIM() {return IM;}

    public void setDetectability(float detectability) {this.detectability = detectability;}

    public void setCounter(int counts) {this.counter = counts;}

    public void setFragmentIonTypes(String[] ions) { this.fragmentIonTypes = ions; }
    public void setFragmentIonTypes() {
        this.fragmentIonTypes = new String[flags.length];
        for (int i = 0; i < flags.length; i++) {
            this.fragmentIonTypes[i] = Constants.flagTOion.get(flags[i]);
        }
    }
    public String[] getFragmentIonTypes() {return fragmentIonTypes;}
}
