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

import java.io.IOException;
import java.util.*;
import java.util.stream.IntStream;

public class IMFunctions {
    public static int numCharges = 7;

    public static double[][][] getIMarrays(MzmlReader mzml, int IMregressionSize) {
        double[][][] finalIMs = new double[numCharges][2][];

        //ArrayList<ArrayList<Float>> expIMs = new ArrayList<ArrayList<Float>>(numCharges);
        ArrayList<ArrayList<Float>> predIMs = new ArrayList<ArrayList<Float>>(numCharges);
        ArrayList<ArrayList<Float>> eScores = new ArrayList<ArrayList<Float>>(numCharges); //for sorting
        ArrayList<ArrayList<String>> peptides = new ArrayList<>(numCharges);
        HashMap<String, ArrayList<Float>> IMoccurences = new HashMap<>();
        for (int i = 0; i < numCharges; i++) {
            //expIMs.add(new ArrayList<>());
            predIMs.add(new ArrayList<>());
            eScores.add(new ArrayList<>());
            peptides.add(new ArrayList<>());
        }

        //collect RTs and escores

        for (int scanNum : new TreeSet<Integer>(mzml.scanNumberObjects.keySet())) {
            MzmlScanNumber scanNumObj = mzml.getScanNumObject(scanNum);
            float im = scanNumObj.IM; //experimental RT for this scan

            for (int i = 1; i < scanNumObj.peptideObjects.length + 1; i++) {
                PeptideObj pep = scanNumObj.getPeptideObject(i);
                if (pep == null) {
                    break;
                }
                int charge = pep.charge - 1;
                float e = Float.parseFloat(pep.escore);
                if (e > Constants.IMescoreCutoff) {
                    break;
                }
                //use this if using experimental IM straight up
                //expIMs.get(charge).add(im);

                //use this if using some summary of all occurences of peptide
                if (! IMoccurences.containsKey(pep.name)) {
                    ArrayList<Float> arrayList = new ArrayList<>();
                    arrayList.add(im);
                    IMoccurences.put(pep.name, arrayList);
                } else {
                    IMoccurences.get(pep.name).add(im);
                }

                predIMs.get(charge).add(pep.IM);
                eScores.get(charge).add(e);
                peptides.get(charge).add(pep.name);
            }
        }

        //get top peptides based on eScore
        //https://stackoverflow.com/questions/4859261/get-the-indices-of-an-array-after-sorting
        //also consider taking them all, if want more samples for regression/z scoring
        //then divide into bins with constant size (higher precursor density in middle of RT)

        //if negative, use all
        //can consider e score cutoff in constants
        for (int c = 0; c < numCharges; c++) {
            int sizeLimit = predIMs.get(c).size();
            if (sizeLimit == 0) { //no PSMs with this charge
                continue;
            }

            int imSize = Math.min(sizeLimit, IMregressionSize);
            //get top e score PSMs
            int[] sortedIndices = IntStream.range(0, eScores.get(c).size())
                    .boxed().sorted(Comparator.comparing(eScores.get(c)::get))
                    .mapToInt(ele -> ele).toArray();

            double[][] IMs = new double[2][imSize];
            for (int i = 0; i < imSize; i++) {
                int idx = sortedIndices[i];
                //use this if using experimental IM straight up
                //IMs[0][i] = expIMs.get(c).get(idx);

                //use this if using some summary of all occurences of peptide
                IMs[0][i] = StatMethods.median(IMoccurences.get(peptides.get(c).get(idx)));

                IMs[1][i] = predIMs.get(c).get(idx);
            }

            //sort by experimental IM
            int[] sortedIndices2 = IntStream.range(0, imSize)
                    .boxed().sorted(Comparator.comparingDouble(k -> IMs[0][k]))
                    .mapToInt(ele -> ele).toArray();

            double[][] IMs2 = new double[2][imSize];
            for (int i = 0; i < imSize; i++) {
                int idx = sortedIndices2[i];
                IMs2[0][i] = IMs[0][idx];
                IMs2[1][i] = IMs[1][idx];
            }

             finalIMs[c] = IMs2;
        }

        return finalIMs;
    }

    public static ArrayList<Float>[][] IMbins(MzmlReader mzml) throws IOException {
        //hard coded as 2, but if there are higher IM values, this can change
        int numBins = 2 * Constants.IMbinMultiplier;

        ArrayList<Float>[][] predIMround = new ArrayList[numCharges][numBins + 1];
        for (int c = 0; c < numCharges; c++) {
            for (int col = 0; col < numBins + 1; col++) {
                predIMround[c][col] = new ArrayList<Float>();
            }
        }

        //iterate through scanNumbers
        for (int scanNum : mzml.scanNumberObjects.keySet()) {
            MzmlScanNumber scanNumObj = mzml.getScanNumObject(scanNum);
            int round = (int) (scanNumObj.IM * Constants.IMbinMultiplier); //experimental RT for this scan, assume in minutes

            //iterate through PSMs
            for (int i = 1; i < scanNumObj.peptideObjects.length + 1; i++) {
                PeptideObj pep = scanNumObj.getPeptideObject(i);
                if (pep == null) {
                    break;
                }
                int charge = pep.charge - 1;

                int instances = Math.max(1, -1 * (int) Math.ceil(Math.log10(Double.parseDouble(pep.escore)))); //this version avoids empty bins
                for (int j = 0; j < instances; j++) {
                    predIMround[charge][round].add(pep.IM);
                }
            }
        }
        return predIMround;
    }
}
