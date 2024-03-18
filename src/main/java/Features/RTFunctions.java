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

import static utils.Print.printInfo;

import umich.ms.fileio.exceptions.FileParsingException;
import umontreal.ssj.probdist.EmpiricalDist;

import java.io.IOException;
import java.util.*;
import java.util.stream.IntStream;

public class RTFunctions {
    //given experimental RT, what should it be on predicted scale?
    public static float normalizeRT(float[] betas, float expRT) {
        return betas[0] + betas[1] * expRT;
    }

//    public static float[] normalizeRT(float[] betas, float[] expRTs) {
//        float[] newRTs = new float[expRTs.length];
//        for (int i = 0; i < expRTs.length; i++) {
//            newRTs[i] = normalizeRT(betas, expRTs[i]);
//        }
//        return newRTs;
//    }
//
//    public static float[] calculateDeltaRT(float[] newRTs, float[] predRTs) {
//        //newRTs can also be array of mean predRTs from bins
//        float[] deltaRTs = new float[predRTs.length];
//        for (int i = 0; i < predRTs.length; i++) {
//            deltaRTs[i] = Math.abs(newRTs[i] - predRTs[i]);
//        }
//        return deltaRTs;
//    }

    //given mzmlreader, get all peptide objects and their RTs
    //get predicted RTs
    //assumes peptide objects already set
    //TODO: only supporting getting PSMs for supported PTMs
    public static float[] getBetas(MzmlReader mzml, int RTregressionSize) throws FileParsingException {
        double[][] RTs = getRTarrays(mzml, RTregressionSize).get("");
        return StatMethods.linearRegression(RTs[0], RTs[1]);
    }

    public static float[] getBetas(double[][] RTs) {
        return StatMethods.linearRegression(RTs[0], RTs[1]);
    }

    private static void getPSMs(MzmlReader mzml,
                         ArrayList<Float> expRTs, ArrayList<Float> predRTs, ArrayList<Float> eScores,
                         ArrayList<String> peptides, HashMap<String, ArrayList<Integer>> pepIdx,
                         float escoreThreshold) throws FileParsingException {
        int added = 0;
        for (int scanNum : mzml.getScanNums()) {
            MzmlScanNumber scanNumObj = mzml.getScanNumObject(scanNum);
            float rt = scanNumObj.RT; //experimental RT for this scan
            if (Float.isNaN(rt)) {
                continue;
            }
            for (PeptideObj pep : scanNumObj.peptideObjects) {
                if (pep == null) {
                    continue;
                }

                float e = Float.parseFloat(pep.escore);
                if (e > escoreThreshold) {
                    continue;
                }
                expRTs.add(rt);
                predRTs.add(pep.RT);
                eScores.add(e);
                peptides.add(pep.name);
                if (pepIdx.containsKey(pep.name)) {
                    ArrayList<Integer> tmpList = pepIdx.get(pep.name);
                    tmpList.add(added);
                    pepIdx.put(pep.name, tmpList);
                } else {
                    ArrayList<Integer> tmpList = new ArrayList<>();
                    tmpList.add(added);
                    pepIdx.put(pep.name, tmpList);
                }
                added += 1;
            }
        }
    }

    //utility for getBetas and LOESS
    //returns exp and pred RT arrays
    //if speed is issue, don't need extra steps for regular regression
    public static HashMap<String, double[][]> getRTarrays(MzmlReader mzml, int RTregressionSize) throws FileParsingException {
        ArrayList<Float> expRTs = new ArrayList<>();
        ArrayList<Float> predRTs = new ArrayList<>();
        ArrayList<Float> eScores = new ArrayList<>(); //for sorting
        ArrayList<String> peptides = new ArrayList<>();
        HashMap<String, ArrayList<Integer>> pepIdx = new HashMap<>();
        getPSMs(mzml, expRTs, predRTs, eScores, peptides, pepIdx, Constants.RTescoreCutoff);

        if (expRTs.size() < Constants.minRTregressionSize) { //no more e score threshold
            printInfo("Not enough high quality PSMs for RT regression with escore cutoff of "
                    + Constants.RTescoreCutoff + ". Relaxing escore cutoff to 0.01");
            RTregressionSize = Constants.minRTregressionSize;

            expRTs = new ArrayList<>();
            predRTs = new ArrayList<>();
            eScores = new ArrayList<>(); //for sorting
            peptides = new ArrayList<>();
            pepIdx = new HashMap<>();
            getPSMs(mzml, expRTs, predRTs, eScores, peptides, pepIdx, 0.01f);
        }

        //remove duplicate PSMs from same peptide
        ArrayList<Integer> PSMtoRemove = new ArrayList<>();
        for (Map.Entry<String, ArrayList<Integer>> entry : pepIdx.entrySet()) {
            if (entry.getValue().size() == 1) {
                continue;
            }
            ArrayList<Float> escores = new ArrayList<>(); //escores stores e values of PSMs
            for (int escoreIdx : entry.getValue()) {
                escores.add(eScores.get(escoreIdx));
            }
            int bestScore = 0;
            ArrayList<Integer> indicesBestScores = new ArrayList<>();
            indicesBestScores.add(0);
            for (int i = 1; i < escores.size(); i++) {
                if (escores.get(i) < escores.get(bestScore)) {
                    bestScore = i;
                    indicesBestScores.clear();
                    indicesBestScores.add(i);
                } else if (Objects.equals(escores.get(i), escores.get(bestScore))) {
                    indicesBestScores.add(i);
                }
            }
            if (indicesBestScores.size() > 1) {
                bestScore = indicesBestScores.get(indicesBestScores.size() / 2);
            }
            for (int i = 0; i < escores.size(); i++) {
                if (i != bestScore) {
                    PSMtoRemove.add(entry.getValue().get(i));
                }
            }
        }
        Collections.sort(PSMtoRemove);
        Collections.reverse(PSMtoRemove);
        for (int removed : PSMtoRemove) {
            expRTs.remove(removed);
            predRTs.remove(removed);
            eScores.remove(removed);
            peptides.remove(removed);
        }

        //get masses that need to be calibrated
        HashMap<String, double[][]> RTs = new HashMap<>();
        if (!Constants.RTmassesForCalibration.isEmpty()) {
            String[] masses = Constants.RTmassesForCalibration.split(",");
            for (String mass : masses) {
                RTs.put(mass, new double[2][]);
            }
            RTs.put("others", new double[2][]);
        } else {
            RTs.put("", new double[2][]);
        }

        for (String mass : RTs.keySet()) {
            ArrayList<Float> thisExpRTs = new ArrayList<>();
            ArrayList<Float> thisPredRTs = new ArrayList<>();
            ArrayList<Float> thisEscores = new ArrayList<>();
            //get PSMs specific to this mass
            if (mass.isEmpty()) {
                thisExpRTs = expRTs;
                thisPredRTs = predRTs;
                thisEscores = eScores;
            } else if (mass.equals("others")) {
                for (int i = 0; i < peptides.size(); i++) {
                    boolean peptideContains = false;
                    for (String m : RTs.keySet()) {
                        String[] masses = m.split("/");
                        for (String minimass : masses) {
                            if (peptides.get(i).contains(minimass)) {
                                peptideContains = true;
                                break;
                            }
                        }
                    }

                    if (!peptideContains) {
                        thisExpRTs.add(expRTs.get(i));
                        thisPredRTs.add(predRTs.get(i));
                        thisEscores.add(eScores.get(i));
                    }
                }
            } else {
                String[] masses = mass.split("/");
                for (int i = 0; i < peptides.size(); i++) {
                    for (String minimass : masses) {
                        if (peptides.get(i).contains(minimass)) {
                            thisExpRTs.add(expRTs.get(i));
                            thisPredRTs.add(predRTs.get(i));
                            thisEscores.add(eScores.get(i));
                            break;
                        }
                    }
                }
            }

            //get top peptides based on eScore
            //https://stackoverflow.com/questions/4859261/get-the-indices-of-an-array-after-sorting
            //also consider taking them all, if want more samples for regression/z scoring
            //then divide into bins with constant size (higher precursor density in middle of RT)

            //if negative, use all
            //can consider e score cutoff in constants
            int sizeLimit = thisExpRTs.size();
            if (sizeLimit < 2) { //hard coded
                if (mass.isEmpty()) {
                    printInfo("Warning: not enough target PSMs (" + sizeLimit + ") are available for regression" +
                            ", setting RT scores equal to 0");
                } else {
                    printInfo("Warning: not enough target PSMs (" + sizeLimit + ") are available for regression" +
                            " for mass " + mass + ", will use RT calibration curve for regular peptides if available");
                }

                //just so that there's an output
                double[][] thisRTs = new double[2][2];
                thisRTs[0] = new double[2];
                thisRTs[1] = new double[2];
                RTs.put(mass, thisRTs);
            } else if (RTregressionSize > 0 && RTregressionSize <= sizeLimit) {
                if (mass.isEmpty()) {
                    printInfo("RT regression using " + RTregressionSize + " PSMs");
                } else {
                    printInfo("RT regression for mass " + mass + " using " + RTregressionSize + " PSMs");
                }
                int[] sortedIndices = IntStream.range(0, thisEscores.size())
                        .boxed().sorted(Comparator.comparing(thisEscores::get))
                        .mapToInt(ele -> ele).toArray();

                int[] sortedIndices2 = Arrays.copyOfRange(sortedIndices, 0, RTregressionSize);

                double[][] thisRTs = new double[2][RTregressionSize];
                for (int i = 0; i < RTregressionSize; i++) {
                    int idx = sortedIndices2[i];
                    thisRTs[0][i] = thisExpRTs.get(idx);
                    thisRTs[1][i] = thisPredRTs.get(idx);
                }
                RTs.put(mass, thisRTs);
            } else {
                if (mass.isEmpty()) {
                    printInfo("RT regression using " + sizeLimit + " PSMs");
                } else {
                    printInfo("RT regression for mass " + mass + " using " + sizeLimit + " PSMs");
                }
                double[][] thisRTs = new double[2][];
                thisRTs[0] = thisExpRTs.stream().mapToDouble(i -> i).toArray();
                thisRTs[1] = thisPredRTs.stream().mapToDouble(i -> i).toArray();
                RTs.put(mass, thisRTs);
            }
        }
        return RTs;
    }

    //function that returns double[] of bin boundaries, with mean and var od each
    //exp on x axis, pred on y
    //to do: test for normal dist nature of each bin (Shapiro Wilk test)
    public static ArrayList<Float>[] RTbins(MzmlReader mzml) throws IOException, FileParsingException {
        //get max index
        float maxRT = -1;
        for (int num : mzml.getScanNums()) {
            MzmlScanNumber msn = mzml.getScanNumObject(num);
            if (msn.RT > maxRT) {
                maxRT = msn.RT;
            }
        }
        int numBins = (int) (maxRT * Constants.RTbinMultiplier) + 1;

        ArrayList<Float>[] predRTround = new ArrayList[numBins + 1];
        for (int col = 0; col < numBins + 1; col++) {
            predRTround[col] = new ArrayList<Float>();
        }

        //iterate through scanNumbers
        for (int scanNum : mzml.getScanNums()) {
            MzmlScanNumber scanNumObj = mzml.getScanNumObject(scanNum);
            int round = (int) (scanNumObj.RT *  Constants.RTbinMultiplier); //experimental RT for this scan, assume in minutes

            //iterate through PSMs
            for (PeptideObj pep : scanNumObj.peptideObjects) {
                if (pep == null) {
                    continue;
                }

                int instances = Math.max(1, -1 * (int) Math.ceil(Math.log10(Double.parseDouble(pep.escore)))); //this version avoids empty bins
                for (int j = 0; j < instances; j++) {
                    predRTround[round].add(pep.RT);
                }
            }
        }

        return predRTround;
    }

    private static EmpiricalDist generateEmpiricalDist(ArrayList<Float> bin) { //could also use umontreal version
        int binSize = bin.size();

        if (binSize > 1) {
            //check that this is sorted
            float maxFloat = bin.get(0);
            boolean needToSort = false;
            for (float f : bin) {
                if (f >= maxFloat) {
                    maxFloat = f;
                } else {
                    needToSort = true;
                    break;
                }
            }
            if (needToSort) {
                Collections.sort(bin);
            }
            return new EmpiricalDist(FloatUtils.floatToDouble(bin));
        } else if (binSize == 1) {
            return new EmpiricalDist(new double[]{(double) bin.get(0), (double) bin.get(0)});
        } else {
            return null; //empty bin never accessed
        }
        //get probability at point with KernelDensity.p();
    }

    public static EmpiricalDist[] generateEmpiricalDist(ArrayList<Float>[] bins) {
        EmpiricalDist[] empDists = new EmpiricalDist[bins.length];

        int i = 0;
        for (ArrayList<Float> bin : bins) {
            empDists[i] = generateEmpiricalDist(bin);
            i++;
        }

        return empDists;
    }

//    public static double RTprobability(float expRT, float predRT, KernelDensity[] bins) {
//        //get right bin to search
//        KernelDensity kd = bins[Math.round(expRT)];
//
//        //check probability at point
//        try {
//            return kd.p(predRT);
//        } catch (Exception e) { //nothing in bin
//            return 0;
//        }
//    }
//
//    public static float RTprobabilityWithUniformPrior(int unifPriorSize, float unifProb,
//                                                       int binSize, float empiricalProb) {
//        float w1 = (float) unifPriorSize / (float) (unifPriorSize + binSize);
//        float w2 = (float) binSize / (float) (unifPriorSize + binSize);
//
//        return w1 * unifProb + w2 * empiricalProb;
//    }
}
