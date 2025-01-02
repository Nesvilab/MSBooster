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

package features.rtandim;

import allconstants.Constants;
import com.github.sanity.pav.PairAdjacentViolators;
import com.github.sanity.pav.Point;
import kotlin.jvm.functions.Function1;
import mainsteps.MzmlScanNumber;
import mainsteps.PeptideObj;
import org.apache.commons.math3.analysis.interpolation.LoessInterpolator;
import readers.datareaders.MzmlReader;
import umich.ms.fileio.exceptions.FileParsingException;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.HashSet;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;
import java.util.Objects;
import java.util.stream.IntStream;

import static allconstants.Constants.minLinearRegressionSize;
import static allconstants.Constants.minLoessRegressionSize;
import static utils.Print.printError;
import static utils.Print.printInfo;
import static utils.StatMethods.*;

public class LoessUtilities {
    //getPSMs returns number of PSMs fitting all criteria except e value cutoff
    //private function only used by getArrays
    private static int getPSMs(MzmlReader mzml,
                               ArrayList<Float> expValues, ArrayList<Float> predValues, ArrayList<Float> eScores,
                               ArrayList<String> peptides, HashMap<String, ArrayList<Integer>> pepIdx,
                               float escoreThreshold, String mode, int charge) throws FileParsingException {
        int added = 0;
        int numPSMsIgnoreEvalue = 0;
        for (int scanNum : mzml.getScanNums()) {
            MzmlScanNumber scanNumObj = mzml.getScanNumObject(scanNum);
            float value = Float.NaN;
            if (mode.equals("RT")) {
                value = scanNumObj.RT;
            } else if (mode.equals("IM")) {
                value = scanNumObj.IM;
            } else {
                printError("Mode must be RT or IM. Exiting.");
                System.exit(1);
            }
            if (Float.isNaN(value)) {
                continue;
            }
            for (PeptideObj pep : scanNumObj.peptideObjects) {
                if (pep == null) {
                    continue;
                }

                if (mode.equals("IM") && pep.charge != charge) {
                    continue;
                }
                numPSMsIgnoreEvalue++;

                float e = Float.parseFloat(pep.escore);
                if (e > escoreThreshold) {
                    continue;
                }

                //add values once all criteria is met
                if (mode.equals("RT")) {
                    if (pep.RT == 0) {
                        continue;
                    }
                    predValues.add(pep.RT);
                } else {
                    if (pep.IM == 0) {
                        continue;
                    }
                    predValues.add(pep.IM);
                }
                expValues.add(value);
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
                added++;
            }
        }
        return numPSMsIgnoreEvalue;
    }

    //utility for getBetas and LOESS
    //returns exp and pred RT arrays
    public static Object[] getArrays(MzmlReader mzml, int regressionSize, String mode, int charge)
            throws FileParsingException {
        //returns HashMap<String, double[][]> for arrays, and HashMap<String, ArrayList<String>> for peptides
        ArrayList<Float> expValues = new ArrayList<>();
        ArrayList<Float> predValues = new ArrayList<>();
        ArrayList<Float> eScores = new ArrayList<>(); //for sorting
        ArrayList<String> peptides = new ArrayList<>();
        HashMap<String, ArrayList<Integer>> pepIdx = new HashMap<>();
        int numPSMsIgnoreEvalue = getPSMs(mzml, expValues, predValues, eScores, peptides, pepIdx,
                Constants.loessEscoreCutoff, mode, charge);

        if (expValues.size() < minLinearRegressionSize &&
                Constants.loessEscoreCutoff < 0.01 && numPSMsIgnoreEvalue > 0) { //no more e score threshold
            printInfo("Not enough high quality PSMs for " + mode + " regression with escore cutoff of "
                    + Constants.loessEscoreCutoff + ". Relaxing escore cutoff to 0.01");
            regressionSize = minLinearRegressionSize;

            expValues = new ArrayList<>();
            predValues = new ArrayList<>();
            eScores = new ArrayList<>(); //for sorting
            peptides = new ArrayList<>();
            pepIdx = new HashMap<>();
            getPSMs(mzml, expValues, predValues, eScores, peptides, pepIdx, 0.01f, mode, charge);
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
            expValues.remove(removed);
            predValues.remove(removed);
            eScores.remove(removed);
            peptides.remove(removed);
        }

        //get masses that need to be calibrated
        HashMap<String, double[][]> massToDataMap = new HashMap<>();
        ArrayList<String> massesList = new ArrayList<>();
        if (!Constants.massesForLoessCalibration.isEmpty()) {
            String[] masses = Constants.massesForLoessCalibration.split(",");
            massesList.addAll(Arrays.asList(masses));
            massesList.add("others");
        } else {
            massesList.add("");
        }

        HashMap<String, ArrayList<String>> peptideMap = new HashMap<>();

        //create ignorable mass offsets
        HashSet<String> ignorableMassOffsets = new HashSet<>();
        for (String mass : massesList) {
            if (mass.contains("&")) {
                String[] masses = mass.split("&");
                ignorableMassOffsets.addAll(List.of(masses));
            }
        }

        for (String mass : massesList) {
            ArrayList<Float> thisExpValues = new ArrayList<>();
            ArrayList<Float> thisPredValues = new ArrayList<>();
            ArrayList<Float> thisEscores = new ArrayList<>();
            ArrayList<String> finalPeptides = new ArrayList<>();
            //get PSMs specific to this mass
            if (mass.isEmpty()) {
                thisExpValues = expValues;
                thisPredValues = predValues;
                thisEscores = eScores;
                finalPeptides = peptides;
            } else if (mass.equals("others")) {
                for (int i = 0; i < peptides.size(); i++) {
                    boolean peptideContains = false;
                    for (String m : massesList) {
                        if (m.equals("others")) {
                            continue;
                        }
                        String[] masses = m.split("&");
                        for (String minimass : masses) {
                            if (peptides.get(i).contains(minimass)) {
                                peptideContains = true;
                                break;
                            }
                        }
                    }

                    if (!peptideContains) {
                        thisExpValues.add(expValues.get(i));
                        thisPredValues.add(predValues.get(i));
                        thisEscores.add(eScores.get(i));
                        finalPeptides.add(peptides.get(i));
                    }
                }
            } else {
                //if it's a regular variable mod, don't include mass offsets
                if (mass.contains("&")) {
                    String[] masses = mass.split("&");
                    for (int i = 0; i < peptides.size(); i++) {
                        for (String minimass : masses) {
                            if (peptides.get(i).contains(minimass)) {
                                thisExpValues.add(expValues.get(i));
                                thisPredValues.add(predValues.get(i));
                                thisEscores.add(eScores.get(i));
                                finalPeptides.add(peptides.get(i));
                                break;
                            }
                        }
                    }
                } else {
                    for (int i = 0; i < peptides.size(); i++) {
                        boolean offsetContinue = false;
                        for (String massOffset : ignorableMassOffsets) {
                            if (peptides.get(i).contains(massOffset)) {
                                offsetContinue = true;
                                break;
                            }
                        }
                        if (offsetContinue) {
                            continue;
                        }
                        if (peptides.get(i).contains(mass)) {
                            thisExpValues.add(expValues.get(i));
                            thisPredValues.add(predValues.get(i));
                            thisEscores.add(eScores.get(i));
                            finalPeptides.add(peptides.get(i));
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
            int sizeLimit = thisExpValues.size();
            StringBuilder messageEnding = new StringBuilder(" PSMs");
            String defaultValue = "0";
            if (mode.equals("IM")) {
                messageEnding.append(" for charge ").append(charge);
                defaultValue = "500";
            }
            if (sizeLimit < minLinearRegressionSize) { //hard coded
                if (numPSMsIgnoreEvalue != 0) {
                    if (mass.isEmpty() || mass.equals("others")) {
                        printInfo("Warning: not enough target PSMs (" + sizeLimit + ") are available for regression" +
                                ", setting " + mode + " scores equal to " + defaultValue);
                        //just so that there's an output
                        massToDataMap.put(mass, null);
                    } else {
                        printInfo("Warning: not enough target PSMs (" + sizeLimit + ") are available for regression" +
                                " for mass " + mass + ", will use " + mode + " calibration curve for regular peptides if available");
                    }
                }
            } else if (regressionSize > 0 && regressionSize <= sizeLimit) {

                if (mass.isEmpty()) {
                    printInfo(mode + " regression using " + regressionSize + messageEnding);
                } else {
                    printInfo(mode + " regression for mass " + mass + " using " + regressionSize + messageEnding);
                }
                int[] sortedIndices = IntStream.range(0, thisEscores.size())
                        .boxed().sorted(Comparator.comparing(thisEscores::get))
                        .mapToInt(ele -> ele).toArray();

                int[] sortedIndices2 = Arrays.copyOfRange(sortedIndices, 0, regressionSize);

                double[][] thisValues = new double[2][regressionSize];
                ArrayList<String> newFinalPeptides = new ArrayList<>();
                for (int i = 0; i < regressionSize; i++) {
                    int idx = sortedIndices2[i];
                    thisValues[0][i] = thisExpValues.get(idx);
                    thisValues[1][i] = thisPredValues.get(idx);
                    newFinalPeptides.add(finalPeptides.get(idx));
                }
                finalPeptides = newFinalPeptides;
                massToDataMap.put(mass, thisValues);
            } else {
                if (mass.isEmpty()) {
                    printInfo(mode + " regression using " + sizeLimit + messageEnding);
                } else {
                    printInfo(mode + " regression for mass " + mass + " using " + sizeLimit + messageEnding);
                }
                double[][] thisValues = new double[2][];
                thisValues[0] = thisExpValues.stream().mapToDouble(i -> i).toArray();
                thisValues[1] = thisPredValues.stream().mapToDouble(i -> i).toArray();
                massToDataMap.put(mass, thisValues);
            }

            peptideMap.put(mass, finalPeptides);
        }
        return new Object[] {massToDataMap, peptideMap};
    }

    private static double getBandwidth(double bandwidth, int numDatapoints) {
        if (numDatapoints < minLoessRegressionSize || bandwidth > 1) {
            bandwidth = 1d; //linear regression
        } else if (bandwidth < (double) 2 / numDatapoints) {
            //the bandwidth must be larger than 2/n
            //https://commons.apache.org/proper/commons-math/javadocs/api-3.6.1/org/apache/commons/math3/analysis/interpolation/LoessInterpolator.html
            bandwidth = (double) 3 / numDatapoints;
        }
        return bandwidth;
    }
    public static Function1<Double, Double> LOESS(double[][] bins, double bandwidth, int robustIters) {
        //need to sort arrays (DIA-U mzml not in order)
        int[] sortedIndices = IntStream.range(0, bins[0].length)
                .boxed().sorted(Comparator.comparingDouble(k -> bins[0][k])).mapToInt(ele -> ele).toArray();
        double[] newX = new double[sortedIndices.length];
        double[] newY = new double[sortedIndices.length];
        for (int i = 0; i < sortedIndices.length; i++) {
            newX[i] = bins[0][sortedIndices[i]];
        }
        for (int i = 0; i < sortedIndices.length; i++) {
            newY[i] = bins[1][sortedIndices[i]];
        }

        //solve monotonicity issue
        double compare = -1;
        for (int i = 0; i < newX.length; i++) {
            double d = newX[i];
            if (d == compare) {
                newX[i] = newX[i - 1] + 0.00000001; //arbitrary increment to handle smooth method
            } else {
                compare = d;
            }
        }

        //check bandwidth
        bandwidth = getBandwidth(bandwidth, newX.length);

        //fit loess
        LoessInterpolator loessInterpolator = new LoessInterpolator(bandwidth, robustIters);
        double[][] results = fittingRound(loessInterpolator, newX, newY, null, false);
        double[] fittedY = results[0];
        newX = results[1];
        newY = results[2];

        //remove Nan
        ArrayList<Integer> nanIdx = new ArrayList<>();
        int i = 0;
        for (double yval : fittedY) {
            if (Double.isNaN(yval)) {
                nanIdx.add(i);
            }
            i += 1;
        }

        fittedY = removeIdx(fittedY, nanIdx);
        newX = removeIdx(newX, nanIdx);
        newY = removeIdx(newY, nanIdx);

        Function1<Double, Double> ir = isotonicRegressor(newX, fittedY);
        for (i = 0; i < fittedY.length; i++) {
            fittedY[i] = ir.invoke(newX[i]);
        }
        results = fittingRound(loessInterpolator, newX, newY, fittedY, true);
        fittedY = results[0];
        newX = results[1];

        return isotonicRegressor(newX, fittedY);
    }

    private static Function1<Double, Double> isotonicRegressor(double[] x, double[] y) {
        List<Point> points = new LinkedList<>();
        for (int i = 0; i < y.length; i++) {
            points.add(new Point(x[i], y[i]));
        }
        PairAdjacentViolators pav = new PairAdjacentViolators(points);

        return pav.interpolator(); //extrapolationStrategy can use flat or tangent (default tangent)
    }

    private static double[][] fittingRound(LoessInterpolator lr, double[] newX, double[] newY, double[] smoothedY,
                                           boolean extra) {
        double[] y;
        if (smoothedY != null) {
            y = smoothedY;
        } else {
            y = lr.smooth(newX, newY);
        }
        if (y.length > 100) {
            if (extra) {
                for (int i = 0; i < y.length; i++) {
                    y[i] = y[i] - newY[i]; //y now represents difference
                }

                //remove outliers
                float meanDiff = mean(y);
                float stdDiff = (float) Math.sqrt(variance(y));
                ArrayList<Integer> outliersIdx = new ArrayList<>();

                for (int i = 0; i < y.length; i++) {
                    if (Math.abs(zscore((float) y[i], meanDiff, stdDiff)) > 2) {
                        outliersIdx.add(i);
                    }
                }
                y = removeIdx(y, outliersIdx);
                newX = removeIdx(newX, outliersIdx);
                newY = removeIdx(newY, outliersIdx);
            }

            double[] weights = new double[y.length];
            for (int i = 0; i < y.length; i++) {
                int start = Math.max(i - y.length / 100, 0);
                int end = Math.min(i + y.length / 100, y.length);
                if (start != end) {
                    weights[i] = Math.abs(median(Arrays.copyOfRange(y, start, end)));
                } else {
                    weights[i] = y[i];
                }
            }

            //redo loess with weights
            y = lr.smooth(newX, newY, weights); //weighting so ends become better fit
        }
        double[][] results = new double[3][];
        results[0] = y;
        results[1] = newX;
        results[2] = newY;
        return results;
    }

    private static double[] removeIdx(double[] array, ArrayList<Integer> idx) {
        if (!idx.isEmpty()) {
            ArrayList<Double> newX = new ArrayList<>();
            for (int i = 0; i < array.length; i++) {
                if (! idx.contains(i)) {
                    newX.add(array[i]);
                }
            }
            array = new double[newX.size()];
            for (int i = 0; i < newX.size(); i++) {
                array[i] = newX.get(i);
            }
        }
        return array;
    }

    public static ArrayList<double[][][]> trainTestSplit(double[][] expAndPred) {
        //arraylist of N splits
        //train-test,exp-pred,data
        double[] exp = expAndPred[0];
        double[] pred = expAndPred[1];

        ArrayList<double[][][]> splits = new ArrayList<>();

        for (int rep = 0; rep < Constants.regressionSplits; rep++) {
            //sort into train and test
            List<Double> trainExpList = new ArrayList<>();
            List<Double> trainPredList = new ArrayList<>();
            List<Double> testExpList = new ArrayList<>();
            List<Double> testPredList = new ArrayList<>();

            for (int i = 0; i < exp.length; i++) {
                if (i % Constants.regressionSplits != rep) {
                    trainExpList.add(exp[i]);
                    trainPredList.add(pred[i]);
                } else {
                    testExpList.add(exp[i]);
                    testPredList.add(pred[i]);
                }
            }

            double[][][] split = new double[2][2][];
            split[0][0] = trainExpList.stream().mapToDouble(Double::doubleValue).toArray();
            split[0][1] = trainPredList.stream().mapToDouble(Double::doubleValue).toArray();
            split[1][0] = testExpList.stream().mapToDouble(Double::doubleValue).toArray();
            split[1][1] = testPredList.stream().mapToDouble(Double::doubleValue).toArray();
            splits.add(split);
        }
        return splits;
    }

    //finds the best bandwidth and absolute error between all points
    public static Object[] gridSearchCV(double[][] rts, float[] bandwidths) {
        float[] bestBandwidths = new float[Constants.regressionSplits];

        //divide into train and test sets
        ArrayList<double[][][]> splits = trainTestSplit(rts);

        System.out.print("Iteration ");
        double bestMSE;
        for (int Nsplit = 0; Nsplit < splits.size(); Nsplit++) {
            bestMSE = Double.MAX_VALUE;
            System.out.print(Nsplit + 1 + "...");
            float bestBandwidth = 1f;

            double[][][] split = splits.get(Nsplit);
            double[][] train = split[0];
            double[][] test = split[1];

            for (float floatb : bandwidths) { //for bandwidth in grid search
                //get the loess model
                try {
                    Function1<Double, Double> loess = LOESS(train, floatb, Constants.robustIters);

                    //calculate MSE by comparing calibrated expRT to predRT
                    double[] calibratedRTs = new double[test[0].length];
                    for (int i = 0; i < calibratedRTs.length; i++) {
                        double rt = test[0][i];
                        calibratedRTs[i] = loess.invoke(rt);
                    }
                    double mse = meanSquaredError(calibratedRTs, test[1]);

                    //choose best model
                    if (mse < bestMSE) {
                        bestMSE = mse;
                        bestBandwidth = (float) getBandwidth(floatb, train[0].length);
                    }
                } catch (Exception ignored) {} //bandwidth too small?
            }
            bestBandwidths[Nsplit] = bestBandwidth;
        }
        System.out.println();
        float finalBandwidth = Float.parseFloat(String.format("%.4f", mean(bestBandwidths)));

        //train one more loess and calculate mse
        //or just have all indexes past 1 as rtDiffs
        //final model trained on all data
        bestMSE = Double.MAX_VALUE;
        Function1<Double, Double> loess = null;
        while (true) {
            try {
                loess = LOESS(rts, finalBandwidth, Constants.robustIters);
                float[] rtDiffs = new float[rts[0].length];
                for (int i = 0; i < rtDiffs.length - 1; i++) {
                    double rt = rts[0][i];
                    rtDiffs[i] = (float) Math.abs(loess.invoke(rt) - rts[1][i]);
                }
                return new Object[]{finalBandwidth, loess, rtDiffs};
            } catch (Exception e) {
                //e.printStackTrace();
                if (finalBandwidth == 1) {
                    return new Object[]{finalBandwidth, loess, (float) bestMSE};
                }
                finalBandwidth = Math.min(finalBandwidth * 2, 1);
                //printInfo("Regression failed, retrying with double the bandwidth: " + finalBandwidth);
            }
        }
    }
}
