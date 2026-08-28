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

import allconstants.Constants;
import features.rtandim.fragalign.FragAlignRegression;
import features.rtandim.fragalign.HermiteSpline;
import mainsteps.MzmlScanNumber;
import mainsteps.PeptideObj;
import readers.datareaders.MzmlReader;
import umich.ms.fileio.exceptions.FileParsingException;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.Map;
import java.util.Objects;
import java.util.stream.IntStream;

import static allconstants.Constants.minLinearRegressionSize;
import static utils.Print.printError;
import static utils.Print.printInfo;

/**
 * Assembles the (experimental, predicted) anchor arrays used for RT and IM
 * calibration, split by mass offset group.
 *
 * <p>Curve fitting itself lives in {@code features.rtandim.fragalign}.
 */
public class LoessUtilities {

    /**
     * Residual reported for every anchor when no calibration curve can be fit,
     * large enough that the model ranks last but small enough that squaring it
     * for an RMSE stays well inside float range.
     */
    private static final float UNFITTABLE_RESIDUAL = 1e6f;

    /**
     * Fits a calibration curve to the given anchors and returns each anchor's
     * absolute residual, which is the per-model score {@code BestModelSearcher}
     * ranks candidate prediction models by.
     *
     * <p>This replaces the bandwidth grid search that used to serve the same
     * purpose. The fitter chooses its own smoothing, so a single fit per model
     * is enough where the grid search needed one per bandwidth per split.
     *
     * @param xy {@code [0]} independent values, {@code [1]} dependent values
     * @return per-anchor absolute residuals; all {@link #UNFITTABLE_RESIDUAL}
     *         when the fit fails, so an unfittable model ranks last instead of
     *         appearing perfect
     */
    public static float[] calibrationResiduals(double[][] xy) {
        float[] residuals = new float[xy[0].length];
        HermiteSpline spline = FragAlignRegression.fit(xy[0], xy[1], null);
        if (spline == null) {
            Arrays.fill(residuals, UNFITTABLE_RESIDUAL);
            return residuals;
        }
        for (int i = 0; i < residuals.length; i++) {
            residuals[i] = (float) Math.abs(spline.evaluate(xy[0][i]) - xy[1][i]);
        }
        return residuals;
    }

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

        //remove duplicate PSMs from same peptide //TODO run this before previous step
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

        //parse each group once, and each peptide's delta masses once, so the assignment below is
        //a numeric comparison rather than a substring search (see MassOffsetGroup)
        HashMap<String, MassOffsetGroup> groups = new HashMap<>();
        for (String mass : massesList) {
            groups.put(mass, MassOffsetGroup.of(mass));
        }
        double[][] peptideDeltaMasses = new double[peptides.size()][];
        for (int i = 0; i < peptides.size(); i++) {
            peptideDeltaMasses[i] = MassOffsetGroup.deltaMasses(peptides.get(i));
        }

        //mass offsets that already share a curve must not be pulled into a plain single-mass group,
        //which stands for a regular variable mod
        ArrayList<MassOffsetGroup> ignorableMassOffsets = new ArrayList<>();
        for (String mass : massesList) {
            if (mass.contains("&")) {
                ignorableMassOffsets.add(groups.get(mass));
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
            } else if (mass.equals(MassOffsetGroup.OTHERS)) {
                for (int i = 0; i < peptides.size(); i++) {
                    boolean peptideContains = false;
                    for (String m : massesList) {
                        if (m.equals(MassOffsetGroup.OTHERS)) {
                            continue;
                        }
                        if (groups.get(m).matches(peptides.get(i), peptideDeltaMasses[i])) {
                            peptideContains = true;
                            break;
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
                MassOffsetGroup group = groups.get(mass);
                //if it's a regular variable mod, don't include mass offsets
                boolean isRegularVariableMod = !mass.contains("&");
                for (int i = 0; i < peptides.size(); i++) {
                    if (isRegularVariableMod) {
                        boolean offsetContinue = false;
                        for (MassOffsetGroup massOffset : ignorableMassOffsets) {
                            if (massOffset.matches(peptides.get(i), peptideDeltaMasses[i])) {
                                offsetContinue = true;
                                break;
                            }
                        }
                        if (offsetContinue) {
                            continue;
                        }
                    }
                    if (group.matches(peptides.get(i), peptideDeltaMasses[i])) {
                        thisExpValues.add(expValues.get(i));
                        thisPredValues.add(predValues.get(i));
                        thisEscores.add(eScores.get(i));
                        finalPeptides.add(peptides.get(i));
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
            StringBuilder messageEnding = new StringBuilder(" precursors");
            String defaultValue = "0";
            if (mode.equals("IM")) {
                messageEnding.append(" for charge ").append(charge);
                defaultValue = "500";
            }
            if (sizeLimit < minLinearRegressionSize) { //hard coded
                if (numPSMsIgnoreEvalue != 0) {
                    if (mass.isEmpty() || mass.equals("others")) {
                        printInfo("Warning: not enough target precursors (" + sizeLimit + ") are available for regression" +
                                ", setting " + mode + " scores equal to " + defaultValue);
                        //just so that there's an output
                        massToDataMap.put(mass, null);
                    } else {
                        printInfo("Warning: not enough target precursors (" + sizeLimit + ") are available for regression" +
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

}
