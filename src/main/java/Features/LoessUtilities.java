package Features;

import umich.ms.fileio.exceptions.FileParsingException;

import java.util.*;
import java.util.stream.IntStream;

import static Features.Constants.minLinearRegressionSize;
import static utils.Print.printError;
import static utils.Print.printInfo;

public class LoessUtilities {
    //getPSMs returns number of PSMs fitting all criteria except e value cutoff
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

                expValues.add(value);
                if (mode.equals("RT")) {
                    predValues.add(pep.RT);
                } else {
                    predValues.add(pep.IM);
                }
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
}
