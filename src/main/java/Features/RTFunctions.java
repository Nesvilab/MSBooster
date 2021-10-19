package Features;

import smile.stat.distribution.KernelDensity;

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
    public static float[] getBetas(mzMLReader mzml, int RTregressionSize) {
        double[][] RTs = getRTarrays(mzml, RTregressionSize);
        return StatMethods.linearRegression(RTs[0], RTs[1]);
    }

    public static float[] getBetas(double[][] RTs) {
        return StatMethods.linearRegression(RTs[0], RTs[1]);
    }

    //utility for getBetas and LOESS
    //returns exp and pred RT arrays
    //if speed is issue, don't need extra steps for regular regression
    public static double[][] getRTarrays(mzMLReader mzml, int RTregressionSize) {
        ArrayList<Float> expRTs = new ArrayList<>();
        ArrayList<Float> predRTs = new ArrayList<>();
        ArrayList<Float> eScores = new ArrayList<>(); //for sorting
        //collect RTs and escores

        for (int scanNum : new TreeSet<Integer>(mzml.scanNumberObjects.keySet())) {
            mzmlScanNumber scanNumObj = mzml.getScanNumObject(scanNum);
            float rt = scanNumObj.RT; //experimental RT for this scan

            //add RT until you reach decoy
            //when doing good regression, decoys don't appear, so once decoy appears, expectation score is already too low

            for (int i = 1; i < scanNumObj.peptideObjects.size() + 1; i++) {
                peptideObj pep = scanNumObj.getPeptideObject(i);
//                if (pep.targetORdecoy == 0) {
//                    break;
//                }
                float e = Float.parseFloat(pep.escore);
                if (e > Constants.RTescoreCutoff) {
                    break;
                }
                expRTs.add(rt);
                predRTs.add(pep.RT);
                eScores.add(e);
            }
        }

        //get top peptides based on eScore
        //https://stackoverflow.com/questions/4859261/get-the-indices-of-an-array-after-sorting
        //also consider taking them all, if want more samples for regression/z scoring
        //then divide into bins with constant size (higher precursor density in middle of RT)

        //if negative, use all
        //can consider e score cutoff in constants
        int sizeLimit = expRTs.size();

        if (RTregressionSize > 0 && RTregressionSize <= sizeLimit) {
            int[] sortedIndices = IntStream.range(0, eScores.size())
                    .boxed().sorted(Comparator.comparing(eScores::get))
                    .mapToInt(ele -> ele).toArray();

            int[] sortedIndices2 = Arrays.copyOfRange(sortedIndices, 0, RTregressionSize);
            Arrays.sort(sortedIndices2); //this ensures increasing RT

            double[][] RTs = new double[2][RTregressionSize];
            for (int i = 0; i < RTregressionSize; i++) {
                int idx = sortedIndices2[i];
                RTs[0][i] = expRTs.get(idx);
                RTs[1][i] = predRTs.get(idx);
            }
            return RTs;

        } else {
            double[][] RTs = new double[2][];
            RTs[0] = expRTs.stream().mapToDouble(i -> i).toArray();
            RTs[1] = predRTs.stream().mapToDouble(i -> i).toArray();
            return RTs;
        }
    }

//    public static Function1<Double, Double> LOESS(mzMLReader mzml, int RTregressionSize,
//                             double bandwidth, int robustIters) {
//        double[][] RTs = getRTarrays(mzml, RTregressionSize);
//
//        return LOESS(RTs, bandwidth, robustIters);
//    }

    //function that returns double[] of bin boundaries, with mean and var od each
    //exp on x axis, pred on y
    //to do: test for normal dist nature of each bin (Shapiro Wilk test)
    public static ArrayList<Float>[] RTbins(mzMLReader mzml) throws IOException {
        //get max index
        float maxRT = -1;
        for (mzmlScanNumber msn : mzml.scanNumberObjects.values()) {
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
        for (int scanNum : mzml.scanNumberObjects.keySet()) {
            mzmlScanNumber scanNumObj = mzml.getScanNumObject(scanNum);
            int round = (int) (scanNumObj.RT *  Constants.RTbinMultiplier); //experimental RT for this scan, assume in minutes

            //iterate through PSMs
            for (int i = 1; i < scanNumObj.peptideObjects.size() + 1; i++) {
                peptideObj pep = scanNumObj.getPeptideObject(i);

                int instances = Math.max(1, -1 * (int) Math.ceil(Math.log10(Double.parseDouble(pep.escore)))); //this version avoids empty bins
                for (int j = 0; j < instances; j++) {
                    predRTround[round].add(pep.RT);
                }
            }
        }

        return predRTround;
    }

    private static KernelDensity generateEmpiricalDist(ArrayList<Float> bin) { //could also use umontreal version
        int binSize = bin.size();
        if (binSize > 1) {
            return new KernelDensity(floatUtils.floatToDouble(bin));
        } else if (binSize == 1) {
            bin.add(bin.get(0)); //just add another entry to make it work
            return new KernelDensity(floatUtils.floatToDouble(bin));
        } else {
            return null; //empty bin never accessed
        }
        //get probability at point with KernelDensity.p();
    }

    public static KernelDensity[] generateEmpiricalDist(ArrayList<Float>[] bins) {
        KernelDensity[] empDists = new KernelDensity[bins.length];

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

    public static void main(String[] args) throws Exception {

    }
}
