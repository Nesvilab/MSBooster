package Features;

import com.github.sanity.pav.PairAdjacentViolators;
import com.github.sanity.pav.Point;
import kotlin.jvm.functions.Function1;
import org.apache.commons.math3.analysis.interpolation.LoessInterpolator;
import smile.stat.distribution.KernelDensity;

import java.io.IOException;
import java.util.*;
import java.util.stream.IntStream;

public class RTFunctions {
    //given experimental RT, what should it be on predicted scale?
    public static float normalizeRT(float[] betas, float expRT) {
        return betas[0] + betas[1] * expRT;
    }

    public static float[] normalizeRT(float[] betas, float[] expRTs) {
        float[] newRTs = new float[expRTs.length];
        for (int i = 0; i < expRTs.length; i++) {
            newRTs[i] = normalizeRT(betas, expRTs[i]);
        }
        return newRTs;
    }

    public static float[] calculateDeltaRT(float[] newRTs, float[] predRTs) {
        //newRTs can also be array of mean predRTs from bins
        float[] deltaRTs = new float[predRTs.length];
        for (int i = 0; i < predRTs.length; i++) {
            deltaRTs[i] = Math.abs(newRTs[i] - predRTs[i]);
        }
        return deltaRTs;
    }

    //given mzmlreader, get all peptide objects and their RTs
    //get predicted RTs
    //assumes peptide objects already set
    public static float[] getBetas(mzMLReader mzml, SpectralPredictionMapper preds, int RTregressionSize) {
        double[][] RTs = getRTarrays(mzml, preds, RTregressionSize);
        return StatMethods.linearRegression(RTs[0], RTs[1]);
    }

    public static float[] getBetas(double[][] RTs) {
        return StatMethods.linearRegression(RTs[0], RTs[1]);
    }

    //utility for getBetas and LOESS
    //returns exp and pred RT arrays
    //if speed is issue, don't need extra steps for regular regression
    public static double[][] getRTarrays(mzMLReader mzml, SpectralPredictionMapper preds, int RTregressionSize) {
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
            Arrays.sort(sortedIndices2);

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

    public static Function1<Double, Double> LOESS(mzMLReader mzml, SpectralPredictionMapper preds, int RTregressionSize,
                             double bandwidth, int robustIters) {
        double[][] RTs = getRTarrays(mzml, preds, RTregressionSize);

        return LOESS(RTs, bandwidth, robustIters);
    }

    public static Function1<Double, Double> LOESS(double[][] RTs, double bandwidth, int robustIters) {
        //solve monotonicity issue
        double compare = -1;
        for (int i = 0; i < RTs[0].length; i++) {
            double d = RTs[0][i];
            if (d == compare) {
                RTs[0][i] = RTs[0][i - 1] + 0.00000001; //arbitrary increment to handle smooth method
            } else {
                compare = d;
            }
        }

        //fit loess
        //may not need if sanity version uses spline
        LoessInterpolator loessInterpolator = new LoessInterpolator(bandwidth, robustIters);
        double[] y = loessInterpolator.smooth(RTs[0], RTs[1]);

        //isotonic regression
        List<Point> points = new LinkedList<>();
        for (int i = 0; i < y.length; i++) {
            points.add(new Point(RTs[0][i], y[i]));
        }
        PairAdjacentViolators pav = new PairAdjacentViolators(points);

        return pav.interpolator(); //extrapolationStrategy can use flat or tangent (default tangent)
    }

    //function that returns double[] of bin boundaries, with mean and var od each
    //exp on x axis, pred on y
    //to do: test for normal dist nature of each bin (Shapiro Wilk test)
    public static ArrayList<Float>[] RTbins(mzMLReader mzml, SpectralPredictionMapper preds) throws IOException {
        //get max index
        int maxKey = Collections.max(mzml.scanNumberObjects.keySet());
        int numBins = (int) mzml.scanNumberObjects.get(maxKey).RT + 1;

        ArrayList<Float>[] predRTround = new ArrayList[numBins + 1];
        for (int col = 0; col < numBins + 1; col++) {
            predRTround[col] = new ArrayList<Float>();
        }

        //iterate through scanNumbers
        for (int scanNum : mzml.scanNumberObjects.keySet()) {
            mzmlScanNumber scanNumObj = mzml.getScanNumObject(scanNum);
            float rt = scanNumObj.RT; //experimental RT for this scan, assume in minutes
            //int floor = (int) rt;
            int round = Math.round(rt);

            //iterate through PSMs
            for (int i = 1; i < scanNumObj.peptideObjects.size() + 1; i++) {
                peptideObj pep = scanNumObj.getPeptideObject(i);

                //decide instances to add based on eScore (subject to change)
//                int instances = -1 * (int) Math.ceil(Math.log10(Double.parseDouble(pep.escore)));
//                if (instances > 0) { //add directly to distribution now
//                    for (int j = 0; j < instances; j++) {
//                        //predRTfloor[td][floor].add(pep.RT);
//                        predRTround[round].add(pep.RT);
//                    }
//                } else { //expect score already too high, no more instances to add
//                    break;
//                }
                int instances = Math.max(1, -1 * (int) Math.ceil(Math.log10(Double.parseDouble(pep.escore)))); //this version avoids empty bins
                for (int j = 0; j < instances; j++) {
                    predRTround[round].add(pep.RT);
                }
            }
        }

        return predRTround;
    }

    //get means and standard devs of each RTbin
    public static float[][] characterizeRTbins(ArrayList<Float>[] bins) {
        float[][] binStats = new float[bins.length][2]; //index by expRT, then mean or standard deviation
        for (int i = 0; i < bins.length; i++) {
            float m = StatMethods.mean(bins[i]);
            float sd = (float) Math.sqrt(StatMethods.variance(bins[i], m));
            //deals with NaN sd, basically discard it
            //adding small constant if one entry, lower expectation score
            if (Float.isNaN(sd) || sd == 0f) {
                sd = Float.POSITIVE_INFINITY;
            }

            binStats[i][0] = m;
            binStats[i][1] = sd;
        }
        return binStats;
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

    public static double RTprobability(float expRT, float predRT, KernelDensity[] bins) {
        //get right bin to search
        KernelDensity kd = bins[Math.round(expRT)];

        //check probability at point
        try {
            return kd.p(predRT);
        } catch (Exception e) { //nothing in bin
            return 0;
        }
    }

    public static float RTprobabilityWithUniformPrior(int unifPriorSize, float unifProb,
                                                       int binSize, float empiricalProb) {
        float w1 = (float) unifPriorSize / (float) (unifPriorSize + binSize);
        float w2 = (float) binSize / (float) (unifPriorSize + binSize);

        return w1 * unifProb + w2 * empiricalProb;
    }

    public static void main(String[] args) throws Exception {
        int window = 6;
        int maxRank = 3;

        //mgfFileReader preds = new mgfFileReader("C:/Users/kevin/Downloads/proteomics/pDeep3preds_wide1.mgf");
        mgfFileReader preds = new mgfFileReader("C:/Users/kevin/OneDriveUmich/proteomics/preds/narrowPDeep3.mgf");

        System.out.println("loading mzml");
        mzMLReader m = new mzMLReader("C:/Users/kevin/OneDriveUmich/proteomics/mzml/narrowWindow/" +
                "23aug2017_hela_serum_timecourse_4mz_narrow_" + window + ".mzML");

        //fill with pepxml files
        System.out.println("loading pepxml");
        for (int rank = 1; rank < maxRank + 1; rank++) {
            pepXMLReader p = new pepXMLReader("C:/Users/kevin/Downloads/proteomics/narrow/"
                    + "23aug2017_hela_serum_timecourse_4mz_narrow_" + window + "_rank" + rank + ".pepXML");
            m.setPepxmlEntries(p, rank, preds);
        }

        Function1<Double, Double> irm = LOESS(m, preds, Constants.RTregressionSize, 0.3, 2); //default according to documentation
        System.out.println(irm.invoke(10.0));
        System.out.println(irm.invoke(80.0));
        System.out.println(irm.invoke(120.0));
        float[] betas = RTFunctions.getBetas(m, preds, Constants.RTregressionSize);
        System.out.println(Arrays.toString(betas));
    }
}
