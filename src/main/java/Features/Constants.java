package Features;

import java.util.Arrays;
import java.util.HashMap;
import java.util.List;

public class Constants {
    //these two constants are for weighted spectral similarity features, not currently supported
    public static final int binwidth = 1;
    public static final int mzFreqWindow = 1;

//    public static final double highScoringProp = 0.1;

    public static float ppmTolerance = 20f; //ppm tolerance of MS2 scans

    public static final boolean basePeakNormalization = true;

    //these two constants for limiting number of fragments used (need to adapt for DIANN)
    public static boolean useTopFragments = true;
    public static int topFragments = 12;

    public static final int fineTuneSize = 100; //for generating a finetune file for pDeep3

    //these two constants for RT features
    public static int RTregressionSize = 5000;
    public static double uniformPriorPercentile = 10;

    public static String outputFolder = "C:/Users/kevin/Downloads/proteomics/"; //where to write all intermediate and final files

    public static final long uid = 0L; //for when I tried to serialize classes, not used now

    public static boolean removeRankPeaks = true; //whether to remove peaks from higher ranks

    //LOESS
    public static double bandwidth = 0.1;
    public static int robustIters = 2;

    public static String predFileFormat = "mgf"; //or mgf
    public static String predFile = null; //use this if predFile already made

    public static String[] usedFeatures = new String[] {"brayCurtis", "euclideanDistance", "cosineSimilarity",
            "spectralContrastAngle", "pearsonCorr", "dotProduct", "deltaRTlinear", "deltaRTbins", "RTzscore",
            "RTprobability", "RTprobabilityUnifPrior"};

    //don't currently support weighted similarity features
    public static final List<String> allowedFeatures = Arrays.asList("cosineSimilarity", "weightedCosineSimilarity",
            "spectralContrastAngle", "weightedSpectralContrastAngle", "euclideanDistance", "weightedEuclideanDistance",
            "brayCurtis", "weightedBrayCurtis", "pearsonCorr", "weightedPearsonCorr", "dotProduct", "weightedDotProduct",
            "detectability", "deltaRTlinear", "deltaRTbins", "deltaRTLOESS", "RTzscore", "RTprobability", "RTprobabilityUnifPrior");
    private static HashMap<Double, Integer> makeModAAToUnimod() {
        HashMap<Double, Integer> map = new HashMap<>();
        map.put(160.03065, 4);
        map.put(147.0354, 35);
        return map;
    }
    public static final HashMap<Double, Integer> modAAmassToUnimod = makeModAAToUnimod();
    private static HashMap<Double, Integer> makeModNtermToUnimod() {
        HashMap<Double, Integer> map = new HashMap<>();
        map.put(43.018425, 1);
        return map;
    }
    public static final HashMap<Double, Integer> modNtermToUnimod = makeModNtermToUnimod();

    public Constants() {

    }
}
