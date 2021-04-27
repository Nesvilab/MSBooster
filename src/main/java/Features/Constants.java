package Features;

import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;

public class Constants {
    //file input
    public static String paramsList = null;

    //these two constants are for weighted spectral similarity features, not currently supported
    public static final Integer binwidth = 1;
    public static final Integer mzFreqWindow = 1;

//    public static final double highScoringProp = 0.1;

    public static Float ppmTolerance = 20f; //ppm tolerance of MS2 scans

    public static final Boolean basePeakNormalization = true;

    //these two constants for limiting number of fragments used (need to adapt for DIANN)
    public static Boolean useTopFragments = true;
    public static Integer topFragments = 12;

    public static final Integer fineTuneSize = 100; //for generating a finetune file for pDeep3

    //these constants for RT features
    public static Integer RTregressionSize = 5000;
    public static Double uniformPriorPercentile = 10d;
    public static Float RTescoreCutoff = 1f; //at this point, these won't make it into regression modeling

    //for detectability filtering
    public static Float detectThreshold = 0.0000002f;

    public static String fasta = "C:/Users/kevin/OneDriveUmich/proteomics/fasta/2019-09-30-td-rev-UP000005640.fas";
    public static String decoyPrefix = ">rev_";
    public static String cutAfter = "KR";
    public static String butNotAfter = "P";
    public static Integer digestMinLength = 7;
    public static Integer digestMaxLength = 50;
    public static Integer digestMinMass = 500; //Da
    public static Integer digestMaxMass = 5000;

    public static String outputFolder = "C:/Users/kevin/Downloads/proteomics/"; //where to write all intermediate and final files

    public static final long uid = 0L; //for when I tried to serialize classes, not used now

    public static Boolean removeRankPeaks = true; //whether to remove peaks from higher ranks

    //LOESS
    public static Double bandwidth = 0.1;
    public static Integer robustIters = 2;

    public static String spectraRTPredFileFormat = "bin"; //mgf, bin, msp
                                                          //DIANN by default
    public static String spectraRTPredFile = null; //use this if predFile already made
    public static String detectPredFile = null;

    //use single string sep by delimiter
    //by default include everything
    //should include parameter to calculate correlation and then choose
    //default auto? Or something a combination I figure out empirically
    public static String features = "brayCurtis,euclideanDistance,cosineSimilarity," +
            "spectralContrastAngle,pearsonCorr,dotProduct," +
            "deltaRTlinear,deltaRTbins,RTzscore,RTprobability,RTprobabilityUnifPrior," +
            "detectablity,detectFractionGreater";
    //public static String features = "auto";

    //don't currently support weighted similarity features
    public static final HashSet<String> allowedFeatures = new HashSet<>(Arrays.asList("all", "auto",
            "cosineSimilarity", "weightedCosineSimilarity", "spectralContrastAngle", "weightedSpectralContrastAngle",
            "euclideanDistance", "weightedEuclideanDistance", "brayCurtis", "weightedBrayCurtis",
            "pearsonCorr", "weightedPearsonCorr", "dotProduct", "weightedDotProduct",
            "deltaRTlinear", "deltaRTbins", "deltaRTLOESS", "RTzscore", "RTprobability", "RTprobabilityUnifPrior",
            "detectFractionGreater", "detectability"));
    public static final HashSet<String> detectFeatures = new HashSet<>(Arrays.asList("detectFractionGreater", "detectability"));
    public static final HashSet<String> spectraRTFeatures = new HashSet<>(Arrays.asList(
            "cosineSimilarity", "weightedCosineSimilarity", "spectralContrastAngle", "weightedSpectralContrastAngle",
            "euclideanDistance", "weightedEuclideanDistance", "brayCurtis", "weightedBrayCurtis",
            "pearsonCorr", "weightedPearsonCorr", "dotProduct", "weightedDotProduct",
            "deltaRTlinear", "deltaRTbins", "deltaRTLOESS", "RTzscore", "RTprobability", "RTprobabilityUnifPrior"));

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
