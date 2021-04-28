package Features;

import java.io.File;
import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;

public class Constants {
    //file input
    public static String paramsList = null;

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    //file locations
    public static String pinPepXMLDirectory = null; //C:/Users/kevin/Downloads/proteomics/cptac/2021-2-21/
    public static String mzmlDirectory = null; //C:/Users/kevin/OneDriveUmich/proteomics/mzml/cptac/
    public static String outputDirectory = null; //where to write all intermediate and final files
    public static String editedPin = null;
    public static String spectraRTPredInput = null;
    public static String detectPredInput = null;
    public static String spectraRTPredFile = null; //use this if predFile already made
    public static String detectPredFile = null;

    //optional file locations and parameters
    //if calculating detectFractionGreater, these are used for FastaReader class
    public static String fasta = null; //C:/Users/kevin/OneDriveUmich/proteomics/fasta/2019-09-30-td-rev-UP000005640.fas
    public static String decoyPrefix = ">rev_";
    public static String cutAfter = "KR";
    public static String butNotAfter = "P";
    public static Integer digestMinLength = 7;
    public static Integer digestMaxLength = 50;
    public static Integer digestMinMass = 500; //Da
    public static Integer digestMaxMass = 5000;

    //locations of executables and other models
    public static Integer numThreads = 11;
    public static String DiaNN = null; //C:/DIA-NN/1.7.15beta1/DiaNN.exe
    public static String spectraRTPredModel = "DIA-NN"; //mgf, bin, msp
                                                        //pDeep3, DIA-NN, Prosit
                                                        //DIANN by default

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //miscellaneous inner workings

    //these two constants are for weighted spectral similarity features, not currently supported
    public static Integer binwidth = 1;
    public static Integer mzFreqWindow = 1;

    public static Float ppmTolerance = 20f; //ppm tolerance of MS2 scans

    //for limiting number of fragments used
    //TODO need to adapt for DIANN
    public static Boolean useTopFragments = true;
    public static Integer topFragments = 12;
    public static Boolean removeRankPeaks = true; //whether to remove peaks from higher ranks

    public static Integer fineTuneSize = 100; //for generating a finetune file for pDeep3

    //these constants for RT features
    public static Integer RTregressionSize = 5000;
    public static Double uniformPriorPercentile = 10d;
    public static Float RTescoreCutoff = 1f; //PSMs with e score higher than this won't make it into RT linear regression modeling

    //LOESS
    public static Double bandwidth = 0.1;
    public static Integer robustIters = 2;

    public static Float detectThreshold = 0.0000002f; //for detectability filtering

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    //use single string sep by comma delimiter
    //should include parameter to calculate correlation and then choose
    //default auto, everything, or all? Or a combination I figure out empirically
    public static String features = "brayCurtis,euclideanDistance,cosineSimilarity," +
            "spectralContrastAngle,pearsonCorr,dotProduct," +
            "deltaRTlinear,deltaRTbins,RTzscore,RTprobability,RTprobabilityUnifPrior,deltaRTLOESS";
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

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //Handling PTMs

    //TODO: better handling of PTMs, all in one location
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

    //methods
    public void updatePaths() {
        if (outputDirectory == null) {
            outputDirectory = pinPepXMLDirectory;
        }
        if (editedPin == null) {
            editedPin = outputDirectory + File.separator + "edited.pin";
        }
        if (spectraRTPredInput == null) {
            spectraRTPredInput = outputDirectory + File.separator + "spectraRT.tsv";
        }
        if (detectPredInput == null) {
            detectPredInput = outputDirectory + File.separator + "detect.tsv";
        }
    }
}
