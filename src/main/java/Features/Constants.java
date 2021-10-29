package Features;

import java.io.File;
import java.util.*;

public class Constants {
    //file input
    public static String paramsList = null;
    public static String fragger = null;

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    //file locations
    //check that pinPepXMLDirectory, mzmlDirectory, spectraRTPredInput, detectPredInput,
    //      spectraRTPredFile, detectPredFile, fasta exist
    public static String pinPepXMLDirectory = null; //C:/Users/kevin/Downloads/proteomics/cptac/2021-2-21/
    public static String mzmlDirectory = null; //C:/Users/kevin/OneDriveUmich/proteomics/mzml/cptac/
    public static String outputDirectory = null; //where to write all intermediate and final files
    public static String editedPin = null; //a prefix. Default is edited_
    public static Integer renamePin = 1;
    public static String spectraRTPredInput = null;
    public static String detectPredInput = null;
    public static String spectraRTPredFile = null; //use this if predFile already made
    public static String detectPredFile = null;

    //optional file locations and parameters
    //if calculating detectFractionGreater, these are used for FastaReader class
    //C:/Users/kevin/OneDriveUmich/proteomics/fasta/2020-12-07-decoys-reviewed-contam-UP000005640.fas
    public static String fasta = "C:/Users/kevin/OneDriveUmich/proteomics/fasta/2020-12-07-decoys-reviewed-contam-UP000005640.fas";
    public static String decoyPrefix = ">rev_";
    public static String cutAfter = "KR";
    public static String butNotAfter = "P";
    public static Integer digestMinLength = 7;
    public static Integer digestMaxLength = 50;
    public static Float digestMinMass = 500f; //Da
    public static Float digestMaxMass = 5000f;
    //public static Boolean includeDecoy = false;
    private static FastaReader fastaReader = null;
    public static void setFastaReader(FastaReader f) {
        fastaReader = f;
    }
    public static FastaReader getFastaReader() {
        return fastaReader;
    }

    //locations of executables and other models
    public static Integer numThreads = 0;
    public static String DiaNN = null; //C:/DIA-NN/1.7.15beta1/DiaNN.exe
    public static final String spectraRTPredModel = "DIA-NN"; //mgf, bin, msp
                                                        //pDeep3, DIA-NN, Prosit
                                                        //DIANN by default
                                                        //currently don't support changing this

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //miscellaneous inner workings

    //these two constants are for weighted spectral similarity features, not currently supported
    public static final Integer binwidth = 1;
    public static final Integer mzFreqWindow = 1;

    public static Float ppmTolerance = 20f; //ppm tolerance of MS2 scans
    public static Float lowResppmTolerance = 300f;
    public static Float highResppmTolerance = 20f;

    //for limiting number of fragments used
    public static Boolean useSpectra = null;
    public static Boolean useTopFragments = true;
    public static Integer topFragments = 12;
    public static Boolean removeRankPeaks = true; //whether to remove peaks from higher ranks

    public static final Integer fineTuneSize = 100; //for generating a finetune file for pDeep3

    //these constants for RT features
    public static Boolean useRT = null;
    public static Integer RTregressionSize = 5000;
    public static Double uniformPriorPercentile = 10d;
    public static Float RTescoreCutoff = (float) Math.pow(10, -3.5); //PSMs with e score higher than this won't make it into RT regression modeling
    public static Integer RTbinMultiplier = 1;
    public static Float RTIQR = 50f;

    //LOESS
    public static Double bandwidth = 0.25;
    public static final Integer robustIters = 2;

    //detect
    public static Boolean useDetect = false;
    public static final Float detectThreshold = 0.0000002f; //for detectability filtering
    public static final Float detectFractionGreaterNumerator = 1f;
    public static final Float detectFractionGreaterDenominator = 2f; //prior

    //ion mobility
    public static Boolean useIM = null;
    public static Integer IMregressionSize = 5000;
    public static Float IMescoreCutoff = (float) Math.pow(10, -3.5);
    public static Integer IMbinMultiplier = 100;
    public static final Float IMIQR = 50f;

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    //use single string sep by comma delimiter
    //should include parameter to calculate correlation and then choose
    //default auto, everything, or all? Or a combination I figure out empirically
    public static String features = "cosineSimilarity,spectralContrastAngle,euclideanDistance,brayCurtis,pearsonCorr,dotProduct," +
            "deltaRTLOESS,deltaRTLOESSnormalized,RTprobabilityUnifPrior," +
            "detectProtSpearmanDiff";
    //public static String features = "auto";

    //don't currently support weighted similarity features
    public static final HashSet<String> allowedFeatures = new HashSet<>(Arrays.asList(
            "cosineSimilarity", "spectralContrastAngle",
            "euclideanDistance", "brayCurtis",
            "pearsonCorr", "dotProduct",
            "deltaRTlinear", "deltaRTbins", "deltaRTLOESS", "RTzscore", "RTprobability", "RTprobabilityUnifPrior", "deltaRTLOESSnormalized",
            "detectFractionGreater", "detectability", "detectSubtractMissing", "detectProtSpearmanDiff",
            "deltaIMLOESS", "deltaIMLOESSnormalized", "IMprobabilityUnifPrior",
            "maxConsecutiveFragments"));
    public static final HashSet<String> detectFeatures =
            new HashSet<>(Arrays.asList("detectFractionGreater", "detectability", "detectSubtractMissing", "detectProtSpearmanDiff"));
    public static final HashSet<String> spectraRTFeatures = new HashSet<>(Arrays.asList(
            "cosineSimilarity", "weightedCosineSimilarity", "spectralContrastAngle", "weightedSpectralContrastAngle",
            "euclideanDistance", "weightedEuclideanDistance", "brayCurtis", "weightedBrayCurtis",
            "pearsonCorr", "weightedPearsonCorr", "dotProduct", "weightedDotProduct",
            "deltaRTlinear", "deltaRTbins", "deltaRTLOESS", "RTzscore", "RTprobability", "RTprobabilityUnifPrior"));
    public static final HashSet<String> spectraFeatures = new HashSet<>(Arrays.asList(
            "cosineSimilarity", "weightedCosineSimilarity", "spectralContrastAngle", "weightedSpectralContrastAngle",
            "euclideanDistance", "weightedEuclideanDistance", "brayCurtis", "weightedBrayCurtis",
            "pearsonCorr", "weightedPearsonCorr", "dotProduct", "weightedDotProduct"));
    public static final HashSet<String> rtFeatures = new HashSet<>(Arrays.asList(
            "deltaRTlinear", "deltaRTbins", "deltaRTLOESS", "RTzscore", "RTprobability", "RTprobabilityUnifPrior",
            "deltaRTLOESSnormalized"));
    public static final HashSet<String> imFeatures =
            new HashSet<>(Arrays.asList("deltaIMLOESS", "deltaIMLOESSnormalized", "IMprobabilityUnifPrior"));

    public static final HashMap<String, String> camelToUnderscore = makeCamelToUnderscore();
    private static HashMap<String, String> makeCamelToUnderscore() {
        HashMap<String, String> map = new HashMap<>();
        map.put("cosineSimilarity", "cosine_similarity");
        map.put("spectralContrastAngle", "spectral_contrast_angle");
        map.put("euclideanDistance", "euclidean_distance");
        map.put("brayCurtis", "bray_curtis");
        map.put("pearsonCorr", "pearson_corr");
        map.put("dotProduct", "dot_product");
        map.put("deltaRTLOESS", "delta_RT_loess");
        map.put("deltaRTLOESSnormalized", "delta_RT_loess_normalized");
        map.put("RTprobabilityUnifPrior", "RT_probability_unif_prior");
        map.put("deltaIMLOESS", "delta_IM_loess");
        map.put("deltaIMLOESSnormalized", "delta_IM_loess_normalized");
        map.put("IMprobabilityUnifPrior", "IM_probability_unif_prior");
        map.put("detectProtSpearmanDiff", "detect_prot_spearman_diff");
        return map;
    }
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //Handling PTMs

    //TODO: better handling of PTMs, all in one location
    private static HashMap<Double, Integer> makeModAAToUnimod() {
        HashMap<Double, Integer> map = new HashMap<>();
        map.put(57.0215, 4);
        map.put(15.9949, 35);
        map.put(42.0106, 1);
        map.put(79.96633, 21);
        map.put(114.042927, 121);
        return map;
    }
    public static final HashMap<Double, Integer> modAAmassToUnimod = makeModAAToUnimod();
    private static HashMap<String, Double> makeUnimodtoModAA() {
        HashMap<String, Double> map = new HashMap<>();
        for (Map.Entry<Double, Integer> entry : modAAmassToUnimod.entrySet()) {
            map.put(String.valueOf(entry.getValue()), entry.getKey());
        }
        return map;
    }
    public static final HashMap<String, Double> AAmassToUnimod = makeUnimodtoModAA();

    //methods
    public void updatePaths() {
        if (outputDirectory == null) {
            String firstFile = pinPepXMLDirectory.split(" ")[0];
            File newFile = new File(firstFile);
            if (newFile.isDirectory()) {
                outputDirectory = firstFile;
            } else { //file
                outputDirectory = newFile.getAbsoluteFile().getParent();
            }
        }
        if (editedPin == null || renamePin == 0) { //if 0, replace at end
            editedPin = "edited";
        }
        if (spectraRTPredInput == null) {
            spectraRTPredInput = outputDirectory + File.separator + "spectraRT.tsv";
        }
        if (detectPredInput == null) {
            detectPredInput = outputDirectory + File.separator + "detect.tsv";
        }
    }
}
