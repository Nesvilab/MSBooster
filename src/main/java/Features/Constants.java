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
    public static Boolean deletePreds = false;

    //optional file locations and parameters
    //if calculating detectFractionGreater, these are used for FastaReader class
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
    public static String DiaNN = null;
    public static String spectraRTPredModel = "DIA-NN";
    public static Boolean replaceYBintensities = true;
    public static Integer splitPredInputFile = 1;

    //additional modifications for alphapeptdeep
    public static String additionalMods = ""; //this used in python script for common/user_defined_modifications
    public static Boolean predict = true;
    public static Boolean transfer = false;
    public static String AlphaPeptDeep = null; //path to exe
    public static String yaml = ""; //this used in python script for parameters
    public static double MS2escore = -3.5;

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //miscellaneous inner workings
    //convert int flag to fragment ion type
    private static HashMap<Integer, String> makeFlagTOion() {
        HashMap<Integer, String> map = new HashMap<>();
        map.put(0, "b");
        map.put(1, "y");
        map.put(2, "b-NL");
        map.put(3, "y-NL");
        return map;
    }
    public static HashMap<Integer, String> flagTOion = makeFlagTOion();

    private static HashMap<String, Integer> makeIonToFlag() {
        HashMap<String, Integer> map = new HashMap<>();
        map.put("b", 0);
        map.put("y", 1);
        map.put("b-NL", 2);
        map.put("y-NL", 3);
        map.put("c", 4);
        map.put("z", 5);
        return map;
    }
    public static HashMap<String, Integer> ionTOflag = makeIonToFlag();

    private static HashMap<String, Integer> makeModelMaxIntensity() {
        HashMap<String, Integer> map = new HashMap<>();
        map.put("DIA-NN", 60000);
        map.put("Prosit", 1);
        map.put("alphapeptdeep", 1);
        map.put("PredFull", 1000);
        return map;
    }
    public static HashMap<String, Integer> modelMaxIntensity = makeModelMaxIntensity();

    //these two constants are for weighted spectral similarity features, not currently supported
    public static final Integer binwidth = 1;
    public static final Integer mzFreqWindow = 1;

    public static Float ppmTolerance = 20f; //ppm tolerance of MS2 scans
    public static Float lowResppmTolerance = 300f;
    public static Float highResppmTolerance = 20f;
    public static Boolean matchWithDaltons = false;
    public static Float DaTolerance = 0.05f;

    //for limiting number of fragments used
    public static Boolean useSpectra = null;
    public static Boolean useTopFragments = true;
    public static Integer topFragments = 20;
    public static Boolean removeRankPeaks = true; //whether to remove peaks from higher ranks
    public static Boolean useBasePeak = true;
    public static Float percentBasePeak = 1f;
    //public static Boolean sqrtPredIntensities = false;
    //public static Float percentBasePeakExperimental = 1f;

    public static final Integer fineTuneSize = 100; //for generating a finetune file for pDeep3

    //these constants for RT features
    public static Boolean useRT = null;
    public static Integer RTregressionSize = 5000;
    public static Double uniformPriorPercentile = 10d;
    public static Float RTescoreCutoff = -3.5f; //PSMs with e score higher than this won't make it into RT regression modeling
    public static Integer RTbinMultiplier = 1;
    public static Float RTIQR = 50f;
    public static boolean noRTscores = false; //TODO: better handling of this

    //LOESS
    public static Double bandwidth = 0.05;
    public static Integer robustIters = 2;

    //detect
    public static Boolean useDetect = false;
    public static final Float detectThreshold = 0.0000002f; //for detectability filtering
    public static final Float detectFractionGreaterNumerator = 1f;
    public static final Float detectFractionGreaterDenominator = 2f; //prior

    //ion mobility
    public static Boolean useIM = null;
    public static Integer IMregressionSize = 5000;
    public static Float IMescoreCutoff = -3.5f;
    public static Integer IMbinMultiplier = 100;
    public static final Float IMIQR = 50f;

    //support for PredFull and Prosit
    public static String FragmentationType = "";
    public static String NCE = "30";
    public static String instrument = "";
    public static Integer maxPredictedFragmentCharge = 100;
    public static Integer minPredictedFragmentNum = 0;
    public static Boolean createPredFileOnly = false;
    public static String ignoredFragmentIonTypes = ""; //split with commas
    public static String onlyFragmentIonTypes = ""; //split with commas
    public static Set<String> makeIgnoredFragmentIonTypes() {
        Set<String> ignoredFragmentIonTypes = new HashSet<>();
        Set<String> onlyFragmentIonTypes = new HashSet<>();
        if (! Constants.onlyFragmentIonTypes.equals("")) {
            String[] commaSplit = Constants.onlyFragmentIonTypes.split(",");
            for (int i = 0; i < commaSplit.length; i++) {
                String fragmentIonType = commaSplit[i].trim();
                if (MassCalculator.allowedFragmentIonTypes.contains(fragmentIonType)) {
                    onlyFragmentIonTypes.add(fragmentIonType);
                } else {
                    System.out.println(fragmentIonType + " is not a supported fragment ion type to include. " +
                            "Please choose from " + MassCalculator.allowedFragmentIonTypes);
                    System.exit(-1);
                }
            }
            for (String fragment : MassCalculator.allowedFragmentIonTypes) {
                if (! onlyFragmentIonTypes.contains(fragment)) {
                    ignoredFragmentIonTypes.add(fragment);
                }
            }
        } else if (! Constants.ignoredFragmentIonTypes.equals("")) {
            //only filter if not excluding certain fragment ion types
            //check that this is allowed
            String[] commaSplit = Constants.ignoredFragmentIonTypes.split(",");
            for (int i = 0; i < commaSplit.length; i++) {
                String fragmentIonType = commaSplit[i].trim();
                if (MassCalculator.allowedFragmentIonTypes.contains(fragmentIonType)) {
                    ignoredFragmentIonTypes.add(fragmentIonType);
                } else {
                    System.out.println(fragmentIonType + " is not a supported fragment ion type to exclude. " +
                            "Please choose from " + MassCalculator.allowedFragmentIonTypes);
                    System.exit(-1);
                }
            }
        }
        return ignoredFragmentIonTypes;
    }
    public static String[] fragmentIonHierarchy = makeFragmentIonHierarchy();
    public static String[] makeFragmentIonHierarchy() {
        switch (Constants.FragmentationType) {
            case "HCD":
                return new String[]{"immonium", "y", "b", "a",
                        "y-NL", "b-NL", "a-NL", "internal", "internal-NL", "unknown"};
            case "ETD":
                return new String[]{"zdot", "c", "z", "y", "unknown"};
            case "ETHCD":
                return new String[]{"immonium", "y", "b", "a", "zdot", "c", "z", "cdot",
                        "y-NL", "b-NL", "a-NL", "internal", "internal-NL", "unknown"};
            default:  //everything else, like CID
                return new String[]{"immonium", "y", "b", "a",
                        "y-NL", "b-NL", "a-NL", "internal", "internal-NL", "unknown"};
        }
    }
    public static Set<String> lowestFragmentIonType = makeLowestFragmentIonType();
    public static Set<String> makeLowestFragmentIonType() {
        Set<String> ignoredFragmentIonTypesSet = makeIgnoredFragmentIonTypes();
        int index = 0;
        for (int i = fragmentIonHierarchy.length - 1; i > -1; i--) {
            String ion = fragmentIonHierarchy[i];
            if (! ignoredFragmentIonTypesSet.contains(ion)) {
                index = i;
                break;
            }
        }
        return new HashSet<>(Arrays.asList(fragmentIonHierarchy).subList(0, index + 1));
    }
    public static String divideFragments = "0";

    //PredFull fragment ion annotation
    //TODO: can only use for PredFull
    public static Boolean useMatchedIntensities = false;
    public static Boolean usePredIntensities = false;
    public static Boolean usePeakCounts = false;
    public static Boolean useIndividualSpectralSimilarities = false;
    public static Boolean useIntensitiesDifference = false;
    public static Boolean useIntensityDistributionSimilarity = false;

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    //use single string sep by comma delimiter
    //should include parameter to calculate correlation and then choose
    //default auto, everything, or all? Or a combination I figure out empirically
    public static String features = "unweightedSpectralEntropy,deltaRTLOESS";
    public static Boolean useMultipleCorrelatedFeatures = false;
    //public static String features = "auto";

    //don't currently support weighted similarity features
    public static final HashSet<String> detectFeatures =
            new HashSet<>(Arrays.asList("detectFractionGreater", "detectability", "detectSubtractMissing", "detectProtSpearmanDiff"));
    public static final HashSet<String> spectraRTFeatures = new HashSet<>(Arrays.asList(
            "cosineSimilarity", "weightedCosineSimilarity", "spectralContrastAngle", "weightedSpectralContrastAngle",
            "euclideanDistance", "weightedEuclideanDistance", "brayCurtis", "weightedBrayCurtis",
            "pearsonCorr", "weightedPearsonCorr", "dotProduct", "weightedDotProduct", "unweightedSpectralEntropy",
            "deltaRTlinear", "deltaRTbins", "deltaRTLOESS", "RTzscore", "RTprobability", "RTprobabilityUnifPrior",
            "deltaRTLOESSnormalized", "calibratedRT", "predictedRT", "numMatchedFragments"));
    public static final HashSet<String> spectraFeatures = new HashSet<>(Arrays.asList(
            "cosineSimilarity", "weightedCosineSimilarity", "spectralContrastAngle", "weightedSpectralContrastAngle",
            "euclideanDistance", "weightedEuclideanDistance", "brayCurtis", "weightedBrayCurtis", "unweightedSpectralEntropy",
            "pearsonCorr", "weightedPearsonCorr", "dotProduct", "weightedDotProduct", "numMatchedFragments"));
    public static final HashSet<String> rtFeatures = new HashSet<>(Arrays.asList(
            "deltaRTlinear", "deltaRTbins", "deltaRTLOESS", "RTzscore", "RTprobability", "RTprobabilityUnifPrior",
            "deltaRTLOESSnormalized", "calibratedRT", "predictedRT"));
    public static final HashSet<String> imFeatures =
            new HashSet<>(Arrays.asList("deltaIMLOESS", "deltaIMLOESSnormalized", "IMprobabilityUnifPrior", "predictedIM",
                    "ionmobility"));
    //TODO: add to features list
    public static HashSet<String> matchedIntensitiesFeatures = null;
    public static HashSet<String> makeMatchedIntensitiesFeatures() {
        HashSet<String> set = new HashSet<>();
        for (String s : Constants.fragmentIonHierarchy) {
            set.add(s + "_matched_intensity");
        }
        return set;
    }
    public static HashSet<String> peakCountsFeatures = null;
    public static HashSet<String> makePeakCountsFeatures() {
        HashSet<String> set = new HashSet<>();
        for (String s : Constants.fragmentIonHierarchy) {
            set.add(s + "_peak_counts");
        }
        return set;
    }
    public static HashSet<String> predIntensitiesFeatures = null;
    public static HashSet<String> makePredIntensitiesFeatures() {
        HashSet<String> set = new HashSet<>();
        for (String s : Constants.fragmentIonHierarchy) {
            set.add(s + "_pred_intensity");
        }
        return set;
    }
    public static HashSet<String> individualSpectralSimilaritiesFeatures = null;
    public static HashSet<String> makeIndividualSpectralSimilarities() {
        HashSet<String> set = new HashSet<>();
        for (String s : Constants.fragmentIonHierarchy) {
            set.add(s + "_spectral_similarity");
        }
        return set;
    }
    public static HashSet<String> intensitiesDifferenceFeatures = null;
    public static HashSet<String> makeintensitiesDifference() {
        HashSet<String> set = new HashSet<>();
        for (String s : Constants.fragmentIonHierarchy) {
            set.add(s + "_intensities_difference");
        }
        return set;
    }
    public static HashSet<String> allowedFeatures = null;
    public static HashSet<String> makeAllowedFeatures() {
        HashSet<String> hs = new HashSet<>();
        hs.add("");
        hs.addAll(spectraRTFeatures);
        hs.addAll(imFeatures);
        hs.addAll(matchedIntensitiesFeatures);
        hs.addAll(predIntensitiesFeatures);
        hs.addAll(peakCountsFeatures);
        hs.addAll(individualSpectralSimilaritiesFeatures);
        hs.addAll(intensitiesDifferenceFeatures);
        hs.add("intensity_distribution_similarity");
        return hs;
    }


    //TODO: got lazy with naming, remove this and adjust code
    public static final HashMap<String, String> camelToUnderscore = makeCamelToUnderscore();
    private static HashMap<String, String> makeCamelToUnderscore() {
        HashMap<String, String> map = new HashMap<>();
        map.put("cosineSimilarity", "cosine_similarity");
        map.put("spectralContrastAngle", "spectral_contrast_angle");
        map.put("euclideanDistance", "euclidean_distance");
        map.put("brayCurtis", "bray_curtis");
        map.put("pearsonCorr", "pearson_corr");
        map.put("dotProduct", "dot_product");
        map.put("unweightedSpectralEntropy", "unweighted_spectral_entropy");
        map.put("numMatchedFragments", "num_matched_fragments");
        map.put("deltaRTLOESS", "delta_RT_loess");
        map.put("deltaRTLOESSnormalized", "delta_RT_loess_normalized");
        map.put("RTprobabilityUnifPrior", "RT_probability_unif_prior");
        map.put("calibratedRT", "calibrated_RT");
        map.put("predictedRT", "predicted_RT");
        map.put("deltaIMLOESS", "delta_IM_loess");
        map.put("deltaIMLOESSnormalized", "delta_IM_loess_normalized");
        map.put("IMprobabilityUnifPrior", "IM_probability_unif_prior");
        map.put("predictedIM", "predicted_IM");
        map.put("detectProtSpearmanDiff", "detect_prot_spearman_diff");
        map.put("detectSubtractMissing", "detect_subtract_missing");
        return map;
    }
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

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
        if (editedPin == null) { //if 0, replace at end
            editedPin = "edited";
        }
        if (spectraRTPredInput == null) {
            if (Constants.spectraRTPredModel.contains("Prosit") ||
                    Constants.spectraRTPredModel.contains("alphapeptdeep")) {
                spectraRTPredInput = outputDirectory + File.separator + "spectraRT.csv";
            } else {
                spectraRTPredInput = outputDirectory + File.separator + "spectraRT.tsv";
            }
        }
        if (detectPredInput == null) {
            detectPredInput = outputDirectory + File.separator + "detect.tsv";
        }
    }
}
