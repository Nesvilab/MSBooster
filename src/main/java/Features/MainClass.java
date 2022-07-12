package Features;

import External.ExternalModelCaller;

import java.io.*;
import java.lang.reflect.Field;
import java.lang.reflect.Modifier;
import java.util.*;

//this is what I use in the java jar file
public class MainClass {
    public static void main(String[] args) throws Exception {
        System.out.println("MSBooster v1.1.5");
        try {
            //accept command line inputs
            HashSet<String> fields = new HashSet<>();
            for (Field f : Constants.class.getDeclaredFields()) {
                fields.add(f.getName());
            }

            HashMap<String, String> params = new HashMap<String, String>();

            //setting new values
            for (int i = 0; i < args.length; i++) {
                String key = args[i].substring(2); //remove --
                if (key.equals("help")) { //help message
                    System.out.println("Usage: java -jar MSFraggerDIA_postprocess-1.0-SNAPSHOT-jar-with-dependencies.jar [flags]");
                    System.out.println("Usage: java -jar MSFraggerDIA_postprocess-1.0-SNAPSHOT-jar-with-dependencies.jar " +
                            "--paramsList {*.txt}");
                    System.out.println("A tool for annotating PSM pin files with deep learning-based features. " +
                            "This has been tested on DIA and DDA data as processed by MSFragger. " +
                            "A variety of features can be added to the provided pin file, " +
                            "including MS2 spectral, retention time (RT), ion monbility (IM), and detectability prediction-based features. " +
                            "At a minimum, an edited pin file is produced. " +
                            "This tool can also run Dia-NN for spectral and RT prediction (path to a copy of DiaNN.exe must be provided), " +
                            "and DeepMSPeptide for detectability prediction (provided by DeepMSPeptideRevised.exe). " +
                            "Predictions for all peptides from the PSMs, both target and decoy, are also saved.");

                    System.out.println("");
                    System.out.println("General flags:");
                    System.out.println("\t--paramsList: Text file containing all parameters to use. Specifying this will override any other parameters provided via command line");
                    System.out.println("\t--pinPepXMLDirectory: One or more directories where pin/pepXML files are located, " +
                            "or one or more file names. These should be separated by a space.");
                    System.out.println("\t--mzmlDirectory: One or more directories where mzML files are stored, " +
                            "or one or more file names. These shoudl be separated by a space. " +
                            "It is assumed that pin and mzML files have the same prefix (ex. sample1.pin and sample1.mzML)");
                    System.out.println("\t--outputDirectory: Directory to store all intermediate and final files (default: pinPepXMLDirectory. " +
                            "If multiple directories/files are provided, the first directory or folder of the file will be used.)");
                    System.out.println("\t--editedPin: prefix for edited pin file (default: {outputDirectory}/edited_). " +
                            "If {renamePin} is not 1, this is ignored, as the edited pin file will replace the old files.");
                    System.out.println("\t--renamePin: If 1, new pin files will be produced using the {editedPin} prefix, " +
                            "and old pin files will be kept. Otherwise, the new pin files will replace the old pin files using the same name " +
                            "(Default = 1).");
                    System.out.println("\t--spectraRTPredInput: path to the prediction input file for DIA-NN. " +
                            "If not yet produced, it will be generated as saved with this name (default: {outputDirectory}/spectraRT.tsv)");
                    System.out.println("\t--detectPredInput: path to the detectability prediction input file (default: {outputDirectory}/detect.tsv)");
                    System.out.println("\t--spectraRTPredFile: path to the spectral/RT prediction file, if produced from a previous run. " +
                            "This parameter must be provided to skip regenerating the predictions");
                    System.out.println("\t--detectPredFile: path to the detectability prediction file, if produced from a previous run.");
                    System.out.println("\t--features: features to add to the edited pin file, " +
                            "separated by commas with no white space (ex. brayCurtis,deltaRTlinear). " +
                            "Default = " + Constants.features);
                    System.out.println("\t--numThreads: number of threads available. " +
                            "numThreads <= 0 indicate for all available processors to be used (Default: 0)");
                    System.out.println("\t--fragger: path to fragger.params used for MSFragger. Parameters such as the " +
                            "fragment error tolerance are read in for matching predicted and experimental fragments. " +
                            "If the detectability feature is used, the in-silico digestion rules are reused.");
                    System.out.println("\t--DiaNN: path to DiaNN.exe. Only needed if spectral/RT/IM predictions need to be produced");
                    System.out.println("\t--useSpectra: whether or not to use spectral features (default: true)");
                    System.out.println("\t--useRT: whether or not to use RT features (default: true)");
                    System.out.println("\t--useDetect: whether or not to use detectability feature (default: false)");
                    System.out.println("\t--useIM: whether or not to use IM features (default: false). This needs to be turned to true manually" +
                            "if analyzing ion mobility data.");

                    System.out.println("");
                    System.out.println("Flags that are only used if calculating detectability features (these can all be provided using" +
                            "the --fragger flag):");
                    System.out.println("\t--fasta: path to fasta file");
                    System.out.println("\t--decoyPrefix: prefix for decoys in fasta file (default: >rev_)");
                    System.out.println("\t--cutAfter: amino acids after which the enzyme digested the peptide (default: KR)");
                    System.out.println("\t--butNotAfter: amino acids that are an exception to the enzyme digestion rules (default: P)");
                    System.out.println("\t--digestMinLength: minimum length peptide searched in datbase search tool (default: 7)");
                    System.out.println("\t--digestMaxLength: maximum length peptide searched in database search tool (default: 50)");
                    System.out.println("\t--digestMinMass: minimum mass peptide searched in database search tool in Daltons (default: 500)");
                    System.out.println("\t--digestMaxMass: maxmimum mass peptide searched in database search tool in Daltons (default: 5000)");

                    System.out.println("");
                    System.out.println("Flags that are only used if adding spectral and/or RT features to the edited pin file:");
                    System.out.println("\t--ppmTolerance: fragment error tolerance (in ppm) of the database search. " +
                            "Used for matching fragment peaks from experimental spectra to predicted spectra (default: 20)");
                    System.out.println("\t--useTopFragments: rather than using all predicted fragments, " +
                            "use only top N predicted intensity peaks (default: true)");
                    System.out.println("\t--topFragments: how many top intensity fragments to use, if useTopFragments is true (default: 12)");
                    System.out.println("\t--removeRankPeaks: if true, fragments within the ppmTolerance of matched predicted fragments" +
                            "are filtered from the experimental spectra of lower ranking peaks from the same scan number (default: true)");
                    System.out.println("\t--RTregressionSize: how many PSMs used for loess regression, sorted by lowest e score, " +
                            "are used for the linear regression of predicted and experimental RTs (default: 5000)");
                    System.out.println("\t--RTescoreCutoff: PSMs with e score above this cutoff are not included in RT linear regression modeling " +
                            "(default: 10e-3.5)");
                    System.out.println("\t--RTbinMultiplier: for generating empirical distributions of predicted RT, how much should the " +
                            "experimental RT bins be split. {RTbinMultiplier} bins per whatever time unit is used in the mzML file. " +
                            "(default: 1)");
                    System.out.println("\t--RTIQR: The interquartile range used when normalizing the deltaRTLOESS for deltaRTLOESSnormalized. " +
                            "This value is between 0 and 100, where 100 would be the range of (max predicted RT - min predicted RT) for each " +
                            "experimental RT bin " +
                            "(default: 50)");

                    System.out.println("");
                    System.out.println("Flags that are only used if adding IM features to the edited pin file:");
                    System.out.println("\t--IMregressionSize: how many PSMs used for loess regression, sorted by lowest e score, " +
                            "are used for the linear regression of predicted and experimental RTs (default: 5000)");
                    System.out.println("\t--IMescoreCutoff: PSMs with e score above this cutoff are not included in IM linear regression modeling " +
                            "(default: 10e-3.5)");
                    System.out.println("\t--IMbinMultiplier: for generating empirical distributions of predicted IM, how much should the " +
                            "experimental IM bins be split. {IMbinMultiplier} bins per 1/K0 unit. " +
                            "(default: 100)");
                    System.out.println("\t--IMIQR: The interquartile range used when normalizing the deltaIMLOESS for deltaIMLOESSnormalized. " +
                            "This value is between 0 and 100, where 100 would be the range of (max predicted 1/K0 - min predicted 1/K0) for each " +
                            "experimental IM bin " +
                            "(default: 50)");

                    System.out.println("");
                    System.out.println("Flags shared by RT and IM features");
                    System.out.println("\t--uniformPriorPercentile: for RTprobabilityUnifPrior and IMprobabilityUnifPrior, " +
                            "how much weight to give to the uniform prior. " +
                            "Range 0 to 100. If 0, this is equal to RT probability. If 100, all PSMs' scores for this feature will be the same " +
                            "(default: 10)");
                    System.out.println("\t--bandwidth: how many of the local points to use for LOESS regression. " +
                            "IM regression is more erratic, so it will use double this value (default = 0.25)");

                    System.out.println("printing example paramsList to exampleParams.txt");
                    printParams(".");
                    System.exit(0);
                }
                i++;
                StringBuilder sb = new StringBuilder(args[i]);
                if (i + 1 >= args.length) {
                    params.put(key, sb.toString());
                } else {
                    while (!args[i + 1].startsWith("--")) {
                        i++;
                        sb.append(" ");
                        sb.append(args[i]);
                        if (i + 1 >= args.length) {
                            break;
                        }
                    }
                    params.put(key, sb.toString());
                }
            }

            //adding to constants class
            if (params.containsKey("paramsList")) { //override previous params input
                //params with nulls are left out
                String line;
                BufferedReader reader = new BufferedReader(new FileReader(params.get("paramsList")));
                while ((line = reader.readLine()) != null) {
                    if (!line.contains("=")) { //maybe empty line or comment line with #
                        continue;
                    }
                    String[] lineSplit = line.split("=", 2);

                    //check if null here
                    if (!lineSplit[1].trim().equals("null")) {
                        params.put(lineSplit[0].trim(), lineSplit[1].trim());
                    }
                }
                reader.close();
            }

            if (params.containsKey("fragger")) { //upload fasta digestion params from fragger file. Does not use PTM info. Will override paramsList
                if (!params.get("fragger").equals("null")) {
                    String line;
                    boolean DaToPPM = false;
                    BufferedReader reader = new BufferedReader(new FileReader(params.get("fragger")));
                    while ((line = reader.readLine()) != null) {
                        String[] lineSplit = line.split("#")[0].split("=");
                        if (lineSplit.length != 2) {
                            continue;
                        }
                        String key = lineSplit[0].trim();
                        String val = lineSplit[1].trim();
                        switch (key) {
                            case "fragment_mass_tolerance":
                                params.put("ppmTolerance", val);
                                break;
                            case "fragment_mass_units":
                                if (val.equals("0")) {
                                    DaToPPM = true;
                                }
                                break;
                            case "decoy_prefix":
                                params.put("decoyPrefix", ">" + val);
                                break;
                            case "search_enzyme_cutafter":
                                params.put("cutAfter", val);
                                break;
                            case "search_enzyme_butnotafter":
                                params.put("butNotAfter", val);
                                break;
                            case "digest_min_length":
                                params.put("digestMinLength", val);
                                break;
                            case "digest_max_length":
                                params.put("digestMaxLength", val);
                                break;
                            case "digest_mass_range":
                                String[] vals = val.split(" ");
                                params.put("digestMinMass", vals[0]);
                                params.put("digestMaxMass", vals[1]);
                                break;
                            case "database_name":
                                params.put("fasta", val);
                                break;
                        }
                    }

                    if (DaToPPM) {
                        float tol = (float) Math.ceil(Float.parseFloat(params.get("ppmTolerance")) * 1000f);
                        if (tol >= 100f) {
                            params.put("lowResppmTolerance", String.valueOf(tol));
                        } else {
                            params.put("highResppmTolerance", String.valueOf(tol));
                        }
                    } else {
                        float tol = Float.parseFloat(params.get("ppmTolerance"));
                        if (tol >= 100f) {
                            params.put("lowResppmTolerance", String.valueOf(tol));
                        } else {
                            params.put("highResppmTolerance", String.valueOf(tol));
                        }
                    }
                }
            }

            Constants c = new Constants();
            for (Map.Entry<String, String> entry : params.entrySet()) {
                String key = entry.getKey();
                if (!fields.contains(key)) {
                    throw new Exception(entry.getKey() + " is not a valid parameter");
                } else {
                    //get class of field
                    Field field = Constants.class.getField(key);
                    Class<?> myClass = field.getType();

                    //do not parse use[something] if null
                    if (myClass.getTypeName().equals("java.lang.Boolean")) {
                        if (entry.getValue().equals("null")) {
                            continue;
                        }
                    }

                    //parse to appropriate type
                    field.set(c, myClass.getConstructor(String.class).newInstance(entry.getValue()));
                }
            }

            if (Constants.divideFragments.equals("1")) { //standard setting of yb vs others
                Constants.divideFragments = "y_b;immonium_a_y-NL_b-NL_a-NL_internal_internal-NL_unknown";
                Constants.topFragments = 12;
            } else if (Constants.divideFragments.equals("2")) {
                Constants.divideFragments = "y;b;immonium;a;y-NL;b-NL;a-NL;internal;internal-NL;unknown";
                Constants.topFragments = 6;
            }

            //defining num threads
            Runtime run = Runtime.getRuntime();
            if (Constants.numThreads <= 0) {
                Constants.numThreads = run.availableProcessors();
            } //otherwise use user-defined
            System.out.println("Using " + Constants.numThreads + " threads");

            //check that at least pinPepXMLDirectory and mzmlDirectory are provided
            if (Constants.pinPepXMLDirectory == null) {
                throw new IllegalArgumentException("pinPepXMLDirectory must be provided");
            }
            if (Constants.mzmlDirectory == null) {
                throw new IllegalArgumentException("mzmlDirectory must be provided");
            }

            //check that features are allowed
            if (Constants.useMultipleCorrelatedFeatures) {
                Constants.features = "brayCurtis,pearsonCorr,dotProduct,unweightedSpectralEntropy," +
                        "deltaRTLOESS,deltaRTLOESSnormalized,RTprobabilityUnifPrior";
            }
            String[] featuresArray = Constants.features.split(",");
            for (String f : featuresArray) {
                if (!Constants.allowedFeatures.contains(f.trim())) {
                    throw new IllegalArgumentException(f + " is not an allowed feature. " +
                            "Please choose from the following: " + Constants.allowedFeatures);
                }
            }

            LinkedList<String> featureLL = new LinkedList<>(Arrays.asList(featuresArray));

            //use "use" variables to update
            try {
                if (Constants.useSpectra) {
                    Set<String> intersection = new HashSet<>(featureLL);
                    intersection.retainAll(Constants.spectraFeatures);
                    if (intersection.size() == 0) {
                        featureLL.add("brayCurtis");
                        featureLL.add("pearsonCorr");
                        featureLL.add("dotProduct");
                    }
                } else {
                    featureLL.removeIf(Constants.spectraFeatures::contains);
                }
            } catch (Exception ignored) {}
            try {
                if (Constants.useRT) {
                    Set<String> intersection = new HashSet<>(featureLL);
                    intersection.retainAll(Constants.rtFeatures);
                    if (intersection.size() == 0) {
                        featureLL.add("deltaRTLOESS");
                        featureLL.add("deltaRTLOESSnormalized");
                        featureLL.add("RTprobabilityUnifPrior");
                    }
                } else {
                    featureLL.removeIf(Constants.rtFeatures::contains);
                }
            } catch (Exception ignored) {}
            try {
                if (Constants.useDetect) {
                    System.out.println("Detect features not fully tested");
                    Set<String> intersection = new HashSet<>(featureLL);
                    intersection.retainAll(Constants.detectFeatures);
                    if (intersection.size() == 0) {
                        featureLL.add("detectFractionGreater");
                        featureLL.add("detectSubtractMissing");
                    }
                } else {
                    featureLL.removeIf(Constants.detectFeatures::contains);
                }
            } catch (Exception ignored) {}
            try {
                if (Constants.useIM) {
                    Set<String> intersection = new HashSet<>(featureLL);
                    intersection.retainAll(Constants.imFeatures);
                    if (intersection.size() == 0) {
                        featureLL.add("deltaIMLOESS");
                        featureLL.add("deltaIMLOESSnormalized");
                        featureLL.add("IMprobabilityUnifPrior");
                    }
                } else {
                    for (String feature : Constants.imFeatures) {
                        featureLL.remove(feature);
                    }
                }
            } catch (Exception ignored) {}
            try {
                if (Constants.useMatchedIntensities) {
                    Set<String> intersection = new HashSet<>(featureLL);
                    intersection.retainAll(Constants.matchedIntensitiesFeatures);
                    if (intersection.size() == 0) {
                        featureLL.addAll(Constants.matchedIntensitiesFeatures);
                    }
                } else {
                    for (String feature : Constants.matchedIntensitiesFeatures) {
                        featureLL.remove(feature);
                    }
                }
            } catch (Exception ignored) {}
//            try {
//                if (Constants.usePredIntensities) {
//                    Set<String> intersection = new HashSet<>(featureLL);
//                    intersection.retainAll(Constants.predIntensitiesFeatures);
//                    if (intersection.size() == 0) {
//                        featureLL.addAll(Constants.predIntensitiesFeatures);
//                    }
//                } else {
//                    for (String feature : Constants.predIntensitiesFeatures) {
//                        featureLL.remove(feature);
//                    }
//                }
//            } catch (Exception ignored) {}
            try {
                if (Constants.usePeakCounts) {
                    Set<String> intersection = new HashSet<>(featureLL);
                    intersection.retainAll(Constants.peakCountsFeatures);
                    if (intersection.size() == 0) {
                        featureLL.addAll(Constants.peakCountsFeatures);
                    }
                } else {
                    for (String feature : Constants.peakCountsFeatures) {
                        featureLL.remove(feature);
                    }
                }
            } catch (Exception ignored) {}
            try {
                if (Constants.useIndividualSpectralSimilarities) {
                    Set<String> intersection = new HashSet<>(featureLL);
                    intersection.retainAll(Constants.individualSpectralSimilaritiesFeatures);
                    if (intersection.size() == 0) {
                        featureLL.addAll(Constants.individualSpectralSimilaritiesFeatures);
                    }
                } else {
                    for (String feature : Constants.individualSpectralSimilaritiesFeatures) {
                        featureLL.remove(feature);
                    }
                }
            } catch (Exception ignored) {}
            try {
                if (Constants.useIntensitiesDifference) {
                    Set<String> intersection = new HashSet<>(featureLL);
                    intersection.retainAll(Constants.intensitiesDifferenceFeatures);
                    if (intersection.size() == 0) {
                        featureLL.addAll(Constants.intensitiesDifferenceFeatures);
                    }
                } else {
                    for (String feature : Constants.intensitiesDifferenceFeatures) {
                        featureLL.remove(feature);
                    }
                }
            } catch (Exception ignored) {}
            try {
                if (Constants.useIntensityDistributionSimilarity) {
                    if (! featureLL.contains("intensity_distribution_similarity")) {
                        featureLL.add("intensity_distribution_similarity");
                    }
                } else {
                    featureLL.remove("intensity_distribution_similarity");
                }
            } catch (Exception ignored) {}

            //update features representation
            featuresArray = new String[featureLL.size()];
            int i = 0;
            for (String feature : featureLL) {
                featuresArray[i] = feature;
                i++;
            }
            Constants.features = String.join(",", featuresArray);

            //if detectFractionGreater, need fasta
            if (featureLL.contains("detectFractionGreater") || featureLL.contains("detectSubtractMissing") ||
                    featureLL.contains("detectProtSpearmanDiff")) {
                if (Constants.fasta == null) {
                    throw new IllegalArgumentException("Using current combination of features, " +
                            "detectFractionGreater is calculated and needs a fasta provided using " +
                            "--fasta <fasta file location>");
                }
            }

            //create file for spectral and RT prediction
            //ignore if files already created
            boolean createSpectraRTPredFile = false;
            boolean createDetectPredFile = false;
            boolean createSpectraRTPredFile2 = false;
            boolean createDetectPredFile2 = false;

            //check which ones we need
            if (featureLL.size() > 0) {
                createSpectraRTPredFile = true;
                createSpectraRTPredFile2 = true;
            }

            featureLL = new LinkedList<>(Arrays.asList(featuresArray));
            featureLL.retainAll(Constants.detectFeatures);
            if (featureLL.size() > 0) {
                createDetectPredFile = true;
                createDetectPredFile2 = true;
            }

            //overriding if intermediate files already made
            if (Constants.spectraRTPredInput != null || Constants.spectraRTPredFile != null) {
                createSpectraRTPredFile = false;
            }
            if ((Constants.spectraRTPredModel.contains("DIA-NN")) && !(Constants.spectraRTPredModel.equals("DIA-NN"))) {
                createSpectraRTPredFile = true;
            }
            if (Constants.detectPredInput != null || Constants.detectPredFile != null) {
                createDetectPredFile = false;
            }

            c.updatePaths(); //setting null paths

            //generate files for prediction models

            //get matched pin files for mzML files
            PinMzmlMatcher pmMatcher = new PinMzmlMatcher(Constants.mzmlDirectory, Constants.pinPepXMLDirectory);
            //TODO: switch case
            if (createSpectraRTPredFile || Constants.createPredFileOnly) {
                if (Constants.spectraRTPredModel.equals("DIA-NN")) {
                    if (Constants.DiaNN == null) {
                        throw new IllegalArgumentException("path to DIA-NN executable must be provided");
                    }
                    System.out.println("Generating input file for DIA-NN");
                    peptideFileCreator.createPeptideFile(pmMatcher.pinFiles, Constants.spectraRTPredInput, "Diann");
                    peptideFileCreator.createPeptideFile(pmMatcher.pinFiles,
                            Constants.spectraRTPredInput.substring(0, Constants.spectraRTPredInput.length() - 4) + "_full.tsv",
                            "createFull");
                } else if (Constants.spectraRTPredModel.equals("pDeep3")) {
                    System.out.println("Generating input file for pDeep3");
                    peptideFileCreator.createPeptideFile(pmMatcher.pinFiles, Constants.spectraRTPredInput, "pDeep3");
                } else if (Constants.spectraRTPredModel.equals("PredFull")) {
                    System.out.println("Generating input file for PredFull");
                    peptideFileCreator.createPeptideFile(pmMatcher.pinFiles, Constants.spectraRTPredInput, "PredFull");
                    peptideFileCreator.createPeptideFile(pmMatcher.pinFiles,
                            Constants.spectraRTPredInput.substring(0, Constants.spectraRTPredInput.length() - 4) + "_full.tsv",
                            "createFull");
                } else if (Constants.spectraRTPredModel.equals("Prosit")) {
                    System.out.println("Generating input file for Prosit");
                    peptideFileCreator.createPeptideFile(pmMatcher.pinFiles, Constants.spectraRTPredInput, "Prosit");
                    peptideFileCreator.createPeptideFile(pmMatcher.pinFiles,
                            Constants.spectraRTPredInput.substring(0, Constants.spectraRTPredInput.length() - 4) + "_full.tsv",
                            "createFull");
                } else if (Constants.spectraRTPredModel.equals("PrositTMT")) {
                    System.out.println("Generating input file for PrositTMT");
                    peptideFileCreator.createPeptideFile(pmMatcher.pinFiles, Constants.spectraRTPredInput, "PrositTMT");
                    peptideFileCreator.createPeptideFile(pmMatcher.pinFiles,
                            Constants.spectraRTPredInput.substring(0, Constants.spectraRTPredInput.length() - 4) + "_full.tsv",
                            "createFull");
                } else if (Constants.spectraRTPredModel.equals("DIA-NN,PredFull")) {
                    if (Constants.DiaNN == null) {
                        throw new IllegalArgumentException("path to DIA-NN executable must be provided");
                    }
                    System.out.println("Generating input file for DIA-NN");
                    peptideFileCreator.createPeptideFile(pmMatcher.pinFiles, Constants.spectraRTPredInput, "Diann");
                } else {
                    System.out.println("spectraRTPredModel must be one of DIA-NN, Prosit, or PredFull");
                    System.exit(-1);
                }

                if (Constants.createPredFileOnly) {
                    System.out.println("Successfully created input file for prediction model. Stopping here");
                    System.exit(0);
                }
            }
//            if (createDetectAllPredFile) {
//                System.out.println("Generating input file for DeepMSPeptide");
//                //long startTime = System.nanoTime();
//                //Constants.setFastaReader(peptideFileCreator.createPeptideFile(pmMatcher.pinFiles, Constants.detectPredInput, "DeepMSPeptideAll", "pin"));
//                //long endTime = System.nanoTime();
//                //long duration = (endTime - startTime);
//            } else if (createDetectPredFile) {
//                System.out.println("Generating input file for DeepMSPeptide");
//                peptideFileCreator.createPeptideFile(pmMatcher.pinFiles, Constants.detectPredInput, "DeepMSPeptide", "pin");
//            }

            //generate predictions
            if ((Constants.spectraRTPredFile == null) && (createSpectraRTPredFile2)) {
                ExternalModelCaller.callModel(run, "DIA-NN");
            }

            //create new pin file with features
            System.out.println("Generating edited pin with following features: " + Arrays.toString(featuresArray));
            long start = System.nanoTime();
            if (Constants.spectraRTPredModel.contains("PredFull")) {
                Constants.matchWithDaltons = true; //they report predictions in bins
            }
            percolatorFormatter.editPin(pmMatcher, Constants.spectraRTPredFile, Constants.detectPredFile, featuresArray, Constants.editedPin);

            //print parameters to ps
            printParamsPS();

            //delete pred files
            if (Constants.deletePreds) {
                File predFile = new File(Constants.outputDirectory + File.separator + "spectraRT.tsv");
                predFile.delete();
                predFile = new File(Constants.outputDirectory + File.separator + "spectraRT_full.tsv");
                predFile.delete();
                predFile = new File(Constants.outputDirectory + File.separator + "spectraRT.predicted.bin");
                predFile.delete();
//                predFile = new File(Constants.outputDirectory + File.separator + "detect.tsv");
//                predFile.delete();
//                predFile = new File(Constants.outputDirectory + File.separator + "detect_Predictions.txt");
//                predFile.delete();
            }

            long end = System.nanoTime();
            long duration = (end - start);
            System.out.println("Done in " + duration / 1000000 + " ms");
        } catch (Exception e) {
            e.printStackTrace();
            System.exit(1);
        }
    }

    static private void printParams(String directory) {
        try {
            Constants c = new Constants();
            BufferedWriter buffer = new BufferedWriter(new FileWriter(directory + File.separator + "exampleParams.txt"));

            Field[] f = Constants.class.getFields();
            for (Field field : f) {
                if ((field.getModifiers() & Modifier.FINAL) != Modifier.FINAL) {
                    if (!field.getName().equals("paramsList")) {
                        buffer.write(field.getName() + " = " + field.get(c) + "\n");
                    }
                }
            }
            buffer.close();
        } catch (Exception e) {
            System.out.println("could not write final params");
            e.getStackTrace();
            System.exit(1);
        }
    }

    static private void printParamsPS() {
        try {
            System.out.println("Final parameters used for feature annotation:");
            Constants c = new Constants();

            Field[] f = Constants.class.getFields();
            for (Field field : f) {
                if ((field.getModifiers() & Modifier.FINAL) != Modifier.FINAL) {
                    if (!field.getName().equals("paramsList")) {
                        if (field.getName().equals("fragmentIonHierarchy")) {
                            System.out.println("\t" + field.getName() + " = " + Arrays.toString((String[]) field.get(c)));
                        } else {
                            System.out.println("\t" + field.getName() + " = " + field.get(c));
                        }
                    }
                }
            }
        } catch (Exception e) {
            System.out.println("could not write final params");
            e.getStackTrace();
            System.exit(1);
        }
    }
}
