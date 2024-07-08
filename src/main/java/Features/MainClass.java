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

import static utils.Print.printError;
import static utils.Print.printInfo;

import External.DiannModelCaller;
import External.KoinaMethods;
import External.KoinaModelCaller;
import External.NCEcalibrator;
import kotlin.jvm.functions.Function1;

import java.io.*;
import java.lang.reflect.Field;
import java.lang.reflect.Modifier;
import java.nio.file.Paths;
import java.nio.file.Files;
import java.util.concurrent.atomic.AtomicBoolean;
import java.util.stream.Collectors;
import java.util.*;
import java.util.concurrent.ScheduledThreadPoolExecutor;

//this is what I use in the java jar file
public class MainClass {
    public static ScheduledThreadPoolExecutor executorService;
    public static void main(String[] args) throws Exception {
        Locale.setDefault(Locale.US);
        printInfo("MSBooster v1.2.35");

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

                    System.out.println();
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

                    System.out.println();
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

                    System.out.println();
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

                    System.out.println();
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

                    System.out.println();
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
                            case "precursor_charge":
                                vals = val.split(" ");
                                params.put("minPrecursorCharge",
                                        String.valueOf(Math.min(Integer.parseInt(vals[0]),
                                                Constants.minPrecursorCharge)));
                                params.put("maxPrecursorCharge",
                                        String.valueOf(Math.max(Integer.parseInt(vals[1]),
                                                Constants.maxPrecursorCharge)));
                                break;
                            case "mass_offsets":
                                if (!val.isEmpty()) {
                                    vals = val.split("/");
                                    String final_vals = "";
                                    for (String v : vals) {
                                        float fv = Float.parseFloat(v);
                                        if (fv != 0) {
                                            if (!final_vals.isEmpty()) {
                                                final_vals += "/";
                                            }
                                            final_vals += String.format("%.4f", fv);
                                        }
                                    }
                                    params.put("massOffsets", final_vals);
                                }
                                break;
                            case "mass_offsets_detailed":
                                if (!val.isEmpty()) {
                                    vals = val.split(";");
                                    String final_vals = "";
                                    for (String v : vals) {
                                        v = v.split("\\(")[0];
                                        float fv = Float.parseFloat(v);
                                        if (fv != 0) {
                                            if (!final_vals.isEmpty()) {
                                                final_vals += "/";
                                            }
                                            final_vals += v;
                                        }
                                    }
                                    params.put("massOffsetsDetailed", final_vals);
                                }
                                break;
                            case "mass_diff_to_variable_mod":
                                params.put("massDiffToVariableMod", val);
                                break;
                        }
                    }

                    float tol;
                    if (DaToPPM) {
                        tol = (float) Math.ceil(Float.parseFloat(params.get("ppmTolerance")) * 1000f);
                    } else {
                        tol = Float.parseFloat(params.get("ppmTolerance"));
                    }
                    if (tol >= 100f) {
                        params.put("lowResppmTolerance", String.valueOf(tol));
                    } else {
                        params.put("highResppmTolerance", String.valueOf(tol));
                    }
                    Constants.ppmTolerance = Constants.highResppmTolerance;
                }
            }

            Constants c = new Constants();
            for (Map.Entry<String, String> entry : params.entrySet()) {
                String key = entry.getKey();
                if (key.charAt(0) == '#') { //comment
                    continue;
                }
                if (key.charAt(0) == '/') { //comment
                    continue;
                }
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
            c.updateOutputDirectory();

            //checking constants are appropriate
            //robust to if url has / or not
            if (!Constants.KoinaURL.isEmpty() &&
                    Constants.KoinaURL.charAt(Constants.KoinaURL.length() - 1) != '/') {
                Constants.KoinaURL += "/";
            }
            //checking plot output is fine extension
            HashSet<String> extensions = new HashSet<>(Arrays.asList("png", "pdf"));
            if (!extensions.contains(Constants.plotExtension.toLowerCase())) {
                printError("plotExtension must be one of png or pdf, not " + Constants.plotExtension);
                printError("Exiting");
                System.exit(1);
            }

            //defining num threads
            if (Constants.numThreads <= 0) {
                Runtime run = Runtime.getRuntime();
                Constants.numThreads = run.availableProcessors() - 1;
            } //otherwise use user-defined
            printInfo("Using " + Constants.numThreads + " threads");
            executorService = new ScheduledThreadPoolExecutor(Constants.numThreads);

            if (!Objects.equals(Constants.massDiffToVariableMod, "0")) {
                if (!Constants.massesForLoessCalibration.isEmpty()) {
                    Constants.massesForLoessCalibration += ",";
                }
                if (!Constants.massOffsetsDetailed.isEmpty()) {
                    Constants.massesForLoessCalibration += Constants.massOffsetsDetailed;
                } else if (!Constants.massOffsets.isEmpty()) {
                    Constants.massesForLoessCalibration += Constants.massOffsets;
                }
            }

            //update fragment ion types based on fragmentation type
            //update Constants to be null initially
            Constants.fragmentIonHierarchy = Constants.makeFragmentIonHierarchy();
            Constants.lowestFragmentIonType = Constants.makeLowestFragmentIonType();
            Constants.matchedIntensitiesFeatures = Constants.makeMatchedIntensitiesFeatures();
            Constants.peakCountsFeatures = Constants.makePeakCountsFeatures();
            Constants.predIntensitiesFeatures = Constants.makePredIntensitiesFeatures();
            Constants.individualSpectralSimilaritiesFeatures = Constants.makeIndividualSpectralSimilarities();
            Constants.intensitiesDifferenceFeatures = Constants.makeintensitiesDifference();

            //get matched pin files for mzML files
            PinMzmlMatcher pmMatcher = new PinMzmlMatcher(Constants.mzmlDirectory, Constants.pinPepXMLDirectory);
            //if no pin files, continue
            if (pmMatcher.pinFiles.length == 0) {
                printInfo("No pin files to process. Continuing without MSBooster.");
                executorService.shutdown();
                System.exit(0);
            }

            //make models properly uppercased, or throw error if not right
            HashMap<String, String> modelMapper = LowercaseModelMapper.lowercaseToModel;
            if (modelMapper.containsKey(Constants.spectraModel.toLowerCase())) {
                Constants.spectraModel = modelMapper.get(Constants.spectraModel.toLowerCase());
            } else {
                printError("No spectra model called " + Constants.spectraModel + ". Exiting.");
                System.exit(0);
            }
            if (modelMapper.containsKey(Constants.rtModel.toLowerCase())) {
                Constants.rtModel = modelMapper.get(Constants.rtModel.toLowerCase());
            } else {
                printError("No RT model called " + Constants.rtModel + ". Exiting.");
                System.exit(0);
            }

            //needed for nce calibration and best model search
            //TODO: could also make a new koinamethods object for decoys
            KoinaMethods km = new KoinaMethods(pmMatcher);
            boolean TMT = false;
            if (Constants.findBestRtModel || Constants.findBestSpectraModel ||
            Constants.KoinaMS2models.contains(Constants.spectraModel) ||
                    Constants.KoinaRTmodels.contains(Constants.rtModel)) {
                km.getTopPeptides();
                //km.getDecoyPeptides();
                for (String pep : km.peptideSet) {
                    if (pep.contains(String.valueOf(PTMhandler.tmtMass))) {
                        TMT = true;
                        break;
                    }
                }
            }

            //if different RT and spectra models
            Constants.spectraRTPredModel = "";
            HashMap<String, float[]> datapointsRT = new HashMap<>();
            HashMap<String, float[]> datapointsRTDecoys = new HashMap<>();
            if (Constants.useRT) {
                //here, look for best rt model
                if (Constants.findBestRtModel) {
                    printInfo("Searching for best RT model for your data");
                    ArrayList<String> consideredModels = new ArrayList<>();
                    if (TMT) {
                        consideredModels.addAll(Constants.rtSearchModelsTMT);
                    } else {
                        String[] rtSearchModels = Constants.rtSearchModelsString.split(",");
                        for (String model : rtSearchModels) {
                            if (modelMapper.containsKey(model.toLowerCase())) {
                                consideredModels.add(modelMapper.get(model.toLowerCase()));
                            } else {
                                printError("No RT model called " + model + ". Exiting.");
                                System.exit(0);
                            }
                        }
                    }
                    printInfo("Searching the following models: ");
                    printInfo(String.valueOf(consideredModels));

                    String jsonOutFolder = Constants.outputDirectory + File.separator + "best_model";
                    MyFileUtils.createWholeDirectory(jsonOutFolder);

                    for (String model : consideredModels) {
                        PredictionEntryHashMap allPreds = null;
                        if (model.equals("DIA-NN")) { //mode for DIA-NN
                            MyFileUtils.createWholeDirectory(jsonOutFolder + File.separator + model);

                            km.writeFullPeptideFile(jsonOutFolder + File.separator + model + File.separator +
                                    "spectraRT_full.tsv", model, km.peptideSet);
                            String inputFile = jsonOutFolder + File.separator + model + File.separator + "spectraRT.tsv";
                            FileWriter myWriter = new FileWriter(inputFile);
                            myWriter.write("peptide" + "\t" + "charge\n");
                            HashSet<String> peptides = new HashSet<>();
                            for (String pep : km.peptideSet) {
                                String[] pepSplit = pep.split(",");
                                String line = pepSplit[1] + "\t" + pepSplit[0].split("\\|")[1] + "\n";
                                if (!peptides.contains(line)) {
                                    myWriter.write(line);
                                    peptides.add(line);
                                }
                            }
                            myWriter.close();

                            DiannModelCaller.callModel(inputFile, false);
                            allPreds = SpectralPredictionMapper.createSpectralPredictionMapper(
                                    Constants.spectraRTPredFile, "DIA-NN", executorService).getPreds();
                            Constants.spectraRTPredFile = null;
                        } else { //mode for koina
                            HashSet<String> allHits = km.writeFullPeptideFile(
                                    jsonOutFolder + File.separator + model + "_full.tsv", model,
                                    km.peptideSet);
                            allPreds = km.getKoinaPredictions(allHits, model, 30,
                                    jsonOutFolder + File.separator + model,
                                    jsonOutFolder + File.separator + model + "_full.tsv");
                        }

                        //collect the data for testing
                        ArrayList<Float> expRTs = new ArrayList<>();
                        ArrayList<Float> predRTs = new ArrayList<>();
                        for (int j = 0; j < pmMatcher.mzmlReaders.length; j++) {
                            MzmlReader mzmlReader = pmMatcher.mzmlReaders[j];
                            LinkedList<Integer> thisScanNums = km.scanNums.get(pmMatcher.mzmlFiles[j].getName());
                            LinkedList<String> thisPeptides = km.peptides.get(pmMatcher.mzmlFiles[j].getName());

                            for (int k = 0; k < thisScanNums.size(); k++) {
                                try {
                                    int scanNum = thisScanNums.get(k);
                                    String peptide = thisPeptides.get(k).split(",")[0];
                                    MzmlScanNumber msn = mzmlReader.getScanNumObject(scanNum);

                                    float predRT = allPreds.get(peptide).RT;
                                    float expRT = msn.RT;
                                    predRTs.add(predRT);
                                    expRTs.add(expRT);
                                } catch (Exception e) {}
                            }
                        }
                        double[][] rts = new double[2][expRTs.size()];
                        for (int i = 0; i < expRTs.size(); i++) {
                            //put pred first so it can be normalized to exp scale
                            rts[0][i] = Math.round(predRTs.get(i) * 100d) / 100d;
                            rts[1][i] = Math.round(expRTs.get(i) * 100d) / 100d;
                        }

                        //testing
                        String[] bandwidths = Constants.loessBandwidth.split(",");
                        float[] floatBandwidths = new float[bandwidths.length];
                        for (int i = 0; i < bandwidths.length; i++) {
                            floatBandwidths[i] = Float.parseFloat(bandwidths[i]);
                        }
                        Object[] bandwidth_loess_rtdiffs = StatMethods.gridSearchCV(rts, floatBandwidths);
                        Function1<Double, Double> loess = (Function1<Double, Double>) bandwidth_loess_rtdiffs[1];
                        float[] rtdiffs = (float[]) bandwidth_loess_rtdiffs[2];
                        float mse = StatMethods.meanSquaredError(rtdiffs);
                        printInfo(model + " has root mean squared error of " +
                                String.format("%.4f", Math.sqrt(mse)));
                        Arrays.sort(rtdiffs);
                        datapointsRT.put(model, rtdiffs);

//                        //decoy
//                        if (model.equals("DIA-NN")) { //mode for DIA-NN
//                            MyFileUtils.createWholeDirectory(jsonOutFolder + File.separator + model);
//
//                            km.writeFullPeptideFile(jsonOutFolder + File.separator + model + File.separator +
//                                    "spectraRT_full.tsv", model, km.peptideSetDecoys);
//                            String inputFile = jsonOutFolder + File.separator + model + File.separator + "spectraRT.tsv";
//                            FileWriter myWriter = new FileWriter(inputFile);
//                            myWriter.write("peptide" + "\t" + "charge\n");
//                            HashSet<String> peptides = new HashSet<>();
//                            for (String pep : km.peptideSetDecoys) {
//                                String[] pepSplit = pep.split(",");
//                                String line = pepSplit[1] + "\t" + pepSplit[0].split("\\|")[1] + "\n";
//                                if (!peptides.contains(line)) {
//                                    myWriter.write(line);
//                                    peptides.add(line);
//                                }
//                            }
//                            myWriter.close();
//
//                            DiannModelCaller.callModel(inputFile, false);
//                            allPreds = SpectralPredictionMapper.createSpectralPredictionMapper(
//                                    Constants.spectraRTPredFile, "DIA-NN", executorService).getPreds();
//                            Constants.spectraRTPredFile = null;
//                        } else { //mode for koina
//                            HashSet<String> allHits = km.writeFullPeptideFile(
//                                    jsonOutFolder + File.separator + model + "_full.tsv", model,
//                                    km.peptideSetDecoys);
//                            allPreds = km.getKoinaPredictions(allHits, model, 30,
//                                    jsonOutFolder + File.separator + model,
//                                    jsonOutFolder + File.separator + model + "_full.tsv");
//                        }
//
//                        //collect the data for testing
//                        expRTs = new ArrayList<>();
//                        predRTs = new ArrayList<>();
//                        for (int j = 0; j < pmMatcher.mzmlReaders.length; j++) {
//                            MzmlReader mzmlReader = pmMatcher.mzmlReaders[j];
//                            LinkedList<Integer> thisScanNums = km.scanNumsDecoys.get(pmMatcher.mzmlFiles[j].getName());
//                            LinkedList<String> thisPeptides = km.peptidesDecoys.get(pmMatcher.mzmlFiles[j].getName());
//
//                            for (int k = 0; k < thisScanNums.size(); k++) {
//                                try {
//                                    int scanNum = thisScanNums.get(k);
//                                    String peptide = thisPeptides.get(k).split(",")[0];
//                                    MzmlScanNumber msn = mzmlReader.getScanNumObject(scanNum);
//
//                                    float predRT = allPreds.get(peptide).RT;
//                                    float expRT = msn.RT;
//                                    predRTs.add(predRT);
//                                    expRTs.add(expRT);
//                                } catch (Exception e) {}
//                            }
//                        }
//                        rts = new double[2][expRTs.size()];
//                        for (int i = 0; i < expRTs.size(); i++) {
//                            //put pred first so it can be normalized to exp scale
//                            rts[0][i] = Math.round(predRTs.get(i) * 100d) / 100d;
//                            rts[1][i] = Math.round(expRTs.get(i) * 100d) / 100d;
//                        }
//
//                        //testing
//                        float[] rtDiffsDecoys = new float[rts[0].length];
//                        for (int i = 0; i < rtDiffsDecoys.length - 1; i++) {
//                            double rt = rts[0][i];
//                            rtDiffsDecoys[i] = (float) Math.abs(loess.invoke(rt) - rts[1][i]);
//                        }
//                        Arrays.sort(rtDiffsDecoys);
//                        datapointsRTDecoys.put(model, rtDiffsDecoys);
                    }

                    //get best model
                    HashMap<String, Integer> votes = new HashMap<>();
                    for (String model : datapointsRT.keySet()) {
                        votes.put(model, 0);
                    }
                    for (int i = 0; i < 10; i++) {
                        String bestModel = "";
                        float bestRtDiff = Float.MAX_VALUE;

                        for (String model : datapointsRT.keySet()) {
                            float[] rtDiffs = datapointsRT.get(model);
                            float deltaRT = rtDiffs[rtDiffs.length - 1 - i];
                            if (deltaRT < bestRtDiff) {
                                bestRtDiff = deltaRT;
                                bestModel = model;
                            }
                        }

                        votes.put(bestModel, votes.get(bestModel) + 1);
                    }

                    int mostVotes = 0;
                    for (String model : votes.keySet()) {
                        int vote = votes.get(model);
                        if (vote > mostVotes) {
                            mostVotes = vote;
                            Constants.rtModel = model;
                        }
                    }

                    printInfo("RT model chosen is " + Constants.rtModel);

                    //write file that has all values for each model
                    BufferedWriter bw = new BufferedWriter(new FileWriter(
                            Constants.outputDirectory + File.separator + "bestrtmodel.tsv"));
                    StringBuilder sb = new StringBuilder();
                    ArrayList<String> models = new ArrayList<>();
                    for (String model : datapointsRT.keySet()) {
                        sb.append(model).append("\t");
                        models.add(model);
                    }
                    sb.deleteCharAt(sb.length() - 1);
                    bw.write(sb.toString() + "\n");
                    for (int i = 0; i < datapointsRT.get(models.get(0)).length; i++) {
                        for (String model : models) {
                            bw.write(datapointsRT.get(model)[i] + "\t");
                        }
                        bw.write("\n");
                    }
                    bw.close();

//                    //decoy
//                    bw = new BufferedWriter(new FileWriter(
//                            Constants.outputDirectory + File.separator + "bestrtmodelDecoy.tsv"));
//                    sb = new StringBuilder();
//                    models = new ArrayList<>();
//                    for (String model : datapointsRTDecoys.keySet()) {
//                        sb.append(model).append("\t");
//                        models.add(model);
//                    }
//                    sb.deleteCharAt(sb.length() - 1);
//                    bw.write(sb.toString() + "\n");
//                    for (int i = 0; i < datapointsRTDecoys.get(models.get(0)).length; i++) {
//                        for (String model : models) {
//                            bw.write(datapointsRTDecoys.get(model)[i] + "\t");
//                        }
//                        bw.write("\n");
//                    }
//                    bw.close();
                }

                if (Constants.rtModel.isEmpty()) {
                    Constants.rtModel = "DIA-NN";
                }
                for (String model : Constants.KoinaRTmodels) {
                    if (model.equalsIgnoreCase(Constants.rtModel)) {
                        Constants.rtModel = model;
                    }
                }
                Constants.spectraRTPredModel += Constants.rtModel;
            } else {
                Constants.rtModel = "";
            }
            if (Constants.useSpectra) {
                //here, look for best spectra model
                if (Constants.findBestSpectraModel) {
                    printInfo("Searching for best spectra model for your data");
                    ArrayList<String> consideredModels = new ArrayList<>();
                    if (TMT) {
                        consideredModels.addAll(Constants.ms2SearchModelsTMT);
                    } else {
                        String[] ms2SearchModels = Constants.ms2SearchModelsString.split(",");
                        for (String model : ms2SearchModels) {
                            if (modelMapper.containsKey(model.toLowerCase())) {
                                consideredModels.add(modelMapper.get(model.toLowerCase()));
                            } else {
                                printError("No spectra model called " + model + ". Exiting.");
                                System.exit(0);
                            }
                        }
                    }
                    printInfo("Searching the following models: ");
                    printInfo(String.valueOf(consideredModels));

                    String jsonOutFolder = Constants.outputDirectory + File.separator + "best_model";
                    MyFileUtils.createWholeDirectory(jsonOutFolder);

                    HashMap<String, Double> medianSimilarities = new HashMap<>();
                    HashMap<String, TreeMap<Integer, ArrayList<Double>>> similarities = new HashMap<>();
                    HashMap<String, ArrayList<Double>> datapoints = new HashMap<>();
                    HashMap<String, ArrayList<Double>> datapointsDecoys = new HashMap<>();
                    for (String model : consideredModels) {
                        MyFileUtils.createWholeDirectory(jsonOutFolder + File.separator + model);
                        if (model.equals("DIA-NN")) { //mode for DIA-NN
                            km.writeFullPeptideFile(jsonOutFolder + File.separator + model + File.separator +
                                    "spectraRT_full.tsv", model, km.peptideSet);
                            String inputFile = jsonOutFolder + File.separator + model + File.separator + "spectraRT.tsv";
                            FileWriter myWriter = new FileWriter(inputFile);
                            myWriter.write("peptide" + "\t" + "charge\n");
                            HashSet<String> peptides = new HashSet<>();
                            for (String pep : km.peptideSet) {
                                String[] pepSplit = pep.split(",");
                                String line = pepSplit[1] + "\t" + pepSplit[0].split("\\|")[1] + "\n";
                                if (!peptides.contains(line)) {
                                    myWriter.write(line);
                                    peptides.add(line);
                                }
                            }
                            myWriter.close();

                            DiannModelCaller.callModel(inputFile, false);
                            PredictionEntryHashMap allPreds =
                                    SpectralPredictionMapper.createSpectralPredictionMapper(
                                    Constants.spectraRTPredFile, "DIA-NN", executorService).getPreds();
                            Constants.spectraRTPredFile = null;

                            //get median similarity
                            ArrayList<Double> similarity = new ArrayList<>();
                            for (PeptideObj peptideObj : km.getPeptideObjects(allPreds, km.scanNums, km.peptides)) {
                                similarity.add(peptideObj.spectralSimObj.unweightedSpectralEntropy());
                            }
                            try {
                                similarity.sort(Comparator.naturalOrder());
                                datapoints.put(model, similarity);
                            } catch (Exception e) {}
                            double median = StatMethods.medianDouble(similarity);
                            medianSimilarities.put(model, median);
                            printInfo("Median similarity for DIA-NN is " +
                                    String.format("%.4f", median));

//                            //decoy
//                            km.writeFullPeptideFile(jsonOutFolder + File.separator + model + File.separator +
//                                    "spectraRT_full.tsv", model, km.peptideSetDecoys);
//                            inputFile = jsonOutFolder + File.separator + model + File.separator + "spectraRT.tsv";
//                            myWriter = new FileWriter(inputFile);
//                            myWriter.write("peptide" + "\t" + "charge\n");
//                            peptides = new HashSet<>();
//                            for (String pep : km.peptideSetDecoys) {
//                                String[] pepSplit = pep.split(",");
//                                String line = pepSplit[1] + "\t" + pepSplit[0].split("\\|")[1] + "\n";
//                                if (!peptides.contains(line)) {
//                                    myWriter.write(line);
//                                    peptides.add(line);
//                                }
//                            }
//                            myWriter.close();
//
//                            DiannModelCaller.callModel(inputFile, false);
//                            allPreds = SpectralPredictionMapper.createSpectralPredictionMapper(
//                                            Constants.spectraRTPredFile, "DIA-NN", executorService).getPreds();
//                            Constants.spectraRTPredFile = null;
//
//                            //get median similarity
//                            similarity = new ArrayList<>();
//                            for (PeptideObj peptideObj : km.getPeptideObjects(allPreds, km.scanNumsDecoys, km.peptidesDecoys)) {
//                                similarity.add(peptideObj.spectralSimObj.unweightedSpectralEntropy());
//                            }
//                            try {
//                                similarity.sort(Comparator.naturalOrder());
//                                datapointsDecoys.put(model, similarity);
//                            } catch (Exception e) {}
                        } else { //mode for koina
                            if (Constants.nceModels.contains(model)) {
                                Object[] results = NCEcalibrator.calibrateNCE(model, new ArrayList<>(), km,
                                        Constants.outputDirectory + File.separator + "best_model" +
                                                File.separator + model, km.peptideSet,
                                        km.scanNums, km.peptides);
                                medianSimilarities.put(model + "&" + results[4], (Double) results[3]);
                                int bestNCE = (int) results[4];
                                similarities.put(model + "&" + bestNCE,
                                        (TreeMap<Integer, ArrayList<Double>>) results[2]);
                                try {
                                    ArrayList<Double> similarity = similarities.get(model + "&" + bestNCE).get(bestNCE);
                                    similarity.sort(Comparator.naturalOrder());
                                    datapoints.put(model, similarity);
                                } catch (Exception e) {}

//                                //decoy
//                                //just get similarities for decoys at same nce as targets
//                                results = NCEcalibrator.calibrateNCE(model, new ArrayList<>(), km,
//                                        Constants.outputDirectory + File.separator + "best_model" +
//                                                File.separator + model, km.peptideSetDecoys,
//                                        km.scanNumsDecoys, km.peptidesDecoys);
//                                TreeMap<Integer, ArrayList<Double>> decoysSims =
//                                        (TreeMap<Integer, ArrayList<Double>>) results[2];
//                                try {
//                                    ArrayList<Double> similarity = decoysSims.get(bestNCE);
//                                    similarity.sort(Comparator.naturalOrder());
//                                    datapointsDecoys.put(model, similarity);
//                                } catch (Exception e) {}
                            } else {
                                HashSet<String> allHits = km.writeFullPeptideFile(jsonOutFolder +
                                        File.separator + model + "_full.tsv", model, km.peptideSet);
                                PredictionEntryHashMap allPreds =
                                        km.getKoinaPredictions(allHits, model, 30,
                                                jsonOutFolder + File.separator + model,
                                                jsonOutFolder + File.separator + model + "_full.tsv");

                                //get median similarity
                                ArrayList<Double> similarity = new ArrayList<>();
                                for (PeptideObj peptideObj : km.getPeptideObjects(allPreds, km.scanNums, km.peptides)) {
                                    similarity.add(peptideObj.spectralSimObj.unweightedSpectralEntropy());
                                }
                                try {
                                    similarity.sort(Comparator.naturalOrder());
                                    datapoints.put(model, similarity);
                                } catch (Exception e) {}
                                double median = StatMethods.medianDouble(similarity);
                                medianSimilarities.put(model, median);
                                printInfo("Median similarity for " + model + " is " +
                                        String.format("%.4f", median));

//                                //decoy
//                                allHits = km.writeFullPeptideFile(jsonOutFolder +
//                                        File.separator + model + "_full.tsv", model, km.peptideSetDecoys);
//                                allPreds =
//                                        km.getKoinaPredictions(allHits, model, 30,
//                                                jsonOutFolder + File.separator + model,
//                                                jsonOutFolder + File.separator + model + "_full.tsv");
//
//                                //get median similarity
//                                similarity = new ArrayList<>();
//                                for (PeptideObj peptideObj : km.getPeptideObjects(allPreds, km.scanNumsDecoys, km.peptidesDecoys)) {
//                                    similarity.add(peptideObj.spectralSimObj.unweightedSpectralEntropy());
//                                }
//                                try {
//                                    similarity.sort(Comparator.naturalOrder());
//                                    datapointsDecoys.put(model, similarity);
//                                } catch (Exception e) {}
                            }
                        }
                    }

                    //find best median similarity
                    double bestMedian = 0d;
                    String bestModel = "";
                    for (Map.Entry<String, Double> entry : medianSimilarities.entrySet()) {
                        double median = entry.getValue();
                        if (median > bestMedian) {
                            bestMedian = median;
                            bestModel = entry.getKey();
                        }
                    }
                    if (bestModel.contains("&")) { //model that uses NCE
                        String[] modelSplit = bestModel.split("&");
                        Constants.spectraModel = modelSplit[0];
                        Constants.NCE = modelSplit[1];
                        NCEcalibrator.plotNCEchart(similarities.get(bestModel));
                    } else {
                        Constants.spectraModel = bestModel;
                    }
                    Constants.calibrateNCE = false;
                    Constants.autoSwitchFragmentation = false;
                    printInfo("Spectra model chosen is " + Constants.spectraModel);

                    //write file that has all values for each model
                    BufferedWriter bw = new BufferedWriter(new FileWriter(
                            Constants.outputDirectory + File.separator + "bestms2model.tsv"));
                    StringBuilder sb = new StringBuilder();
                    ArrayList<String> models = new ArrayList<>();
                    for (String model : datapoints.keySet()) {
                        sb.append(model).append("\t");
                        models.add(model);
                    }
                    sb.deleteCharAt(sb.length() - 1);
                    bw.write(sb.toString() + "\n");
                    for (int i = 0; i < datapoints.get(models.get(0)).size(); i++) {
                        for (String model : models) {
                            bw.write(datapoints.get(model).get(i) + "\t");
                        }
                        bw.write("\n");
                    }
                    bw.close();

//                    //decoy
//                    bw = new BufferedWriter(new FileWriter(
//                            Constants.outputDirectory + File.separator + "bestms2modelDecoy.tsv"));
//                    sb = new StringBuilder();
//                    models = new ArrayList<>();
//                    for (String model : datapointsDecoys.keySet()) {
//                        sb.append(model).append("\t");
//                        models.add(model);
//                    }
//                    sb.deleteCharAt(sb.length() - 1);
//                    bw.write(sb.toString() + "\n");
//                    for (int i = 0; i < datapointsDecoys.get(models.get(0)).size(); i++) {
//                        for (String model : models) {
//                            bw.write(datapointsDecoys.get(model).get(i) + "\t");
//                        }
//                        bw.write("\n");
//                    }
//                    bw.close();
                }

                if (Constants.spectraModel.isEmpty()) {
                    Constants.spectraModel = "DIA-NN";
                }
                for (String model : Constants.KoinaMS2models) {
                    if (model.equalsIgnoreCase(Constants.spectraModel)) {
                        Constants.spectraModel = model;
                    }
                }
                if (!Constants.spectraRTPredModel.isEmpty()) {
                    Constants.spectraRTPredModel += ",";
                }
                Constants.spectraRTPredModel += Constants.spectraModel;
            } else {
                Constants.spectraModel = "";
            }
            Constants.foundBest = true;

            //set useKoina based on model
            if (Constants.KoinaRTmodels.contains(Constants.rtModel) ||
            Constants.KoinaMS2models.contains(Constants.spectraModel)) {
                Constants.useKoina = true;
            }

            if (Constants.adaptiveFragmentNum) {
                Constants.topFragments = 36; //TODO think of better way than hardcoding
            } else if (Constants.divideFragments.equals("1")) { //standard setting of yb vs others
                Constants.divideFragments = "y_b;immonium_a_y-NL_b-NL_a-NL_internal_internal-NL_unknown";
                Constants.topFragments = 12;
            } else if (Constants.divideFragments.equals("2")) {
                Constants.divideFragments = "y;b;immonium;a;y-NL;b-NL;a-NL;internal;internal-NL;unknown";
                Constants.topFragments = 6;
            } else if (Constants.divideFragments.equals("3")) { //standard setting of yb vs others
                Constants.divideFragments = "y_b_y-NL_b-NL;immonium_a_a-NL_internal_internal-NL_unknown";
                Constants.topFragments = 12;
            } else if (Constants.divideFragments.equals("4")) { //etd
                Constants.divideFragments = "c_z;zdot_y_unknown";
                Constants.topFragments = 12;
            } else if (Constants.divideFragments.equals("5")) { //ethcd
                Constants.divideFragments = "b_y_c_z;immonium_a_cdot_zdot_y-NL_b-NL_a-NL_internal_internal-NL_unknown";
                Constants.topFragments = 12;
            } else if (Constants.divideFragments.equals("0") && Constants.spectraRTPredModel.equals("DIA-NN")) {
                //Constants.topFragments = 12; //may update in future
            }

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

            Constants.allowedFeatures = Constants.makeAllowedFeatures();
            String[] featuresArray = Constants.features.replaceAll("\\s", "").split(",");
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
                        featureLL.add("unweightedSpectralEntropy");
                    }
                } else {
                    featureLL.removeIf(Constants.spectraFeatures::contains);
                }
            } catch (Exception ignored) {
            }
            try {
                if (Constants.useRT) {
                    Set<String> intersection = new HashSet<>(featureLL);
                    intersection.retainAll(Constants.rtFeatures);
                    if (intersection.size() == 0) {
                        featureLL.add("deltaRTLOESS");
                    }
                } else {
                    featureLL.removeIf(Constants.rtFeatures::contains);
                }
            } catch (Exception ignored) {
            }
            try {
                if (Constants.useDetect) {
                    printInfo("Detect features not fully tested");
                    Set<String> intersection = new HashSet<>(featureLL);
                    intersection.retainAll(Constants.detectFeatures);
                    if (intersection.size() == 0) {
                        featureLL.add("detectFractionGreater");
                        featureLL.add("detectSubtractMissing");
                    }
                } else {
                    featureLL.removeIf(Constants.detectFeatures::contains);
                }
            } catch (Exception ignored) {
            }
            try {
                if (Constants.useIM) {
                    Set<String> intersection = new HashSet<>(featureLL);
                    intersection.retainAll(Constants.imFeatures);
                    if (intersection.isEmpty()) {
                        featureLL.add("deltaIMLOESS");
                    }
                } else {
                    for (String feature : Constants.imFeatures) {
                        featureLL.remove(feature);
                    }
                }
            } catch (Exception ignored) {
            }
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
            } catch (Exception ignored) {
            }
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
            } catch (Exception ignored) {
            }
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
            } catch (Exception ignored) {
            }
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
            } catch (Exception ignored) {
            }
            try {
                if (Constants.useIntensityDistributionSimilarity) {
                    if (!featureLL.contains("intensity_distribution_similarity")) {
                        featureLL.add("intensity_distribution_similarity");
                    }
                } else {
                    featureLL.remove("intensity_distribution_similarity");
                }
            } catch (Exception ignored) {
            }

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

            //don't need koina if pred file ready
            if (Constants.spectraRTPredFile != null) {
                Constants.useKoina = false;
            }
            //overriding if intermediate files already made
            if (Constants.spectraRTPrefix != null || Constants.spectraRTPredFile != null) {
                createSpectraRTPredFile = false;
            }
            if (Constants.detectPredInput != null || Constants.detectPredFile != null) {
                createDetectPredFile = false;
            }
            c.updateInputPaths(); //setting null paths

            //generate files for prediction models
            List<String> modelsList = Arrays.asList(Constants.spectraRTPredModel.split(","));
            ArrayList<String> models = new ArrayList<>(modelsList);
            if (models.size() != 1) {
                if (models.get(0).equals(models.get(1))) {
                    models.remove(1);
                    Constants.spectraRTPredModel = models.get(0);
                    if (Constants.spectraRTPredModel.equals("DIA-NN")) {
                        Constants.useKoina = false;
                    }
                }
            }
            if (createSpectraRTPredFile || Constants.createPredFileOnly) {
                //createfull is needed for everything
                boolean revertToKoina = Constants.useKoina;
                Constants.useKoina = false;
                PeptideFileCreator.createPeptideFile(pmMatcher,
                        Constants.spectraRTPrefix + "_full.tsv",
                        "createFull");
                if (revertToKoina) {
                    Constants.useKoina = true;
                }
                for (String currentModel : Constants.spectraRTPredModel.split(",")) {
                    if (Constants.useKoina && !currentModel.equals("DIA-NN")) {
                        //make sure using right model
                        boolean changed = KoinaMethods.switchModel();
                        if (changed) {
                            for (int m = 0; m < models.size(); m++) {
                                if (models.get(m).equals(currentModel)) {
                                    models.set(m, Constants.spectraModel);
                                }
                            }
                            currentModel = Constants.spectraModel;
                        }

                        if (Constants.KoinaMS2models.contains(currentModel) && Constants.calibrateNCE &&
                        Constants.nceModels.contains(currentModel)) {
                            Object[] modelInfo = NCEcalibrator.calibrateNCE(currentModel, models, km,
                                    Constants.outputDirectory + File.separator + "NCE_calibration",
                                    km.peptideSet, km.scanNums, km.peptides);
                            currentModel = (String) modelInfo[0];
                            models = (ArrayList<String>) modelInfo[1];
                            Constants.NCE = String.valueOf((int) modelInfo[4]);

                            //model may have changed
                            if (Constants.nceModels.contains(currentModel)) {
                                NCEcalibrator.plotNCEchart((TreeMap<Integer, ArrayList<Double>>) modelInfo[2]);
                            }
                        }

                        PeptideFileCreator.createPeptideFile(pmMatcher,
                                Constants.spectraRTPrefix + "_" + currentModel + ".json", currentModel);
                    } else {
                        switch (currentModel) {
                            case "DIA-NN":
                                if (Constants.DiaNN == null) {
                                    throw new IllegalArgumentException("path to DIA-NN executable must be provided");
                                }
                                printInfo("Generating input file for DIA-NN");
                                PeptideFileCreator.createPeptideFile(pmMatcher, Constants.spectraRTPrefix + ".tsv", "Diann");
                                break;
                            case "pDeep2":
                                printInfo("Generating input file for pDeep2");
                                PeptideFileCreator.createPeptideFile(pmMatcher, Constants.spectraRTPrefix + ".tsv", "pDeep2");
                                break;
                            case "pDeep3":
                                printInfo("Generating input file for pDeep3");
                                PeptideFileCreator.createPeptideFile(pmMatcher, Constants.spectraRTPrefix + ".tsv", "pDeep3");
                                break;
                            case "PredFull":
                                printInfo("Generating input file for PredFull");
                                PeptideFileCreator.createPeptideFile(pmMatcher, Constants.spectraRTPrefix + ".tsv", "PredFull");
                                break;
                            case "Prosit":
                                printInfo("Generating input file for Prosit");
                                PeptideFileCreator.createPeptideFile(pmMatcher, Constants.spectraRTPrefix + ".csv", "Prosit");
                                break;
                            case "PrositTMT":
                                printInfo("Generating input file for PrositTMT");
                                PeptideFileCreator.createPeptideFile(pmMatcher, Constants.spectraRTPrefix + ".csv", "PrositTMT");
                                break;
                            case "alphapeptdeep":
                                printInfo("Generating input file for alphapeptdeep");
                                PeptideFileCreator.createPeptideFile(pmMatcher, Constants.spectraRTPrefix + ".csv", "alphapeptdeep");
                                break;
                            default:
                                printError("spectraRTPredModel must be one of DIA-NN, Prosit, PrositTMT, " +
                                        "PredFull, pDeep2, pDeep3, or alphapeptdeep");
                                System.exit(-1);
                        }
                    }
                }

                if (Constants.createPredFileOnly) {
                    printInfo("Successfully created input file for prediction model. Stopping here");
                    System.exit(0);
                }
            }
//            if (createDetectAllPredFile) {
//                printInfo("Generating input file for DeepMSPeptide");
//                //long startTime = System.nanoTime();
//                //Constants.setFastaReader(peptideFileCreator.createPeptideFile(pmMatcher.pinFiles, Constants.detectPredInput, "DeepMSPeptideAll", "pin"));
//                //long endTime = System.nanoTime();
//                //long duration = (endTime - startTime);
//            } else if (createDetectPredFile) {
//                printInfo("Generating input file for DeepMSPeptide");
//                peptideFileCreator.createPeptideFile(pmMatcher.pinFiles, Constants.detectPredInput, "DeepMSPeptide", "pin");
//            }


            //generate predictions
            KoinaLibReader klr = new KoinaLibReader();
            KoinaModelCaller kmc = new KoinaModelCaller();
            //this is just so that ms2 is predicted first, and which works for Koina
            Collections.reverse(models);
            boolean onlyUsedKoina = true;
            String spectraRTPredFile = "";
            if ((Constants.spectraRTPredFile == null) && (createSpectraRTPredFile2)) {
                for (String currentModel : models) {
                    if (Constants.useKoina && !currentModel.equals("DIA-NN")) {
                        kmc.callModel(currentModel, klr, Constants.JsonDirectory, executorService, true, true);
                        spectraRTPredFile = Constants.outputDirectory + File.separator + "spectraRT_koina.mgf" +
                                spectraRTPredFile;
                    } else {
                        DiannModelCaller.callModel(Constants.spectraRTPrefix + ".tsv", true);
                        onlyUsedKoina = false;
                        spectraRTPredFile = Constants.spectraRTPrefix + ".predicted.bin" + spectraRTPredFile;
                    }
                    spectraRTPredFile = "," + spectraRTPredFile;
                }
            }
            if (Constants.useKoina) {
                kmc.assignMissingPeptidePredictions(klr,
                        Constants.spectraRTPrefix + "_full.tsv");
                MgfFileWriter mfw = new MgfFileWriter(klr);
                mfw.write(Constants.outputDirectory + File.separator + "spectraRT_koina.mgf");
            }
            if (onlyUsedKoina) {
                Constants.spectralPredictionMapper = klr;
            } else {
                Constants.spectraRTPredFile = spectraRTPredFile.substring(1);
            }

            //create new pin file with features
            printInfo("Generating edited pin with following features: " + Arrays.toString(featuresArray));
            long start = System.nanoTime();
            if (Constants.spectraRTPredModel.contains("PredFull")) {
                Constants.matchWithDaltons = true; //they report predictions in bins
            }
            PercolatorFormatter.editPin(pmMatcher, Constants.spectraRTPredFile, Constants.detectPredFile,
                    featuresArray, Constants.editedPin, executorService);
            executorService.shutdown();

            //print parameters to ps
            //printParamsPS();

            //delete pred files
            if (Constants.deletePreds) {
                File predFile = new File(Constants.outputDirectory + File.separator + "spectraRT.tsv");
                predFile.delete();
                predFile = new File(Constants.outputDirectory + File.separator + "spectraRT_full.tsv");
                predFile.delete();
                predFile = new File(Constants.outputDirectory + File.separator + "spectraRT.predicted.bin");
                predFile.delete();
            }
            //clean directories
            MyFileUtils.deleteWholeDirectory(Constants.outputDirectory + File.separator + "best_model");
            MyFileUtils.deleteWholeDirectory(Constants.outputDirectory + File.separator + "NCE_calibration");
            MyFileUtils.deleteWholeDirectory(Constants.JsonDirectory);

            //file with used models
            try {
                // Read all lines from the file into a list
                List<String> lines = Files.readAllLines(Paths.get(Constants.paramsList));
                AtomicBoolean spectrafound = new AtomicBoolean(false);
                AtomicBoolean rtfound = new AtomicBoolean(false);

                // Stream the list, and if a line starts with the specified prefix,
                // replace everything after the prefix with newSuffix
                // Otherwise keep the line as is
                List<String> modifiedLines = lines.stream()
                        .map(line -> {
                            if (line.trim().startsWith("spectraModel")) {
                                spectrafound.set(true);
                                return "spectraModel=" + Constants.spectraModel;
                            } else if (line.trim().startsWith("rtModel")) {
                                rtfound.set(true);
                                return "rtModel=" + Constants.rtModel;
                            }
                            return line;
                        })
                        .collect(Collectors.toList());

                // If no line was found with the specified prefix, add the new line
                if (!spectrafound.get()) {
                    modifiedLines.add("spectraModel=" + Constants.spectraModel);
                }
                if (!rtfound.get()) {
                    modifiedLines.add("rtModel=" + Constants.rtModel);
                }

                // Write the modified lines back to the file
                Files.write(Paths.get(Constants.paramsList), modifiedLines);
            } catch (IOException e) {
                e.printStackTrace();
            }

            long end = System.nanoTime();
            long duration = (end - start);
            printInfo("Feature calculation, edited pin writing, and QC plot generation done in " +
                    duration / 1000000 + " ms");
            System.exit(0);
        } catch (Exception e) {
            e.printStackTrace();
            executorService.shutdown();
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
            printError("could not write final params");
            e.getStackTrace();
            System.exit(1);
        }
    }

    static private void printParamsPS() {
        try {
            printInfo("Final parameters used for feature annotation:");
            //Constants c = new Constants();

            Field[] f = Constants.class.getFields();
            for (Field field : f) {
                if ((field.getModifiers() & Modifier.FINAL) != Modifier.FINAL) {
                    if (!field.getName().equals("paramsList")) {
                        if (field.getName().equals("fragmentIonHierarchy")) {
                            //printInfo("\t" + field.getName() + " = " + Arrays.toString((String[]) field.get(c)));
                            printInfo("\t" + field.getName() + " = " + Arrays.toString((String[]) field.get(Constants.class)));
                        } else {
                            //printInfo("\t" + field.getName() + " = " + field.get(c));
                            printInfo("\t" + field.getName() + " = " + field.get(Constants.class));
                        }
                    }
                }
            }
        } catch (Exception e) {
            printError("could not write final params");
            e.getStackTrace();
            System.exit(1);
        }
    }
}
