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

package mainsteps;

import allconstants.Constants;
import allconstants.ConstantsInterface;
import allconstants.FragmentIonConstants;
import allconstants.LowercaseModelMapper;
import allconstants.NceConstants;
import bestmodelsearch.BestModelSearcher;
import features.spectra.MassCalculator;
import koinaclasses.KoinaMethods;
import koinaclasses.NCEcalibrator;
import modelcallers.DiannModelCaller;
import modelcallers.KoinaModelCaller;
import peptideptmformatting.PTMhandler;
import peptideptmformatting.PeptideFormatter;
import predictions.PredictionEntry;
import predictions.PredictionEntryHashMap;
import readers.datareaders.MzmlReader;
import readers.predictionreaders.KoinaLibReader;
import readers.predictionreaders.LibraryPredictionMapper;
import utils.MyFileUtils;
import utils.NumericUtils;
import writers.MgfFileWriter;
import writers.PeptideFileCreator;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.lang.reflect.Field;
import java.lang.reflect.InvocationTargetException;
import java.lang.reflect.Modifier;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;
import java.util.LinkedList;
import java.util.List;
import java.util.Locale;
import java.util.Map;
import java.util.Objects;
import java.util.Set;
import java.util.TreeMap;
import java.util.concurrent.ScheduledThreadPoolExecutor;
import java.util.concurrent.atomic.AtomicBoolean;
import java.util.stream.Collectors;

import static features.rtandim.LoessUtilities.gridSearchCV;
import static peptideptmformatting.PTMhandler.tmtMasses;
import static utils.Print.*;

//this is what I use in the java jar file
public class MainClass {
    public static ScheduledThreadPoolExecutor executorService;
    public static void main(String[] args) throws Exception {
        Locale.setDefault(Locale.US);
        printInfo("MSBooster v1.3.12");

        try {
            //accept command line inputs
            HashMap<String, String> params = new HashMap<String, String>();

            //setting new values
            for (int i = 0; i < args.length; i++) {
                String key = args[i].substring(2); //remove --
                if (key.equals("help")) { //help message
                    System.out.println("Usage: java -jar MSBooster-1.2.57.jar [flags]");
                    System.out.println("Usage: java -jar MSBooster-1.2.57.jar --paramsList {*.txt}");
                    System.out.println("A tool for annotating PSM pin files with deep learning-based features. " +
                            "This has been tested on MSFragger pin files and may generalize to pin files from any tool," +
                            " although this has not been specifically tested. ");
                    System.out.println("For more information on required and optional parameters, please refer to the " +
                            "MSBooster Github at https://github.com/Nesvilab/MSBooster");
                    System.out.println("Printing example paramsList to exampleParams.txt. It will be a bit verbose, so " +
                            "please delete any parameters not relevant to your analysis.");
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
                    boolean ppmToDa = false; //set Da tolernace once all parameters read in. No guarantee which param is read first
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
                                    ppmToDa = true;
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
                                                final_vals += "&";
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
                                                final_vals += "&";
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
                            case "add_B_user_amino_acid":
                                if (! val.equals("0.0")) {
                                    MassCalculator.AAmap.put('B', Float.valueOf(val));
                                }
                                break;
                            case "add_J_user_amino_acid":
                                if (! val.equals("0.0")) {
                                    MassCalculator.AAmap.put('J', Float.valueOf(val));
                                }
                                break;
                            case "add_O_user_amino_acid":
                                if (! val.equals("0.0")) {
                                    MassCalculator.AAmap.put('O', Float.valueOf(val));
                                }
                                break;
                            case "add_U_user_amino_acid":
                                if (! val.equals("0.0")) {
                                    MassCalculator.AAmap.put('U', Float.valueOf(val));
                                }
                                break;
                            case "add_Z_user_amino_acid":
                                if (! val.equals("0.0")) {
                                    MassCalculator.AAmap.put('Z', Float.valueOf(val));
                                }
                                break;
                            case "add_X_user_amino_acid":
                                if (! val.equals("0.0")) {
                                    MassCalculator.AAmap.put('X', Float.valueOf(val));
                                }
                                break;
                            default:
                                if (key.startsWith("variable_mod_")) { //reading in TMT/iTraq values when variable
                                    Double varmod = Double.valueOf(val.trim().split(" ")[0]);
                                    for (double potentialTmtMass : tmtMasses) {
                                        if (NumericUtils.massesCloseEnough(varmod, potentialTmtMass)) {
                                            printInfo("TMT/iTRAQ mass detected in fragger.params as variable modification: " +
                                                    potentialTmtMass);
                                            PTMhandler.setTmtMass(potentialTmtMass);
                                            Constants.searchTMTmodels = true;
                                            break;
                                        }
                                    }
                                } else if (key.startsWith("add_")) { //reading in TMT/iTraq values when fixed mode
                                    for (double potentialTmtMass : tmtMasses) {
                                        double fixedMod = Double.parseDouble(val);
                                        if (NumericUtils.massesCloseEnough(fixedMod, potentialTmtMass)) {
                                            printInfo("TMT/iTRAQ mass detected in fragger.params as fixed modification: " +
                                                    potentialTmtMass);
                                            PTMhandler.setTmtMass(potentialTmtMass);
                                            Constants.searchTMTmodels = true;
                                            break;
                                        }
                                    }
                                }
                        }
                    }

                    float tol = Float.parseFloat(params.get("ppmTolerance"));
                    if (ppmToDa) { //read in from msfragger params. Low res tolerance used for all matching
                        Constants.matchWithDaltons = true;
                        Constants.matchWithDaltonsAux = true;
                        Constants.matchWithDaltonsDefault = true;
                        Constants.DaTolerance = tol;
                        printInfo("Using Dalton tolerance of " + Constants.DaTolerance + " Da");
                    } else {
                        if (tol >= 100f) {
                            params.put("lowResppmTolerance", String.valueOf(tol));
                        } else {
                            params.put("highResppmTolerance", String.valueOf(tol));
                        }
                        Constants.ppmTolerance = Constants.highResppmTolerance;
                    }
                }
            }

            HashSet<String> fields = new HashSet<>();
            for (Field f : Constants.class.getDeclaredFields()) {
                fields.add(f.getName());
            }
            Constants c = new Constants();

            HashSet<String> fieldsFragmentIon = new HashSet<>();
            for (Field f : FragmentIonConstants.class.getDeclaredFields()) {
                fieldsFragmentIon.add(f.getName());
            }
            FragmentIonConstants fic = new FragmentIonConstants();

            HashSet<String> fieldsNCE = new HashSet<>();
            for (Field f : NceConstants.class.getDeclaredFields()) {
                fieldsNCE.add(f.getName());
            }
            NceConstants nc = new NceConstants();

            for (Map.Entry<String, String> entry : params.entrySet()) {
                String key = entry.getKey();
                if (key.charAt(0) == '#') { //comment
                    continue;
                }
                if (key.charAt(0) == '/') { //comment
                    continue;
                }
                if (!fields.contains(key) && !fieldsFragmentIon.contains(key) && !fieldsNCE.contains(key)) {
                    throw new Exception(entry.getKey() + " is not a valid parameter");
                } else if (fields.contains(key)) {
                    Field field = Constants.class.getField(key);
                    updateField(field, entry.getValue(), c);
                    if (key.equals("ppmTolerance")) {
                        updateField(field, String.valueOf(Float.parseFloat(entry.getValue()) / 1000000f), c);
                    }
                } else if (fieldsFragmentIon.contains(key)) {
                    Field field = FragmentIonConstants.class.getField(key);
                    updateField(field, entry.getValue(), fic);
                } else {
                    Field field = NceConstants.class.getField(key);
                    updateField(field, entry.getValue(), nc);
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

            //collecting mass offsets
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

            //get matched pin files for mzML files
            PinMzmlMatcher pmMatcher = new PinMzmlMatcher(Constants.mzmlDirectory, Constants.pinPepXMLDirectory);

            //exit if no models used
            if (!Constants.useSpectra && !Constants.useRT && !Constants.useIM) {
                printInfo("useSpectra, useRT, and useIM are all false. Exiting.");
                System.exit(0);
            }

            //update fragment ion types based on fragmentation type
            FragmentIonConstants.makeFragmentIonHierarchy();
            Constants.matchedIntensitiesFeatures = Constants.makeMatchedIntensitiesFeatures();
            Constants.peakCountsFeatures = Constants.makePeakCountsFeatures();
            Constants.predIntensitiesFeatures = Constants.makePredIntensitiesFeatures();
            Constants.individualSpectralSimilaritiesFeatures = Constants.makeIndividualSpectralSimilarities();
            Constants.intensitiesDifferenceFeatures = Constants.makeintensitiesDifference();

            //if best model search, initially set models to nothing
            if (Constants.findBestSpectraModel) {
                Constants.spectraModel = "";
            }
            if (Constants.findBestRtModel) {
                Constants.rtModel = "";
            }
            if (Constants.findBestImModel) {
                Constants.imModel = "";
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
            if (modelMapper.containsKey(Constants.imModel.toLowerCase())) {
                Constants.imModel = modelMapper.get(Constants.imModel.toLowerCase());
            } else {
                printError("No IM model called " + Constants.imModel + ". Exiting.");
                System.exit(0);
            }
            if (modelMapper.containsKey(Constants.auxSpectraModel.toLowerCase())) {
                Constants.auxSpectraModel = modelMapper.get(Constants.auxSpectraModel.toLowerCase());
            } else {
                printError("No auxiliary spectra model called " + Constants.auxSpectraModel + ". Exiting.");
                System.exit(0);
            }

            //needed for nce calibration and best model search
            //TODO: could also make a new koinamethods object for decoys
            KoinaMethods km = new KoinaMethods(pmMatcher);
            if (Constants.findBestRtModel || Constants.findBestSpectraModel || Constants.findBestImModel ||
                    Constants.KoinaMS2models.contains(Constants.spectraModel) ||
                    Constants.KoinaRTmodels.contains(Constants.rtModel) ||
                    Constants.KoinaIMmodels.contains(Constants.imModel) ||
                    ! Constants.auxSpectraModel.isEmpty()
            ) {
                km.getTopPeptides();
            }

            //if different RT and spectra models
            HashMap<String, float[]> datapointsRT = new HashMap<>();
            if (Constants.useRT) {
                //here, look for best rt model
                if (Constants.findBestRtModel) {
                    printInfo("Searching for best RT model for your data");
                    ArrayList<String> consideredModels = new ArrayList<>();
                    if (Constants.searchTMTmodels) {
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
                    printInfo("Searching the following models: " + consideredModels);

                    String jsonOutFolder = Constants.outputDirectory + File.separator + "best_model";
                    MyFileUtils.createWholeDirectory(jsonOutFolder);

                    for (String model : consideredModels) {
                        Constants.rtModel = model;
                        PredictionEntryHashMap rtPreds = null;
                        if (model.equals("DIA-NN")) { //mode for DIA-NN
                            MyFileUtils.createWholeDirectory(jsonOutFolder + File.separator + model);

                            PeptideFileCreator.createPartialFile(
                                    jsonOutFolder + File.separator + model + File.separator + "spectraRT_full.tsv",
                                    model, km.peptideArraylist);
                            String inputFile = jsonOutFolder + File.separator + model + File.separator + "spectraRT.tsv";
                            FileWriter myWriter = new FileWriter(inputFile);
                            myWriter.write("peptide" + "\t" + "charge\n");
                            HashSet<String> peptides = new HashSet<>();
                            for (PeptideFormatter pf : km.peptideArraylist) {
                                String line = pf.getDiann() + "\t" + pf.getCharge() + "\n";
                                if (!peptides.contains(line)) {
                                    myWriter.write(line);
                                    peptides.add(line);
                                }
                            }
                            for (PeptideFormatter pf : km.peptideArraylist) {
                                pf.foundUnimods.clear();
                            }
                            myWriter.close();

                            String predFileString = DiannModelCaller.callModel(inputFile, false);
                            rtPreds = LibraryPredictionMapper.createLibraryPredictionMapper(
                                    predFileString, "DIA-NN", executorService).getPreds();
                        } else { //mode for koina
                            PeptideFileCreator.createPartialFile(
                                    jsonOutFolder + File.separator + model + "_full.tsv",
                                    model, km.peptideArraylist);
                            HashSet<String> allHits = KoinaMethods.createPartialKoinaSet(model, km.peptideArraylist);
                            for (PeptideFormatter pf : km.peptideArraylist) {
                                pf.foundUnimods.clear();
                            }
                            rtPreds = km.getKoinaPredictions(allHits, model, 30,
                                    jsonOutFolder + File.separator + model,
                                    jsonOutFolder + File.separator + model + "_full.tsv");
                        }

                        //collect the data for testing
                        ArrayList<Float> expRTs = new ArrayList<>();
                        ArrayList<Float> predRTs = new ArrayList<>();
                        for (int j = 0; j < pmMatcher.mzmlReaders.length; j++) {
                            MzmlReader mzmlReader = pmMatcher.mzmlReaders[j];
                            LinkedList<Integer> thisScanNums = km.scanNums.get(pmMatcher.mzmlFiles[j].getName());
                            LinkedList<PeptideFormatter> thisPeptides = km.peptides.get(pmMatcher.mzmlFiles[j].getName());

                            for (int k = 0; k < thisScanNums.size(); k++) {
                                try {
                                    int scanNum = thisScanNums.get(k);
                                    String peptide = thisPeptides.get(k).getBaseCharge();
                                    MzmlScanNumber msn = mzmlReader.getScanNumObject(scanNum);

                                    float predRT = rtPreds.get(peptide).RT;
                                    float expRT = msn.RT;
                                    predRTs.add(predRT);
                                    expRTs.add(expRT);
                                } catch (Exception ignored) {}
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
                        Object[] bandwidth_loess_rtdiffs = gridSearchCV(rts, floatBandwidths, false);
                        float[] rtdiffs = (float[]) bandwidth_loess_rtdiffs[2];
                        datapointsRT.put(model, rtdiffs);
                    }

                    //get best model
                    BestModelSearcher bms = new BestModelSearcher();
                    if (Constants.rtBestModelSearchMetric.equals("top")) {
                        printInfo("Choosing RT model based on top 10 consensus method");
                        Constants.rtModel = bms.consensus(datapointsRT, 10, false, true);
                    } else if (Constants.rtBestModelSearchMetric.equals("RMSE")) {
                        printInfo("Choosing RT model based on RMSE method");
                        Constants.rtModel = bms.RMSE(datapointsRT, false);
                    } else {
                        printError("Unknown RT model search metric: " + Constants.rtBestModelSearchMetric);
                        printError("Exiting");
                        System.exit(1);
                    }
                    printInfoHighlight("RT model chosen is " + Constants.rtModel);

                    //write file that has all values for each model
                    bms.writeScores(Constants.outputDirectory + File.separator + "bestrtmodel.tsv");
                }

                if (Constants.rtModel.isEmpty()) {
                    Constants.rtModel = "DIA-NN";
                }
                for (String model : Constants.KoinaRTmodels) {
                    if (model.equalsIgnoreCase(Constants.rtModel)) {
                        Constants.rtModel = model;
                    }
                }
            } else {
                Constants.rtModel = "";
            }

            HashMap<String, float[]> datapointsIM = new HashMap<>();
            if (Constants.useIM) {
                if (Constants.findBestImModel) {
                    printInfo("Searching for best IM model for your data");
                    ArrayList<String> consideredModels = new ArrayList<>();
                    String[] imSearchModels = Constants.imSearchModelsString.split(",");
                    for (String model : imSearchModels) {
                        if (modelMapper.containsKey(model.toLowerCase())) {
                            consideredModels.add(modelMapper.get(model.toLowerCase()));
                        } else {
                            printError("No IM model called " + model + ". Exiting.");
                            System.exit(0);
                        }
                    }
                    printInfo("Searching the following models: " + consideredModels);

                    String jsonOutFolder = Constants.outputDirectory + File.separator + "best_model";
                    MyFileUtils.createWholeDirectory(jsonOutFolder);

                    for (String model : consideredModels) {
                        Constants.imModel = model;
                        PredictionEntryHashMap imPreds = null;
                        if (model.equals("DIA-NN")) { //mode for DIA-NN
                            MyFileUtils.createWholeDirectory(jsonOutFolder + File.separator + model);

                            PeptideFileCreator.createPartialFile(
                                    jsonOutFolder + File.separator + model + File.separator + "spectraRT_full.tsv",
                                    model, km.peptideArrayListIM);
                            KoinaMethods.createPartialKoinaSet(model, km.peptideArrayListIM);
                            String inputFile = jsonOutFolder + File.separator + model + File.separator + "spectraRT.tsv";
                            FileWriter myWriter = new FileWriter(inputFile);
                            myWriter.write("peptide" + "\t" + "charge\n");
                            HashSet<String> peptides = new HashSet<>();
                            for (PeptideFormatter pf : km.peptideArrayListIM) {
                                String line = pf.getDiann() + "\t" + pf.getCharge() + "\n";
                                if (!peptides.contains(line)) {
                                    myWriter.write(line);
                                    peptides.add(line);
                                }
                            }
                            for (PeptideFormatter pf : km.peptideArrayListIM) {
                                pf.foundUnimods.clear();
                            }
                            myWriter.close();

                            String predFileString = DiannModelCaller.callModel(inputFile, false);
                            imPreds = LibraryPredictionMapper.createLibraryPredictionMapper(
                                    predFileString, "DIA-NN", executorService).getPreds();
                        } else { //mode for koina
                            PeptideFileCreator.createPartialFile(
                                    jsonOutFolder + File.separator + model + "_full.tsv",
                                    model, km.peptideArrayListIM);
                            HashSet<String> allHits = KoinaMethods.createPartialKoinaSet(model, km.peptideArrayListIM);
                            for (PeptideFormatter pf : km.peptideArrayListIM) {
                                pf.foundUnimods.clear();
                            }
                            imPreds = km.getKoinaPredictions(allHits, model, 30,
                                    jsonOutFolder + File.separator + model,
                                    jsonOutFolder + File.separator + model + "_full.tsv");
                        }

                        //collect the data for testing
                        ArrayList<Float> expIMs = new ArrayList<>();
                        ArrayList<Float> predIMs = new ArrayList<>();
                        for (int j = 0; j < pmMatcher.mzmlReaders.length; j++) {
                            MzmlReader mzmlReader = pmMatcher.mzmlReaders[j];
                            LinkedList<Integer> thisScanNums = km.scanNumsIM.get(pmMatcher.mzmlFiles[j].getName());
                            LinkedList<PeptideFormatter> thisPeptides = km.peptidesIM.get(pmMatcher.mzmlFiles[j].getName());

                            for (int k = 0; k < thisScanNums.size(); k++) {
                                try {
                                    int scanNum = thisScanNums.get(k);
                                    String peptide = thisPeptides.get(k).getBaseCharge();
                                    MzmlScanNumber msn = mzmlReader.getScanNumObject(scanNum);

                                    float predIM = imPreds.get(peptide).IM;
                                    float expIM = msn.IM;
                                    predIMs.add(predIM);
                                    expIMs.add(expIM);
                                } catch (Exception ignored) {}
                            }
                        }
                        double[][] ims = new double[2][expIMs.size()];
                        for (int i = 0; i < expIMs.size(); i++) {
                            //put pred first so it can be normalized to exp scale
                            ims[0][i] = Math.round(predIMs.get(i) * 100d) / 100d;
                            ims[1][i] = Math.round(expIMs.get(i) * 100d) / 100d;
                        }

                        //testing
                        String[] bandwidths = Constants.loessBandwidth.split(",");
                        float[] floatBandwidths = new float[bandwidths.length];
                        for (int i = 0; i < bandwidths.length; i++) {
                            floatBandwidths[i] = Float.parseFloat(bandwidths[i]);
                        }
                        Object[] bandwidth_loess_imdiffs = gridSearchCV(ims, floatBandwidths, false);
                        float[] imdiffs = (float[]) bandwidth_loess_imdiffs[2];
                        datapointsIM.put(model, imdiffs);
                    }

                    //get best model //TODO: test which method is more accurate
                    BestModelSearcher bms = new BestModelSearcher();
                    if (Constants.imBestModelSearchMetric.equals("top")) {
                        printInfo("Choosing IM model based on top 10 consensus method");
                        Constants.imModel = bms.consensus(datapointsIM, 10, false, true);
                    } else if (Constants.imBestModelSearchMetric.equals("RMSE")) {
                        printInfo("Choosing IM model based on RMSE method");
                        Constants.imModel = bms.RMSE(datapointsIM, false);
                    } else {
                        printError("Unknown IM model search metric: " + Constants.imBestModelSearchMetric);
                        printError("Exiting");
                        System.exit(1);
                    }
                    printInfoHighlight("IM model chosen is " + Constants.imModel);

                    //write file that has all values for each model
                    bms.writeScores(Constants.outputDirectory + File.separator + "bestimmodel.tsv");
                }

                if (Constants.imModel.isEmpty()) {
                    Constants.imModel = "DIA-NN";
                }
                for (String model : Constants.KoinaIMmodels) {
                    if (model.equalsIgnoreCase(Constants.imModel)) {
                        Constants.imModel = model;
                    }
                }
            } else {
                Constants.imModel = "";
            }

            if (Constants.useSpectra) {
                //here, look for best spectra model
                if (Constants.findBestSpectraModel) {
                    printInfo("Searching for best spectra model for your data");
                    ArrayList<String> consideredModels = new ArrayList<>();
                    if (Constants.searchTMTmodels) {
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
                    printInfo("Searching the following models: " + consideredModels);

                    String jsonOutFolder = Constants.outputDirectory + File.separator + "best_model";
                    MyFileUtils.createWholeDirectory(jsonOutFolder);

                    HashMap<String, TreeMap<Integer, ArrayList<Double>>> similarities = new HashMap<>();
                    HashMap<String, float[]> datapointsSpectra = new HashMap<>();
                    for (String model : consideredModels) {
                        Constants.spectraModel = model;

                        //set matching with Da or not
                        //if matchWithDaltons are true, also accept that (e.g. from reading it from fragger.params)
                        if (model.equalsIgnoreCase("predfull")) {
                            Constants.matchWithDaltonsDefault = true;
                        } else {
                            if (Constants.matchWithDaltons == null) {
                                Constants.matchWithDaltonsDefault = false;
                            } else {
                                Constants.matchWithDaltonsDefault = Constants.matchWithDaltons;
                            }
                        }

                        MyFileUtils.createWholeDirectory(jsonOutFolder + File.separator + model);
                        if (model.equals("DIA-NN")) { //mode for DIA-NN
                            PeptideFileCreator.createPartialFile(
                                    jsonOutFolder + File.separator + model + File.separator + "spectraRT_full.tsv",
                                    model, km.peptideArraylist);
                            KoinaMethods.createPartialKoinaSet(model, km.peptideArraylist);
                            String inputFile = jsonOutFolder + File.separator + model + File.separator + "spectraRT.tsv";
                            FileWriter myWriter = new FileWriter(inputFile);
                            myWriter.write("peptide" + "\t" + "charge\n");
                            HashSet<String> peptides = new HashSet<>();
                            for (PeptideFormatter pf : km.peptideArraylist) {
                                String line = pf.getDiann() + "\t" + pf.getCharge() + "\n";
                                if (!peptides.contains(line)) {
                                    myWriter.write(line);
                                    peptides.add(line);
                                }
                            }
                            for (PeptideFormatter pf : km.peptideArraylist) {
                                pf.foundUnimods.clear();
                            }
                            myWriter.close();

                            String predFileString = DiannModelCaller.callModel(inputFile, false);
                            PredictionEntryHashMap allPreds =
                                    LibraryPredictionMapper.createLibraryPredictionMapper(
                                            predFileString, "DIA-NN", executorService).getPreds();

                            ArrayList<Double> similarity = new ArrayList<>();
                            for (PeptideObj peptideObj : km.getPeptideObjects(allPreds, km.scanNums, km.peptides)) {
                                similarity.add(peptideObj.spectralSimObj.unweightedSpectralEntropy());
                            }
                            try {
                                float[] floatSimilarities = new float[similarity.size()];
                                for (int i = 0; i < similarity.size(); i++) {
                                    floatSimilarities[i] = similarity.get(i).floatValue();
                                }
                                datapointsSpectra.put(model, floatSimilarities);
                            } catch (Exception e) {}
                        } else { //mode for koina
                            if (NceConstants.nceModels.contains(model) && NceConstants.calibrateNCE) {
                                Object[] results = NCEcalibrator.calibrateNCE(model, km,
                                        Constants.outputDirectory + File.separator + "best_model" +
                                                File.separator + model, false);
                                for (PeptideFormatter pf : km.peptideArraylist) {
                                    pf.foundUnimods.clear();
                                }
                                int bestNCE = (int) results[2];
                                NceConstants.calibratedModels.put(model, String.valueOf(bestNCE));
                                similarities.put(model + "&" + bestNCE,
                                        (TreeMap<Integer, ArrayList<Double>>) results[0]);
                                try {
                                    ArrayList<Double> similarity = similarities.get(model + "&" + bestNCE).get(bestNCE);
                                    float[] floatSimilarities = new float[similarity.size()];
                                    for (int i = 0; i < similarity.size(); i++) {
                                        floatSimilarities[i] = similarity.get(i).floatValue();
                                    }
                                    datapointsSpectra.put(model, floatSimilarities);
                                } catch (Exception e) {}
                            } else {
                                PeptideFileCreator.createPartialFile(
                                        jsonOutFolder + File.separator + model + "_full.tsv",
                                        model, km.peptideArraylist);
                                HashSet<String> allHits = KoinaMethods.createPartialKoinaSet(model, km.peptideArraylist);
                                for (PeptideFormatter pf : km.peptideArraylist) {
                                    pf.foundUnimods.clear();
                                }
                                PredictionEntryHashMap allPreds =
                                        km.getKoinaPredictions(allHits, model, (int) Float.parseFloat(NceConstants.getNCE()),
                                                jsonOutFolder + File.separator + model,
                                                jsonOutFolder + File.separator + model + "_full.tsv");

                                ArrayList<Double> similarity = new ArrayList<>();
                                for (PeptideObj peptideObj : km.getPeptideObjects(allPreds, km.scanNums, km.peptides)) {
                                    if (peptideObj != null) {
                                        similarity.add(peptideObj.spectralSimObj.unweightedSpectralEntropy());
                                    }
                                }
                                try {
                                    float[] floatSimilarities = new float[similarity.size()];
                                    for (int i = 0; i < similarity.size(); i++) {
                                        floatSimilarities[i] = similarity.get(i).floatValue();
                                    }
                                    datapointsSpectra.put(model, floatSimilarities);
                                } catch (Exception e) {}
                            }
                        }
                    }

                    //get best model
                    BestModelSearcher bms = new BestModelSearcher();
                    String bestModel = "";
                    if (Constants.spectraBestModelSearchMetric.equals("median")) {
                        printInfo("Choosing spectra model based on median method");
                        bestModel = bms.medianMethod(datapointsSpectra, true);
                    } else if (Constants.spectraBestModelSearchMetric.equals("bottom")) {
                        printInfo("Choosing spectra model based on bottom 100 consensus method");
                        bestModel = bms.consensus(datapointsSpectra, 100, true, false);
                    } else {
                        printError("Unknown spectra model search metric: " + Constants.spectraBestModelSearchMetric);
                        printError("Exiting");
                        System.exit(1);
                    }

                    if (bestModel.contains("&")) { //model that uses NCE
                        String[] modelSplit = bestModel.split("&");
                        Constants.spectraModel = modelSplit[0];
                        NCEcalibrator.plotNCEchart(Constants.spectraModel, similarities.get(bestModel));
                    } else {
                        Constants.spectraModel = bestModel;
                    }
                    Constants.autoSwitchFragmentation = false;
                    printInfoHighlight("Spectra model chosen is " + Constants.spectraModel);

                    //write file that has all values for each model
                    bms.writeScores(Constants.outputDirectory + File.separator + "bestms2model.tsv");
                }

                if (Constants.spectraModel.isEmpty()) {
                    Constants.spectraModel = "DIA-NN";
                }
                for (String model : Constants.KoinaMS2models) {
                    if (model.equalsIgnoreCase(Constants.spectraModel)) {
                        Constants.spectraModel = model;
                    }
                }

                if (Constants.spectraModel.equals(Constants.auxSpectraModel)) {
                    printInfo("Spectra and aux spectra model are the same. Only predicting once.");
                    Constants.auxSpectraModel = "";
                }
            } else {
                Constants.spectraModel = "";
                Constants.auxSpectraModel = "";
            }
            Constants.foundBest = true;

            //check that at least pinPepXMLDirectory and mzmlDirectory are provided
            if (Constants.pinPepXMLDirectory == null) {
                throw new IllegalArgumentException("pinPepXMLDirectory must be provided");
            }
            if (Constants.mzmlDirectory == null) {
                throw new IllegalArgumentException("mzmlDirectory must be provided");
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
                    if (intersection.isEmpty()) {
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
                    if (intersection.isEmpty()) {
                        featureLL.add("deltaRTLOESS");
                        featureLL.add("predRTrealUnits");
                    }
                } else {
                    featureLL.removeIf(Constants.rtFeatures::contains);
                }
            } catch (Exception ignored) {
            }

            try {
                if (Constants.useIM) {
                    Set<String> intersection = new HashSet<>(featureLL);
                    intersection.retainAll(Constants.imFeatures);
                    if (intersection.isEmpty()) {
                        featureLL.add("deltaIMLOESS");
                        featureLL.add("ionmobility");
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
                    if (intersection.isEmpty()) {
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
                    if (intersection.isEmpty()) {
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
                    if (intersection.isEmpty()) {
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
                    if (intersection.isEmpty()) {
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

            //create file for spectral and RT prediction
            //ignore if files already created
            boolean createSpectraRTPredFile = false;

            //check which ones we need
            if (!featureLL.isEmpty()) {
                createSpectraRTPredFile = true;
            }

            //if pred file ready. Assumes that either all or none of pred files are ready
            if (Constants.spectraPredFile != null || Constants.RTPredFile != null || Constants.IMPredFile != null) {
                createSpectraRTPredFile = false;
            }

            c.updateInputPaths(); //setting null paths

            //set useKoina based on model
            if (Constants.KoinaRTmodels.contains(Constants.rtModel) ||
                    Constants.KoinaMS2models.contains(Constants.spectraModel) ||
                    Constants.KoinaIMmodels.contains(Constants.imModel)) {
                Constants.useKoina = true;
            }

            //go through all models and do any preprocessing needed
            //models needs to be in this order, as other code relies on this order, like spectra before aux spectra
            ArrayList<String> models = new ArrayList<>();
            ArrayList<String> modelTypes = new ArrayList<>();
            if (! Constants.spectraModel.isEmpty()) {
                KoinaMethods.switchModel(); //make sure using right model for fragmentation type
                models.add(Constants.spectraModel);
                modelTypes.add("spectra");
            }
            if (! Constants.rtModel.isEmpty()) {
                models.add(Constants.rtModel);
                modelTypes.add("RT");
            }
            if (! Constants.imModel.isEmpty()) {
                models.add(Constants.imModel);
                modelTypes.add("IM");
            }
            if (! Constants.auxSpectraModel.isEmpty()) {
                models.add(Constants.auxSpectraModel);
                modelTypes.add("auxSpectra");
            }

            //generate input files for prediction models
            MyFileUtils.deleteWholeDirectory(Constants.outputDirectory + File.separator + "NCE_calibration");
            if (createSpectraRTPredFile || Constants.createPredFileOnly) {
                //createfull is needed for everything
                PeptideFileCreator.createPeptideFile(pmMatcher,
                        Constants.spectraRTPrefix + "_full.tsv",
                        "createFull");

                HashSet<String> modelsRan = new HashSet<>();
                for (String currentModel : models) {
                    if (! modelsRan.contains(currentModel)) { //do not rerun (e.g. spectra and RT same model)
                        if (Constants.KoinaModels.contains(currentModel)) {
                            if (NceConstants.nceModels.contains(currentModel) && NceConstants.calibrateNCE &&
                                    ! NceConstants.calibratedModels.containsKey(currentModel)) {
                                //set matching with Da or not
                                //if matchWithDaltons are true, also accept that (e.g. from reading it from fragger.params)
                                if (currentModel.equalsIgnoreCase("predfull")) {
                                    Constants.matchWithDaltonsDefault = true;
                                } else {
                                    if (Constants.matchWithDaltons == null) {
                                        Constants.matchWithDaltonsDefault = false;
                                    } else {
                                        Constants.matchWithDaltonsDefault = Constants.matchWithDaltons;
                                    }
                                }

                                Object[] modelInfo = NCEcalibrator.calibrateNCE(currentModel, km,
                                        Constants.outputDirectory + File.separator + "NCE_calibration", true);
                                for (PeptideFormatter pf : km.peptideArraylist) {
                                    pf.foundUnimods.clear();
                                }

                                String bestNCE = String.valueOf((int) modelInfo[2]);
                                NceConstants.calibratedModels.put(currentModel, bestNCE);
                                NCEcalibrator.plotNCEchart(currentModel, (TreeMap<Integer, ArrayList<Double>>) modelInfo[0]);
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
                                    System.exit(1);
                            }
                        }
                        modelsRan.add(currentModel);
                    }
                }

                if (Constants.createPredFileOnly) {
                    printInfo("Successfully created input file for prediction model. Stopping here");
                    System.exit(0);
                }
            }

            //generate predictions
            //send input files to prediction models
            ArrayList<KoinaLibReader> klrs = new ArrayList<>();
            String DiannPredFilePath = "";
            ArrayList<String> predFilePaths = new ArrayList<>(); //replace "koina" with final name later
            if (createSpectraRTPredFile) {
                boolean ranKoina = false;
                for (String currentModel : models) {
                    KoinaModelCaller kmc = new KoinaModelCaller();
                    if (Constants.KoinaModels.contains(currentModel)) {
                        KoinaLibReader klr = new KoinaLibReader(currentModel);
                        kmc.callModel(currentModel, klr, Constants.JsonDirectory, executorService,
                                true, true);
                        ranKoina = true;
                        predFilePaths.add("koina" + currentModel);
                        klrs.add(klr);
                    } else {
                        if (DiannPredFilePath.isEmpty()) {
                            DiannPredFilePath = DiannModelCaller.callModel(
                                    Constants.spectraRTPrefix + ".tsv", true);
                        }
                        predFilePaths.add(DiannPredFilePath);
                    }
                }
                PredictionEntryHashMap koinaPreds = new PredictionEntryHashMap();
                if (ranKoina) {
                    koinaPreds.transferKoinaPreds(klrs, Constants.spectraRTPrefix + "_full.tsv");
                }

                //filter out fragment ion types not considered by fragmentation type before writing predictions
                koinaPreds.preprocessPredictedSpectra(executorService,
                        FragmentIonConstants.fragmentIonHierarchySet, FragmentIonConstants.fragmentIonHierarchySet);

                StringBuilder koinaPredFilePath = new StringBuilder("koina.mgf");
                StringBuilder auxPredFilePath = new StringBuilder("koina.mgf");
                for (int j = 0; j < models.size(); j++) {
                    switch (modelTypes.get(j)) {
                        case "spectra":
                            if (predFilePaths.get(j).startsWith("koina")) {
                                koinaPredFilePath.insert(0, "spectra-" +
                                        predFilePaths.get(j).substring(5) + ".");
                            } else { //DIANN
                                Constants.spectraPredFile = DiannPredFilePath;
                            }
                            break;
                        case "RT":
                            if (predFilePaths.get(j).startsWith("koina")) {
                                koinaPredFilePath.insert(0, "RT-" +
                                        predFilePaths.get(j).substring(5) + ".");
                            } else {
                                Constants.RTPredFile = DiannPredFilePath;
                            }
                            break;
                        case "IM":
                            if (predFilePaths.get(j).startsWith("koina")) {
                                koinaPredFilePath.insert(0, "IM-" +
                                        predFilePaths.get(j).substring(5) + ".");
                            } else {
                                Constants.IMPredFile = DiannPredFilePath;
                            }
                            break;
                        case "auxSpectra": //DIANN cannot predict this
                            auxPredFilePath.insert(0, "auxSpectra-" +
                                    predFilePaths.get(j).substring(5) + ".");
                            break;
                    }
                }

                if (!koinaPredFilePath.toString().equals("koina.mgf")) { //koina was used
                    MgfFileWriter mfw = new MgfFileWriter(koinaPreds);
                    koinaPredFilePath.insert(0, Constants.outputDirectory + File.separator);
                    mfw.write(koinaPredFilePath.toString());
                    for (int j = 0; j < models.size(); j++) {
                        if (predFilePaths.get(j).startsWith("koina")) {
                            switch (modelTypes.get(j)) {
                                case "spectra":
                                    Constants.spectraPredFile = koinaPredFilePath.toString();
                                    break;
                                case "RT":
                                    Constants.RTPredFile = koinaPredFilePath.toString();
                                    break;
                                case "IM":
                                    Constants.IMPredFile = koinaPredFilePath.toString();
                                    break;
                            }
                        }
                    }
                }
                if (!auxPredFilePath.toString().equals("koina.mgf")) { //koina was used
                    //make mini prediction entry hashmap
                    PredictionEntryHashMap miniMap = new PredictionEntryHashMap();
                    for (Map.Entry<String, PredictionEntry> pe : koinaPreds.entrySet()) {
                        miniMap.put(pe.getKey(), pe.getValue().auxSpectra);
                    }

                    MgfFileWriter mfwAux = new MgfFileWriter(miniMap);
                    auxPredFilePath.insert(0, Constants.outputDirectory + File.separator);
                    mfwAux.write(auxPredFilePath.toString());
                    Constants.auxSpectraPredFile = auxPredFilePath.toString();
                }
            }

            //create new pin file with features
            printInfo("Generating edited pin with following features: " + Arrays.toString(featuresArray));
            long start = System.nanoTime();
            PercolatorFormatter.editPin(pmMatcher, featuresArray, Constants.editedPinSuffix, executorService);
            executorService.shutdown();

            //print parameters to ps
            //printParamsPS();

            //delete pred files
            if (Constants.deletePreds) {
                //DIA-NN input files
                File predFile = new File(Constants.outputDirectory + File.separator + "spectraRT.tsv");
                predFile.delete();
                predFile = new File(Constants.outputDirectory + File.separator + "spectraRT_full.tsv");
                predFile.delete();

                //prediction files
                if (Constants.spectraPredFile != null) {
                    predFile = new File(Constants.spectraPredFile);
                    predFile.delete();
                }
                if (Constants.RTPredFile != null) {
                    predFile = new File(Constants.RTPredFile);
                    predFile.delete();
                }
                if (Constants.IMPredFile != null) {
                    predFile = new File(Constants.IMPredFile);
                    predFile.delete();
                }
                if (Constants.auxSpectraPredFile != null) {
                    predFile = new File(Constants.auxSpectraPredFile);
                    predFile.delete();
                }
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
                AtomicBoolean imfound = new AtomicBoolean(false);
                AtomicBoolean spectraPredFilePDVFound = new AtomicBoolean(false);

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
                            } else if (line.trim().startsWith("imModel")) {
                                imfound.set(true);
                                return "imModel=" + Constants.imModel;
                            } else if (line.trim().startsWith("spectraPredFilePDV")) {
                                spectraPredFilePDVFound.set(true);
                                return "spectraPredFilePDV=" + Constants.spectraPredFile;
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
                if (!imfound.get()) {
                    modifiedLines.add("imModel=" + Constants.imModel);
                }
                //FragPipe PDV needs to know which file has spectral predictions
                if (Constants.useSpectra && !spectraPredFilePDVFound.get()) {
                    modifiedLines.add("spectraPredFilePDV=" + Constants.spectraPredFile);
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

    //for updating constants read from parameter file
    static private void updateField(Field field, String value, ConstantsInterface c)
            throws NoSuchMethodException, InvocationTargetException, InstantiationException, IllegalAccessException {
        //get class of field
        Class<?> myClass = field.getType();

        //do not parse use[something] if null
        if (myClass.getTypeName().equals("java.lang.Boolean")) {
            if (value.equals("null")) {
                return;
            }
        }

        //parse to appropriate type
        field.set(c, myClass.getConstructor(String.class).newInstance(value));
    }
}
