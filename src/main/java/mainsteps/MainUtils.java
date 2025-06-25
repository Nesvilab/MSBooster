package mainsteps;

import allconstants.*;
import bestmodelsearch.BestModelSearcher;
import koinaclasses.KoinaMethods;
import koinaclasses.NCEcalibrator;
import modelcallers.DiannModelCaller;
import modelcallers.KoinaModelCaller;
import peptideptmformatting.PeptideFormatter;
import predictions.PredictionEntry;
import predictions.PredictionEntryHashMap;
import readers.datareaders.MzmlReader;
import readers.predictionreaders.KoinaLibReader;
import readers.predictionreaders.LibraryPredictionMapper;
import umich.ms.fileio.exceptions.FileParsingException;
import utils.Model;
import utils.MyFileUtils;
import writers.MgfFileWriter;
import writers.PeptideFileCreator;

import java.io.*;
import java.util.*;
import java.util.concurrent.ExecutionException;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.ScheduledThreadPoolExecutor;

import static bestmodelsearch.ModelCollectionDecider.decideCollection;
import static features.rtandim.LoessUtilities.gridSearchCV;
import static utils.Print.*;

public class MainUtils {
    static Model setRTmodel(KoinaMethods km, PinMzmlMatcher pmMatcher, ExecutorService executorService) throws Exception {
        if (Constants.useRT) {
            //here, look for best rt model
            if (Constants.findBestRtModel) {
                HashMap<String, float[]> datapointsRT = new HashMap<>();
                ArrayList<String> consideredModels = decideCollection("RT");

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
            for (String model : ModelCollections.KoinaRTmodels) {
                if (model.equalsIgnoreCase(Constants.rtModel)) {
                    Constants.rtModel = model;
                }
            }
        } else {
            Constants.rtModel = "";
        }
        return new Model(Constants.rtModel, "RT");
    }

    static Model setIMmodel(KoinaMethods km, PinMzmlMatcher pmMatcher, ExecutorService executorService) throws Exception {
        if (Constants.useIM) {
            if (Constants.findBestImModel) {
                HashMap<String, float[]> datapointsIM = new HashMap<>();
                ArrayList<String> consideredModels = decideCollection("IM");

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
            for (String model : ModelCollections.KoinaIMmodels) {
                if (model.equalsIgnoreCase(Constants.imModel)) {
                    Constants.imModel = model;
                }
            }
        } else {
            Constants.imModel = "";
        }
        return new Model(Constants.imModel, "IM");
    }

    static Model setMS2model(KoinaMethods km, ExecutorService executorService) throws Exception {
        if (Constants.useSpectra) {
            //here, look for best spectra model
            if (Constants.findBestSpectraModel) {
                ArrayList<String> consideredModels = decideCollection("MS2");

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
            for (String model : ModelCollections.KoinaMS2models) {
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
        return new Model(Constants.spectraModel, "spectra");
    }

    //if more full spectrum models in future, can adapt
    static Model setAuxMS2model() {
        return new  Model(Constants.auxSpectraModel, "auxSpectra");
    }

    static String[] updateFeatures() {
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
        if (Constants.useSpectra) {
            Set<String> intersection = new HashSet<>(featureLL);
            intersection.retainAll(Constants.spectraFeatures);
            if (intersection.isEmpty()) {
                featureLL.add("unweightedSpectralEntropy");
                featureLL.add("weightedSpectralEntropy");
                featureLL.add("hypergeometricProbability");
                featureLL.add("intersection");
            }
        } else {
            featureLL.removeIf(Constants.spectraFeatures::contains);
        }
        if (Constants.useRT) {
            Set<String> intersection = new HashSet<>(featureLL);
            intersection.retainAll(Constants.rtFeatures);
            if (intersection.isEmpty()) {
                featureLL.add("deltaRTLOESS");
                featureLL.add("deltaRTLOESSreal");
            }
        } else {
            featureLL.removeIf(Constants.rtFeatures::contains);
        }
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
        if (Constants.useIntensityDistributionSimilarity) {
            if (!featureLL.contains("intensity_distribution_similarity")) {
                featureLL.add("intensity_distribution_similarity");
            }
        } else {
            featureLL.remove("intensity_distribution_similarity");
        }

        //update features representation
        featuresArray = new String[featureLL.size()];
        int i = 0;
        for (String feature : featureLL) {
            featuresArray[i] = feature;
            i++;
        }
        Constants.features = String.join(",", featuresArray);

        return featuresArray;
    }

    static void makeInputFiles(PinMzmlMatcher pmMatcher, ArrayList<Model> models, KoinaMethods km)
            throws FileParsingException, IOException, ExecutionException, InterruptedException {
        MyFileUtils.deleteWholeDirectory(Constants.outputDirectory + File.separator + "NCE_calibration");
        //createfull is needed for everything
        PeptideFileCreator.createPeptideFile(pmMatcher,
                Constants.spectraRTPrefix + "_full.tsv",
                "createFull");

        HashSet<String> modelsRan = new HashSet<>();
        for (Model model : models) {
            String currentModel = model.name;
            if (! modelsRan.contains(currentModel)) { //do not rerun (e.g. spectra and RT same model)
                if (ModelCollections.KoinaModels.contains(currentModel)) {
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
    }

    static void getPredictionFiles(ArrayList<Model> models, ScheduledThreadPoolExecutor executorService) throws Exception {
        ArrayList<KoinaLibReader> klrs = new ArrayList<>();
        String DiannPredFilePath = "";
        ArrayList<String> predFilePaths = new ArrayList<>(); //replace "koina" with final name later

        boolean ranKoina = false;
        for (Model model : models) {
            String currentModel = model.name;
            KoinaModelCaller kmc = new KoinaModelCaller();
            if (ModelCollections.KoinaModels.contains(currentModel)) {
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
            switch (models.get(j).modelType) {
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
                    switch (models.get(j).modelType) {
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

    static void deletePredFiles() {
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
}
