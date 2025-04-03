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

package koinaclasses;

import allconstants.Constants;
import com.google.common.collect.ImmutableMap;
import mainsteps.MzmlScanNumber;
import mainsteps.PeptideObj;
import mainsteps.PinMzmlMatcher;
import modelcallers.KoinaModelCaller;
import peptideptmformatting.PeptideFormatter;
import peptideptmformatting.PeptideSkipper;
import predictions.PredictionEntryHashMap;
import readers.datareaders.MzmlReader;
import readers.datareaders.PinReader;
import readers.predictionreaders.KoinaLibReader;
import umich.ms.fileio.exceptions.FileParsingException;
import writers.JSONWriter;

import java.io.File;
import java.io.IOException;
import java.net.URISyntaxException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.LinkedList;
import java.util.Map;
import java.util.concurrent.ExecutionException;
import java.util.concurrent.ScheduledThreadPoolExecutor;

import static utils.InstrumentUtils.mapInstrumentToModelSpecific;
import static utils.Print.printInfo;

public class KoinaMethods {
    //these fields are shared regardless of which model is called
    public PinMzmlMatcher pmMatcher;
    public ArrayList<PeptideFormatter> peptideArraylist = new ArrayList<>();
    public HashMap<String, LinkedList<Integer>> scanNums = new HashMap<>();
    public HashMap<String, LinkedList<PeptideFormatter>> peptides = new HashMap<>();
    public ArrayList<PeptideFormatter> peptideArrayListIM = new ArrayList<>();
    public HashMap<String, LinkedList<Integer>> scanNumsIM = new HashMap<>();
    public HashMap<String, LinkedList<PeptideFormatter>> peptidesIM = new HashMap<>();

    //decoys
//    public HashSet<String> peptideSetDecoys = new HashSet<>();
//    public HashMap<String, LinkedList<Integer>> scanNumsDecoys = new HashMap<>();
//    public HashMap<String, LinkedList<String>> peptidesDecoys = new HashMap<>();

    private static final Map<String, String> generalModels = ImmutableMap.of(
            "Prosit_2020_intensity_HCD", "Prosit",
            "Prosit_2020_intensity_CID", "Prosit");
    private static final Map<String, String> modelConversion = ImmutableMap.of(
            "HCD.Prosit", "Prosit_2020_intensity_HCD",
            "CID.Prosit", "Prosit_2020_intensity_CID");

    public KoinaMethods(PinMzmlMatcher pmMatcher) {
        this.pmMatcher = pmMatcher;
    }

    public void getTopPeptides() throws IOException {
        //need to collect top 1000 peptides for calibration
        //approximate by doing a subset per pin
        int numTopPSMs = (int) Math.ceil((float) Constants.numPSMsToCalibrate /
                (float) pmMatcher.pinFiles.length);

        for (int j = 0; j < pmMatcher.pinFiles.length; j++) {
            File pinFile = pmMatcher.pinFiles[j];
            PinReader pinReader = new PinReader(pinFile.getAbsolutePath());
            LinkedList[] topPSMs = pinReader.getTopPSMs(numTopPSMs, false);
            peptideArraylist.addAll(topPSMs[0]);
            scanNums.put(pmMatcher.mzmlFiles[j].getName(), topPSMs[1]);
            peptides.put(pmMatcher.mzmlFiles[j].getName(), topPSMs[0]);
            if (Constants.useIM) {
                LinkedList[] topPSMsIM = pinReader.getTopPSMs(numTopPSMs, true);
                peptideArrayListIM.addAll(topPSMsIM[0]);
                scanNumsIM.put(pmMatcher.mzmlFiles[j].getName(), topPSMsIM[1]);
                peptidesIM.put(pmMatcher.mzmlFiles[j].getName(), topPSMsIM[0]);
            }
        }
    }

//    public void getDecoyPeptides() throws IOException {
//        //need to collect top 1000 peptides for calibration
//        //approximate by doing a subset per pin
//        int numTopPSMs = (int) Math.ceil((float) Constants.numPSMsToCalibrate /
//                (float) pmMatcher.pinFiles.length);
//
//        for (int j = 0; j < pmMatcher.pinFiles.length; j++) {
//            File pinFile = pmMatcher.pinFiles[j];
//            PinReader pinReader = new PinReader(pinFile.getAbsolutePath());
//            LinkedList[] decoyPSMs = pinReader.getDecoyPSMs(numTopPSMs);
//            peptideSetDecoys.addAll(decoyPSMs[0]);
//            scanNumsDecoys.put(pmMatcher.mzmlFiles[j].getName(), decoyPSMs[1]);
//            peptidesDecoys.put(pmMatcher.mzmlFiles[j].getName(), decoyPSMs[0]);
//        }
//    }

    //used for partial prediction. allHits is only a subset of all the hits in a pin file
    public PredictionEntryHashMap getKoinaPredictions(
            HashSet<String> allHits, String model, int NCE, String folder, String fulltsv) {
        ScheduledThreadPoolExecutor executorService = new ScheduledThreadPoolExecutor(1);

        HashSet<String> hits = new HashSet<>();
        for (String s : allHits) { //s is peptide,charge
            String instrument = mapInstrumentToModelSpecific(model, s.split(",")[1]);
            hits.add(s + "," + NCE + "," + instrument + "," + Constants.FragmentationType);
        }

        JSONWriter jw = new JSONWriter(model, hits, false);

        String jsonFolder = "";
        try {
            jsonFolder = jw.write(true, folder, executorService);
        } catch (Exception e) {
            e.printStackTrace();
            System.exit(1);
        }

        //send predictions to Koina
        KoinaLibReader klr = new KoinaLibReader(model);
        KoinaModelCaller kmc = new KoinaModelCaller();
        kmc.callModel(model, klr, jsonFolder, executorService, false, false);
        executorService.shutdown();
        try {
            PredictionEntryHashMap assignedPreds = new PredictionEntryHashMap();
            ArrayList<KoinaLibReader> klrs = new ArrayList<>();
            klrs.add(klr);
            assignedPreds.transferKoinaPreds(klrs, fulltsv);
            klr.setPreds(assignedPreds);
        } catch (Exception e) {
            e.printStackTrace();
        }
        return klr.allPreds;
    }

    public PeptideObj[] getPeptideObjects(PredictionEntryHashMap allPreds,
                                          HashMap<String, LinkedList<Integer>> scanNums,
                                          HashMap<String, LinkedList<PeptideFormatter>> peptides)
            throws FileParsingException, ExecutionException, InterruptedException, IOException, URISyntaxException {
        //pass empty hashsets since no fragment filtering at this step
        allPreds.preprocessPredictedSpectra(new ScheduledThreadPoolExecutor(Constants.numThreads), new HashSet<>(), new HashSet<>());

        int arrayLength = 0;
        for (LinkedList<Integer> scanNum : scanNums.values()) {
            arrayLength += scanNum.size();
        }
        PeptideObj[] peptideObjs = new PeptideObj[arrayLength];

        int index = 0;
        for (int j = 0; j < pmMatcher.mzmlReaders.length; j++) {
            MzmlReader mzmlReader = pmMatcher.mzmlReaders[j];
            LinkedList<Integer> thisScanNums = scanNums.get(pmMatcher.mzmlFiles[j].getName());
            LinkedList<PeptideFormatter> thisPeptides = peptides.get(pmMatcher.mzmlFiles[j].getName());

            for (int k = 0; k < thisScanNums.size(); k++) {
                int scanNum = thisScanNums.get(k);
                MzmlScanNumber msn = mzmlReader.getScanNumObject(scanNum);

                PeptideFormatter pf = thisPeptides.get(k);
                if (allPreds.containsKey(pf.getBaseCharge())) {
                    PeptideObj pobj = msn.setPeptideObject(
                            new PeptideFormatter(pf.getBase(), pf.getCharge(), "base"),
                            1, 1, "0", allPreds, false);
                    peptideObjs[index] = pobj;
                }
                index++;
            }
        }
        return peptideObjs;
    }

    //use nested hashmap to make sure correct version of model is assigned for given fragmentation type
    //uses Constants information (FragmentationType and current spectraModel) to correct spectraModel
    //solution: use composite key fragmentationtype + "." + modelType, where modelType is a HashMap of
    //exact model name and more general type (i.e. Prosit models, APD models, etc)
    //if hashset contains exact model name, return that general type
    //returns true if model changed
    public static boolean switchModel() {
        //correct to CID model
        //TODO do in opposite direction. Or more general method to get right fragmentation model
        if (Constants.autoSwitchFragmentation && generalModels.containsKey(Constants.spectraModel)) {
            String genModel = generalModels.get(Constants.spectraModel);
            String newModel = modelConversion.get(Constants.FragmentationType + "." + genModel);
            if (!newModel.equals(Constants.spectraModel)) {
                printInfo("Switching from " + Constants.spectraModel + " to " + newModel);
                Constants.spectraModel = newModel;
                return true;
            }
        }
        return false;
    }

    public static HashSet<String> createPartialKoinaSet(String currentModel, ArrayList<PeptideFormatter> peptideFormatters) {
        HashSet<String> allHits = new HashSet<>();
        for (PeptideFormatter pf : peptideFormatters) {
            if (PeptideSkipper.skipPeptide(pf.getStripped(), pf.getCharge(), currentModel)) {
                continue;
            }

            //Koina requires all uppercase UNIMOD
            String pep = pf.getModel(currentModel);

            allHits.add(pep + "," + pf.getCharge());
        }

        return allHits;
    }
}
