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

import allconstants.*;

import koinaclasses.KoinaMethods;

import utils.Model;
import utils.MyFileUtils;

import java.io.*;

import java.util.*;

import java.util.concurrent.ScheduledThreadPoolExecutor;

import static mainsteps.MainUtils.*;
import static mainsteps.ParameterUtils.*;
import static utils.Print.*;

//this is what I use in the java jar file
public class MainClass {
    public static ScheduledThreadPoolExecutor executorService;

    public static void main(String[] args) throws Exception {
        Locale.setDefault(Locale.US);
        printInfo("MSBooster v1.3.23");

        try {
            HashMap<String, String> params = new HashMap<>();

            //accept command line inputs
            processCommandLineInputs(args, params);

            //adding parameter list
            if (params.containsKey("paramsList")) { //override previous params input
                processParameterList(params);
            }

            //read in MSFragger parameters
            if (params.containsKey("fragger")) { //Will override paramsList
                if (!params.get("fragger").equals("null")) {
                    readFraggerParams(params);
                }
            }

            //adding to different constants classes, and setting input/output paths
            updatesConstants(params);

            //defining num threads
            if (Constants.numThreads <= 0) {
                Runtime run = Runtime.getRuntime();
                Constants.numThreads = run.availableProcessors() - 1;
            } //otherwise use user-defined
            printInfo("Using " + Constants.numThreads + " threads");
            executorService = new ScheduledThreadPoolExecutor(Constants.numThreads);

            //check that at least pinPepXMLDirectory and mzmlDirectory are provided
            if (Constants.pinPepXMLDirectory == null) {
                throw new IllegalArgumentException("pinPepXMLDirectory must be provided");
            }
            if (Constants.mzmlDirectory == null) {
                throw new IllegalArgumentException("mzmlDirectory must be provided");
            }
            //get matched pin files for mzML files
            PinMzmlMatcher pmMatcher = new PinMzmlMatcher(Constants.mzmlDirectory, Constants.pinPepXMLDirectory);

            //checking plot output is fine extension
            HashSet<String> extensions = new HashSet<>(Arrays.asList("png", "pdf"));
            if (!extensions.contains(Constants.plotExtension.toLowerCase())) {
                printError("plotExtension must be one of png or pdf, not " + Constants.plotExtension);
                printError("Exiting");
                System.exit(1);
            }

            //calculate actual parts per million
            if (Constants.ppmTolerance <= 0f) { //not manually set
                Constants.ppmTolerance = Constants.highResppmTolerance;
            }
            Constants.lowResppmTolerance = Constants.lowResppmTolerance / 1000000f;
            Constants.highResppmTolerance = Constants.highResppmTolerance / 1000000f;
            Constants.ppmTolerance = Constants.ppmTolerance / 1000000f;

            //checking constants are appropriate
            //robust to if url has / or not
            if (!Constants.KoinaURL.isEmpty() &&
                    Constants.KoinaURL.charAt(Constants.KoinaURL.length() - 1) != '/') {
                Constants.KoinaURL += "/";
            }

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

            //update fragment ion types based on fragmentation type
            FragmentIonConstants.makeFragmentIonHierarchy();
            Constants.matchedIntensitiesFeatures = Constants.makeMatchedIntensitiesFeatures();
            Constants.peakCountsFeatures = Constants.makePeakCountsFeatures();
            Constants.predIntensitiesFeatures = Constants.makePredIntensitiesFeatures();
            Constants.individualSpectralSimilaritiesFeatures = Constants.makeIndividualSpectralSimilarities();
            Constants.intensitiesDifferenceFeatures = Constants.makeintensitiesDifference();

            //exit if no models used
            if (!Constants.useSpectra && !Constants.useRT && !Constants.useIM) {
                printInfo("useSpectra, useRT, and useIM are all false. Exiting.");
                System.exit(0);
            }

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
            LowercaseModelMapper lmm = new LowercaseModelMapper();
            HashMap<String, String> modelMapper = lmm.getLowercaseToModel();
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

            //set models
            KoinaMethods km = new KoinaMethods(pmMatcher);
            if (Constants.findBestRtModel || Constants.findBestSpectraModel || Constants.findBestImModel ||
                    ModelCollections.KoinaMS2models.contains(Constants.spectraModel) ||
                    ModelCollections.KoinaRTmodels.contains(Constants.rtModel) ||
                    ModelCollections.KoinaIMmodels.contains(Constants.imModel) ||
                    ! Constants.auxSpectraModel.isEmpty()
            ) {
                km.getTopPeptides();
            }
            ArrayList<Model> models = new ArrayList<>();
            Model spectraModel = setMS2model(km, executorService);
            if (!spectraModel.name.isEmpty()) {
                models.add(spectraModel);
            }
            Model rtModel = setRTmodel(km, pmMatcher, executorService);
            if (!rtModel.name.isEmpty()) {
                models.add(rtModel);
            }
            Model imModel = setIMmodel(km, pmMatcher, executorService);
            if (!imModel.name.isEmpty()) {
                models.add(imModel);
            }
            Model auxSpectraModel = setAuxMS2model();
            if (!auxSpectraModel.name.isEmpty()) {
                models.add(auxSpectraModel);
            }
            Constants.foundBest = true;
            if (ModelCollections.KoinaRTmodels.contains(Constants.rtModel) ||
                    ModelCollections.KoinaMS2models.contains(Constants.spectraModel) ||
                    ModelCollections.KoinaIMmodels.contains(Constants.imModel)) {
                Constants.useKoina = true;
            }

            //update features
            String[] featuresArray = updateFeatures();

            //create file for spectral and RT prediction
            //ignore if files already created
            boolean createSpectraRTPredFile = false;
            if (featuresArray.length > 0) {
                createSpectraRTPredFile = true;
            }
            //if pred file ready, skip pred file creation. Assumes that either all or none of pred files are ready
            if (Constants.spectraPredFile != null || Constants.RTPredFile != null || Constants.IMPredFile != null) {
                createSpectraRTPredFile = false;
            }

            //generate input files for prediction models
            if (createSpectraRTPredFile || Constants.createPredFileOnly) {
                makeInputFiles(pmMatcher, models, km);
                if (Constants.createPredFileOnly) {
                    printInfo("Successfully created input file for prediction model. Stopping here");
                    System.exit(0);
                }
            }

            //generate predictions
            if (createSpectraRTPredFile) {
                getPredictionFiles(models, executorService);
            }

            //create new pin file with features
            printInfo("Generating edited pin with following features: " + Arrays.toString(featuresArray));
            long start = System.nanoTime();
            PercolatorFormatter.editPin(pmMatcher, featuresArray, Constants.editedPinSuffix, executorService);
            long end = System.nanoTime();
            long duration = (end - start);
            printInfo("Feature calculation, edited pin writing, and QC plot generation done in " +
                    duration / 1000000 + " ms");
            executorService.shutdown();

            //delete pred files
            if (Constants.deletePreds) {
                deletePredFiles();
            }

            //clean directories
            MyFileUtils.deleteWholeDirectory(Constants.outputDirectory + File.separator + "best_model");
            MyFileUtils.deleteWholeDirectory(Constants.outputDirectory + File.separator + "NCE_calibration");
            MyFileUtils.deleteWholeDirectory(Constants.JsonDirectory);

            //finalize parameter file with final selected models, citations, and spectral file for PDV
            finalizeParameterFile(models);

            System.exit(0);
        } catch (Exception e) {
            e.printStackTrace();
            executorService.shutdown();
            System.exit(1);
        }
    }
}
