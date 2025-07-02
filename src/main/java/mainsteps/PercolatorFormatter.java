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
import allconstants.FragmentIonConstants;
import features.FeatureCalculator;
import features.rtandim.IMFunctions;
import figures.*;
import kotlin.jvm.functions.Function1;
import predictions.PredictionEntry;
import predictions.PredictionEntryHashMap;
import readers.MgfFileReader;
import readers.datareaders.MzmlReader;
import readers.datareaders.PinReader;
import readers.predictionreaders.LibraryPredictionMapper;
import utils.MyFileUtils;
import writers.PinWriter;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.*;
import java.util.concurrent.ExecutorService;

import static allconstants.Constants.figureDirectory;
import static allconstants.FragmentIonConstants.createFragmentGroups;
import static utils.Print.printError;
import static utils.Print.printInfo;

public class PercolatorFormatter {

    public static PredictionEntryHashMap allPreds = new PredictionEntryHashMap();
    static HashMap<String, TreeMap<String, List<Float>[]>> RTregressionCurves = new HashMap<>(); //<mass, <mzml name, [x, y]>>
    public static HashMap<String, List<Double>> scoreDescriptors = new HashMap<>(); //<mzml, scores for PSMs below evalue>
    public static HashMap<String, HashMap<String, double[][]>> absDeltaRTs = new HashMap<>();

    public static void editPin(PinMzmlMatcher pmMatcher, String[] features, String outfile,
                               ExecutorService executorService) throws Exception {

        ArrayList<String> featuresList = new ArrayList<>(Arrays.asList(features));

        File[] pinFiles = pmMatcher.pinFiles;
        File[] mzmlFiles = pmMatcher.mzmlFiles;
        MzmlReader[] mzmlReaders = pmMatcher.mzmlReaders;

        //load predicted spectra
        LibraryPredictionMapper predictedSpectra;
        LibraryPredictionMapper predictedRT;
        LibraryPredictionMapper predictedIM;
        LibraryPredictionMapper predictedAuxSpectra;

        ArrayList<String> libraryFilePaths = new ArrayList<>();
        ArrayList<LibraryPredictionMapper> libraries = new ArrayList<>();
        HashMap<String, String> allProperties = new HashMap<>(); //key: property, value: library file path

        //check that all needed files are here
        if (Constants.useSpectra && Constants.spectraPredFile == null) {
            printError("Spectral prediction file is missing. " +
                    "Please specify its path with the --spectraPredFile parameter. Exiting");
            System.exit(1);
        }
        if (Constants.useRT && Constants.RTPredFile == null) {
            printError("RT prediction file is missing. " +
                    "Please specify its path with the --RTPredFile parameter. Exiting");
            System.exit(1);
        }
        if (Constants.useIM && Constants.IMPredFile == null) {
            printError("IM prediction file is missing. " +
                    "Please specify its path with the --IMPredFile parameter. Exiting");
            System.exit(1);
        }

        //could use aux spectra if primary spectra missing
        if (Constants.spectraPredFile != null) {
            printInfo("Loading predicted spectra: " + Constants.spectraPredFile);
            predictedSpectra = LibraryPredictionMapper.createLibraryPredictionMapper(
                    Constants.spectraPredFile, Constants.spectraModel, executorService);
            libraryFilePaths.add(Constants.spectraPredFile);
            libraries.add(predictedSpectra);
            allProperties.put("spectra", Constants.spectraPredFile);

            //deciding which fragments to use for primary and aux models, as well as setting which fragments to use for multiple scores
            for (PredictionEntry pe : predictedSpectra.getPreds().values()) {
                FragmentIonConstants.setPrimaryAndAuxFragmentIonTypes(pe.fragmentIonTypes); //this can be approximated if too slow
            }
            FragmentIonConstants.fragmentGroups = createFragmentGroups();

            if (Constants.spectraModel.equals("PredFull")) {
                Constants.matchWithDaltons = true; //they report predictions in bins
                Constants.matchWithDaltonsDefault = true;
            } else if (Constants.matchWithDaltons == null) {
                Constants.matchWithDaltons = false;
                Constants.matchWithDaltonsDefault = false;
            }
        }

        if (Constants.RTPredFile != null) {
            printInfo("Loading predicted retention times: " + Constants.RTPredFile);
            if (! libraryFilePaths.contains(Constants.RTPredFile)) {
                predictedRT = LibraryPredictionMapper.createLibraryPredictionMapper(
                        Constants.RTPredFile, Constants.rtModel, executorService);
                libraryFilePaths.add(Constants.RTPredFile);
                libraries.add(predictedRT);
            }
            allProperties.put("RT", Constants.RTPredFile);
        }

        if (Constants.IMPredFile != null) {
            printInfo("Loading predicted ion mobilities: " + Constants.IMPredFile);
            if (! libraryFilePaths.contains(Constants.IMPredFile)) {
                predictedIM = LibraryPredictionMapper.createLibraryPredictionMapper(
                        Constants.IMPredFile, Constants.imModel, executorService);
                libraryFilePaths.add(Constants.IMPredFile);
                libraries.add(predictedIM);
            }
            allProperties.put("IM", Constants.IMPredFile);
        }

        if (Constants.auxSpectraPredFile != null) {
            printInfo("Loading predicted auxiliary spectra: " + Constants.auxSpectraPredFile);
            if (! libraryFilePaths.contains(Constants.auxSpectraPredFile)) {
                predictedAuxSpectra = LibraryPredictionMapper.createLibraryPredictionMapper(
                        Constants.auxSpectraPredFile, Constants.auxSpectraModel, executorService);
                libraryFilePaths.add(Constants.auxSpectraPredFile);
                libraries.add(predictedAuxSpectra);
            }
            allProperties.put("auxSpectra", Constants.auxSpectraPredFile);
            if (Constants.auxSpectraModel.equals("PredFull")) {
                Constants.matchWithDaltonsAux = true; //they report predictions in bins
            } else if (Constants.matchWithDaltonsAux == null) {
                Constants.matchWithDaltonsAux = false;
            }
        }

        //merging libraries
        if (libraries.size() > 1) {
            printInfo("Merging libraries");
            for (int i = 0; i < libraries.size(); i++) {
                String libraryPath = libraryFilePaths.get(i);
                LibraryPredictionMapper library = libraries.get(i);
                for (Map.Entry<String, String> prop : allProperties.entrySet()) {
                    if (prop.getValue().equals(libraryPath)) {
                        library.getPreds().mergeIntoLibrary(allPreds, prop.getKey());
                    }
                }
            }
        } else {
            for (LibraryPredictionMapper lpm : libraries) {
                allPreds = lpm.getPreds();
            }
        }

        try {
            //////////////////////////////iterate through pin and mzml files//////////////////////////////////////////
            for (int i = 0; i < pinFiles.length; i++) {
                String newOutfile = pinFiles[i].getAbsolutePath().replaceAll("\\.pin$", "_" + outfile + ".pin");

                //load mzml file
                MzmlReader mzml;
                if (mzmlFiles[i].getName().substring( mzmlFiles[i].getName().length() - 3).equalsIgnoreCase("mgf")) {
                    mzml = new MzmlReader(new MgfFileReader(mzmlFiles[i].getCanonicalPath(),
                            true, executorService, ""));
                } else {
                    if (mzmlReaders[i] == null) {
                        mzml = new MzmlReader(mzmlFiles[i].getCanonicalPath());
                    } else {
                        mzml = mzmlReaders[i];
                    }
                }

                //load pin file, which already includes all ranks
                PinReader pin = new PinReader(pinFiles[i].getCanonicalPath());

                //handling for empty pin files
                if (pin.getLength() < 2) {
                    printInfo(pinFiles[i].getCanonicalPath() + " is empty. Writing empty edited file");
                    PinWriter pw = new PinWriter(newOutfile, pin, featuresList, mzml);
                    pw.write();
                    continue;
                }

                printInfo("Processing pin " + pin.name);
                if (Constants.removeWashGradient) {
                    if (Constants.rtCutoff.isNaN()) {
                        pin.attachMzml(mzml);
                        pin.findWashGradientCutoff();
                    } else {
                        pin.rtCutoff = Constants.rtCutoff;
                    }
                }

                //Special preparations dependent on features we require
                mzml.setPinEntries(pin, allPreds, executorService);

                if (featuresList.contains("deltaRTLOESS") || featuresList.contains("deltaRTLOESSnormalized") ||
                        featuresList.contains("deltaRTLOESSreal")) {
                    mzml.setLOESS(Constants.rtLoessRegressionSize, Constants.loessBandwidth, Constants.robustIters, "RT");
                    mzml.setInverseLoess(executorService);
                    mzml.predictRTLOESS(executorService); //potentially only invoke once if normalized included
                    mzml.getBestCalibratedRTs();

                    //generate calibration figure, need mzml and loess
                    boolean plot = false;
                    for (Map.Entry<String, double[][]> entry : mzml.expAndPredRTs.entrySet()) {
                        if (entry.getValue() != null) {
                            plot = true;
                            break;
                        }
                    }
                    if (plot) {
                        CalibrationFigure cf = new RTCalibrationFigure(mzml, pinFiles[i].getCanonicalPath(),
                                Constants.loessScatterOpacity, mzml.expAndPredRTs, mzml.RTLOESS);

                        //save curve for final plot
                        for (Map.Entry<String, List<Float>[]> entry : cf.curves.entrySet()) {
                            TreeMap<String, List<Float>[]> newMap = new TreeMap<>();
                            String name = new File(mzml.pathStr).getName();
                            name = name.substring(0, name.length() - 5);
                            newMap.put(name, entry.getValue());

                            TreeMap<String, List<Float>[]> totalMap;
                            if (RTregressionCurves.containsKey(entry.getKey())) {
                                totalMap = RTregressionCurves.get(entry.getKey());
                                totalMap.put(name, entry.getValue());
                            } else {
                                totalMap = newMap;
                            }

                            RTregressionCurves.put(entry.getKey(), totalMap);
                        }

                        File f = new File(pinFiles[i].getCanonicalPath());
                        String calibrationPeptideFilePathBase =
                                figureDirectory + File.separator + cf.folderString +
                                        File.separator + f.getName().substring(0, f.getName().length() - 4) +
                                        "_RTcalibrationPeptides";

                        //figure in real units
                        new RTCalibrationFigure(mzml, pinFiles[i].getCanonicalPath(),
                                Constants.loessScatterOpacity, mzml.expAndPredRTsMinutes);

                        //save absolute delta RTs
                        for (Map.Entry<String, double[][]> entry : mzml.expAndPredRTsMinutes.entrySet()) {
                            double[][] absDeltaRT = new double[2][entry.getValue()[0].length];
                            for (int j = 0; j < entry.getValue()[0].length; j++) {
                                absDeltaRT[1][j] = entry.getValue()[1][j] - entry.getValue()[0][j];
                            }
                            absDeltaRT[0] = entry.getValue()[0];

                            String mzmlname = new File(mzml.pathStr).getName();
                            if (absDeltaRTs.containsKey(entry.getKey())) {
                                HashMap<String, double[][]> oldMap = absDeltaRTs.get(entry.getKey());
                                oldMap.put(mzmlname.substring(0, mzmlname.length() - 5), absDeltaRT);
                                absDeltaRTs.put(entry.getKey(), oldMap);
                            } else {
                                HashMap<String, double[][]> newMap = new HashMap<>();
                                newMap.put(mzmlname.substring(0, mzmlname.length() - 5), absDeltaRT);
                                absDeltaRTs.put(entry.getKey(), newMap);
                            }
                        }

                        //individual figures by mass
                        //separate figures
                        for (String mass : mzml.expAndPredRTs.keySet()) {
                            if (!mass.isEmpty()) {
                                //path name
                                String massName = mass.split("&")[0];

                                HashMap<String, double[][]> miniMassToData = new HashMap<>();
                                miniMassToData.put(mass, mzml.expAndPredRTs.get(mass));
                                HashMap<String, Function1<Double, Double>> miniLoessFunctions = new HashMap<>();
                                miniLoessFunctions.put(mass, mzml.RTLOESS.get(mass));
                                new RTCalibrationFigure(mzml,
                                        pinFiles[i].getCanonicalPath().substring(0,
                                                pinFiles[i].getCanonicalPath().length() - 4) + "_" + massName + ".pin",
                                        Constants.loessScatterOpacity, miniMassToData, miniLoessFunctions);

                                miniMassToData = new HashMap<>();
                                miniMassToData.put(mass, mzml.expAndPredRTsMinutes.get(mass));
                                new RTCalibrationFigure(mzml,
                                        pinFiles[i].getCanonicalPath().substring(0,
                                                pinFiles[i].getCanonicalPath().length() - 4) + "_" + massName + ".pin",
                                        Constants.loessScatterOpacity, miniMassToData);

                                //write peptides used
                                String calibrationPeptideFilePath = calibrationPeptideFilePathBase;
                                calibrationPeptideFilePath += "_" + massName;
                                try (BufferedWriter writer = new BufferedWriter(
                                        new FileWriter(calibrationPeptideFilePath + ".txt"))) {
                                    ArrayList<String> RTpeptides = mzml.RTpeptides.get(mass);
                                    for (String peptide : RTpeptides) {
                                        writer.write(peptide);
                                        writer.newLine();
                                    }
                                } catch (IOException e) {
                                    e.printStackTrace();
                                    System.exit(1);
                                }
                            }
                        }
                    }
                }
                if (featuresList.contains("deltaRTlinear")) {
                    if (mzml.expAndPredRTs != null) {
                        mzml.setBetas();
                    } else { mzml.setBetas(Constants.rtLoessRegressionSize);
                    }
                    mzml.normalizeRTs(executorService);
                }
                if (featuresList.contains("deltaRTbins") || featuresList.contains("RTzscore") ||
                        featuresList.contains("RTprobability") || featuresList.contains("RTprobabilityUnifPrior") ||
                        featuresList.contains("deltaRTLOESSnormalized")) {
                    mzml.setRTbins();
                    mzml.calculateBinStats("RT");
                }
                if (featuresList.contains("deltaRTbins") || featuresList.contains("RTzscore")) {
                    mzml.calculateDeltaRTbinAndRTzscore(executorService);
                }
                if (featuresList.contains("deltaRTLOESSnormalized")) {
                    mzml.calculateDeltaRTLOESSnormalized(executorService);
                }

                if (featuresList.contains("RTprobabilityUnifPrior")) {
                    mzml.setRTBinSizes(executorService);

                    //decide how many uniform points to add
                    int[] binSizes = new int[mzml.RTbins.length];
                    for (int bin = 0; bin < mzml.RTbins.length; bin++) {
                        binSizes[bin] = mzml.RTbins[bin].size();
                    }
                    Arrays.sort(binSizes);
                    int cutoff = (int) Math.floor(((double) mzml.RTbins.length / 100.0) * Constants.uniformPriorPercentile);
                    mzml.unifPriorSize = binSizes[cutoff];

                    //also need uniform probability with right bound of max predicted RT
                    mzml.unifProb = 1.0f / allPreds.getMaxPredRT(); //TODO: might need to change, if this is ever revisited. Need negative range
                }
                if (featuresList.contains("RTprobability") || featuresList.contains("RTprobabilityUnifPrior")) {
                    mzml.setKernelDensities(executorService, "RT");
                }
                if (featuresList.contains("deltaIMLOESS") || featuresList.contains("deltaIMLOESSnormalized")) {
                    mzml.setLOESS(Constants.imLoessRegressionSize, Constants.loessBandwidth, Constants.robustIters, "IM");
                    mzml.predictIMLOESS(executorService);

                    CalibrationFigure cf;
                    HashSet<String> writtenMasses = new HashSet<>();
                    for (int charge : mzml.expAndPredIMsHashMap.keySet()) {
                        boolean plot = false;
                        for (Map.Entry<String, double[][]> entry : mzml.expAndPredIMsHashMap.get(charge).entrySet()) {
                            if (entry.getValue() != null) {
                                plot = true;
                                break;
                            }
                        }
                        if (plot) {
                            cf = new IMCalibrationFigure(mzml, pinFiles[i].getCanonicalPath(), Constants.loessScatterOpacity,
                                    mzml.expAndPredIMsHashMap.get(charge), mzml.IMLOESS.get(charge - 1), charge);

                            File f = new File(pinFiles[i].getCanonicalPath());
                            String calibrationPeptideFilePathBase =
                                    figureDirectory + File.separator + cf.folderString +
                                            File.separator + f.getName().substring(0, f.getName().length() - 4) +
                                            "_IMcalibrationPeptides";

                            //individual figures by mass
                            //separate figures
                            for (String mass : mzml.expAndPredIMsHashMap.get(charge).keySet()) {
                                if (!mass.isEmpty()) {
                                    String massName = mass.split("&")[0];

                                    HashMap<String, double[][]> miniMassToData = new HashMap<>();
                                    miniMassToData.put(mass, mzml.expAndPredIMsHashMap.get(charge).get(mass));
                                    HashMap<String, Function1<Double, Double>> miniLoessFunctions = new HashMap<>();
                                    miniLoessFunctions.put(mass, mzml.IMLOESS.get(charge - 1).get(mass));
                                    new IMCalibrationFigure(mzml,
                                            pinFiles[i].getCanonicalPath().substring(0,
                                                    pinFiles[i].getCanonicalPath().length() - 4) + "_" + massName + ".pin",
                                            Constants.loessScatterOpacity, miniMassToData, miniLoessFunctions, charge);

                                    //write peptides used
                                    if (!writtenMasses.contains(mass)) {
                                        String calibrationPeptideFilePath = calibrationPeptideFilePathBase;
                                        calibrationPeptideFilePath += "_" + massName;

                                        try (BufferedWriter writer = new BufferedWriter(
                                                new FileWriter(calibrationPeptideFilePath + ".txt"))) {
                                            ArrayList<String> IMpeptides = mzml.IMpeptides.get(mass);
                                            for (String peptide : IMpeptides) {
                                                writer.write(peptide);
                                                writer.newLine();
                                            }
                                            writtenMasses.add(mass);
                                        } catch (IOException e) {
                                            e.printStackTrace();
                                            System.exit(1);
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
                if (featuresList.contains("deltaIMLOESSnormalized") || featuresList.contains("IMprobabilityUnifPrior")) {
                    mzml.setIMbins();
                    mzml.calculateBinStats("IM");
                }
                if (featuresList.contains("deltaIMLOESSnormalized")) {
                    mzml.calculateDeltaIMLOESSnormalized(executorService);
                }

                if (featuresList.contains("IMprobabilityUnifPrior")) {
                    mzml.setIMBinSizes(executorService);

                    //decide how many uniform points to add
                    for (int charge = 0; charge < IMFunctions.numCharges; charge++) {
                        int[] binSizes = new int[mzml.IMbins[charge].length];
                        for (int bin = 0; bin < mzml.IMbins[charge].length; bin++) {
                            binSizes[bin] = mzml.IMbins[charge][bin].size();
                        }
                        Arrays.sort(binSizes);
                        int cutoff = (int) Math.floor(((double) mzml.IMbins[charge].length / 100.0) * Constants.uniformPriorPercentile);
                        mzml.unifPriorSizeIM[charge] = binSizes[cutoff];

                        //also need uniform probability with right bound of max predicted RT
                        mzml.unifProbIM[charge] = 1.0f / (2 * Constants.IMbinMultiplier);

                        mzml.setKernelDensities(executorService, "IM");
                    }
                }

                printInfo("Calculating features");
                FeatureCalculator fc = new FeatureCalculator(pin, featuresList, mzml);
                fc.calculate(executorService);
                if (!fc.ms2Scores.isEmpty()) {
                    String name = new File(mzml.pathStr).getName();
                    scoreDescriptors.put(name.substring(0, name.length() - 5), fc.ms2Scores);
                }

                printInfo("Writing features");
                PinWriter pw = new PinWriter(newOutfile, pin, featuresList, mzml);
                pw.write();
                String histFile;
                if (Constants.renamePin == 1) {
                    printInfo("Edited pin file at " + newOutfile);
                    histFile = newOutfile;
                } else { //really should be 0
                    //move file at newOutfile to pinFiles[i] canonical name
                    File movedFile = new File(newOutfile);
                    pinFiles[i].delete();
                    movedFile.renameTo(pinFiles[i]);
                    histFile = pinFiles[i].getCanonicalPath();
                }

                //plot hist
                new ScoreHistogram(new PinReader(histFile), pw.newColumnNames);
                mzml.clear();
            }
        } catch (Exception e) {
            executorService.shutdown();
            e.printStackTrace();
            System.exit(1);
        }

        //cumulative plots
        if (pmMatcher.pinFiles.length > 1) {
            MyFileUtils.createWholeDirectory(figureDirectory + File.separator + "cumulativeQC");

            //RT plot
            for (Map.Entry<String, TreeMap<String, List<Float>[]>> entry : RTregressionCurves.entrySet()) {
                new CumulativeCurves(entry.getKey(), entry.getValue(), "RT");
            }

            //delta RT plot
            if (!absDeltaRTs.isEmpty()) {
                new CumulativeAbsDeltaLinePlots(absDeltaRTs);
            }

            //spectra plot
            if (!scoreDescriptors.isEmpty()) {
                new CumulativeMS2LinePlots(scoreDescriptors);
            }
        }
    }
}
