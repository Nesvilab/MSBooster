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
import com.google.common.collect.RangeMap;
import com.google.common.collect.TreeRangeMap;
import features.FeatureCalculator;
import features.rtandim.IMFunctions;
import figures.CalibrationFigure;
import figures.IMCalibrationFigure;
import figures.RTCalibrationFigure;
import figures.ScoreHistogram;
import kotlin.jvm.functions.Function1;
import predictions.PredictionEntry;
import predictions.PredictionEntryHashMap;
import readers.MgfFileReader;
import readers.datareaders.MzmlReader;
import readers.datareaders.PinReader;
import readers.predictionreaders.LibraryPredictionMapper;
import writers.PinWriter;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Map;
import java.util.concurrent.ExecutorService;

import static allconstants.Constants.figureDirectory;
import static allconstants.FragmentIonConstants.createFragmentGroups;
import static utils.Print.printInfo;

public class PercolatorFormatter {

    public static PredictionEntryHashMap allPreds = new PredictionEntryHashMap();
    static RangeMap<Double, ArrayList<Integer>> allMatchedScans = TreeRangeMap.create();

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
            } else if (Constants.matchWithDaltons == null) {
                Constants.matchWithDaltons = false;
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

//                if (featuresList.contains("adjacentSimilarity")) {
//                    printInfo("Calculating adjacent similarity");
//                    //check that scan m/z ranges are in mzml
//                    if (mzml.getScanNumObject(mzml.getScanNums().first()).isolationUpper == 0.0d) {
//                        printError("Adjacent similarity feature requires m/z ranges " +
//                                "for each scan. This mzml is incompatible. Exiting.");
//                        System.exit(1);
//                    }
//
//                    SortedSet<Double> scanNumbers = new TreeSet<>();
//                    for (int num : mzml.getScanNums()) {
//                        MzmlScanNumber msn = mzml.getScanNumObject(num);
//                        scanNumbers.add(msn.isolationLower);
//                        scanNumbers.add(msn.isolationUpper);
//                    }
//                    Double[] scanNumbersList = new Double[scanNumbers.size()];
//                    int scanNumberIdx = 0;
//                    for (Double d : scanNumbers) {
//                        scanNumbersList[scanNumberIdx] = d;
//                        scanNumberIdx += 1;
//                    }
//
//                    //initialize
//                    for (int j = 0; j < scanNumbersList.length - 1; j++) {
//                        allMatchedScans.put(Range.open(scanNumbersList[j], scanNumbersList[j + 1]),
//                                new ArrayList<>());
//                    }
//
//                    //add scan numbers to ranges with sub range maps?
//                    for (int num : mzml.getScanNums()) {
//                        MzmlScanNumber msn = mzml.getScanNumObject(num);
//                        RangeMap<Double, ArrayList<Integer>> rangemap =
//                                allMatchedScans.subRangeMap(Range.open(msn.isolationLower, msn.isolationUpper));
//                        for (ArrayList<Integer> entry :
//                                rangemap.asMapOfRanges().values()) {
//                            entry.add(msn.scanNum);
//                        }
//                    }
//
//                    String[] precursors = new String[allPreds.size()];
//                    int precursorNum = 0;
//                    for (String s : allPreds.keySet()) {
//                        precursors[precursorNum] = s;
//                        precursorNum += 1;
//                    }
//
//                    //TODO: make less clunky
//                    for (int num : mzml.getScanNums()) {
//                        MzmlScanNumber msn = mzml.getScanNumObject(num);
//                        for (PeptideObj pobj : msn.peptideObjects) {
//                            if (pobj == null) {
//                                break;
//                            }
//                            PredictionEntry pe = allPreds.get(pobj.name);
//                            double mz = pobj.precursorMz;
//                            if (pe.precursorMz != 0d) {
//                                mz = pe.precursorMz;
//                            }
//
//                            ArrayList<Integer> scans = allMatchedScans.get(mz);
//                            if (scans.size() == 0) {
//                                mz -= 0.0001;
//                                scans = allMatchedScans.get(mz);
//                            }
//
//                            //adjacent windows where one has more scans than other
//                            if (allMatchedScans.get(mz + 0.0002).size() > scans.size() &&
//                                    allMatchedScans.get(mz + 0.0002).contains(scans.get(0))) {
//                                mz += 0.0002;
//                                scans = allMatchedScans.get(mz);
//                            }
//
//                            int scanIdx = scans.indexOf(msn.scanNum);
//                            if (scanIdx == -1) {
//                                mz += 0.0002;
//                                scans = allMatchedScans.get(mz);
//                                scanIdx = scans.indexOf(msn.scanNum);
//
//                                if (scanIdx == -1) {
//                                    mz -= 0.0004;
//                                    scans = allMatchedScans.get(mz);
//                                    scanIdx = scans.indexOf(msn.scanNum);
//                                }
//                            }
//
//                            if (pe.times.size() == 0) {
//                                pobj.chromatogramWindowQuery = Math.min(Constants.chromatogramWindow, scanIdx);
//                            } else {
//                                pobj.chromatogramWindowQuery = scanIdx - pe.times.get(0) +
//                                        Math.min(Constants.chromatogramWindow, pe.times.get(0));
//                            }
//                            pe.times.add(scanIdx);
//                            pe.precursorMz = mz;
//                            allPreds.put(pobj.name, pe);
//                        }
//                    }
//
//                    mzml.futureList.clear();
//                    Multithreader mt = new Multithreader(allPreds.size(), Constants.numThreads);
//                    for (int j = 0; j < Constants.numThreads; j++) {
//                        int finalI = j;
//                        mzml.futureList.add(executorService.submit(() -> {
//                            ProgressReporter pr = new ProgressReporter(mt.indices[finalI + 1] - mt.indices[finalI]);
//                            PredictionEntry predictionEntry;
////                                String[] entrySplit;
////                                MassCalculator mc;
//                            ArrayList<Integer> scans;
//                            MzmlScanNumber msn;
//                            Float[] scores;
//                            for (int k = mt.indices[finalI]; k < mt.indices[finalI + 1]; k++) {
//                                pr.progress();
//
//                                SpectrumComparison sc = null;
//
//                                String precursor = precursors[k];
//                                precursors[k] = null;
//                                predictionEntry = allPreds.get(precursor);
//
////                                    entrySplit = precursor.split("\\|");
////                                    int charge = Integer.parseInt(entrySplit[1]);
////                                    mc = new MassCalculator(entrySplit[0], entrySplit[1]);
////                                    double precursorMz = (mc.mass + charge * mc.proton) / charge;
//
//                                scans = allMatchedScans.get(predictionEntry.precursorMz);
////                                    if (scans.size() == 0) {
////                                        scans = allMatchedScans.get(precursorMz - 0.0001);
////                                        scans.addAll(allMatchedScans.get(precursorMz + 0.0001));
////                                        scans.sort(Comparator.naturalOrder());
////                                    }
//                                List<Integer> scansList = scans.subList(Math.max(0, predictionEntry.times.get(0) - Constants.chromatogramWindow),
//                                            Math.min(scans.size(), predictionEntry.times.get(predictionEntry.times.size() - 1) + Constants.chromatogramWindow));
//
////                                    double maxScore = 0d;
////                                    int bestScan = 0;
//                                predictionEntry.scores.put("entropy", new Float[scansList.size()]);
//                                predictionEntry.scores.put("hypergeom", new Float[scansList.size()]);
//                                predictionEntry.scores.put("spearman", new Float[scansList.size()]);
//                                int scoreIdx = 0;
////                                    double minDeltaRT = Double.MAX_VALUE;
////
////                                    int RTidx = 0;
//                                for (int scan : scansList) {
//                                    try {
//                                        msn = mzml.getScanNumObject(scan);
//                                    } catch (FileParsingException e) {
//                                        throw new RuntimeException(e);
//                                    }
////                                        double deltaRT = Math.abs(msn.calibratedRT - predictionEntry.RT);
////                                        if (deltaRT < minDeltaRT) {
////                                            minDeltaRT = deltaRT;
////                                            predictionEntry.bestScanIdx = RTidx;
////                                        }
//                                    PeptideObj pobj = new PeptideObj();
//                                    pobj.name = precursor;
//                                    pobj.charge = Integer.parseInt(precursor.split("\\|")[1]);
//                                    pobj.scanNumObj = msn;
//                                    pobj.length = 0;
//                                    for (int l = 0; l < pobj.name.length() - 2; l++) {
//                                        if (Character.isAlphabetic(pobj.name.charAt(l))) {
//                                            pobj.length += 1;
//                                        }
//                                    }
//
//                                    if (sc != null) {
//                                        sc.reload(pobj, msn.getExpMZs(), msn.getExpIntensities());
//                                    } else {
//                                        sc = new SpectrumComparison(pobj,
//                                                msn.getExpMZs(), msn.getExpIntensities(),
//                                                predictionEntry.mzs, predictionEntry.intensities,
//                                                predictionEntry.fragmentIonTypes, true);
//                                    }
//                                    float score = (float) sc.unweightedSpectralEntropy();
////                                        if (score > maxScore) {
////                                            maxScore = score;
////                                            bestScan = scan;
////                                        }
//                                    scores = predictionEntry.scores.get("entropy");
//                                    scores[scoreIdx] = score;
//                                    predictionEntry.scores.put("entropy", scores);
//
//                                    try {
//                                        score = (float) sc.hypergeometricProbability();
//                                    } catch (IOException | URISyntaxException e) {
//                                        throw new RuntimeException(e);
//                                    }
//                                    scores = predictionEntry.scores.get("hypergeom");
//                                    scores[scoreIdx] = score;
//                                    predictionEntry.scores.put("hypergeom", scores);
//
//                                    score = (float) sc.spearmanCorr();
//                                    scores = predictionEntry.scores.get("spearman");
//                                    scores[scoreIdx] = score;
//                                    predictionEntry.scores.put("spearman", scores);
//
//                                    scoreIdx += 1;
////                                        RTidx += 1;
//                                }
//                                //predictionEntry.bestScan = bestScan;
//                                //TODO: what if repeated scans? Longer than other arrays
////                                    scores = predictionEntry.scores.get("entropy");
////
////                                    int startShift = 0;
////                                    int endShift = 0;
////                                    int arrayStart = predictionEntry.bestScanIdx - window;
////                                    int arrayEnd = predictionEntry.bestScanIdx + window;
////                                    if (arrayStart < 0) {
////                                        endShift = window - predictionEntry.bestScanIdx;
////                                        arrayStart = 0;
////                                    }
////                                    if (arrayEnd >= scores.length) {
////                                        startShift = arrayEnd - scores.length + 1;
////                                        arrayEnd = scores.length - 1;
////                                    }
////
////                                    Float[] newScores1 = Arrays.copyOfRange(scores,
////                                            arrayStart - startShift,
////                                            arrayEnd + endShift + 1);
////                                    predictionEntry.scores.put("entropy", newScores1);
////
////                                    scores = predictionEntry.scores.get("hypergeom");
////                                    Float[] newScores2 = Arrays.copyOfRange(scores,
////                                            arrayStart - startShift,
////                                            arrayEnd + endShift + 1);
////                                    predictionEntry.scores.put("hypergeom", newScores2);
////                                    allPreds.put(precursor, predictionEntry);
//                            }
//                        }));
//                    }
//                    for (Future future : mzml.futureList) {
//                        future.get();
//                    }
//                }
                if (featuresList.contains("deltaRTLOESS") || featuresList.contains("deltaRTLOESSnormalized")) {
                    mzml.setLOESS(Constants.rtLoessRegressionSize, Constants.loessBandwidth, Constants.robustIters, "RT");
                    mzml.predictRTLOESS(executorService); //potentially only invoke once if normalized included

                    //generate calibration figure, need mzml and loess
                    boolean plot = false;
                    for (Map.Entry<String, double[][]> entry : mzml.expAndPredRTs.entrySet()) {
                        if (entry.getValue() != null) {
                            plot = true;
                            break;
                        }
                    }
                    if (plot) {
                        CalibrationFigure cf = new RTCalibrationFigure(mzml, pinFiles[i].getCanonicalPath(), Constants.loessScatterOpacity,
                                mzml.expAndPredRTs, mzml.RTLOESS);

                        File f = new File(pinFiles[i].getCanonicalPath());
                        String calibrationPeptideFilePathBase =
                                figureDirectory + File.separator + cf.folderString +
                                        File.separator + f.getName().substring(0, f.getName().length() - 4) +
                                        "_RTcalibrationPeptides";

                        //individual figures by mass
                        //separate figures
                        int massI = 1;
                        for (String mass : mzml.expAndPredRTs.keySet()) {
                            if (!mass.isEmpty()) {
                                HashMap<String, double[][]> miniMassToData = new HashMap<>();
                                miniMassToData.put(mass, mzml.expAndPredRTs.get(mass));
                                HashMap<String, Function1<Double, Double>> miniLoessFunctions = new HashMap<>();
                                miniLoessFunctions.put(mass, mzml.RTLOESS.get(mass));
                                new RTCalibrationFigure(mzml,
                                        pinFiles[i].getCanonicalPath().substring(0,
                                                pinFiles[i].getCanonicalPath().length() - 4) + "_" + massI + ".pin",
                                        Constants.loessScatterOpacity, miniMassToData, miniLoessFunctions);

                                //write peptides used
                                String calibrationPeptideFilePath = calibrationPeptideFilePathBase;
                                calibrationPeptideFilePath += "_" + massI;
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

                                massI++;
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
                            int massI = 1;
                            for (String mass : mzml.expAndPredIMsHashMap.get(charge).keySet()) {
                                if (!mass.isEmpty()) {
                                    HashMap<String, double[][]> miniMassToData = new HashMap<>();
                                    miniMassToData.put(mass, mzml.expAndPredIMsHashMap.get(charge).get(mass));
                                    HashMap<String, Function1<Double, Double>> miniLoessFunctions = new HashMap<>();
                                    miniLoessFunctions.put(mass, mzml.IMLOESS.get(charge - 1).get(mass));
                                    new IMCalibrationFigure(mzml,
                                            pinFiles[i].getCanonicalPath().substring(0,
                                                    pinFiles[i].getCanonicalPath().length() - 4) + "_" + massI + ".pin",
                                            Constants.loessScatterOpacity, miniMassToData, miniLoessFunctions, charge);

                                    //write peptides used
                                    if (!writtenMasses.contains(mass)) {
                                        String calibrationPeptideFilePath = calibrationPeptideFilePathBase;
                                        calibrationPeptideFilePath += "_" + massI;

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

                                    massI++;
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

                //TODO: multithread?
                printInfo("Calculating features");
                FeatureCalculator fc = new FeatureCalculator(pin, featuresList, mzml);
                fc.calculate(executorService);

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
    }
}
