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

import static utils.Print.printError;
import static utils.Print.printInfo;

import allconstants.Constants;
import features.FeatureCalculator;
import features.detectability.FastaReader;
import features.rtandim.IMFunctions;
import features.spectra.SpectrumComparison;
import figures.CalibrationFigure;
import figures.IMCalibrationFigure;
import figures.RTCalibrationFigure;
import figures.ScoreHistogram;
import peptideptmformatting.PeptideFormatter;
import predictions.PredictionEntry;
import predictions.PredictionEntryHashMap;
import readers.datareaders.MzmlReader;
import readers.datareaders.PinReader;
import readers.MgfFileReader;
import readers.predictionreaders.LibraryPredictionMapper;
import writers.PinWriter;
import com.google.common.collect.Range;
import com.google.common.collect.RangeMap;
import com.google.common.collect.TreeRangeMap;
import kotlin.jvm.functions.Function1;
import org.apache.commons.lang3.ArrayUtils;
import umich.ms.fileio.exceptions.FileParsingException;
import utils.ProgressReporter;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.sql.SQLException;
import java.util.*;
import java.util.concurrent.ExecutionException;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Future;

public class PercolatorFormatter {

    public static PredictionEntryHashMap allPreds = new PredictionEntryHashMap();
    static RangeMap<Double, ArrayList<Integer>> allMatchedScans = TreeRangeMap.create();

    public static void editPin(PinMzmlMatcher pmMatcher, String[] features, String outfile,
                               ExecutorService executorService)
            throws IOException, InterruptedException, ExecutionException, FileParsingException, SQLException {

        ArrayList<String> featuresList = new ArrayList<>(Arrays.asList(features));
        //remove features from array for multiple protein formatting
        if (featuresList.contains("detectability")) {
            int idx = ArrayUtils.indexOf(features, "detectability");
            features = ArrayUtils.remove(features, idx);
        }

        File[] pinFiles = pmMatcher.pinFiles;
        File[] mzmlFiles = pmMatcher.mzmlFiles;
        MzmlReader[] mzmlReaders = pmMatcher.mzmlReaders;

        //load predicted spectra
        LibraryPredictionMapper predictedSpectra;
        LibraryPredictionMapper predictedRT;
        LibraryPredictionMapper predictedIM;
        LibraryPredictionMapper predictedAuxSpectra;
        //SpectralPredictionMapper predictedSpectra2 = null; //first is prosit/diann, second predfull

        if (Constants.predictedLibrary != null) { //library ready from koina predictions
            allPreds = Constants.predictedLibrary.getPreds();
        } else {
            HashMap<String, LibraryPredictionMapper> allLibraries = new HashMap<>(); //key: library file path, value: library
            HashMap<String, String> allProperties = new HashMap<>(); //key: property, value: library file path

            //could use aux spectra if primary spectra missing
            if (Constants.spectraPredFile != null) {
                printInfo("Loading predicted spectra");
                if (Constants.spectraModel.equals("PredFull")) {
                    predictedSpectra = LibraryPredictionMapper.createLibraryPredictionMapper(
                            Constants.spectraPredFile, pinFiles, executorService);
                } else {
                    predictedSpectra = LibraryPredictionMapper.createLibraryPredictionMapper(
                            Constants.spectraPredFile, Constants.spectraModel, executorService);
                }
                allLibraries.put(Constants.spectraPredFile, predictedSpectra);
                allProperties.put("spectra", Constants.spectraPredFile);
            }

            if (Constants.RTPredFile != null) {
                if (! allLibraries.containsKey(Constants.RTPredFile)) {
                    printInfo("Loading predicted retention times");
                    predictedRT = LibraryPredictionMapper.createLibraryPredictionMapper(
                            Constants.RTPredFile, Constants.rtModel, executorService);
                    allLibraries.put(Constants.RTPredFile, predictedRT);
                }
                allProperties.put("RT", Constants.RTPredFile);
            }

            if (Constants.IMPredFile != null) {
                if (! allLibraries.containsKey(Constants.IMPredFile)) {
                    printInfo("Loading predicted ion mobilities");
                    predictedIM = LibraryPredictionMapper.createLibraryPredictionMapper(
                            Constants.IMPredFile, Constants.imModel, executorService);
                    allLibraries.put(Constants.IMPredFile, predictedIM);
                }
                allProperties.put("IM", Constants.IMPredFile);
            }

            if (Constants.auxSpectraPredFile != null) {
                if (! allLibraries.containsKey(Constants.auxSpectraPredFile)) {
                    printInfo("Loading predicted auxiliary spectra");
                    predictedAuxSpectra = LibraryPredictionMapper.createLibraryPredictionMapper(
                            Constants.auxSpectraPredFile, Constants.auxSpectraModel, executorService);
                    allLibraries.put(Constants.auxSpectraPredFile, predictedAuxSpectra);
                }
                allProperties.put("auxSpectra", Constants.auxSpectraPredFile);
            }

            //merging libraries
            printInfo("Merging libraries");
            for (Map.Entry<String, LibraryPredictionMapper> entry : allLibraries.entrySet()) {
                String libraryPath = entry.getKey();
                LibraryPredictionMapper library = entry.getValue();
                for (Map.Entry<String, String> prop : allProperties.entrySet()) {
                    if (prop.getValue().equals(libraryPath)) {
                        library.mergeLibraries(allPreds, prop.getKey());
                    }
                }

            }
        }
        //TODO test to make sure predicted spectra and allpreds are same

        //create detectMap to store detectabilities for base sequence peptides
        //store peptide detectabilities in PredictionEntry
//        DetectMap dm = null;
//        ArrayList<String> dFeatures = new ArrayList<String>(Constants.detectFeatures);
//        dFeatures.retainAll(featuresList);
//        //long startTime = System.nanoTime();
//        if (dFeatures.size() > 0) {
//            dm = new DetectMap(detectFile);
//            for (Map.Entry<String, PredictionEntry> e : allPreds.entrySet()) {
//                e.getValue().setDetectability(dm.getDetectability(
//                        new PeptideFormatter(e.getKey().split("\\|")[0], e.getKey().split("\\|")[1], "pin").stripped));
//            }
//        }

        FastaReader fasta = null;

        if (featuresList.contains("detectFractionGreater") || featuresList.contains("detectSubtractMissing")
                || featuresList.contains("detectProtSpearmanDiff") || featuresList.contains("peptideCounts")) {
            if (Constants.peptideCounter.isEmpty()) {
                //get all peptides present in pin
                for (File pinFile : pmMatcher.pinFiles) {
                    PinReader pin = new PinReader(pinFile.getCanonicalPath());

                    //add to counter
                    while (pin.next(true)) {
                        PeptideFormatter pf = pin.getPep();
                        if (Float.valueOf(pin.getColumn("hyperscore")) > 10) {
                            if (Constants.peptideCounter.containsKey(pf.getStripped())) {
                                HashSet<String> peptideSet = Constants.peptideCounter.get(pf.getStripped());
                                peptideSet.add(pin.name);
                                Constants.peptideCounter.put(pf.getStripped(), peptideSet);
                            } else {
                                HashSet<String> peptideSet = new HashSet<>();
                                peptideSet.add(pin.name);
                                Constants.peptideCounter.put(pf.getStripped(), peptideSet);
                            }
                        } else {
                            if (! Constants.peptideCounter.containsKey(pf.getStripped())) {
                                HashSet<String> peptideSet = new HashSet<>();
                                Constants.peptideCounter.put(pf.getStripped(), peptideSet);
                            }
                        }
                    }
                    pin.close();
                }
            }

//            //load fasta
//            if (Constants.getFastaReader() == null) {
//                printInfo("Creating fasta object");
//                fasta = new FastaReader(Constants.fasta);
//            } else {
//                printInfo("Loading fasta");
//                fasta = Constants.getFastaReader();
//            }
//
//            printInfo("Loading detectabilities for unique peptides from each protein");
//            for (Map.Entry<String, ProteinEntry> e : fasta.protToPep.entrySet()) {
//                ArrayList<String> pepList = e.getValue().peptides;
//                float[] protDetects = new float[pepList.size()]; //for storing initial detect order
//
//                //store detect unsorted
//                for (int pep = 0; pep < pepList.size(); pep++) {
//                    protDetects[pep] = dm.getDetectability(pepList.get(pep));
//                }
//
//                //dual pivot quicksort
//                //sorted indices
//                int[] sortedIndices = IntStream.range(0, protDetects.length)
//                        .boxed().sorted((k, j) -> Float.compare(protDetects[k], protDetects[j]))
//                        .mapToInt(ele -> ele).toArray();
//
//                float[] sortedDetect = new float[protDetects.length];
//                for (int j = 0; j < protDetects.length; j++) {
//                    sortedDetect[j] = protDetects[sortedIndices[j]];
//                }
//                e.getValue().detects = sortedDetect;
//
//                //check which peptides present, and get spectral counts
//                float[] protPresence = new float[protDetects.length];
//                float[] pepCounts = new float[protDetects.length];
//                float numPresent = 1f;
//                for (int j = protDetects.length - 1; j > -1; j--) {
//                    String currentPep = pepList.get(sortedIndices[j]);
//
//                    if (Constants.peptideCounter.containsKey(currentPep)) {
//                        protPresence[j] = numPresent;
//                        numPresent += 1f;
//                        pepCounts[j] = Constants.peptideCounter.get(currentPep);
//                    }
//                }
//                fasta.protToPep.get(e.getKey()).presence = protPresence;
//                fasta.protToPep.get(e.getKey()).spectralCounts = pepCounts;
//            }
//
//            for (Map.Entry<String, PredictionEntry> e : allPreds.entrySet()) {
//                try {
//                    e.getValue().setCounter(Constants.peptideCounter.get(e.getKey().split("\\|")[0]));
//                } catch (Exception ee) { //peptide was in a pin file from another run
//                }
//            }
//            dm.clear();
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
                //these require all experimental peaks before removing higher rank peaks
                if (Constants.removeRankPeaks &&
                        (featuresList.contains("hypergeometricProbability") ||
                                featuresList.contains("intersection") ||
                                featuresList.contains("adjacentSimilarity"))) {
                    for (int num : mzml.getScanNums()) {
                        MzmlScanNumber msn = mzml.getScanNumObject(num);
                        msn.expMZs = msn.savedExpMZs;
                        msn.expIntensities = msn.savedExpIntensities;
                        msn.savedExpMZs = null;
                        msn.savedExpIntensities = null;
                    }
                    Constants.removeRankPeaks = false;
                }

                if (featuresList.contains("adjacentSimilarity")) {
                    printInfo("Calculating adjacent similarity");
                    //check that scan m/z ranges are in mzml
                    if (mzml.getScanNumObject(mzml.getScanNums().first()).isolationUpper == 0.0d) {
                        printError("Adjacent similarity feature requires m/z ranges " +
                                "for each scan. This mzml is incompatible. Exiting.");
                        System.exit(1);
                    }

                    SortedSet<Double> scanNumbers = new TreeSet<>();
                    for (int num : mzml.getScanNums()) {
                        MzmlScanNumber msn = mzml.getScanNumObject(num);
                        scanNumbers.add(msn.isolationLower);
                        scanNumbers.add(msn.isolationUpper);
                    }
                    Double[] scanNumbersList = new Double[scanNumbers.size()];
                    int scanNumberIdx = 0;
                    for (Double d : scanNumbers) {
                        scanNumbersList[scanNumberIdx] = d;
                        scanNumberIdx += 1;
                    }

                    //initialize
                    for (int j = 0; j < scanNumbersList.length - 1; j++) {
                        allMatchedScans.put(Range.open(scanNumbersList[j], scanNumbersList[j + 1]),
                                new ArrayList<>());
                    }

                    //add scan numbers to ranges with sub range maps?
                    for (int num : mzml.getScanNums()) {
                        MzmlScanNumber msn = mzml.getScanNumObject(num);
                        RangeMap<Double, ArrayList<Integer>> rangemap =
                                allMatchedScans.subRangeMap(Range.open(msn.isolationLower, msn.isolationUpper));
                        for (ArrayList<Integer> entry :
                                rangemap.asMapOfRanges().values()) {
                            entry.add(msn.scanNum);
                        }
                    }

                    String[] precursors = new String[allPreds.size()];
                    int precursorNum = 0;
                    for (String s : allPreds.keySet()) {
                        precursors[precursorNum] = s;
                        precursorNum += 1;
                    }

                    //TODO: make less clunky
                    for (int num : mzml.getScanNums()) {
                        MzmlScanNumber msn = mzml.getScanNumObject(num);
                        for (PeptideObj pobj : msn.peptideObjects) {
                            if (pobj == null) {
                                break;
                            }
                            PredictionEntry pe = allPreds.get(pobj.name);
                            double mz = pobj.precursorMz;
                            if (pe.precursorMz != 0d) {
                                mz = pe.precursorMz;
                            }

                            ArrayList<Integer> scans = allMatchedScans.get(mz);
                            if (scans.size() == 0) {
                                mz -= 0.0001;
                                scans = allMatchedScans.get(mz);
                            }

                            //adjacent windows where one has more scans than other
                            if (allMatchedScans.get(mz + 0.0002).size() > scans.size() &&
                                    allMatchedScans.get(mz + 0.0002).contains(scans.get(0))) {
                                mz += 0.0002;
                                scans = allMatchedScans.get(mz);
                            }

                            int scanIdx = scans.indexOf(msn.scanNum);
                            if (scanIdx == -1) {
                                mz += 0.0002;
                                scans = allMatchedScans.get(mz);
                                scanIdx = scans.indexOf(msn.scanNum);

                                if (scanIdx == -1) {
                                    mz -= 0.0004;
                                    scans = allMatchedScans.get(mz);
                                    scanIdx = scans.indexOf(msn.scanNum);
                                }
                            }

                            if (pe.times.size() == 0) {
                                pobj.chromatogramWindowQuery = Math.min(Constants.chromatogramWindow, scanIdx);
                            } else {
                                pobj.chromatogramWindowQuery = scanIdx - pe.times.get(0) +
                                        Math.min(Constants.chromatogramWindow, pe.times.get(0));
                            }
                            pe.times.add(scanIdx);
                            pe.precursorMz = mz;
                            allPreds.put(pobj.name, pe);
                        }
                    }

                    mzml.futureList.clear();
                    for (int j = 0; j < Constants.numThreads; j++) {
                        int start = (int) (allPreds.size() * (long) j) / Constants.numThreads;
                        int end = (int) (allPreds.size() * (long) (j + 1)) / Constants.numThreads;
                        mzml.futureList.add(executorService.submit(() -> {
                            ProgressReporter pr = new ProgressReporter(end - start);
                            PredictionEntry predictionEntry;
//                                String[] entrySplit;
//                                MassCalculator mc;
                            ArrayList<Integer> scans;
                            MzmlScanNumber msn;
                            Float[] scores;
                            for (int k = start; k < end; k++) {
                                pr.progress();

                                SpectrumComparison sc = null;

                                String precursor = precursors[k];
                                precursors[k] = null;
                                predictionEntry = allPreds.get(precursor);

//                                    entrySplit = precursor.split("\\|");
//                                    int charge = Integer.parseInt(entrySplit[1]);
//                                    mc = new MassCalculator(entrySplit[0], entrySplit[1]);
//                                    double precursorMz = (mc.mass + charge * mc.proton) / charge;

                                scans = allMatchedScans.get(predictionEntry.precursorMz);
//                                    if (scans.size() == 0) {
//                                        scans = allMatchedScans.get(precursorMz - 0.0001);
//                                        scans.addAll(allMatchedScans.get(precursorMz + 0.0001));
//                                        scans.sort(Comparator.naturalOrder());
//                                    }
                                List<Integer> scansList = scans.subList(Math.max(0, predictionEntry.times.get(0) - Constants.chromatogramWindow),
                                            Math.min(scans.size(), predictionEntry.times.get(predictionEntry.times.size() - 1) + Constants.chromatogramWindow));

//                                    double maxScore = 0d;
//                                    int bestScan = 0;
                                predictionEntry.scores.put("entropy", new Float[scansList.size()]);
                                predictionEntry.scores.put("hypergeom", new Float[scansList.size()]);
                                predictionEntry.scores.put("spearman", new Float[scansList.size()]);
                                int scoreIdx = 0;
//                                    double minDeltaRT = Double.MAX_VALUE;
//
//                                    int RTidx = 0;
                                for (int scan : scansList) {
                                    try {
                                        msn = mzml.getScanNumObject(scan);
                                    } catch (FileParsingException e) {
                                        throw new RuntimeException(e);
                                    }
//                                        double deltaRT = Math.abs(msn.calibratedRT - predictionEntry.RT);
//                                        if (deltaRT < minDeltaRT) {
//                                            minDeltaRT = deltaRT;
//                                            predictionEntry.bestScanIdx = RTidx;
//                                        }
                                    PeptideObj pobj = new PeptideObj();
                                    pobj.name = precursor;
                                    pobj.charge = Integer.parseInt(precursor.split("\\|")[1]);
                                    pobj.scanNumObj = msn;
                                    pobj.length = 0;
                                    for (int l = 0; l < pobj.name.length() - 2; l++) {
                                        if (Character.isAlphabetic(pobj.name.charAt(l))) {
                                            pobj.length += 1;
                                        }
                                    }

                                    if (sc != null) {
                                        sc.reload(pobj, msn.getExpMZs(), msn.getExpIntensities());
                                    } else {
                                        sc = new SpectrumComparison(pobj,
                                                msn.getExpMZs(), msn.getExpIntensities(),
                                                predictionEntry.mzs, predictionEntry.intensities, pobj.length,
                                                true);
                                    }
                                    float score = (float) sc.unweightedSpectralEntropy();
//                                        if (score > maxScore) {
//                                            maxScore = score;
//                                            bestScan = scan;
//                                        }
                                    scores = predictionEntry.scores.get("entropy");
                                    scores[scoreIdx] = score;
                                    predictionEntry.scores.put("entropy", scores);

                                    score = (float) sc.hyperGeometricProbability();
                                    scores = predictionEntry.scores.get("hypergeom");
                                    scores[scoreIdx] = score;
                                    predictionEntry.scores.put("hypergeom", scores);

                                    score = (float) sc.spearmanCorr();
                                    scores = predictionEntry.scores.get("spearman");
                                    scores[scoreIdx] = score;
                                    predictionEntry.scores.put("spearman", scores);

                                    scoreIdx += 1;
//                                        RTidx += 1;
                                }
                                //predictionEntry.bestScan = bestScan;
                                //TODO: what if repeated scans? Longer than other arrays
//                                    scores = predictionEntry.scores.get("entropy");
//
//                                    int startShift = 0;
//                                    int endShift = 0;
//                                    int arrayStart = predictionEntry.bestScanIdx - window;
//                                    int arrayEnd = predictionEntry.bestScanIdx + window;
//                                    if (arrayStart < 0) {
//                                        endShift = window - predictionEntry.bestScanIdx;
//                                        arrayStart = 0;
//                                    }
//                                    if (arrayEnd >= scores.length) {
//                                        startShift = arrayEnd - scores.length + 1;
//                                        arrayEnd = scores.length - 1;
//                                    }
//
//                                    Float[] newScores1 = Arrays.copyOfRange(scores,
//                                            arrayStart - startShift,
//                                            arrayEnd + endShift + 1);
//                                    predictionEntry.scores.put("entropy", newScores1);
//
//                                    scores = predictionEntry.scores.get("hypergeom");
//                                    Float[] newScores2 = Arrays.copyOfRange(scores,
//                                            arrayStart - startShift,
//                                            arrayEnd + endShift + 1);
//                                    predictionEntry.scores.put("hypergeom", newScores2);
//                                    allPreds.put(precursor, predictionEntry);
                            }
                        }));
                    }
                    for (Future future : mzml.futureList) {
                        future.get();
                    }
                }
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

                        //individual figures by mass
                        //separate figures
                        for (String mass : mzml.expAndPredRTs.keySet()) {
                            if (!mass.isEmpty()) {
                                HashMap<String, double[][]> miniMassToData = new HashMap<>();
                                miniMassToData.put(mass, mzml.expAndPredRTs.get(mass));
                                HashMap<String, Function1<Double, Double>> miniLoessFunctions = new HashMap<>();
                                miniLoessFunctions.put(mass, mzml.RTLOESS.get(mass));
                                new RTCalibrationFigure(mzml,
                                        pinFiles[i].getCanonicalPath().substring(0,
                                                pinFiles[i].getCanonicalPath().length() - 4) + "_" + mass + ".pin",
                                        Constants.loessScatterOpacity, miniMassToData, miniLoessFunctions);
                            }
                        }

                        //write peptides used
                        File f = new File(pinFiles[i].getCanonicalPath());
                        String calibrationPeptideFilePathBase =
                                f.getParent() + File.separator + "MSBooster_plots" + File.separator + cf.folderString +
                                        File.separator + f.getName().substring(0, f.getName().length() - 4) +
                                        "_RTcalibrationPeptides";
                        for (Map.Entry<String, ArrayList<String>> entry : mzml.RTpeptides.entrySet()) {
                            String calibrationPeptideFilePath = calibrationPeptideFilePathBase;
                            if (!entry.getKey().isEmpty()) {
                                calibrationPeptideFilePath += "_" + entry.getKey();
                            }
                            try (BufferedWriter writer = new BufferedWriter(
                                    new FileWriter(calibrationPeptideFilePath + ".txt"))) {
                                for (String peptide : entry.getValue()) {
                                    writer.write(peptide);
                                    writer.newLine();
                                }
                            } catch (IOException e) {
                                e.printStackTrace();
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

                    CalibrationFigure cf = null;
                    boolean writePeptides = false;
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
                            writePeptides = true;

                            //individual figures by mass
                            //separate figures
                            for (String mass : mzml.expAndPredIMsHashMap.get(charge).keySet()) {
                                if (!mass.isEmpty()) {
                                    HashMap<String, double[][]> miniMassToData = new HashMap<>();
                                    miniMassToData.put(mass, mzml.expAndPredIMsHashMap.get(charge).get(mass));
                                    HashMap<String, Function1<Double, Double>> miniLoessFunctions = new HashMap<>();
                                    miniLoessFunctions.put(mass, mzml.IMLOESS.get(charge - 1).get(mass));
                                    new IMCalibrationFigure(mzml,
                                            pinFiles[i].getCanonicalPath().substring(0,
                                                    pinFiles[i].getCanonicalPath().length() - 4) + "_" + mass + ".pin",
                                            Constants.loessScatterOpacity, miniMassToData, miniLoessFunctions, charge);
                                }
                            }
                        }
                    }
                    if (writePeptides) {
                        //write peptides used
                        File f = new File(pinFiles[i].getCanonicalPath());
                        String calibrationPeptideFilePathBase =
                                f.getParent() + File.separator + "MSBooster_plots" + File.separator + cf.folderString +
                                        File.separator + f.getName().substring(0, f.getName().length() - 4) +
                                        "_IMcalibrationPeptides";
                        for (Map.Entry<String, ArrayList<String>> entry : mzml.IMpeptides.entrySet()) {
                            String calibrationPeptideFilePath = calibrationPeptideFilePathBase;
                            if (!entry.getKey().isEmpty()) {
                                calibrationPeptideFilePath += "_" + entry.getKey();
                            }
                            try (BufferedWriter writer = new BufferedWriter(
                                    new FileWriter(calibrationPeptideFilePath + ".txt"))) {
                                for (String peptide : entry.getValue()) {
                                    writer.write(peptide);
                                    writer.newLine();
                                }
                            } catch (IOException e) {
                                e.printStackTrace();
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
                PinWriter pw = new PinWriter(newOutfile, pin, featuresList, mzml, fc.featureStats);
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
                new ScoreHistogram(new PinReader(histFile), featuresList);
                mzml.clear();
            }
        } catch (Exception e) {
            executorService.shutdown();
            e.printStackTrace();
            System.exit(1);
        }
    }
}
