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

import com.google.common.collect.Range;
import com.google.common.collect.RangeMap;
import com.google.common.collect.TreeRangeMap;
import org.apache.commons.lang3.ArrayUtils;
import umich.ms.fileio.exceptions.FileParsingException;

import java.io.File;
import java.io.IOException;
import java.sql.SQLException;
import java.util.*;
import java.util.concurrent.*;
import java.util.stream.IntStream;

public class PercolatorFormatter {

    static ConcurrentHashMap<String, PredictionEntry> allPreds;
    static RangeMap<Double, ArrayList<Integer>> allMatchedScans = TreeRangeMap.create();

    public static void editPin(PinMzmlMatcher pmMatcher, String mgf, String detectFile,
                               String[] features, String outfile, ExecutorService executorService)
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
        SpectralPredictionMapper predictedSpectra = null;
        SpectralPredictionMapper predictedSpectra2 = null; //first is prosit/diann, second predfull
        if (Constants.spectralPredictionMapper != null) {
            predictedSpectra = Constants.spectralPredictionMapper;
        }

        //Special preparations dependent on features we require
        //only time this isn't needed is detect
        //could consider an mgf constant
        if (mgf != null) {
            String[] mgfSplit = mgf.split(",");

            if (mgfSplit.length == 1) {
                System.out.println("Loading predicted library");
                if (Constants.spectraRTPredModel.equals("PredFull")) {
                    predictedSpectra = SpectralPredictionMapper.createSpectralPredictionMapper(
                            mgfSplit[1], pinFiles, executorService);
                } else {
                    predictedSpectra = SpectralPredictionMapper.createSpectralPredictionMapper(
                            mgf, Constants.spectraRTPredModel, executorService);
                }
            } else if (mgfSplit.length == 2){
                //if fragment from predfull is not y/b, add.
                //Prosit/diann is first, predfull second
                //can also add two models, the first being for RT, the second for spectra
                String[] modelSplit = Constants.spectraRTPredModel.split(",");

                System.out.println("Loading predicted RT: " + mgfSplit[0]);
                predictedSpectra = SpectralPredictionMapper.createSpectralPredictionMapper(
                        mgfSplit[0], modelSplit[0], executorService);
                System.out.println("Loading predicted spectra: " + mgfSplit[1]);
                if (modelSplit[1].equals("PredFull")) {
                    predictedSpectra2 = SpectralPredictionMapper.createSpectralPredictionMapper(
                            mgfSplit[1], pinFiles, executorService); //get predfull library
                } else {
                    predictedSpectra2 = SpectralPredictionMapper.createSpectralPredictionMapper(
                            mgfSplit[1], modelSplit[1], executorService); //get other library
                    Constants.addNonYb = false; //probably just want to use separate spectra and RT models
                }
                System.out.println("Merging spectral libraries");

                //get all possible keys from both preds1 and preds2
                Set<String> totalKeyset = new HashSet<String>();
                allPreds = predictedSpectra.getPreds();
                totalKeyset.addAll(allPreds.keySet());
                totalKeyset.addAll(predictedSpectra2.getPreds().keySet());

                //check what fragment ion types have been predicted by model 1
                HashSet<String> model1FragmentIonTypes = new HashSet<>();
                for (PredictionEntry pe : allPreds.values()) {
                    if (pe.fragmentIonTypes == null) {
                        pe.setFragmentIonTypes();
                    }
                    model1FragmentIonTypes.addAll(Arrays.asList(pe.fragmentIonTypes));
                }

                for (String key : totalKeyset) {
                    PredictionEntry pe = allPreds.get(key);
                    if (pe == null) { //missing in prosit/diann
                        allPreds.put(key, predictedSpectra2.getPreds().get(key));
                    } else { //add non-y/b ions
                        ArrayList<Float> mzs = new ArrayList<>();
                        ArrayList<Float> intensities = new ArrayList<>();
                        ArrayList<String> fragTypes = new ArrayList<>();

                        if (Constants.addNonYb) {
                            float maxIntensity = Constants.modelMaxIntensity.get(modelSplit[0]);
                            float maxIntensityMZ = Float.NaN;

                            //add original peaks
                            for(int i = 0; i < pe.mzs.length; i++) {
                                float mz = pe.mzs[i];
                                float intensity = pe.intensities[i];
                                String fragType = pe.fragmentIonTypes[i];

                                if (intensity == maxIntensity) {
                                    maxIntensityMZ = mz;
                                }

                                mzs.add(mz);
                                intensities.add(intensity);
                                fragTypes.add(fragType);
                            }

                            float minMZ = maxIntensityMZ - Constants.DaTolerance;
                            float maxMZ = maxIntensityMZ + Constants.DaTolerance;

                            //add new peaks
                            //Scale so that max intensity fragment of diann has same intensity as matched fragment in predfull
                            //TODO: multiply pe2 intensity by (diann max intensity / predfull intensity of matching fragment)
                            PredictionEntry pe2 = predictedSpectra2.getPreds().get(key);
                            //if null, convert to base format

                            if ((!Objects.isNull(pe2)) && (!Objects.isNull(pe2.fragmentIonTypes))) {
                                float matchedFragInt = Constants.modelMaxIntensity.get(modelSplit[1]);
                                for (int i = 0; i < pe2.mzs.length; i++) {
                                    float potentialMZ = pe2.mzs[i];
                                    float potentialInt = pe2.intensities[i];
                                    if ((potentialMZ >= minMZ) & (potentialMZ <= maxMZ) & (potentialInt > matchedFragInt)) {
                                        matchedFragInt = potentialInt;
                                    }
                                }

                                for (int i = 0; i < pe2.fragmentIonTypes.length; i++) {
                                    if (!model1FragmentIonTypes.contains(pe2.fragmentIonTypes[i])) {
                                        mzs.add(pe2.mzs[i]);
                                        intensities.add(pe2.intensities[i] * maxIntensity / matchedFragInt); //putting intensities on same scale
                                        fragTypes.add(pe2.fragmentIonTypes[i]);
                                    }
                                }
                            }

                            //convert back to array
                            float[] mzArray = new float[mzs.size()];
                            float[] intArray = new float[intensities.size()];
                            String[] typeArray = new String[fragTypes.size()];
                            for (int i = 0; i < mzArray.length; i++) {
                                mzArray[i] = mzs.get(i);
                                intArray[i] = intensities.get(i);
                                typeArray[i] = fragTypes.get(i);
                            }

                            PredictionEntry newPe = new PredictionEntry(mzArray, intArray,
                                    pe.getFragNums(), pe.getCharges(), typeArray);
                            allPreds.put(key, newPe);
                        } else { //retain predfull intensities, just add RT from other model
                            //but if predfull has missing entry, use other model instead
                            PredictionEntry pe2 = predictedSpectra2.getPreds().get(key);
                            if (!Objects.isNull(pe2)) {
                                pe2.setRT(pe.getRT());
                                allPreds.put(key, pe2);
                            }
                        }
                    }
                }
                predictedSpectra2 = null; //free up memory
            }
        }
        System.out.println();
        predictedSpectra.setPreds(predictedSpectra.filterTopFragments(executorService));
        allPreds = predictedSpectra.getPreds();
        //TODO test to make sure predicted spectra and allpreds are same

        //create detectMap to store detectabilities for base sequence peptides
        //store peptide detectabilities in PredictionEntry
        DetectMap dm = null;
        ArrayList<String> dFeatures = new ArrayList<String>(Constants.detectFeatures);
        dFeatures.retainAll(featuresList);
        //long startTime = System.nanoTime();
        if (dFeatures.size() > 0) {
            dm = new DetectMap(detectFile);
            for (Map.Entry<String, PredictionEntry> e : allPreds.entrySet()) {
                e.getValue().setDetectability(dm.getDetectability(
                        new PeptideFormatter(e.getKey().split("\\|")[0], e.getKey().split("\\|")[1], "pin").stripped));
            }
        }

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
                        if (Constants.peptideCounter.containsKey(pf.stripped)) {
                            Constants.peptideCounter.put(pf.stripped,
                                    Constants.peptideCounter.get(pf.stripped) + 1);
                        } else {
                            Constants.peptideCounter.put(pf.stripped, 1);
                        }
                    }
                    pin.close();
                }
            }

//            //load fasta
//            if (Constants.getFastaReader() == null) {
//                System.out.println("Creating fasta object");
//                fasta = new FastaReader(Constants.fasta);
//            } else {
//                System.out.println("Loading fasta");
//                fasta = Constants.getFastaReader();
//            }
//
//            System.out.println("Loading detectabilities for unique peptides from each protein");
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
                System.out.println("Processing pin " + pin.name);
                if (Constants.removeWashGradient) {
                    if (Constants.rtCutoff.isNaN()) {
                        pin.attachMzml(mzml);
                        pin.findWashGradientCutoff();
                    } else {
                        pin.rtCutoff = Constants.rtCutoff;
                    }
                }

                //Special preparations dependent on features we require
                mzml.setPinEntries(pin, predictedSpectra, executorService);
                //these require all experimental peaks before removing higher rank peaks
                if (Constants.removeRankPeaks &&
                        (featuresList.contains("hypergeometricProbability") ||
                                featuresList.contains("intersection") ||
                                featuresList.contains("adjacentSimilarity"))) {
                    for (MzmlScanNumber msn : mzml.scanNumberObjects.values()) {
                        msn.expMZs = msn.savedExpMZs;
                        msn.expIntensities = msn.savedExpIntensities;
                        msn.savedExpMZs = null;
                        msn.savedExpIntensities = null;
                    }
                    Constants.removeRankPeaks = false;
                }

                if (featuresList.contains("adjacentSimilarity")) {
                    System.out.println("Calculating adjacent similarity");
                    //check that scan m/z ranges are in mzml
                    if (mzml.scanNumberObjects.firstEntry().getValue().isolationUpper == 0.0d) {
                        System.out.println("Adjacent similarity feature requires m/z ranges " +
                                "for each scan. This mzml is incompatible. Exiting.");
                        System.exit(1);
                    }

                    SortedSet<Double> scanNumbers = new TreeSet<>();
                    for (MzmlScanNumber msn : mzml.scanNumberObjects.values()) {
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
                    for (MzmlScanNumber msn : mzml.scanNumberObjects.values()) {
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
                    for (MzmlScanNumber msn : mzml.scanNumberObjects.values()) {
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
                                    msn = mzml.scanNumberObjects.get(scan);
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
                    mzml.setLOESS(Constants.RTregressionSize, Constants.rtBandwidth, Constants.robustIters, "RT");
                    mzml.predictRTLOESS(executorService); //potentially only invoke once if normalized included

                    //generate calibration figure, need mzml and loess
                    if (! Constants.noRTscores) {
                        new RTCalibrationFigure(mzml, pinFiles[i].getCanonicalPath(), 0.2f);
                    }
                }
                if (featuresList.contains("deltaRTlinear")) {
                    if (mzml.expAndPredRTs != null) {
                        mzml.setBetas();
                    } else { mzml.setBetas(predictedSpectra, Constants.RTregressionSize);
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
                    mzml.unifProb = 1.0f / predictedSpectra.getMaxPredRT(); //TODO: might need to change, if this is ever revisited. Need negative range
                }
                if (featuresList.contains("RTprobability") || featuresList.contains("RTprobabilityUnifPrior")) {
                    mzml.setKernelDensities(executorService, "RT");
                }
                if (featuresList.contains("deltaIMLOESS") || featuresList.contains("deltaIMLOESSnormalized")) {
                    mzml.setLOESS(Constants.IMregressionSize, Constants.imBandwidth, Constants.robustIters, "IM");
                    mzml.predictIMLOESS(executorService);
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
                System.out.println("Calculating features");
                FeatureCalculator fc = new FeatureCalculator(pin, featuresList, mzml);
                fc.calculate(executorService);

                System.out.println("Writing features");
                PinWriter pw = new PinWriter(newOutfile, pin, featuresList, mzml, fc.featureStats);
                pw.write();
                String histFile;
                if (Constants.renamePin == 1) {
                    System.out.println("Edited pin file at " + newOutfile);
                    System.out.println();
                    histFile = newOutfile;
                } else { //really should be 0
                    //move file at newOutfile to pinFiles[i] canonical name
                    File movedFile = new File(newOutfile);
                    pinFiles[i].delete();
                    movedFile.renameTo(pinFiles[i]);
                    histFile = pinFiles[i].getCanonicalPath();
                }

                //plot hist
                PinReader pinReader = new PinReader(histFile);
                new ScoreHistogram(pinReader, featuresList);
            }
        } catch (Exception e) {
            executorService.shutdown();
            e.printStackTrace();
            System.exit(1);
        }
    }
}
