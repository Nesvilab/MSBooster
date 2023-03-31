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

import com.univocity.parsers.tsv.TsvWriter;
import com.univocity.parsers.tsv.TsvWriterSettings;
import org.apache.commons.lang3.ArrayUtils;
import org.apache.commons.math3.stat.correlation.PearsonsCorrelation;
import org.apache.commons.math3.stat.correlation.SpearmansCorrelation;
import umich.ms.fileio.exceptions.FileParsingException;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.sql.SQLException;
import java.util.*;
import java.util.concurrent.ExecutionException;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.stream.IntStream;

import static Features.Constants.camelToUnderscore;

public class percolatorFormatter {

    //set mgf or detectFile as null if not applicable
    //baseNames is the part before mzml or pin extensions
    public static void editPin(String pinDirectory, String mzmlDirectory, String mgf, String detectFile,
                               String[] features, String outfile)
            throws IOException, InterruptedException, ExecutionException, FileParsingException, SQLException {

        PinMzmlMatcher pmMatcher = new PinMzmlMatcher(mzmlDirectory, pinDirectory);

        editPin(pmMatcher, mgf, detectFile, features, outfile);
    }

    public static void editPin(PinMzmlMatcher pmMatcher, String mgf, String detectFile,
                               String[] features, String outfile)
            throws IOException, InterruptedException, ExecutionException, FileParsingException, SQLException {
        ExecutorService executorService = Executors.newFixedThreadPool(Constants.numThreads);

        //defining num threads, in case using this outside of jar file
        Runtime run  = Runtime.getRuntime();
        if (Constants.numThreads <= 0) {
            Constants.numThreads = run.availableProcessors();
        }

        ArrayList<String> featuresList = new ArrayList<>(Arrays.asList(features));
        //remove features from array for multiple protein formatting
        if (featuresList.contains("detectability")) {
            int idx = ArrayUtils.indexOf(features, "detectability");
            features = ArrayUtils.remove(features, idx);
        }

        //booleans for future determination of what to do
        boolean needsMGF = false;

        File[] pinFiles = pmMatcher.pinFiles;
        File[] mzmlFiles = pmMatcher.mzmlFiles;

        //load predicted spectra
        SpectralPredictionMapper predictedSpectra = null;
        SpectralPredictionMapper predictedSpectra2 = null; //first is prosit/diann, second predfull

        //Special preparations dependent on features we require
        //only time this isn't needed is detect
        //could consider an mgf constant
        if (mgf != null) {
            needsMGF = true;
            String[] mgfSplit = mgf.split(",");

            if (mgfSplit.length == 1) {
                System.out.println("Loading predicted spectra");
                if (Constants.spectraRTPredModel.equals("PredFull")) {
                    predictedSpectra = SpectralPredictionMapper.createSpectralPredictionMapper(
                            mgfSplit[1], pinFiles, executorService);
                } else {
                    predictedSpectra = SpectralPredictionMapper.createSpectralPredictionMapper(
                            mgf, Constants.spectraRTPredModel, executorService);
                }
            } else if (mgfSplit.length == 2){ //if fragment from predfull is not y/b, add.
                                              //Prosit/diann is first, predfull second
                String[] modelSplit = Constants.spectraRTPredModel.split(",");

                System.out.println("Loading predicted spectra 1");
                predictedSpectra = SpectralPredictionMapper.createSpectralPredictionMapper(
                        mgfSplit[0], modelSplit[0], executorService);
                System.out.println("Loading predicted spectra 2");
                if (modelSplit[1].equals("PredFull")) {
                    predictedSpectra2 = SpectralPredictionMapper.createSpectralPredictionMapper(
                            mgfSplit[1], pinFiles, executorService);
                } else {
                    predictedSpectra2 = SpectralPredictionMapper.createSpectralPredictionMapper(
                            mgfSplit[1], modelSplit[1], executorService);
                }
                System.out.println("Merging spectral libraries");

                //get all possible keys from both preds1 and preds2
                Set<String> totalKeyset = new HashSet<String>();
                totalKeyset.addAll(predictedSpectra.getPreds().keySet());
                totalKeyset.addAll(predictedSpectra2.getPreds().keySet());
                float maxIntensity = Constants.modelMaxIntensity.get(modelSplit[0]);

                //check what fragment ion types have been predicted by model 1
                HashSet<String> model1FragmentIonTypes = new HashSet<>();
                for (PredictionEntry pe : predictedSpectra.getPreds().values()) {
                    if (pe.fragmentIonTypes == null) {
                        pe.setFragmentIonTypes();
                    }
                    model1FragmentIonTypes.addAll(Arrays.asList(pe.fragmentIonTypes));
                }

                for (String key : totalKeyset) {
//                    if (entry.getValue().fragmentIonTypes == null) { //something with base modifications that is never actually queried
//                        predictedSpectra.getPreds().remove(entry.getKey());
//                        //continue; //DIA-NN may still have the prediction
//                    }
                    PredictionEntry pe = predictedSpectra.getPreds().get(key);
                    if (pe == null) { //missing in prosit/diann
                        predictedSpectra.getPreds().put(key, predictedSpectra2.getPreds().get(key));
                    } else { //add non-y/b ions
                        ArrayList<Float> mzs = new ArrayList<>();
                        ArrayList<Float> intensities = new ArrayList<>();
                        ArrayList<String> fragTypes = new ArrayList<>();

                        if (Constants.replaceYBintensities) {
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

                            pe.setMzs(mzArray);
                            pe.setIntensities(intArray);
                            pe.setFragmentIonTypes(typeArray);
                            predictedSpectra.getPreds().put(key, pe);
                        } else { //retain predfull intensities, just add RT from other model
                            //but if predfull has missing entry, use other model instead
                            PredictionEntry pe2 = predictedSpectra2.getPreds().get(key);
                            if (!Objects.isNull(pe2)) {
                                pe.setMzs(pe2.mzs);
                                pe.setIntensities(pe2.intensities);
                                pe.setFragmentIonTypes(pe2.fragmentIonTypes);
                            }
                            predictedSpectra.getPreds().put(key, pe);
                        }

                        predictedSpectra2.getPreds().put(key, null);
                    }
                }
                predictedSpectra2 = null; //free up memory
            }
        }
        //create detectMap to store detectabilities for base sequence peptides
        //store peptide detectabilities in PredictionEntry
        detectMap dm = null;
        ArrayList<String> dFeatures = new ArrayList<String>(Constants.detectFeatures);
        dFeatures.retainAll(featuresList);
        //long startTime = System.nanoTime();
        if (dFeatures.size() > 0) {
            dm = new detectMap(detectFile);
            HashMap<String, PredictionEntry> allPreds = predictedSpectra.getPreds();
            for (Map.Entry<String, PredictionEntry> e : allPreds.entrySet()) {
                e.getValue().setDetectability(dm.getDetectability(
                        new PeptideFormatter(e.getKey().split("\\|")[0], e.getKey().split("\\|")[1], "pin").stripped));
            }
        }

        FastaReader fasta = null;

        if (featuresList.contains("detectFractionGreater") || featuresList.contains("detectSubtractMissing")
                || featuresList.contains("detectProtSpearmanDiff")) {

            HashMap<String, Integer> pepCounter = new HashMap<>();

            //get all peptides present in pin
            for (File pinFile : pmMatcher.pinFiles) {
                pinReader pin = new pinReader(pinFile.getCanonicalPath());

                //add to counter
                while (pin.next()) {
                    String pep = pin.getPep().base;
                    if (pepCounter.containsKey(pep)) {
                        pepCounter.put(pep, pepCounter.get(pep) + 1);
                    } else {
                        pepCounter.put(pep, 1);
                    }
                }
                pin.close();
            }

            //load fasta
            if (Constants.getFastaReader() == null) {
                System.out.println("Creating fasta object");
                fasta = new FastaReader(Constants.fasta);
            } else {
                System.out.println("Loading fasta");
                fasta = Constants.getFastaReader();
            }

            System.out.println("Loading detectabilities for unique peptides from each protein");
            for (Map.Entry<String, ProteinEntry> e : fasta.protToPep.entrySet()) {
                ArrayList<String> pepList = e.getValue().peptides;
                float[] protDetects = new float[pepList.size()]; //for storing initial detect order

                //store detect unsorted
                for (int pep = 0; pep < pepList.size(); pep++) {
                    protDetects[pep] = dm.getDetectability(pepList.get(pep));
                }

                //dual pivot quicksort
                //sorted indices
                int[] sortedIndices = IntStream.range(0, protDetects.length)
                        .boxed().sorted((k, j) -> Float.compare(protDetects[k], protDetects[j]))
                        .mapToInt(ele -> ele).toArray();

                float[] sortedDetect = new float[protDetects.length];
                for (int j = 0; j < protDetects.length; j++) {
                    sortedDetect[j] = protDetects[sortedIndices[j]];
                }
                e.getValue().detects = sortedDetect;

                //check which peptides present, and get spectral counts
                float[] protPresence = new float[protDetects.length];
                float[] pepCounts = new float[protDetects.length];
                float numPresent = 1f;
                for (int j = protDetects.length - 1; j > -1; j--) {
                    String currentPep = pepList.get(sortedIndices[j]);

                    if (pepCounter.containsKey(currentPep)) {
                        protPresence[j] = numPresent;
                        numPresent += 1f;
                        pepCounts[j] = pepCounter.get(currentPep);
                    }
                }
                fasta.protToPep.get(e.getKey()).presence = protPresence;
                fasta.protToPep.get(e.getKey()).spectralCounts = pepCounts;
            }

            HashMap<String, PredictionEntry> allPreds = predictedSpectra.getPreds();
            for (Map.Entry<String, PredictionEntry> e : allPreds.entrySet()) {
                try {
                    e.getValue().setCounter(pepCounter.get(e.getKey().split("\\|")[0]));
                } catch (Exception ee) { //peptide was in a pin file from another run
                }
            }
            dm.clear();
        }
        //long endTime = System.nanoTime();
        //long duration = (endTime - startTime);
        //System.out.println("Detectability map and formatting loading took " + duration / 1000000 +" milliseconds");

        try {
            //////////////////////////////iterate through pin and mzml files//////////////////////////////////////////
            for (int i = 0; i < pinFiles.length; i++) {
                //startTime = System.nanoTime();
                String newOutfile = pinFiles[i].getAbsolutePath().replaceAll("\\.pin$", "_" + outfile + ".pin");

                TsvWriterSettings tws = new TsvWriterSettings();
                tws.setMaxCharsPerColumn(-1);
                tws.setMaxColumns(50000); //who knows if it needs to be longer?
                TsvWriter writer = new TsvWriter(new File(newOutfile), tws);
                //load mzml file
                System.out.println("Processing " + mzmlFiles[i].getName());

                mzMLReader mzml;
                if (mzmlFiles[i].getName().substring( mzmlFiles[i].getName().length() - 3).toLowerCase().equals("mgf")) {
                    //mzml = new mzMLReader(new mgfFileReader(mzmlFiles[i].getCanonicalPath()));
                    mzml = new mzMLReader(new mgfFileReader(mzmlFiles[i].getCanonicalPath(), true, executorService, ""));
                    //endTime = System.nanoTime();
                    //duration = (endTime - startTime);
                    //System.out.println("mgf loading took " + duration / 1000000 +" milliseconds");
                } else {
                    mzml = new mzMLReader(mzmlFiles[i].getCanonicalPath());
                    //endTime = System.nanoTime();
                    //duration = (endTime - startTime);
                    //System.out.println("mzML loading took " + duration / 1000000 +" milliseconds");
                }

                //load pin file, which already includes all ranks
                pinReader pin = new pinReader(pinFiles[i].getCanonicalPath());

                //add header to written tsv
                ArrayList<String> newHeader = new ArrayList<>();
                newHeader.addAll(Arrays.asList(pin.header));
                    //replace column names
                ArrayList<String> newNames = new ArrayList<>(featuresList.size());
                for (String s : featuresList) {
                    String newName = s;
                    if (camelToUnderscore.containsKey(s)) {
                        newName = camelToUnderscore.get(s);
                    }
                    //add columns for spectral features divided by fragment ion type
                    if (Constants.spectraFeatures.contains(s)) {
                        if (! Constants.divideFragments.equals("0")) {
                            String[] divisions = Constants.divideFragments.split(";");
                            for (String div : divisions) {
                                newNames.add(newName + "_" + div);
                            }
                        } else {
                            newNames.add(newName);
                        }
                    } else {
                            newNames.add(newName);
                        }
                }
                newHeader.addAll(pin.pepIdx, newNames); //add features before Peptide
                newHeader.remove("detectability");
                writer.writeHeaders(newHeader);

                //Special preparations dependent on features we require
                if (needsMGF) {
                    mzml.setPinEntries(pin, predictedSpectra);
                }
                if (featuresList.contains("deltaRTLOESS") || featuresList.contains("deltaRTLOESSnormalized")) {
                    mzml.setLOESS(Constants.RTregressionSize, Constants.bandwidth, Constants.robustIters, "RT");
                    mzml.predictRTLOESS(executorService); //potentially only invoke once if normalized included

                    //generate calibration figure, need mzml and loess
                    if (! Constants.noRTscores) {
                        RTCalibrationFigure rtfc = new RTCalibrationFigure(mzml, pinFiles[i].getCanonicalPath(), 0.2f);
                    }
                }
                if (featuresList.contains("deltaRTlinear")) {
                    //System.out.println("Calculating delta RT linear");
                    if (mzml.expAndPredRTs != null) {
                        mzml.setBetas();
                    } else { mzml.setBetas(predictedSpectra, Constants.RTregressionSize);
                    }
                    mzml.normalizeRTs(executorService);
                }
                if (featuresList.contains("deltaRTbins") || featuresList.contains("RTzscore") ||
                        featuresList.contains("RTprobability") || featuresList.contains("RTprobabilityUnifPrior") ||
                        featuresList.contains("deltaRTLOESSnormalized")) {
                    //System.out.println("Generating RT bins");
                    mzml.setRTbins();
                    mzml.calculateBinStats("RT");
                }
                if (featuresList.contains("deltaRTbins") || featuresList.contains("RTzscore")) {
                    mzml.calculateDeltaRTbinAndRTzscore(executorService);
                }
                if (featuresList.contains("deltaRTLOESSnormalized")) {
                    mzml.calculateDeltaRTLOESSnormalized(executorService);
                }
                int unifPriorSize = 0; //only used if using uniform prior
                float unifProb = 0;
                if (featuresList.contains("RTprobabilityUnifPrior")) {
                    mzml.setRTBinSizes(executorService);

                    //decide how many uniform points to add
                    int[] binSizes = new int[mzml.RTbins.length];
                    for (int bin = 0; bin < mzml.RTbins.length; bin++) {
                        binSizes[bin] = mzml.RTbins[bin].size();
                    }
                    Arrays.sort(binSizes);
                    int cutoff = (int) Math.floor(((double) mzml.RTbins.length / 100.0) * Constants.uniformPriorPercentile);
                    unifPriorSize = binSizes[cutoff];

                    //also need uniform probability with right bound of max predicted RT
                    unifProb = 1.0f / predictedSpectra.getMaxPredRT(); //TODO: might need to change, if this is ever revisited. Need negative range
                }
                if (featuresList.contains("RTprobability") || featuresList.contains("RTprobabilityUnifPrior")) {
                    //System.out.println("Generating RT empirical densities");
                    mzml.setKernelDensities(executorService, "RT");
                }
                if (featuresList.contains("deltaIMLOESS") || featuresList.contains("deltaIMLOESSnormalized")) {
                    //System.out.println("Generating IM LOESS regression");
                    mzml.setLOESS(Constants.IMregressionSize, Constants.bandwidth, Constants.robustIters, "IM");
                    mzml.predictIMLOESS(executorService);
                }
                if (featuresList.contains("deltaIMLOESSnormalized") || featuresList.contains("IMprobabilityUnifPrior")) {
                    //System.out.println("Generating IM bins");
                    mzml.setIMbins();
                    mzml.calculateBinStats("IM");
                }
                if (featuresList.contains("deltaIMLOESSnormalized")) {
                    mzml.calculateDeltaIMLOESSnormalized(executorService);
                }
                float[] unifProbIM = new float[IMFunctions.numCharges];
                int[] unifPriorSizeIM = new int[IMFunctions.numCharges];;
                if (featuresList.contains("IMprobabilityUnifPrior")) {
                    mzml.setIMBinSizes(executorService);

                    //decide how many uniform points to add
                    //System.out.println("Generating IM empirical densities");
                    for (int charge = 0; charge < IMFunctions.numCharges; charge++) {
                        int[] binSizes = new int[mzml.IMbins[charge].length];
                        for (int bin = 0; bin < mzml.IMbins[charge].length; bin++) {
                            binSizes[bin] = mzml.IMbins[charge][bin].size();
                        }
                        Arrays.sort(binSizes);
                        int cutoff = (int) Math.floor(((double) mzml.IMbins[charge].length / 100.0) * Constants.uniformPriorPercentile);
                        unifPriorSizeIM[charge] = binSizes[cutoff];

                        //also need uniform probability with right bound of max predicted RT
                        unifProbIM[charge] = 1.0f / (2 * Constants.IMbinMultiplier);

                        mzml.setKernelDensities(executorService, "IM");
                    }
                }

                //System.out.println("Getting predictions for each row");
                //int totalPSMs = 0;

                featuresList.remove("detectability");
                while (pin.next()) {
                    //totalPSMs += 1;
                    //peptide name
                    String pep = pin.getPep().baseCharge;

                    //trying filtering out low detectability
//                    if (featuresList.contains("detectability")) {
//                        if (dm.getDetectability(pep) < Constants.detectThreshold) {
//                            continue;
//                        }
//                    }

                    //get entry
                    peptideObj pepObj = null;
                    if (needsMGF) {
                        pepObj = mzml.scanNumberObjects.get(pin.getScanNum()).getPeptideObject(pep);
                    }

                    //write everything we already have, not including extra protein columns
                    String[] row = pin.getRow();
                    for (int j = 0; j < pin.header.length; j++) {
                        writer.addValue(pin.header[j], row[j]);
                    }
                    //add extra protein columns
                    if (pin.getRow().length > pin.header.length) {
                        for (int j = pin.header.length; j < pin.getRow().length; j++) {
                            writer.addValue(newHeader.size(), row[j]);
                        }
                    }

                    //switch case
                    for (String feature : featuresList) {
                        switch (feature) {
                            case "detectFractionGreater":
                                float d = predictedSpectra.getPreds().get(pep).detectability;
                                //for each protein, get the position of pep's detect and see how many peptides with greater detect are present
                                //take max (proxy for protein that actually generated peptide)
                                String[] r = pin.getRow();
                                String[] prots = Arrays.copyOfRange(r, pin.pepIdx + 1, r.length);
                                float maxFraction = 0f;
                                for (String prot : prots) { //if more than one, this peptide is shared among proteins
                                    String protAbr;
                                    //skip protein if it is decoy and looking at target peptide
                                    if (prot.startsWith(Constants.decoyPrefix.substring(1))) {
                                        if (r[pin.labelIdx].equals("1")) {
                                            continue;
                                        } else { //decoy peptide compared to target protein
                                            protAbr = prot.substring(Constants.decoyPrefix.length() - 1);
                                        }
                                    } else {
                                        protAbr = prot;
                                    }

                                    float[] arr;
                                    try {
                                        arr = fasta.protToPep.get(protAbr).detects;
                                    } catch (Exception e) { //no peptides qualify from this protein
                                        continue;
                                    }

                                    int idx = Arrays.binarySearch(arr, d);
                                    if (idx < 0) { //not found
                                        idx = (-1 * idx) - 1;
                                    } else {
                                        idx += 1; //don't want to include itself in calculation
                                    }
                                    float[] presenceArr = Arrays.copyOfRange(fasta.protToPep.get(protAbr).presence, idx, arr.length);
                                    float total = 0f;
                                    for (float j : presenceArr) {
                                        if (j != 0f) {
                                            total = j;
                                            break;
                                        }
                                    }
                                    float fraction = (total + Constants.detectFractionGreaterNumerator) /
                                            (presenceArr.length + Constants.detectFractionGreaterDenominator); //customizable prior
                                    if (fraction > maxFraction) {
                                        maxFraction = fraction;
                                    }
                                }
                                writer.addValue("detect_fraction_greater", maxFraction);
                                break;
                            case "detectSubtractMissing":
                                d = predictedSpectra.getPreds().get(pep).detectability;
                                //for each protein, get the position of pep's detect and see how many peptides with greater detect are present
                                //take max (proxy for protein that actually generated peptide)
                                r = pin.getRow();
                                prots = Arrays.copyOfRange(r, pin.pepIdx + 1, r.length);
                                float minDiff = 1f;
                                for (String prot : prots) { //if more than one, this peptide is shared among proteins
                                    String protAbr;
                                    //skip protein if it is decoy and looking at target peptide
                                    if (prot.startsWith(Constants.decoyPrefix.substring(1))) {
                                        if (r[pin.labelIdx].equals("1")) {
                                            continue;
                                        } else { //decoy
                                            protAbr = prot.substring(Constants.decoyPrefix.length() - 1);
                                        }
                                    } else {
                                        protAbr = prot;
                                    }

                                    float[] arr;
                                    try {
                                        arr = fasta.protToPep.get(protAbr).detects;
                                    } catch (Exception e) { //no peptides qualify from this protein
                                        continue;
                                    }

                                    int idx = Arrays.binarySearch(arr, d);
                                    if (idx < 0) { //not found
                                        idx = (-1 * idx) - 1;
                                    } else {
                                        idx += 1; //don't want to include itself in calculation
                                    }
                                    float[] presenceArr = Arrays.copyOfRange(fasta.protToPep.get(protAbr).presence, idx, arr.length);
                                    float[] detectArr = Arrays.copyOfRange(arr, idx, arr.length);
                                    float total = 0f;
                                    for (int k = 0; k < presenceArr.length; k++) {
                                        if (presenceArr[k] == 0) {
                                            total += detectArr[k] - d;
                                        }
                                    }
                                    float diff = total / presenceArr.length;
                                    if (diff < minDiff) {
                                        minDiff = diff;
                                    }
                                    if (minDiff == 0) {
                                        break;
                                    }
                                }
                                writer.addValue("detect_subtract_missing", minDiff);
                                break;
                            case "detectProtSpearmanDiff":
                                SpearmansCorrelation sc = new SpearmansCorrelation();
                                r = pin.getRow();
                                prots = Arrays.copyOfRange(r, pin.pepIdx + 1, r.length);
                                double maxSpearmanDiff = -3;
                                float detect = predictedSpectra.getPreds().get(pep).detectability;;
                                for (String prot : prots) { //if more than one, this peptide is shared among proteins
                                    //skip protein if it is decoy and looking at target peptide
                                    String protAbr;
                                    if (prot.startsWith(Constants.decoyPrefix.substring(1))) {
                                        if (r[pin.labelIdx].equals("1")) {
                                            continue;
                                        } else {
                                            protAbr = prot.substring(Constants.decoyPrefix.length() - 1);
                                        }
                                    } else {
                                        protAbr = prot;
                                    }

                                    float[] arr;
                                    try {
                                        arr = fasta.protToPep.get(protAbr).detects;
                                    } catch (Exception e) { //no peptides qualify from this protein
                                        continue;
                                    }

                                    float[] counts = fasta.protToPep.get(protAbr).spectralCounts;

                                    //only add if not 0 spectral counts (lots of missing ones, need to be more lenient for targets)
                                    ArrayList<Double> newDetects = new ArrayList<>();
                                    ArrayList<Double> newCounts = new ArrayList<>();
                                    for (int k = 0; k < arr.length; k++) {
                                        if (counts[k] != 0) {
                                            if (! (arr[k] == detect)) { //will add later
                                                newDetects.add((double) arr[k]);
                                                newCounts.add((double) counts[k]);
                                            }
                                        }
                                    }
                                    if (newDetects.size() < 2) {
                                        continue;
                                    }
                                    double spear = sc.correlation(newDetects.stream().mapToDouble(dd -> dd).toArray(),
                                            newCounts.stream().mapToDouble(dd -> dd).toArray() );

                                    //add new pep to this calculation
                                    newDetects.add((double) detect);
                                    newCounts.add((double) predictedSpectra.getPreds().get(pep).counter);
                                    double spearDiff = sc.correlation(newDetects.stream().mapToDouble(dd -> dd).toArray(),
                                            newCounts.stream().mapToDouble(dd -> dd).toArray() ) - spear;
                                    if (spearDiff > maxSpearmanDiff) {
                                        maxSpearmanDiff = spearDiff;
                                    }
                                }
                                if (maxSpearmanDiff == -3) {
                                    maxSpearmanDiff = 0;
                                }
                                writer.addValue("detect_prot_spearman_diff", maxSpearmanDiff);
                                break;
                            case "deltaRTlinear":
                                if (Constants.noRTscores) {
                                    pepObj.deltaRT = 0;
                                }
                                writer.addValue("deltaRTlinear", pepObj.deltaRT);
                                break;
                            case "deltaRTbins":
                                if (Constants.noRTscores) {
                                    pepObj.deltaRTbin = 0;
                                }
                                writer.addValue("deltaRTbins", pepObj.deltaRTbin);
                                break;
                            case "deltaRTLOESS":
                                if (Constants.noRTscores) {
                                    pepObj.deltaRTLOESS = 0;
                                }
                                writer.addValue("delta_RT_loess", pepObj.deltaRTLOESS);
                                break;
                            case "deltaRTLOESSnormalized":
                                if (Constants.noRTscores) {
                                    pepObj.deltaRTLOESSnormalized = 0;
                                }
                                writer.addValue("delta_RT_loess_normalized", pepObj.deltaRTLOESSnormalized);
                                break;
                            case "RTzscore":
                                if (Constants.noRTscores) {
                                    pepObj.RTzscore = 0;
                                }
                                writer.addValue("RTzscore", pepObj.RTzscore);
                                break;
                            case "RTprobability":
                                if (Constants.noRTscores) {
                                    pepObj.RTprob = 0;
                                }
                                writer.addValue("RTprobability", pepObj.RTprob);
                                break;
                            case "RTprobabilityUnifPrior":
                                if (Constants.noRTscores) {
                                    writer.addValue("RT_probability_unif_prior", 0);
                                } else {
                                    float prob = StatMethods.probabilityWithUniformPrior(unifPriorSize, unifProb,
                                            pepObj.scanNumObj.RTbinSize, (float) pepObj.RTprob);
                                    writer.addValue("RT_probability_unif_prior", prob);
                                }
                                break;
                            case "calibratedRT":
                                writer.addValue("calibrated_RT", pepObj.calibratedRT);
                                break;
                            case "predictedRT":
                                writer.addValue("predicted_RT", pepObj.RT);
                                break;
                            case "brayCurtis":
                                if (pepObj.spectralSimObj.spectrumComparisons.size() == 0) {
                                    writer.addValue("bray_curtis", pepObj.spectralSimObj.brayCurtis());
                                } else {
                                    String[] dividedFragments = Constants.divideFragments.split(";");
                                    for (int j = 0; j < dividedFragments.length; j++) {
                                        writer.addValue("bray_curtis_" + dividedFragments[j],
                                                pepObj.spectralSimObj.spectrumComparisons.get(j).brayCurtis());
                                    }
                                }
                                break;
                            case "cosineSimilarity":
                                if (pepObj.spectralSimObj.spectrumComparisons.size() == 0) {
                                    writer.addValue("cosine_similarity", pepObj.spectralSimObj.cosineSimilarity());
                                } else {
                                    String[] dividedFragments = Constants.divideFragments.split(";");
                                    for (int j = 0; j < dividedFragments.length; j++) {
                                        writer.addValue("cosine_similarity_" + dividedFragments[j],
                                                pepObj.spectralSimObj.spectrumComparisons.get(j).cosineSimilarity());
                                    }
                                }
                                break;
                            case "spectralContrastAngle":
                                if (pepObj.spectralSimObj.spectrumComparisons.size() == 0) {
                                    writer.addValue("spectra_contrast_angle", pepObj.spectralSimObj.spectralContrastAngle());
                                } else {
                                    String[] dividedFragments = Constants.divideFragments.split(";");
                                    for (int j = 0; j < dividedFragments.length; j++) {
                                        writer.addValue("spectra_contrast_angle_" + dividedFragments[j],
                                                pepObj.spectralSimObj.spectrumComparisons.get(j).spectralContrastAngle());
                                    }
                                }
                                break;
                            case "euclideanDistance":
                                if (pepObj.spectralSimObj.spectrumComparisons.size() == 0) {
                                    writer.addValue("euclidean_distance", pepObj.spectralSimObj.euclideanDistance());
                                } else {
                                    String[] dividedFragments = Constants.divideFragments.split(";");
                                    for (int j = 0; j < dividedFragments.length; j++) {
                                        writer.addValue("euclidean_distance_" + dividedFragments[j],
                                                pepObj.spectralSimObj.spectrumComparisons.get(j).euclideanDistance());
                                    }
                                }
                                break;
                            case "pearsonCorr":
                                if (pepObj.spectralSimObj.spectrumComparisons.size() == 0) {
                                    writer.addValue("pearson_corr", pepObj.spectralSimObj.pearsonCorr());
                                } else {
                                    String[] dividedFragments = Constants.divideFragments.split(";");
                                    for (int j = 0; j < dividedFragments.length; j++) {
                                        writer.addValue("pearson_corr_" + dividedFragments[j],
                                                pepObj.spectralSimObj.spectrumComparisons.get(j).pearsonCorr());
                                    }
                                }
                                break;
                            case "dotProduct":
                                if (pepObj.spectralSimObj.spectrumComparisons.size() == 0) {
                                    writer.addValue("dot_product", pepObj.spectralSimObj.dotProduct());
                                } else {
                                    String[] dividedFragments = Constants.divideFragments.split(";");
                                    for (int j = 0; j < dividedFragments.length; j++) {
                                        writer.addValue("dot_product_" + dividedFragments[j],
                                                pepObj.spectralSimObj.spectrumComparisons.get(j).dotProduct());
                                    }
                                }
                                break;
                            case "unweightedSpectralEntropy":
                                if (Constants.divideFragments.equals("0")) {
                                    writer.addValue("unweighted_spectral_entropy", pepObj.spectralSimObj.unweightedSpectralEntropy());
                                } else if ((pepObj.spectralSimObj.spectrumComparisons.size() > 0) & (predictedSpectra.getPreds().containsKey(pepObj.name))){
                                    String[] dividedFragments = Constants.divideFragments.split(";");
                                    for (int j = 0; j < dividedFragments.length; j++) {
                                        writer.addValue("unweighted_spectral_entropy_" + dividedFragments[j],
                                                pepObj.spectralSimObj.spectrumComparisons.get(j).unweightedSpectralEntropy());
                                    }
                                } else if (pepObj.spectralSimObj.spectrumComparisons.size() == 0) {
                                    String[] dividedFragments = Constants.divideFragments.split(";");
                                    for (int j = 0; j < dividedFragments.length; j++) {
                                        writer.addValue("unweighted_spectral_entropy_" + dividedFragments[j], 0);
                                    }
                                } else {
                                    System.out.println("Something wrong with feature calculation");
                                    System.exit(-1);
                                }
                                break;
                            case "numMatchedFragments":
                                float[] matchedIntensities = pepObj.spectralSimObj.matchedIntensities;
                                int v = matchedIntensities.length;
                                for (float f : matchedIntensities) {
                                    if (f == 0) {
                                        v -= 1;
                                    }
                                }
                                writer.addValue("num_matched_fragments", v);
                                break;
                            case "deltaIMLOESS":
                                writer.addValue("delta_IM_loess", pepObj.deltaIMLOESS);
                                break;
                            case "deltaIMLOESSnormalized":
                                writer.addValue("delta_IM_loess_normalized", pepObj.deltaIMLOESSnormalized);
                                break;
                            case "IMprobabilityUnifPrior":
                                float prob = StatMethods.probabilityWithUniformPrior(unifPriorSizeIM[pepObj.charge - 1], unifProbIM[pepObj.charge - 1],
                                        pepObj.scanNumObj.IMbinSize, (float) pepObj.IMprob);
                                writer.addValue("IM_probability_unif_prior", prob);
                                break;
                            case "predictedIM":
                                writer.addValue("predicted_IM", pepObj.IM);
                                break;
                            case "ionmobility":
                                writer.addValue("ionmobility", pepObj.scanNumObj.IM);
                                break;
                            case "y_matched_intensity":
                                writer.addValue("y_matched_intensity", pepObj.matchedIntensities.get("y"));
                                break;
                            case "b_matched_intensity":
                                writer.addValue("b_matched_intensity", pepObj.matchedIntensities.get("b"));
                                break;
                            case "a_matched_intensity":
                                writer.addValue("a_matched_intensity", pepObj.matchedIntensities.get("a"));
                                break;
                            case "x_matched_intensity":
                                writer.addValue("x_matched_intensity", pepObj.matchedIntensities.get("x"));
                                break;
                            case "c_matched_intensity":
                                writer.addValue("c_matched_intensity", pepObj.matchedIntensities.get("c"));
                                break;
                            case "z_matched_intensity":
                                writer.addValue("z_matched_intensity", pepObj.matchedIntensities.get("z"));
                                break;
                            case "cdot_matched_intensity":
                                writer.addValue("cdot_matched_intensity", pepObj.matchedIntensities.get("cdot"));
                                break;
                            case "zdot_matched_intensity":
                                writer.addValue("zdot_matched_intensity", pepObj.matchedIntensities.get("zdot"));
                                break;
                            case "y-NL_matched_intensity":
                                writer.addValue("y-NL_matched_intensity", pepObj.matchedIntensities.get("y-NL"));
                                break;
                            case "b-NL_matched_intensity":
                                writer.addValue("b-NL_matched_intensity", pepObj.matchedIntensities.get("b-NL"));
                                break;
                            case "a-NL_matched_intensity":
                                writer.addValue("a-NL_matched_intensity", pepObj.matchedIntensities.get("a-NL"));
                                break;
                            case "x-NL_matched_intensity":
                                writer.addValue("x-NL_matched_intensity", pepObj.matchedIntensities.get("x-NL"));
                                break;
                            case "c-NL_matched_intensity":
                                writer.addValue("c-NL_matched_intensity", pepObj.matchedIntensities.get("c-NL"));
                                break;
                            case "z-NL_matched_intensity":
                                writer.addValue("z-NL_matched_intensity", pepObj.matchedIntensities.get("z-NL"));
                                break;
                            case "precursor_matched_intensity":
                                writer.addValue("precursor_matched_intensity", pepObj.matchedIntensities.get("precursor"));
                                break;
                            case "precursor-NL_matched_intensity":
                                writer.addValue("precursor-NL_matched_intensity", pepObj.matchedIntensities.get("precursor-NL"));
                                break;
                            case "internal_matched_intensity":
                                writer.addValue("internal_matched_intensity", pepObj.matchedIntensities.get("internal"));
                                break;
                            case "internal-NL_matched_intensity":
                                writer.addValue("internal-NL_matched_intensity", pepObj.matchedIntensities.get("internal-NL"));
                                break;
                            case "immonium_matched_intensity":
                                writer.addValue("immonium_matched_intensity", pepObj.matchedIntensities.get("immonium"));
                                break;
                            case "unknown_matched_intensity":
                                writer.addValue("unknown_matched_intensity", pepObj.matchedIntensities.get("unknown"));
                                break;
                            case "y_intensities_difference":
                                writer.addValue("y_intensities_difference",
                                        Math.abs(pepObj.matchedIntensities.get("y") - pepObj.predIntensities.get("y")));
                                break;
                            case "b_intensities_difference":
                                writer.addValue("b_intensities_difference",
                                        Math.abs(pepObj.matchedIntensities.get("b") - pepObj.predIntensities.get("b")));
                                break;
                            case "a_intensities_difference":
                                writer.addValue("a_intensities_difference",
                                        Math.abs(pepObj.matchedIntensities.get("a") - pepObj.predIntensities.get("a")));
                                break;
                            case "x_intensities_difference":
                                writer.addValue("x_intensities_difference",
                                        Math.abs(pepObj.matchedIntensities.get("x") - pepObj.predIntensities.get("x")));
                                break;
                            case "c_intensities_difference":
                                writer.addValue("c_intensities_difference",
                                        Math.abs(pepObj.matchedIntensities.get("c") - pepObj.predIntensities.get("c")));
                                break;
                            case "z_intensities_difference":
                                writer.addValue("z_intensities_difference",
                                        Math.abs(pepObj.matchedIntensities.get("z") - pepObj.predIntensities.get("z")));
                                break;
                            case "cdot_intensities_difference":
                                writer.addValue("cdot_intensities_difference",
                                        Math.abs(pepObj.matchedIntensities.get("cdot") - pepObj.predIntensities.get("cdot")));
                                break;
                            case "zdot_intensities_difference":
                                writer.addValue("zdot_intensities_difference",
                                        Math.abs(pepObj.matchedIntensities.get("zdot") - pepObj.predIntensities.get("zdot")));
                                break;
                            case "y-NL_intensities_difference":
                                writer.addValue("y-NL_intensities_difference",
                                        Math.abs(pepObj.matchedIntensities.get("y-NL") - pepObj.predIntensities.get("y-NL")));
                                break;
                            case "b-NL_intensities_difference":
                                writer.addValue("b-NL_intensities_difference",
                                        Math.abs(pepObj.matchedIntensities.get("b-NL") - pepObj.predIntensities.get("b-NL")));
                                break;
                            case "a-NL_intensities_difference":
                                writer.addValue("a-NL_intensities_difference",
                                        Math.abs(pepObj.matchedIntensities.get("a-NL") - pepObj.predIntensities.get("a-NL")));
                                break;
                            case "x-NL_intensities_difference":
                                writer.addValue("x-NL_intensities_difference",
                                        Math.abs(pepObj.matchedIntensities.get("x-NL") - pepObj.predIntensities.get("x-NL")));
                                break;
                            case "c-NL_intensities_difference":
                                writer.addValue("c-NL_intensities_difference",
                                        Math.abs(pepObj.matchedIntensities.get("c-NL") - pepObj.predIntensities.get("c-NL")));
                                break;
                            case "z-NL_intensities_difference":
                                writer.addValue("z-NL_intensities_difference",
                                        Math.abs(pepObj.matchedIntensities.get("z-NL") - pepObj.predIntensities.get("z-NL")));
                                break;
                            case "precursor_intensities_difference":
                                writer.addValue("precursor_intensities_difference",
                                        Math.abs(pepObj.matchedIntensities.get("precursor") - pepObj.predIntensities.get("precursor")));
                                break;
                            case "precursor-NL_intensities_difference":
                                writer.addValue("precursor-NL_intensities_difference",
                                        Math.abs(pepObj.matchedIntensities.get("precursor-NL") - pepObj.predIntensities.get("precursor-NL")));
                                break;
                            case "internal_intensities_difference":
                                writer.addValue("internal_intensities_difference",
                                        Math.abs(pepObj.matchedIntensities.get("internal") - pepObj.predIntensities.get("internal")));
                                break;
                            case "internal-NL_intensities_difference":
                                writer.addValue("internal-NL_intensities_difference",
                                        Math.abs(pepObj.matchedIntensities.get("internal-NL") - pepObj.predIntensities.get("internal-NL")));
                                break;
                            case "immonium_intensities_difference":
                                writer.addValue("immonium_intensities_difference",
                                        Math.abs(pepObj.matchedIntensities.get("immonium") - pepObj.predIntensities.get("immonium")));
                                break;
                            case "unknown_intensities_difference":
                                writer.addValue("unknown_intensities_difference",
                                        Math.abs(pepObj.matchedIntensities.get("unknown") - pepObj.predIntensities.get("unknown")));
                                break;
                            case "y_peak_counts":
                                writer.addValue("y_peak_counts", pepObj.peakCounts.get("y"));
                                break;
                            case "b_peak_counts":
                                writer.addValue("b_peak_counts", pepObj.peakCounts.get("b"));
                                break;
                            case "a_peak_counts":
                                writer.addValue("a_peak_counts", pepObj.peakCounts.get("a"));
                                break;
                            case "x_peak_counts":
                                writer.addValue("x_peak_counts", pepObj.peakCounts.get("x"));
                                break;
                            case "c_peak_counts":
                                writer.addValue("c_peak_counts", pepObj.peakCounts.get("c"));
                                break;
                            case "z_peak_counts":
                                writer.addValue("z_peak_counts", pepObj.peakCounts.get("z"));
                                break;
                            case "cdot_peak_counts":
                                writer.addValue("cdot_peak_counts", pepObj.peakCounts.get("cdot"));
                                break;
                            case "zdot_peak_counts":
                                writer.addValue("zdot_peak_counts", pepObj.peakCounts.get("zdot"));
                                break;
                            case "y-NL_peak_counts":
                                writer.addValue("y-NL_peak_counts", pepObj.peakCounts.get("y-NL"));
                                break;
                            case "b-NL_peak_counts":
                                writer.addValue("b-NL_peak_counts", pepObj.peakCounts.get("b-NL"));
                                break;
                            case "a-NL_peak_counts":
                                writer.addValue("a-NL_peak_counts", pepObj.peakCounts.get("a-NL"));
                                break;
                            case "x-NL_peak_counts":
                                writer.addValue("x-NL_peak_counts", pepObj.peakCounts.get("x-NL"));
                                break;
                            case "c-NL_peak_counts":
                                writer.addValue("c-NL_peak_counts", pepObj.peakCounts.get("c-NL"));
                                break;
                            case "z-NL_peak_counts":
                                writer.addValue("z-NL_peak_counts", pepObj.peakCounts.get("z-NL"));
                                break;
                            case "precursor_peak_counts":
                                writer.addValue("precursor_peak_counts", pepObj.peakCounts.get("precursor"));
                                break;
                            case "precursor-NL_peak_counts":
                                writer.addValue("precursor-NL_peak_counts", pepObj.peakCounts.get("precursor-NL"));
                                break;
                            case "internal_peak_counts":
                                writer.addValue("internal_peak_counts", pepObj.peakCounts.get("internal"));
                                break;
                            case "internal-NL_peak_counts":
                                writer.addValue("internal-NL_peak_counts", pepObj.peakCounts.get("internal-NL"));
                                break;
                            case "immonium_peak_counts":
                                writer.addValue("immonium_peak_counts", pepObj.peakCounts.get("immonium"));
                                break;
                            case "unknown_peak_counts":
                                writer.addValue("unknown_peak_counts", pepObj.peakCounts.get("unknown"));
                                break;
                            case "y_spectral_similarity":
                                writer.addValue("y_spectral_similarity", pepObj.individualSpectralSimilarities.get("y"));
                                break;
                            case "b_spectral_similarity":
                                writer.addValue("b_spectral_similarity", pepObj.individualSpectralSimilarities.get("b"));
                                break;
                            case "a_spectral_similarity":
                                writer.addValue("a_spectral_similarity", pepObj.individualSpectralSimilarities.get("a"));
                                break;
                            case "x_spectral_similarity":
                                writer.addValue("x_spectral_similarity", pepObj.individualSpectralSimilarities.get("x"));
                                break;
                            case "c_spectral_similarity":
                                writer.addValue("c_spectral_similarity", pepObj.individualSpectralSimilarities.get("c"));
                                break;
                            case "z_spectral_similarity":
                                writer.addValue("z_spectral_similarity", pepObj.individualSpectralSimilarities.get("z"));
                                break;
                            case "cdot_spectral_similarity":
                                writer.addValue("cdot_spectral_similarity", pepObj.individualSpectralSimilarities.get("cdot"));
                                break;
                            case "zdot_spectral_similarity":
                                writer.addValue("zdot_spectral_similarity", pepObj.individualSpectralSimilarities.get("zdot"));
                                break;
                            case "y-NL_spectral_similarity":
                                writer.addValue("y-NL_spectral_similarity", pepObj.individualSpectralSimilarities.get("y-NL"));
                                break;
                            case "b-NL_spectral_similarity":
                                writer.addValue("b-NL_spectral_similarity", pepObj.individualSpectralSimilarities.get("b-NL"));
                                break;
                            case "a-NL_spectral_similarity":
                                writer.addValue("a-NL_spectral_similarity", pepObj.individualSpectralSimilarities.get("a-NL"));
                                break;
                            case "x-NL_spectral_similarity":
                                writer.addValue("x-NL_spectral_similarity", pepObj.individualSpectralSimilarities.get("x-NL"));
                                break;
                            case "c-NL_spectral_similarity":
                                writer.addValue("c-NL_spectral_similarity", pepObj.individualSpectralSimilarities.get("c-NL"));
                                break;
                            case "z-NL_spectral_similarity":
                                writer.addValue("z-NL_spectral_similarity", pepObj.individualSpectralSimilarities.get("z-NL"));
                                break;
                            case "precursor_spectral_similarity":
                                writer.addValue("precursor_spectral_similarity", pepObj.individualSpectralSimilarities.get("precursor"));
                                break;
                            case "precursor-NL_spectral_similarity":
                                writer.addValue("precursor-NL_spectral_similarity", pepObj.individualSpectralSimilarities.get("precursor-NL"));
                                break;
                            case "internal_spectral_similarity":
                                writer.addValue("internal_spectral_similarity", pepObj.individualSpectralSimilarities.get("internal"));
                                break;
                            case "internal-NL_spectral_similarity":
                                writer.addValue("internal-NL_spectral_similarity", pepObj.individualSpectralSimilarities.get("internal-NL"));
                                break;
                            case "immonium_spectral_similarity":
                                writer.addValue("immonium_spectral_similarity", pepObj.individualSpectralSimilarities.get("immonium"));
                                break;
                            case "unknown_spectral_similarity":
                                writer.addValue("unknown_spectral_similarity", pepObj.individualSpectralSimilarities.get("unknown"));
                                break;
                            case "intensity_distribution_similarity":
                                //unknown should be at the end
                                float[] predIntensities = new float[Constants.fragmentIonHierarchy.length - 1];
                                float[] expIntensities = new float[Constants.fragmentIonHierarchy.length - 1];
                                for (int j = 0; j < Constants.fragmentIonHierarchy.length - 1; j++) {
                                    predIntensities[j] = pepObj.predIntensities.get(Constants.fragmentIonHierarchy[j]);
                                    expIntensities[j] = pepObj.matchedIntensities.get(Constants.fragmentIonHierarchy[j]);
                                }
                                double value = new PearsonsCorrelation().correlation(Features.floatUtils.floatToDouble(predIntensities),
                                        Features.floatUtils.floatToDouble(expIntensities));
                                if (Double.isNaN(value)) {
                                    value = -1;
                                }
                                writer.addValue("intensity_distribution_similarity", value);
                                break;
                        }
                    }
                    //flush values to output
                    writer.writeValuesToRow();
                }
                pin.close();
                //endTime = System.nanoTime();
                //duration = (endTime - startTime);
                //System.out.println("Pin editing took " + duration / 1000000 +" milliseconds");
                writer.close();
                mzml.clear();
                if (Constants.renamePin == 1) {
                    System.out.println("Edited pin file at " + newOutfile);
                } else { //really should be 0
                    //move file at newOutfile to pinFiles[i] canonical name
                    File movedFile = new File(newOutfile);
                    pinFiles[i].delete();
                    movedFile.renameTo(pinFiles[i]);
                }
            }
        } catch (Exception e) {
            executorService.shutdown();
            e.printStackTrace();
            System.exit(1);
        }
        executorService.shutdown();
    }
}
