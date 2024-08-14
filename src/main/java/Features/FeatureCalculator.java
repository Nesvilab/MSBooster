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

import org.apache.commons.math3.stat.correlation.PearsonsCorrelation;
import umich.ms.fileio.exceptions.FileParsingException;

import java.io.IOException;
import java.util.*;
import java.util.concurrent.ExecutionException;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.Future;

public class FeatureCalculator {

    HashMap<String, HashMap<Integer, StatMethods>> featureStats = new HashMap<>();

    final Set<String> supportedFeatures = new HashSet<>(Arrays.asList("intersection",
            "hypergeometricProbability", "unweightedSpectralEntropy"));
    final Set<String> medianMethods = new HashSet<>(Arrays.asList("intersection", "unweightedSpectralEntropy"));
    final Set<String> zscoreMethods = new HashSet<>(List.of("hypergeometricProbability"));


    PinReader pin;
    ArrayList<String> featuresList;
    MzmlReader mzml;

    public FeatureCalculator(PinReader pin, ArrayList<String> featuresList, MzmlReader mzml) {
        this.pin = pin;
        this.featuresList = featuresList;
        this.mzml = mzml;
    }

    class calcFeat implements Runnable {
        private final String line;
        private final int specIdx;
        private final int pepIdx;
        private final int scanNumIdx;
        private String pep;
        private String stripped;
        private int scanNum;
        private ProgressReporter pr;
        final ExecutorService executorService = Executors.newFixedThreadPool(2);
        public calcFeat(String line, int specIdx, int pepIdx, int scanNumIdx, ProgressReporter pr) {
            this.line = line;
            this.specIdx = specIdx;
            this.pepIdx = pepIdx;
            this.scanNumIdx = scanNumIdx;
            this.pr = pr;
        }

        @Override
        public void run() {
            String[] row = line.split("\t");
            String[] periodSplit = row[specIdx].split("\\.");
            PeptideFormatter pf = new PeptideFormatter(row[pepIdx],
                    periodSplit[periodSplit.length - 1].split("_")[0], "pin");
            this.pep = pf.baseCharge;
            this.scanNum = Integer.parseInt(row[scanNumIdx]);
            this.stripped = pf.getStripped();

            PeptideObj pepObj = null;
            try {
                pepObj = mzml.getScanNumObject(scanNum).getPeptideObject(pep);
            } catch (FileParsingException e) {
                throw new RuntimeException(e);
            }

            //RT filter
            if (pepObj.deltaRTLOESS_real <= Constants.realMinuteFilter) {
                //switch case
                for (String feature : featuresList) {
                    switch (feature) {
//                    case "detectFractionGreater":
//                        float d = predictedSpectra.getPreds().get(pep).detectability;
//                        //for each protein, get the position of pep's detect and see how many peptides with greater detect are present
//                        //take max (proxy for protein that actually generated peptide)
//                        String[] r = pin.getRow();
//                        String[] prots = Arrays.copyOfRange(r, pin.pepIdx + 1, r.length);
//                        float maxFraction = 0f;
//                        for (String prot : prots) { //if more than one, this peptide is shared among proteins
//                            String protAbr;
//                            //skip protein if it is decoy and looking at target peptide
//                            if (prot.startsWith(Constants.decoyPrefix.substring(1))) {
//                                if (r[pin.labelIdx].equals("1")) {
//                                    continue;
//                                } else { //decoy peptide compared to target protein
//                                    protAbr = prot.substring(Constants.decoyPrefix.length() - 1);
//                                }
//                            } else {
//                                protAbr = prot;
//                            }
//
//                            float[] arr;
//                            try {
//                                arr = fasta.protToPep.get(protAbr).detects;
//                            } catch (Exception e) { //no peptides qualify from this protein
//                                continue;
//                            }
//
//                            int idx = Arrays.binarySearch(arr, d);
//                            if (idx < 0) { //not found
//                                idx = (-1 * idx) - 1;
//                            } else {
//                                idx += 1; //don't want to include itself in calculation
//                            }
//                            float[] presenceArr = Arrays.copyOfRange(fasta.protToPep.get(protAbr).presence, idx, arr.length);
//                            float total = 0f;
//                            for (float j : presenceArr) {
//                                if (j != 0f) {
//                                    total = j;
//                                    break;
//                                }
//                            }
//                            float fraction = (total + Constants.detectFractionGreaterNumerator) /
//                                    (presenceArr.length + Constants.detectFractionGreaterDenominator); //customizable prior
//                            if (fraction > maxFraction) {
//                                maxFraction = fraction;
//                            }
//                        }
//                        writer.addValue("detect_fraction_greater", maxFraction);
//                        break;
//                    case "detectSubtractMissing":
//                        d = predictedSpectra.getPreds().get(pep).detectability;
//                        //for each protein, get the position of pep's detect and see how many peptides with greater detect are present
//                        //take max (proxy for protein that actually generated peptide)
//                        r = pin.getRow();
//                        prots = Arrays.copyOfRange(r, pin.pepIdx + 1, r.length);
//                        float minDiff = 1f;
//                        for (String prot : prots) { //if more than one, this peptide is shared among proteins
//                            String protAbr;
//                            //skip protein if it is decoy and looking at target peptide
//                            if (prot.startsWith(Constants.decoyPrefix.substring(1))) {
//                                if (r[pin.labelIdx].equals("1")) {
//                                    continue;
//                                } else { //decoy
//                                    protAbr = prot.substring(Constants.decoyPrefix.length() - 1);
//                                }
//                            } else {
//                                protAbr = prot;
//                            }
//
//                            float[] arr;
//                            try {
//                                arr = fasta.protToPep.get(protAbr).detects;
//                            } catch (Exception e) { //no peptides qualify from this protein
//                                continue;
//                            }
//
//                            int idx = Arrays.binarySearch(arr, d);
//                            if (idx < 0) { //not found
//                                idx = (-1 * idx) - 1;
//                            } else {
//                                idx += 1; //don't want to include itself in calculation
//                            }
//                            float[] presenceArr = Arrays.copyOfRange(fasta.protToPep.get(protAbr).presence, idx, arr.length);
//                            float[] detectArr = Arrays.copyOfRange(arr, idx, arr.length);
//                            float total = 0f;
//                            for (int k = 0; k < presenceArr.length; k++) {
//                                if (presenceArr[k] == 0) {
//                                    total += detectArr[k] - d;
//                                }
//                            }
//                            float diff = total / presenceArr.length;
//                            if (diff < minDiff) {
//                                minDiff = diff;
//                            }
//                            if (minDiff == 0) {
//                                break;
//                            }
//                        }
//                        writer.addValue("detect_subtract_missing", minDiff);
//                        break;
//                    case "detectProtSpearmanDiff":
//                        SpearmansCorrelation sc = new SpearmansCorrelation();
//                        r = pin.getRow();
//                        prots = Arrays.copyOfRange(r, pin.pepIdx + 1, r.length);
//                        double maxSpearmanDiff = -3;
//                        float detect = predictedSpectra.getPreds().get(pep).detectability;;
//                        for (String prot : prots) { //if more than one, this peptide is shared among proteins
//                            //skip protein if it is decoy and looking at target peptide
//                            String protAbr;
//                            if (prot.startsWith(Constants.decoyPrefix.substring(1))) {
//                                if (r[pin.labelIdx].equals("1")) {
//                                    continue;
//                                } else {
//                                    protAbr = prot.substring(Constants.decoyPrefix.length() - 1);
//                                }
//                            } else {
//                                protAbr = prot;
//                            }
//
//                            float[] arr;
//                            try {
//                                arr = fasta.protToPep.get(protAbr).detects;
//                            } catch (Exception e) { //no peptides qualify from this protein
//                                continue;
//                            }
//
//                            float[] counts = fasta.protToPep.get(protAbr).spectralCounts;
//
//                            //only add if not 0 spectral counts (lots of missing ones, need to be more lenient for targets)
//                            ArrayList<Double> newDetects = new ArrayList<>();
//                            ArrayList<Double> newCounts = new ArrayList<>();
//                            for (int k = 0; k < arr.length; k++) {
//                                if (counts[k] != 0) {
//                                    if (! (arr[k] == detect)) { //will add later
//                                        newDetects.add((double) arr[k]);
//                                        newCounts.add((double) counts[k]);
//                                    }
//                                }
//                            }
//                            if (newDetects.size() < 2) {
//                                continue;
//                            }
//                            double spear = sc.correlation(newDetects.stream().mapToDouble(dd -> dd).toArray(),
//                                    newCounts.stream().mapToDouble(dd -> dd).toArray() );
//
//                            //add new pep to this calculation
//                            newDetects.add((double) detect);
//                            newCounts.add((double) predictedSpectra.getPreds().get(pep).counter);
//                            double spearDiff = sc.correlation(newDetects.stream().mapToDouble(dd -> dd).toArray(),
//                                    newCounts.stream().mapToDouble(dd -> dd).toArray() ) - spear;
//                            if (spearDiff > maxSpearmanDiff) {
//                                maxSpearmanDiff = spearDiff;
//                            }
//                        }
//                        if (maxSpearmanDiff == -3) {
//                            maxSpearmanDiff = 0;
//                        }
//                        writer.addValue("detect_prot_spearman_diff", maxSpearmanDiff);
//                        break;
                        case "deltaRTlinear":
                            break;
                        case "deltaRTbins":
                            break;
                        case "deltaRTLOESS":
                            break;
                        case "deltaRTLOESSreal":
                            break;
                        case "deltaRTLOESSnormalized":
                            break;
                        case "RTzscore":
                            break;
                        case "RTprobability":
                            break;
                        case "RTprobabilityUnifPrior":
                            pepObj.RTprobabilityUnifPrior = StatMethods.probabilityWithUniformPrior(mzml.unifPriorSize,
                                    mzml.unifProb, pepObj.scanNumObj.RTbinSize, (float) pepObj.RTprob);
                            break;
                        case "calibratedRT":
                            break;
                        case "predRTrealUnits":
                            break;
                        case "predictedRT":
                            break;
                        case "brayCurtis":
                            if (pepObj.spectralSimObj.spectrumComparisons.size() == 0) {
                                pepObj.spectralSimObj.scores.put(feature, pepObj.spectralSimObj.brayCurtis());
                            } else {
                                String[] dividedFragments = Constants.divideFragments.split(";");
                                for (int j = 0; j < dividedFragments.length; j++) {
                                    pepObj.spectralSimObj.spectrumComparisons.get(j).scores.put(
                                            feature, pepObj.spectralSimObj.spectrumComparisons.get(j).brayCurtis()
                                    );
                                }
                            }
                            break;
                        case "cosineSimilarity":
                            if (pepObj.spectralSimObj.spectrumComparisons.size() == 0) {
                                pepObj.spectralSimObj.scores.put(feature, pepObj.spectralSimObj.cosineSimilarity());
                            } else {
                                String[] dividedFragments = Constants.divideFragments.split(";");
                                for (int j = 0; j < dividedFragments.length; j++) {
                                    pepObj.spectralSimObj.spectrumComparisons.get(j).scores.put(
                                            feature, pepObj.spectralSimObj.spectrumComparisons.get(j).cosineSimilarity()
                                    );
                                }
                            }
                            break;
                        case "spectralContrastAngle":
                            if (pepObj.spectralSimObj.spectrumComparisons.size() == 0) {
                                pepObj.spectralSimObj.scores.put(feature, pepObj.spectralSimObj.spectralContrastAngle());
                            } else {
                                String[] dividedFragments = Constants.divideFragments.split(";");
                                for (int j = 0; j < dividedFragments.length; j++) {
                                    pepObj.spectralSimObj.spectrumComparisons.get(j).scores.put(
                                            feature, pepObj.spectralSimObj.spectrumComparisons.get(j).spectralContrastAngle()
                                    );
                                }
                            }
                            break;
                        case "euclideanDistance":
                            if (pepObj.spectralSimObj.spectrumComparisons.size() == 0) {
                                pepObj.spectralSimObj.scores.put(feature, pepObj.spectralSimObj.euclideanDistance());
                            } else {
                                String[] dividedFragments = Constants.divideFragments.split(";");
                                for (int j = 0; j < dividedFragments.length; j++) {
                                    pepObj.spectralSimObj.spectrumComparisons.get(j).scores.put(
                                            feature, pepObj.spectralSimObj.spectrumComparisons.get(j).euclideanDistance()
                                    );
                                }
                            }
                            break;
                        case "pearsonCorr":
                            if (pepObj.spectralSimObj.spectrumComparisons.size() == 0) {
                                pepObj.spectralSimObj.scores.put(feature, pepObj.spectralSimObj.pearsonCorr());
                            } else {
                                String[] dividedFragments = Constants.divideFragments.split(";");
                                for (int j = 0; j < dividedFragments.length; j++) {
                                    pepObj.spectralSimObj.spectrumComparisons.get(j).scores.put(
                                            feature, pepObj.spectralSimObj.spectrumComparisons.get(j).pearsonCorr()
                                    );
                                }
                            }
                            break;
                        case "spearmanCorr":
                            if (pepObj.spectralSimObj.spectrumComparisons.size() == 0) {
                                pepObj.spectralSimObj.scores.put(feature, pepObj.spectralSimObj.spearmanCorr());
                            } else {
                                String[] dividedFragments = Constants.divideFragments.split(";");
                                for (int j = 0; j < dividedFragments.length; j++) {
                                    pepObj.spectralSimObj.spectrumComparisons.get(j).scores.put(
                                            feature, pepObj.spectralSimObj.spectrumComparisons.get(j).spearmanCorr()
                                    );
                                }
                            }
                            break;
                        case "hypergeometricProbability":
                            if (pepObj.spectralSimObj.spectrumComparisons.size() == 0) {
                                pepObj.spectralSimObj.scores.put(feature, pepObj.spectralSimObj.hyperGeometricProbability());
                            } else {
                                String[] dividedFragments = Constants.divideFragments.split(";");
                                for (int j = 0; j < dividedFragments.length; j++) {
                                    pepObj.spectralSimObj.spectrumComparisons.get(j).scores.put(
                                            feature, pepObj.spectralSimObj.spectrumComparisons.get(j).hyperGeometricProbability()
                                    );
                                }
                            }
                            break;
                        case "intersection":
                            if (pepObj.spectralSimObj.spectrumComparisons.size() == 0) {
                                pepObj.spectralSimObj.scores.put(feature, pepObj.spectralSimObj.intersection());
                            } else {
                                String[] dividedFragments = Constants.divideFragments.split(";");
                                for (int j = 0; j < dividedFragments.length; j++) {
                                    pepObj.spectralSimObj.spectrumComparisons.get(j).scores.put(
                                            feature, pepObj.spectralSimObj.spectrumComparisons.get(j).intersection()
                                    );
                                }
                            }
                            break;
                        case "dotProduct":
                            if (pepObj.spectralSimObj.spectrumComparisons.size() == 0) {
                                pepObj.spectralSimObj.scores.put(feature, pepObj.spectralSimObj.dotProduct());
                            } else {
                                String[] dividedFragments = Constants.divideFragments.split(";");
                                for (int j = 0; j < dividedFragments.length; j++) {
                                    pepObj.spectralSimObj.spectrumComparisons.get(j).scores.put(
                                            feature, pepObj.spectralSimObj.spectrumComparisons.get(j).dotProduct()
                                    );
                                }
                            }
                            break;
                        case "unweightedSpectralEntropy":
                            if (pepObj.spectralSimObj.spectrumComparisons.size() == 0) {
                                pepObj.spectralSimObj.scores.put(feature, pepObj.spectralSimObj.unweightedSpectralEntropy());
                            } else {
                                String[] dividedFragments = Constants.divideFragments.split(";");
                                for (int j = 0; j < dividedFragments.length; j++) {
                                    pepObj.spectralSimObj.spectrumComparisons.get(j).scores.put(
                                            feature, pepObj.spectralSimObj.spectrumComparisons.get(j).unweightedSpectralEntropy()
                                    );
                                }
                            }
                            break;
                        case "bootstrapSimilarity":
                            if (pepObj.spectralSimObj.spectrumComparisons.size() == 0) {
                                mzml.futureList.clear();
                                for (int j = 0; j < Constants.numThreads; j++) {
                                    int start = (int) (Constants.bootstraps * (long) j) / Constants.numThreads;
                                    int end = (int) (Constants.bootstraps * (long) (j + 1)) / Constants.numThreads;
                                    PeptideObj finalPepObj = pepObj;
                                    mzml.futureList.add(executorService.submit(() -> {
                                        SpectrumComparison sc;
                                        ArrayList<Double> scores = new ArrayList<>();
                                        for (int i = start; i < end; i++) {
                                            sc = finalPepObj.spectralSimObj.pickedPredicted();
                                            scores.add(sc.unweightedSpectralEntropy());
                                        }
                                        Collections.sort(scores);
                                        finalPepObj.spectralSimObj.scores.put(feature,
                                                (scores.get(scores.size() / 2) + scores.get(scores.size() / 2 - 1)) / 2);
                                    }));
                                }
                            }
                            break;
                        case "adjacentSimilarity":
//                        double score = 0;
//                        int divisor = 0;
//                        ArrayList<Double> scores = new ArrayList<>();
//
//                        int previous = pepObj.previousScan;
//                        if (previous != 0) {
//                            MzmlScanNumber msn = mzml.scanNumberObjects.get(previous);
//                            PredictionEntry pe = PercolatorFormatter.allPreds.get(pep);
//                            SpectrumComparison sc = new SpectrumComparison(msn.getExpMZs(), msn.getExpIntensities(),
//                                    pe.mzs, pe.intensities, pepObj.length,
//                                    Constants.useTopFragments, Constants.useBasePeak, false);
//                            score += sc.unweightedSpectralEntropy();
//                            scores.add(sc.unweightedSpectralEntropy());
//                            divisor += 1;
//                        }
//
//                        int next = pepObj.nextScan;
//                        if (next != 0) {
//                            MzmlScanNumber msn = mzml.scanNumberObjects.get(next);
//                            PredictionEntry pe = PercolatorFormatter.allPreds.get(pep);
//                            SpectrumComparison sc = new SpectrumComparison(msn.getExpMZs(), msn.getExpIntensities(),
//                                    pe.mzs, pe.intensities, pepObj.length,
//                                    Constants.useTopFragments, Constants.useBasePeak, false);
//                            score += sc.unweightedSpectralEntropy();
//                            scores.add(sc.unweightedSpectralEntropy());
//                            divisor += 1;
//                        }
//
//                        if (divisor == 0) { //TODO: why are they missing so deep into RT?
//                            pepObj.spectralSimObj.scores.put(feature, 0d);
//                        } else {
//                            //pepObj.spectralSimObj.scores.put(feature, score / divisor);
//                            pepObj.spectralSimObj.scores.put(feature, Collections.max(scores));
//                        }

//                        for (PeptideObj pobj : msn.peptideObjects) {
//                            allPreds.get(pobj.name).times.add(allMatchedScans.get(pobj.precursorMz).indexOf(msn.scanNum));
//                        }
                            break;
                        case "bestScan":
//                        pepObj.spectralSimObj.scores.put(feature,
//                                (double) Math.abs(PercolatorFormatter.allPreds.get(pep).bestScan - pepObj.scanNum));
                            break;
                        case "deltaIMLOESS":
                            break;
                        case "deltaIMLOESSnormalized":
                            break;
                        case "IMprobabilityUnifPrior":
                            float prob = StatMethods.probabilityWithUniformPrior(mzml.unifPriorSizeIM[pepObj.charge - 1],
                                    mzml.unifProbIM[pepObj.charge - 1], pepObj.scanNumObj.IMbinSize, (float) pepObj.IMprob);
                            pepObj.IMprobabilityUnifPrior = prob;
                            break;
//                    case "predictedIM":
//                        break;
                    case "ionmobility":
                        break;
//                    case "y_matched_intensity":
//                        break;
//                    case "b_matched_intensity":
//                        break;
//                    case "a_matched_intensity":
//                        break;
//                    case "x_matched_intensity":
//                        break;
//                    case "c_matched_intensity":
//                        break;
//                    case "z_matched_intensity":
//                        break;
//                    case "cdot_matched_intensity":
//                        break;
//                    case "zdot_matched_intensity":
//                        break;
//                    case "y-NL_matched_intensity":
//                        break;
//                    case "b-NL_matched_intensity":
//                        break;
//                    case "a-NL_matched_intensity":
//                        break;
//                    case "x-NL_matched_intensity":
//                        break;
//                    case "c-NL_matched_intensity":
//                        break;
//                    case "z-NL_matched_intensity":
//                        break;
//                    case "precursor_matched_intensity":
//                        break;
//                    case "precursor-NL_matched_intensity":
//                        break;
//                    case "internal_matched_intensity":
//                        break;
//                    case "internal-NL_matched_intensity":
//                        break;
//                    case "immonium_matched_intensity":
//                        break;
//                    case "unknown_matched_intensity":
//                        break;
//                    case "y_intensities_difference":
//                        break;
//                    case "b_intensities_difference":
//                        break;
//                    case "a_intensities_difference":
//                        break;
//                    case "x_intensities_difference":
//                        break;
//                    case "c_intensities_difference":
//                        break;
//                    case "z_intensities_difference":
//                        break;
//                    case "cdot_intensities_difference":
//                        break;
//                    case "zdot_intensities_difference":
//                        break;
//                    case "y-NL_intensities_difference":
//                        break;
//                    case "b-NL_intensities_difference":
//                        break;
//                    case "a-NL_intensities_difference":
//                        break;
//                    case "x-NL_intensities_difference":
//                        break;
//                    case "c-NL_intensities_difference":
//                        break;
//                    case "z-NL_intensities_difference":
//                        break;
//                    case "precursor_intensities_difference":
//                        break;
//                    case "precursor-NL_intensities_difference":
//                        break;
//                    case "internal_intensities_difference":
//                        break;
//                    case "internal-NL_intensities_difference":
//                        break;
//                    case "immonium_intensities_difference":
//                        break;
//                    case "unknown_intensities_difference":
//                        break;
//                    case "y_peak_counts":
//                        break;
//                    case "b_peak_counts":
//                        break;
//                    case "a_peak_counts":
//                        break;
//                    case "x_peak_counts":
//                        break;
//                    case "c_peak_counts":
//                        break;
//                    case "z_peak_counts":
//                        break;
//                    case "cdot_peak_counts":
//                        break;
//                    case "zdot_peak_counts":
//                        break;
//                    case "y-NL_peak_counts":
//                        break;
//                    case "b-NL_peak_counts":
//                        break;
//                    case "a-NL_peak_counts":
//                        break;
//                    case "x-NL_peak_counts":
//                        break;
//                    case "c-NL_peak_counts":
//                        break;
//                    case "z-NL_peak_counts":
//                        break;
//                    case "precursor_peak_counts":
//                        break;
//                    case "precursor-NL_peak_counts":
//                        break;
//                    case "internal_peak_counts":
//                        break;
//                    case "internal-NL_peak_counts":
//                        break;
//                    case "immonium_peak_counts":
//                        break;
//                    case "unknown_peak_counts":
//                        break;
//                    case "y_spectral_similarity":
//                        break;
//                    case "b_spectral_similarity":
//                        break;
//                    case "a_spectral_similarity":
//                        break;
//                    case "x_spectral_similarity":
//                        break;
//                    case "c_spectral_similarity":
//                        break;
//                    case "z_spectral_similarity":
//                        break;
//                    case "cdot_spectral_similarity":
//                        break;
//                    case "zdot_spectral_similarity":
//                        break;
//                    case "y-NL_spectral_similarity":
//                        break;
//                    case "b-NL_spectral_similarity":
//                        break;
//                    case "a-NL_spectral_similarity":
//                        break;
//                    case "x-NL_spectral_similarity":
//                        break;
//                    case "c-NL_spectral_similarity":
//                        break;
//                    case "z-NL_spectral_similarity":
//                        break;
//                    case "precursor_spectral_similarity":
//                        break;
//                    case "precursor-NL_spectral_similarity":
//                        break;
//                    case "internal_spectral_similarity":
//                        break;
//                    case "internal-NL_spectral_similarity":
//                        break;
//                    case "immonium_spectral_similarity":
//                        break;
//                    case "unknown_spectral_similarity":
//                        break;
                        case "intensity_distribution_similarity":
                            //unknown should be at the end
                            float[] predIntensities = new float[Constants.fragmentIonHierarchy.length - 1];
                            float[] expIntensities = new float[Constants.fragmentIonHierarchy.length - 1];
                            for (int j = 0; j < Constants.fragmentIonHierarchy.length - 1; j++) {
                                predIntensities[j] = pepObj.predIntensities.get(Constants.fragmentIonHierarchy[j]);
                                expIntensities[j] = pepObj.matchedIntensities.get(Constants.fragmentIonHierarchy[j]);
                            }
                            double value = new PearsonsCorrelation().correlation(FloatUtils.floatToDouble(predIntensities),
                                    FloatUtils.floatToDouble(expIntensities));
                            if (Double.isNaN(value)) {
                                value = -1;
                            }
                            pepObj.intensity_distribution_similarity = value;
                            break;
                        case "peptideCounts":
                            pepObj.peptideCounts = Constants.peptideCounter.get(stripped).size();
                            break;
                    }
                }
            }
            pr.progress();
        }
    }

    public void calculate(ExecutorService executorService) throws IOException, ExecutionException, InterruptedException {
        ProgressReporter pr = new ProgressReporter(pin.getLength());
        mzml.futureList.clear();
        while (pin.next(false)) {
            calcFeat task = new calcFeat(pin.line, pin.specIdx, pin.pepIdx, pin.scanNumIdx, pr);
            mzml.futureList.add(executorService.submit(task));
        }
        for (Future future : mzml.futureList) {
            future.get();
        }
        pin.reset();
    }
}
