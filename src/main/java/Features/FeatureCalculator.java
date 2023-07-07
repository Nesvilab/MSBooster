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

import java.io.IOException;
import java.util.*;

public class FeatureCalculator {

    HashMap<String, HashMap<Integer, StatMethods>> featureStats = new HashMap<>();

    final Set<String> supportedFeatures = new HashSet<>(Arrays.asList("intersection",
            "hypergeometricProbability", "unweightedSpectralEntropy"));
    final Set<String> medianMethods = new HashSet<>(Arrays.asList("intersection", "unweightedSpectralEntropy"));
    final Set<String> zscoreMethods = new HashSet<>(Arrays.asList("hypergeometricProbability"));


    PinReader pin;
    ArrayList<String> featuresList;
    MzmlReader mzml;

    public FeatureCalculator(PinReader pin, ArrayList<String> featuresList, MzmlReader mzml) {
        this.pin = pin;
        this.featuresList = featuresList;
        this.mzml = mzml;
    }

    public void calculate() throws IOException {
        PeptideObj pepObj = null;
        HashMap<Integer, StatMethods> hm;
        StatMethods sm;

        ProgressReporter pr = new ProgressReporter(pin.length);
        while (pin.next()) {
            pr.progress();

            String pep = pin.getPep().baseCharge;

            pepObj = mzml.scanNumberObjects.get(pin.getScanNum()).getPeptideObject(pep);

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
                        if (Constants.noRTscores) {
                            pepObj.deltaRT = 0;
                        }
                        break;
                    case "deltaRTbins":
                        if (Constants.noRTscores) {
                            pepObj.deltaRTbin = 0;
                        }
                        break;
                    case "deltaRTLOESS":
                        if (Constants.noRTscores) {
                            pepObj.deltaRTLOESS = 0;
                        }
                        break;
                    case "deltaRTLOESSnormalized":
                        if (Constants.noRTscores) {
                            pepObj.deltaRTLOESSnormalized = 0;
                        }
                        break;
                    case "RTzscore":
                        if (Constants.noRTscores) {
                            pepObj.RTzscore = 0;
                        }
                        break;
                    case "RTprobability":
                        if (Constants.noRTscores) {
                            pepObj.RTprob = 0;
                        }
                        break;
                    case "RTprobabilityUnifPrior":
                        if (Constants.noRTscores) {
                            pepObj.RTprobabilityUnifPrior = 0;
                        } else {
                            pepObj.RTprobabilityUnifPrior = StatMethods.probabilityWithUniformPrior(mzml.unifPriorSize,
                                    mzml.unifProb, pepObj.scanNumObj.RTbinSize, (float) pepObj.RTprob);
                        }
                        break;
                    case "calibratedRT":
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
//                    case "ionmobility":
//                        break;
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
                }

                if (Constants.normalizeScoresByPeptideLength && supportedFeatures.contains(feature)) {
                    int peptideLength = pepObj.length;
                    if (featureStats.containsKey(feature)) {
                        hm = featureStats.get(feature);
                        if (hm.containsKey(peptideLength)) {
                            sm = hm.get(peptideLength);
                        } else {
                            sm = new StatMethods();
                        }
                    } else {
                        hm = new HashMap<>();
                        sm = new StatMethods();
                    }

                    double score = pepObj.spectralSimObj.scores.get(feature);
                    if (medianMethods.contains(feature)) {
                        sm.updateMedian(score);
                    } else if (zscoreMethods.contains(feature)) {
                        sm.updateVariance(score);
                    }

                    hm.put(peptideLength, sm);
                    featureStats.put(feature, hm);
                }
            }
            pepObj.clearArrays();
        }
        System.out.println("");
        pin.reset();
    }
}
