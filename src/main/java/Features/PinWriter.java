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

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;

import static Features.Constants.camelToUnderscore;

public class PinWriter {
    String newOutfile;
    PinReader pin;
    ArrayList<String> featuresList;
    MzmlReader mzml;
    TsvWriter writer;
    ArrayList<String> header = new ArrayList<>();
    HashMap<String, HashMap<Integer, StatMethods>> featureStats;

    public PinWriter(String newOutfile, PinReader pin, ArrayList<String> featuresList, MzmlReader mzml,
                     HashMap<String, HashMap<Integer, StatMethods>> featureStats) {
        this.newOutfile = newOutfile;
        this.pin = pin;
        this.featuresList = featuresList;
        this.mzml = mzml;
        this.featureStats = featureStats;

        TsvWriterSettings tws = new TsvWriterSettings();
        tws.setMaxCharsPerColumn(-1);
        tws.setMaxColumns(50000); //who knows if it needs to be longer?
        writer = new TsvWriter(new File(newOutfile), tws);

        //add header to written tsv
        header.addAll(Arrays.asList(pin.header));
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
        header.addAll(pin.pepIdx, newNames); //add features before Peptide
        header.remove("detectability");
        writer.writeHeaders(header);
    }

    public void write() throws IOException {
        PeptideObj pepObj = null;
        ProgressReporter pr = new ProgressReporter(pin.length);
        while (pin.next()) {
            pr.progress();

            //write everything we already have, not including extra protein columns
            String[] row = pin.getRow();
            for (int j = 0; j < pin.header.length; j++) {
                writer.addValue(pin.header[j], row[j]);
            }
            //add extra protein columns
            if (pin.getRow().length > pin.header.length) {
                for (int j = pin.header.length; j < pin.getRow().length; j++) {
                    writer.addValue(header.size(), row[j]);
                }
            }

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
                        writer.addValue("deltaRTlinear", pepObj.deltaRT);
                        break;
                    case "deltaRTbins":
                        writer.addValue("deltaRTbins", pepObj.deltaRTbin);
                        break;
                    case "deltaRTLOESS":
                        writer.addValue("delta_RT_loess", pepObj.deltaRTLOESS);
                        break;
                    case "deltaRTLOESSnormalized":
                        writer.addValue("delta_RT_loess_normalized", pepObj.deltaRTLOESSnormalized);
                        break;
                    case "RTzscore":
                        writer.addValue("RTzscore", pepObj.RTzscore);
                        break;
                    case "RTprobability":
                        writer.addValue("RTprobability", pepObj.RTprob);
                        break;
                    case "RTprobabilityUnifPrior":
                        writer.addValue("RT_probability_unif_prior", pepObj.RTprobabilityUnifPrior);
                        break;
                    case "calibratedRT":
                        writer.addValue("calibrated_RT", pepObj.calibratedRT);
                        break;
                    case "predictedRT":
                        writer.addValue("predicted_RT", pepObj.RT);
                        break;
                    case "brayCurtis":
                        if (pepObj.spectralSimObj.spectrumComparisons.size() == 0) {
                            double score = pepObj.spectralSimObj.scores.get(feature);
                            writer.addValue("bray_curtis", score);
                        } else {
                            String[] dividedFragments = Constants.divideFragments.split(";");
                            for (int j = 0; j < dividedFragments.length; j++) {
                                double score = pepObj.spectralSimObj.spectrumComparisons.get(j).scores.get(feature);
                                writer.addValue("bray_curtis_" + dividedFragments[j], score);
                            }
                        }
                        break;
                    case "cosineSimilarity":
                        if (pepObj.spectralSimObj.spectrumComparisons.size() == 0) {
                            double score = pepObj.spectralSimObj.scores.get(feature);
                            writer.addValue("cosine_similarity", score);
                        } else {
                            String[] dividedFragments = Constants.divideFragments.split(";");
                            for (int j = 0; j < dividedFragments.length; j++) {
                                double score = pepObj.spectralSimObj.spectrumComparisons.get(j).scores.get(feature);
                                writer.addValue("cosine_similarity_" + dividedFragments[j], score);
                            }
                        }
                        break;
                    case "spectralContrastAngle":
                        if (pepObj.spectralSimObj.spectrumComparisons.size() == 0) {
                            double score = pepObj.spectralSimObj.scores.get(feature);
                            writer.addValue("spectra_contrast_angle", score);
                        } else {
                            String[] dividedFragments = Constants.divideFragments.split(";");
                            for (int j = 0; j < dividedFragments.length; j++) {
                                double score = pepObj.spectralSimObj.spectrumComparisons.get(j).scores.get(feature);
                                writer.addValue("spectra_contrast_angle_" + dividedFragments[j], score);
                            }
                        }
                        break;
                    case "euclideanDistance":
                        if (pepObj.spectralSimObj.spectrumComparisons.size() == 0) {
                            double score = pepObj.spectralSimObj.scores.get(feature);
                            writer.addValue("euclidean_distance", score);
                        } else {
                            String[] dividedFragments = Constants.divideFragments.split(";");
                            for (int j = 0; j < dividedFragments.length; j++) {
                                double score = pepObj.spectralSimObj.spectrumComparisons.get(j).scores.get(feature);
                                writer.addValue("euclidean_distance_" + dividedFragments[j], score);
                            }
                        }
                        break;
                    case "pearsonCorr":
                        if (pepObj.spectralSimObj.spectrumComparisons.size() == 0) {
                            double score = pepObj.spectralSimObj.scores.get(feature);
                            writer.addValue("pearson_corr", score);
                        } else {
                            String[] dividedFragments = Constants.divideFragments.split(";");
                            for (int j = 0; j < dividedFragments.length; j++) {
                                double score = pepObj.spectralSimObj.spectrumComparisons.get(j).scores.get(feature);
                                writer.addValue("pearson_corr_" + dividedFragments[j], score);
                            }
                        }
                        break;
                    case "spearmanCorr":
                        if (pepObj.spectralSimObj.spectrumComparisons.size() == 0) {
                            double score = pepObj.spectralSimObj.scores.get(feature);
                            if (Constants.normalizeScoresByPeptideLength) {
                                score *= Math.log(pepObj.length);
                            }
                            writer.addValue("spearman_corr", score);
                        } else {
                            String[] dividedFragments = Constants.divideFragments.split(";");
                            for (int j = 0; j < dividedFragments.length; j++) {
                                double score = pepObj.spectralSimObj.spectrumComparisons.get(j).scores.get(feature);
                                writer.addValue("spearman_corr_" + dividedFragments[j], score);
                            }
                        }
                        break;
                    case "hypergeometricProbability":
                        if (pepObj.spectralSimObj.spectrumComparisons.size() == 0) {
                            double score = pepObj.spectralSimObj.scores.get(feature);
                            if (Constants.normalizeScoresByPeptideLength) {
                                score = (score - featureStats.get(feature).get(pepObj.length).getMean())
                                        / featureStats.get(feature).get(pepObj.length).getStd();
                            }
                            writer.addValue("hypergeometric_probability", score);
                        } else {
                            String[] dividedFragments = Constants.divideFragments.split(";");
                            for (int j = 0; j < dividedFragments.length; j++) {
                                double score = pepObj.spectralSimObj.spectrumComparisons.get(j).scores.get(feature);
                                writer.addValue("hypergeometric_probability_" + dividedFragments[j], score);
                            }
                        }
                        break;
                    case "intersection":
                        if (pepObj.spectralSimObj.spectrumComparisons.size() == 0) {
                            double score = pepObj.spectralSimObj.scores.get(feature);
                            if (Constants.normalizeScoresByPeptideLength) {
                                score -= featureStats.get(feature).get(pepObj.length).getMedian();
                            }
                            writer.addValue("intersection", score);
                        } else {
                            String[] dividedFragments = Constants.divideFragments.split(";");
                            for (int j = 0; j < dividedFragments.length; j++) {
                                double score = pepObj.spectralSimObj.spectrumComparisons.get(j).scores.get(feature);
                                writer.addValue("intersection_" + dividedFragments[j], score);
                            }
                        }
                        break;
                    case "dotProduct":
                        if (pepObj.spectralSimObj.spectrumComparisons.size() == 0) {
                            double score = pepObj.spectralSimObj.scores.get(feature);
                            writer.addValue("dot_product", score);
                        } else {
                            String[] dividedFragments = Constants.divideFragments.split(";");
                            for (int j = 0; j < dividedFragments.length; j++) {
                                double score = pepObj.spectralSimObj.spectrumComparisons.get(j).scores.get(feature);
                                writer.addValue("dot_product_" + dividedFragments[j], score);
                            }
                        }
                        break;
                    case "unweightedSpectralEntropy":
                        if (pepObj.spectralSimObj.spectrumComparisons.size() == 0) {
                            double score = pepObj.spectralSimObj.scores.get(feature);
                            if (Constants.normalizeScoresByPeptideLength) {
                                score -= featureStats.get(feature).get(pepObj.length).getMedian();
                            }
                            writer.addValue("unweighted_spectral_entropy", score);
                        } else {
                            String[] dividedFragments = Constants.divideFragments.split(";");
                            for (int j = 0; j < dividedFragments.length; j++) {
                                double score = pepObj.spectralSimObj.spectrumComparisons.get(j).scores.get(feature);
                                writer.addValue("unweighted_spectral_entropy_" + dividedFragments[j], score);
                            }
                        }
                        break;
                    case "deltaIMLOESS":
                        writer.addValue("delta_IM_loess", pepObj.deltaIMLOESS);
                        break;
                    case "deltaIMLOESSnormalized":
                        writer.addValue("delta_IM_loess_normalized", pepObj.deltaIMLOESSnormalized);
                        break;
                    case "IMprobabilityUnifPrior":
                        writer.addValue("IM_probability_unif_prior", pepObj.IMprobabilityUnifPrior);
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
                        writer.addValue("intensity_distribution_similarity", pepObj.intensity_distribution_similarity);
                        break;
                }
            }
            //flush values to output
            writer.writeValuesToRow();

            //clear old pep obj
            pepObj.spectralSimObj = null;
        }
        System.out.println("");
        pin.close();
        writer.close();
        mzml.clear();
    }
}
