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

package features;

import allconstants.Constants;
import allconstants.FragmentIonConstants;
import features.spectra.SpectrumComparison;
import mainsteps.PeptideObj;
import mainsteps.PercolatorFormatter;
import org.apache.commons.math3.stat.correlation.PearsonsCorrelation;
import peptideptmformatting.PeptideFormatter;
import readers.datareaders.MzmlReader;
import readers.datareaders.PinReader;
import umich.ms.fileio.exceptions.FileParsingException;
import utils.Multithreader;
import utils.NumericUtils;
import utils.ProgressReporter;
import utils.StatMethods;

import java.io.IOException;
import java.net.URISyntaxException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import java.util.concurrent.ExecutionException;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.Future;

import static allconstants.FragmentIonConstants.fragmentGroups;

public class FeatureCalculator {
    PinReader pin;
    ArrayList<String> featuresList;
    MzmlReader mzml;
    public List<Double> ms2Scores = new ArrayList<>();
    public List<Double> rtScores = new ArrayList<>();

    public FeatureCalculator(PinReader pin, ArrayList<String> featuresList, MzmlReader mzml) {
        this.pin = pin;
        this.featuresList = featuresList;
        this.mzml = mzml;
        if (Constants.useSpectra) {
            List<Double> ms2ScoresUnsync = new ArrayList<>();
            ms2Scores = Collections.synchronizedList(ms2ScoresUnsync);
        }
        if (Constants.useRT) {
            List<Double> rtScoresUnsync = new ArrayList<>();
            rtScores = Collections.synchronizedList(rtScoresUnsync);
        }
    }

    class calcFeat implements Runnable {
        private final String line;
        private final int specIdx;
        private final int pepIdx;
        private final int scanNumIdx;
        private final ProgressReporter pr;
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
            String pep = pf.getBaseCharge();
            int scanNum = Integer.parseInt(row[scanNumIdx]);
            String stripped = pf.getStripped();

            PeptideObj pepObj;
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
                        case "deltaRTlinear":
                        case "deltaRTbins":
                        case "deltaRTLOESS":
                        case "deltaRTLOESSnormalized":
                        case "RTzscore":
                        case "RTprobability":
                        case "calibratedRT":
                        case "predRTrealUnits":
                        case "predictedRT":
                        case "deltaIMLOESS":
                        case "deltaIMLOESSnormalized":
                        case "predictedIM":
                        case "ionmobility":
                            break;

                        case "deltaRTLOESSreal":
                            if (Float.parseFloat(pepObj.escore) < Constants.loessEscoreCutoff) {
                                rtScores.add(pepObj.deltaRTLOESS_real);
                            }
                            break;

                        case "RTprobabilityUnifPrior":
                            pepObj.RTprobabilityUnifPrior = StatMethods.probabilityWithUniformPrior(
                                    mzml.unifPriorSize, mzml.unifProb,
                                    pepObj.scanNumObj.RTbinSize, (float) pepObj.RTprob);
                            break;
                        case "IMprobabilityUnifPrior":
                            pepObj.IMprobabilityUnifPrior = StatMethods.probabilityWithUniformPrior(
                                    mzml.unifPriorSizeIM[pepObj.charge - 1], mzml.unifProbIM[pepObj.charge - 1],
                                    pepObj.scanNumObj.IMbinSize, (float) pepObj.IMprob);
                            break;

                        case "peptideCounts":
                            pepObj.peptideCounts = Constants.peptideCounter.get(stripped).size();
                            break;

                        case "unweightedSpectralEntropy":
                            if (pepObj.spectralSimObj.spectrumComparisons.isEmpty()) {
                                try {
                                    double score = pepObj.spectralSimObj.getScore(feature);
                                    pepObj.spectralSimObj.scores.put(feature, score);
                                    if (Float.parseFloat(pepObj.escore) < Constants.loessEscoreCutoff) {
                                        ms2Scores.add(score);
                                    }
                                } catch (IOException | URISyntaxException e) {
                                    throw new RuntimeException(e);
                                }
                            } else {
                                for (int j = 0; j < fragmentGroups.length; j++) {
                                    try {
                                        pepObj.spectralSimObj.spectrumComparisons.get(j).scores.put(feature,
                                                pepObj.spectralSimObj.spectrumComparisons.get(j).getScore(feature));
                                    } catch (IOException | URISyntaxException e) {
                                        throw new RuntimeException(e);
                                    }
                                }
                            }
                            break;

                        case "brayCurtis":
                        case "cosineSimilarity":
                        case "spectralContrastAngle":
                        case "euclideanDistance":
                        case "pearsonCorr":
                        case "spearmanCorr":
                        case "hypergeometricProbability":
                        case "intersection":
                        case "dotProduct":
                        case "weightedSpectralEntropy":
                        case "heuristicSpectralEntropy":
                        case "top6matchedIntensity":
                            if (pepObj.spectralSimObj.spectrumComparisons.isEmpty()) {
                                try {
                                    pepObj.spectralSimObj.scores.put(feature, pepObj.spectralSimObj.getScore(feature));
                                } catch (IOException | URISyntaxException e) {
                                    throw new RuntimeException(e);
                                }
                            } else {
                                for (int j = 0; j < fragmentGroups.length; j++) {
                                    try {
                                        pepObj.spectralSimObj.spectrumComparisons.get(j).scores.put(feature,
                                                pepObj.spectralSimObj.spectrumComparisons.get(j).getScore(feature));
                                    } catch (IOException | URISyntaxException e) {
                                        throw new RuntimeException(e);
                                    }
                                }
                            }
                            break;

                        case "bootstrapSimilarity":
                            if (pepObj.spectralSimObj.spectrumComparisons.isEmpty()) {
                                mzml.futureList.clear();
                                Multithreader mt = new Multithreader(Constants.bootstraps, Constants.numThreads);
                                for (int j = 0; j < Constants.numThreads; j++) {
                                    int finalI = j;
                                    PeptideObj finalPepObj = pepObj;
                                    mzml.futureList.add(executorService.submit(() -> {
                                        SpectrumComparison sc;
                                        ArrayList<Double> scores = new ArrayList<>();
                                        for (int i = mt.indices[finalI]; i < mt.indices[finalI + 1]; i++) {
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
                            break;
                        case "bestScan":
//                        pepObj.spectralSimObj.scores.put(feature,
//                                (double) Math.abs(PercolatorFormatter.allPreds.get(pep).bestScan - pepObj.scanNum));
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
//                    case "p_matched_intensity":
//                        break;
//                    case "p-NL_matched_intensity":
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
//                    case "p_intensities_difference":
//                        break;
//                    case "p-NL_intensities_difference":
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
//                    case "p_peak_counts":
//                        break;
//                    case "p-NL_peak_counts":
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
//                    case "p_spectral_similarity":
//                        break;
//                    case "p-NL_spectral_similarity":
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
                            float[] predIntensities = new float[FragmentIonConstants.fragmentIonHierarchy.length - 1];
                            float[] expIntensities = new float[FragmentIonConstants.fragmentIonHierarchy.length - 1];
                            for (int j = 0; j < FragmentIonConstants.fragmentIonHierarchy.length - 1; j++) {
                                predIntensities[j] = pepObj.predIntensities.get(FragmentIonConstants.fragmentIonHierarchy[j]);
                                expIntensities[j] = pepObj.matchedIntensities.get(FragmentIonConstants.fragmentIonHierarchy[j]);
                            }
                            double value = new PearsonsCorrelation().correlation(NumericUtils.floatToDouble(predIntensities),
                                    NumericUtils.floatToDouble(expIntensities));
                            if (Double.isNaN(value)) {
                                value = -1;
                            }
                            pepObj.intensity_distribution_similarity = value;
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
