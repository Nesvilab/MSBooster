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

package features.spectra;

import allconstants.Constants;
import allconstants.FragmentIonConstants;
import mainsteps.PeptideObj;
import org.apache.commons.math3.distribution.HypergeometricDistribution;
import org.apache.commons.math3.random.Well19937c;
import org.apache.commons.math3.stat.correlation.PearsonsCorrelation;
import org.apache.commons.math3.stat.correlation.SpearmansCorrelation;

import java.io.IOException;
import java.net.URISyntaxException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;
import java.util.LinkedHashSet;
import java.util.Objects;
import java.util.Random;
import java.util.TreeSet;

import static allconstants.FragmentIonConstants.fragmentGroups;
import static utils.NumericUtils.floatToDouble;
import static utils.Print.printError;
import static utils.Print.printInfo;

//TODO: also square root intensities? Squaring intensities may help for single cell data
public class SpectrumComparison {
    float[] predMZs;
    float[] predIntensities;
    String[] predFragmentIonTypes;
    float[] matchedIntensities;
    float[] unitNormMatchedIntensities;
    float[] unitNormPredIntensities;
    float[] sum1MatchedIntensities;
    float[] sum1PredIntensities;
    float[] allMatchedIntensities;
    public HashMap<String, Double> scores = new HashMap<>();
    public LinkedHashSet<Integer> matchedIdx = new LinkedHashSet<Integer>();
    private static final PearsonsCorrelation pc = new PearsonsCorrelation();
    private static final SpearmansCorrelation sc = new SpearmansCorrelation();
    public ArrayList<SpectrumComparison> spectrumComparisons = new ArrayList<>();
    public static Well19937c rng = new Well19937c(123);
    public PeptideObj pepObj;

    public SpectrumComparison(PeptideObj pepObj, float[] eMZs, float[] eIntensities,
                              float[] pMZs, float[] pIntensities, String[] pFragmentIonTypes, boolean separateByFragmentGroup) {
        predMZs = pMZs;
        predIntensities = pIntensities;
        predFragmentIonTypes = pFragmentIonTypes;

        if (fragmentGroups.length == 1 || !separateByFragmentGroup) {
            matchedIntensities = this.getMatchedIntensities(eMZs, eIntensities, predMZs, pFragmentIonTypes, pepObj);
        } else {
            //get fragments that match the allowed
            for (TreeSet<String> allowedTypesSet : fragmentGroups) {
                ArrayList<Integer> acceptedIdx = new ArrayList<>();

                for (int i = 0; i < pFragmentIonTypes.length; i++) {
                    String ion = pFragmentIonTypes[i];
                    if (allowedTypesSet.contains(ion)) {
                        acceptedIdx.add(i);
                    }
                }
                float[] mzs = new float[acceptedIdx.size()];
                float[] ints = new float[acceptedIdx.size()];
                String[] fits = new String[acceptedIdx.size()];
                for (int i = 0; i < acceptedIdx.size(); i++) {
                    mzs[i] = predMZs[acceptedIdx.get(i)];
                    ints[i] = predIntensities[acceptedIdx.get(i)];
                    fits[i] = predFragmentIonTypes[acceptedIdx.get(i)];
                }

                //create new spectrumComparison obj and add to list
                spectrumComparisons.add(new SpectrumComparison(pepObj, eMZs, eIntensities,
                        mzs, ints, fits, false));
            }
        }
        predMZs = null;

        this.pepObj = pepObj;
        if (Constants.features.contains("adjacent") || Constants.features.contains("bestScan")) {
            MassCalculator mc = new MassCalculator(pepObj.name.split("\\|")[0], pepObj.charge);
            pepObj.precursorMz = (mc.mass + pepObj.charge * mc.proton) / pepObj.charge;
        }
    }

    //TODO: if ever reimplement adjacent similarity, think of solution that does not require another constructor
//    public SpectrumComparison(PeptideObj peptideObj, float[] eMZs, float[] eIntensities,
//                              float[] pMZs, float[] pIntensities, String[] pFragmentIonTypes,
//                              boolean willReload) {
//        pepObj = peptideObj;
//        predMZs = pMZs;
//        predIntensities = pIntensities;
//        predFragmentIonTypes = pFragmentIonTypes;
//
//        matchedIntensities = this.getMatchedIntensities(eMZs, eIntensities, predMZs, pFragmentIonTypes);
//        if (! willReload) {
//            predMZs = null;
//        }
//    }
    private SpectrumComparison() {}

    //get new scan read in
//    public void reload(PeptideObj pobj, float[] eMZs, float[] eIntensities) {
//        pepObj = pobj;
//
//        unitNormMatchedIntensities = null;
//        unitNormPredIntensities = null;
//        sum1MatchedIntensities = null;
//        sum1PredIntensities = null;
//
//        allMatchedIntensities = null;
//        matchedIntensities = getMatchedIntensities(eMZs, eIntensities, predMZs);
//    }

    //TODO sort
    public SpectrumComparison pickedPredicted() {
        SpectrumComparison sc = new SpectrumComparison();

        sc.predIntensities = new float[(int) (predIntensities.length * Constants.bootstrapFragmentProportion)];
        sc.matchedIntensities = new float[sc.predIntensities.length];

        int[] rand = new Random().ints(0, predIntensities.length).distinct().
                limit(sc.predIntensities.length).toArray();


        for (int i = 0; i < rand.length; i++) {
            sc.predIntensities[i] = predIntensities[rand[i]];
            sc.matchedIntensities[i] = matchedIntensities[rand[i]];
        }

        return sc;
    }

    private void squareRootPredIntensities() {
        for (int i = 0; i < predIntensities.length; i++) {
            predIntensities[i] = (float) Math.sqrt(predIntensities[i]);
        }
    }

    /* Get best peaks from experimental spectrum that match to predicted peaks.
           Same experimental peak may match to the multiple predicted peaks,
              if they're close enough and experimental peak is strong.
           Unmatched peaks assigned 0
    */
    private float[] getMatchedIntensities(float[] expMZs, float[] expIntensities,
                                          float[] predMZs, String[] predFragmentIonTypes, PeptideObj pepObj) {
        if (predMZs.length == 1) {
            return predMZs; // I think we can return predMZs here instead of intensities array because a length 1 array will just return score of 1
        }
        int startPos = 0;
        int matchedNum = 0;
        float[] matchedInts = new float[predMZs.length];

        //based on fragment ion type, decide if we need to use ppm or Da tolerance
        boolean[] usePPMs = new boolean[predFragmentIonTypes.length];
        boolean mwd;
        for (int i = 0; i < predFragmentIonTypes.length; i++) {
            if (pepObj.daltonMatching) {
                usePPMs[i] = false;
                continue;
            }
            if (FragmentIonConstants.primaryFragmentIonTypes.isEmpty()) {
                mwd = Constants.matchWithDaltonsDefault;
            } else {
                if (FragmentIonConstants.primaryFragmentIonTypes.contains(predFragmentIonTypes[i])) {
                    mwd = Constants.matchWithDaltons;
                } else { //use aux spectrum model
                    mwd = Constants.matchWithDaltonsAux;
                }
            }
            usePPMs[i] = !mwd;
        }

        //matching fragments
        double ppm = Constants.ppmTolerance / 1000000;

        for (int i = 0; i < predMZs.length; i++) {
            double mz = predMZs[i];
            boolean usePPM = usePPMs[i];

            if (usePPM) {
                //see if any experimental peaks in vicinity
                //double fragmentError = ppm * mz;
                double fragmentMin = mz * (1 - ppm);
                double fragmentMax = mz * (1 + ppm);

                float matchedInt = 0;
                int pastStart = 0;

                while (startPos + pastStart < expMZs.length) {
                    double startMass = expMZs[startPos + pastStart];

                    if (startMass < fragmentMin) { //yet to reach peak within fragment tolerance
                        startPos += 1;
                    } else if (startMass <= fragmentMax) { //peak within fragment tolerance
                        //only for use when removing peaks from lower ranks
                        if (Constants.removeRankPeaks) {
                            matchedIdx.add(startPos + pastStart);
                        }

                        float potentialInt = expIntensities[startPos + pastStart];

                        if (potentialInt > matchedInt) { //new maximum intensity
                            matchedInt = potentialInt;
                        }
                        pastStart += 1;
                    } else { //outside of fragment tolerance range again
                        break;
                    }
                }

                matchedInts[matchedNum] = matchedInt;
                matchedNum += 1;
            } else {
                //see if any experimental peaks in vicinity
                //double fragmentError = ppm * mz;
                double fragmentMin = mz - Constants.DaTolerance;
                double fragmentMax = mz + Constants.DaTolerance;

                float matchedInt = 0;
                int pastStart = 0;

                while (startPos + pastStart < expMZs.length) {
                    double startMass = expMZs[startPos + pastStart];

                    if (startMass < fragmentMin) { //yet to reach peak within fragment tolerance
                        startPos += 1;
                    } else if (startMass <= fragmentMax) { //peak within fragment tolerance
                        //only for use when removing peaks from lower ranks
                        if (Constants.removeRankPeaks) {
                            matchedIdx.add(startPos + pastStart);
                        }

                        float potentialInt = expIntensities[startPos + pastStart];

                        if (potentialInt > matchedInt) { //new maximum intensity
                            matchedInt = potentialInt;
                        }
                        pastStart += 1;
                    } else { //outside of fragment tolerance range again
                        break;
                    }
                }

                matchedInts[matchedNum] = matchedInt;
                matchedNum += 1;
            }
        }

        return matchedInts;
    }

    private void getAllMatchedIntensities() {
        if (allMatchedIntensities == null) {
            //what fragment ion types do we need?
            HashSet<String> fiontypes = new HashSet<>(Arrays.asList(this.predFragmentIonTypes));

            //calculate mzs of possible fragments
            MassCalculator mc = new MassCalculator(pepObj.name.split("\\|")[0], pepObj.charge);
            TreeSet<Float> mzsSet = new TreeSet<>();
            HashMap<Float, String> mzToIontype = new HashMap<>();

            //iterate over fragment ion types
            for (String fiontype : fiontypes) {
                if (
                        (FragmentIonConstants.primaryFragmentIonTypes.contains(fiontype) &&
                                Constants.spectraModel.equalsIgnoreCase("unispec")) ||
                                (!FragmentIonConstants.primaryFragmentIonTypes.contains(fiontype) &&
                                        Constants.auxSpectraModel.equalsIgnoreCase("unispec")) ||
                                FragmentIonConstants.annotatePredfullLikeUnispec
                ) {
                    ArrayList<Float> mzsList = mc.possibleUnispecMzs(fiontype);
                    mzsSet.addAll(mzsList);
                    for (float mz : mzsList) {
                        mzToIontype.put(mz, fiontype);
                    }
                } else {
                    ArrayList<Float> mzsList = mc.possibleFragmentIons(fiontype);
                    mzsSet.addAll(mzsList);
                    for (float mz : mzsList) {
                        mzToIontype.put(mz, fiontype);
                    }
                }
            }

            //filter out peaks out of ms2 m/z range
            HashSet<Float> excludedMzs = new HashSet<>();
            for (float mz : mzsSet) {
                if (mz < pepObj.scanNumObj.lowerLimit || mz > pepObj.scanNumObj.upperLimit) {
                    excludedMzs.add(mz);
                }
            }
            for (Float mz : excludedMzs) {
                mzsSet.remove(mz);
            }

            float[] mzs = new float[mzsSet.size()];
            int i = 0;
            for (float f : mzsSet) {
                mzs[i] = f;
                i++;
            }

            String[] fragmentIonTypes = new String[mzs.length];
            for (int j = 0; j < fragmentIonTypes.length; j++) {
                fragmentIonTypes[j] = mzToIontype.get(mzs[j]);
            }

            allMatchedIntensities = getMatchedIntensities(
                    pepObj.scanNumObj.getSavedExpMZs(), pepObj.scanNumObj.getSavedExpIntensities(),
                    mzs, fragmentIonTypes, pepObj);
        }
    }


    public double[] getWeights(double[] freqs) {
        //will need to use predMZs
        int maxIndex = freqs.length;

        double[] weights = new double[predMZs.length];
        for (int i = 0; i < predMZs.length; i++) {
            int binIndex = (int) Math.floor(predMZs[i] / Constants.binwidth);
            if (binIndex < maxIndex) {
                weights[i] = freqs[binIndex];
            } else { //detected too big a fragment
                weights[i] = 0; //arbitrary, just ignore
            }
        }

        return weights;
    }

    private float[][] filterFragments(int top) {
        printError("Reimplement preprocessFragments!");
        System.exit(1);
//        top = Math.min(top, sortedIndicesList.size());
//
//        //calculate
//        float[] newMatched = new float[top];
//        float[] newPred = new float[top];
//
//        if (Constants.adaptiveFragmentNum && predIntensities.length > top) {
//            //function to filter vectors, and methods will work with these filtered vectors
//            for (int i = 0; i < top; i++) {
//                int index = sortedIndicesList.indexOf(i);
//                newMatched[i] = matchedIntensities[index];
//                newPred[i] = predIntensities[index];
//            }
//        } else {
//            newPred = predIntensities;
//            newMatched = matchedIntensities;
//        }
//
//        return new float[][]{newPred, newMatched};
        return new float[][]{};
    }

    private static float[] unitNormalize(float[] vector) {
        //if size 1
        if (vector.length == 1) {
            return vector;
        }

        //if we wish to normalize to unit vector
        double magnitude = 0;
        for (double i : vector) {
            magnitude += i * i;
        }
        magnitude =  Math.sqrt(magnitude);

        float[] newVector = new float[vector.length];

        if (magnitude != 0) { //fixes Bray curtis
            float mag = (float) magnitude;
            for (int i = 0; i < newVector.length; i++) {
                newVector[i] = vector[i] / mag;
            }
        }

        return newVector;
    }

    public void unitNormalizeIntensities() { //only use after filtering
        unitNormPredIntensities = unitNormalize(predIntensities);
        unitNormMatchedIntensities = unitNormalize(matchedIntensities);
    }

    private static float[] oneNormalize(float[] vector) {
        //if size 1
        if (vector.length == 1) {
            return new float[]{1};
        }

        //if we wish to normalize to unit vector
        double total = 0;
        for (double i : vector) {
            total += i;
        }

        float[] newVector = new float[vector.length];

        if (total != 0) { //fixes Bray curtis
            float t = (float) total;
            for (int i = 0; i < newVector.length; i++) {
                newVector[i] = vector[i] / t;
            }
        }

        return newVector;
    }

    public void oneNormalizeIntensities(float[] predIs, float[] matchIs) {
        sum1MatchedIntensities = oneNormalize(matchIs);
        sum1PredIntensities = oneNormalize(predIs);
    }

    //if constants.dividefragments split is 2 or more,
    //create a new spectrumComparison object for each split
    //and calculate metric for both.
    //If percolatorFormatter gets multiple values back, then write them separately
    public double cosineSimilarity() {
        if (predIntensities.length < 2) {
            return 0;
        }
        //numerator
        double num = 0;
        for (int i = 0; i < predIntensities.length; i++) {
            num += predIntensities[i] * matchedIntensities[i];
        }

        //denominator
        double a = 0;
        double b = 0;
        for (int i = 0; i < predIntensities.length; i++) {
            a += predIntensities[i] * predIntensities[i];
            b += matchedIntensities[i] * matchedIntensities[i];
        }
        double den = Math.sqrt(a * b);

        if (den == 0) { //fixes no matched peaks
            return 0;
        } else {
            return num / den;
        }
    }

    //https://stats.stackexchange.com/questions/384419/weighted-cosine-similarity
    public double weightedCosineSimilarity(double[] weights) {
        if (predIntensities.length < 2) {
            return 0;
        }
        //numerator
        double num = 0;
        for (int i = 0; i < predIntensities.length; i++) {
            num += predIntensities[i] * matchedIntensities[i] * weights[i];
        }

        //denominator
        double a = 0;
        double b = 0;
        for (int i = 0; i < predIntensities.length; i++) {
            a += predIntensities[i] * predIntensities[i] * weights[i];
            b += matchedIntensities[i] * matchedIntensities[i] * weights[i];
        }
        double den = Math.sqrt(a * b);

        if (den == 0) {
            return 0;
        } else {
            return num / den;
        }
    }

    public double spectralContrastAngle() {
        if (predIntensities.length < 2) {
            return 0;
        }
        double cosSim = this.cosineSimilarity();
        return 1 - (2 * Math.acos(cosSim) / Math.PI);
    }

    public double weightedSpectralContrastAngle(double[] weights) {
        if (predIntensities.length < 2) {
            return 0;
        }
        double cosSim = this.weightedCosineSimilarity(weights);
        return 1 - (2 * Math.acos(cosSim) / Math.PI);
    }

    public double euclideanDistance() {
        if (predIntensities.length < 2) {
            return 0;
        }
        if (unitNormPredIntensities == null) {
            this.unitNormalizeIntensities();
        }

        //max distance between two points in the positive quadrant with unit vectors is sqrt(2)
        float floatSum = 0.0f;
        for (float f : unitNormMatchedIntensities){
            floatSum += f;
        }
        if (floatSum == 0.0f) {
            return 1 - Math.sqrt(2);
        } else {
            double numSum = 0;
            for (int i = 0; i < predIntensities.length; i++) {
                double diff = unitNormPredIntensities[i] - unitNormMatchedIntensities[i];
                double square = diff * diff;
                numSum += square;
            }
            return 1 - Math.sqrt(numSum);
        }
    }

    public double weightedEuclideanDistance(double[] weights) {
        if (predIntensities.length < 2) {
            return 0;
        }
        if (unitNormPredIntensities == null) {
            this.unitNormalizeIntensities();
        }

        float floatSum = 0.0f;
        for (float f : unitNormMatchedIntensities){
            floatSum += f;
        }
        if (floatSum == 0.0f) {
            return 1 - Math.sqrt(2);
        } else {
            //now we need to normalize again by weights
            float[] newNormPred = new float[unitNormPredIntensities.length];
            float[] newNormMatched = new float[unitNormMatchedIntensities.length];
            for (int i = 0; i < unitNormPredIntensities.length; i++) {
                newNormPred[i] = (float) weights[i] * unitNormPredIntensities[i];
                newNormMatched[i] = (float) weights[i] * unitNormMatchedIntensities[i];
            }

            newNormPred = unitNormalize(newNormPred);
            newNormMatched = unitNormalize(newNormMatched);

            //now just do euclidean distance
            double numSum = 0;
            for (int i = 0; i < predIntensities.length; i++) {
                double diff = newNormPred[i] - newNormMatched[i];
                double square = diff * diff;
                numSum += square;
            }
            return 1 - Math.sqrt(numSum);
        }
    }

    public double brayCurtis() {
        if (predIntensities.length < 2) {
            return 0;
        }
        if (unitNormPredIntensities == null) {
            this.unitNormalizeIntensities();
        }

        //check if no matched peaks
        float floatSum = 0.0f;
        for (float f : unitNormMatchedIntensities){
            floatSum += f;
        }
        if (floatSum == 0.0f) {
            return 0;
        } else {
            double num = 0;
            double den = 0;
            for (int i = 0; i < predIntensities.length; i++) {
                double exp = unitNormMatchedIntensities[i];
                double pred = unitNormPredIntensities[i];

                num += Math.abs(exp - pred);
                den += (exp + pred);
            }
            return 1 - (num / den);
        }
    }

    public double weightedBrayCurtis(double[] weights) {
        if (predIntensities.length < 2) {
            return 0;
        }
        if (unitNormPredIntensities == null) {
            this.unitNormalizeIntensities();
        }

        //check if no matched peaks
        float floatSum = 0.0f;
        for (float f : unitNormMatchedIntensities){
            floatSum += f;
        }
        if (floatSum == 0.0f) {
            return 0;
        } else {
            double num = 0;
            double den = 0;
            for (int i = 0; i < predIntensities.length; i++) {
                double exp = unitNormMatchedIntensities[i];
                double pred = unitNormPredIntensities[i];

                num += weights[i] * Math.abs(exp - pred);
                den += weights[i] * (exp + pred);
            }
            return 1 - (num / den);
        }
    }

    public double pearsonCorr() {
        if (predIntensities.length < 2) {
            return -1;
        }
        if (Arrays.stream(floatToDouble(matchedIntensities)).sum() == 0 || matchedIntensities.length == 1) {
            return -1; //minimum pearson correlation
        } else {
            //uses Apache
            return pc.correlation(floatToDouble(matchedIntensities), floatToDouble(predIntensities));
        }
    }

    //top 36
    public double spearmanCorr() {
        if (predIntensities.length < 2) {
            return -1;
        }

        if (Arrays.stream(floatToDouble(matchedIntensities)).sum() == 0 || matchedIntensities.length == 1) {
            return -1;
        } else {
            //uses Apache
            return sc.correlation(floatToDouble(matchedIntensities), floatToDouble(predIntensities));
        }
    }

    public double weightedPearsonCorr(double[] weights) {
        if (predIntensities.length < 2) {
            return 0;
        }
        if (Arrays.stream(floatToDouble(matchedIntensities)).sum() == 0) {
            return -1; //minimum pearson correlation
        } else {

            double[] newPred = new double[predIntensities.length];
            double[] newMatched = new double[matchedIntensities.length];
            for (int i = 0; i < predIntensities.length; i++) {
                newPred[i] = weights[i] * predIntensities[i];
                newMatched[i] = weights[i] * matchedIntensities[i];
            }
            return pc.correlation(newMatched, newPred);
        }
    }

    public double dotProduct() {
        if (predIntensities.length < 2) {
            return 0;
        }
        if (unitNormPredIntensities == null) {
            this.unitNormalizeIntensities();
        }

        boolean nonzero = false;
        for (float f : matchedIntensities){
            if (f > 0) {
                nonzero = true;
                break;
            }
        }
        if (nonzero) {
            double num = 0;
            for (int i = 0; i < predIntensities.length; i++) {
                num += unitNormPredIntensities[i] * unitNormMatchedIntensities[i];
            }

            return num;
        } else {
            return 0;
        }
    }

    public double weightedDotProduct(double[] weights) {
        if (predIntensities.length < 2) {
            return 0;
        }
        float floatSum = 0.0f;
        for (float f : matchedIntensities){
            floatSum += f;
        }
        if (floatSum == 0.0f) {
            return 0;
        } else {
            double[] newPred = new double[predIntensities.length];
            double[] newMatched = new double[matchedIntensities.length];
            for (int i = 0; i < predIntensities.length; i++) {
                newPred[i] = weights[i] * predIntensities[i];
                newMatched[i] = weights[i] * matchedIntensities[i];
            }

            double predMax = 1 / Arrays.stream(newPred).max().getAsDouble();
            double matchedMax = 1 / Arrays.stream(newMatched).max().getAsDouble();
            double multiplier = predMax * matchedMax;

            double num = 0;
            for (int i = 0; i < predIntensities.length; i++) {
                num += newPred[i] * newMatched[i] * multiplier;
            }

            return num;
        }
    }

    private double spectralEntropy(float[] vector) {
        double entropy = 0;
        for (float f : vector) {
            if (f != 0) { //log(0) problematic
                entropy += (f * Math.log(f));
            }
        }
        return -1 * entropy;
    }

    //top 12
    public double unweightedSpectralEntropy() { //from https://www.nature.com/articles/s41592-021-01331-z
        if (predIntensities.length < 2) {
            return 0;
        }

        if (sum1PredIntensities == null) {
            oneNormalizeIntensities(predIntensities, matchedIntensities);
        }

        float[] SabVector = new float[sum1PredIntensities.length];
        int numFrags = 0;
        for (float j : sum1MatchedIntensities) {
            if (j != 0) {
                numFrags += 1;
            }
        }

        if (numFrags < 2) {
            return 0;
        } else {
            for (int i = 0; i < SabVector.length; i++) {
                SabVector[i] = (sum1PredIntensities[i] + sum1MatchedIntensities[i]) / 2;
            }
        }

        return 1 - ( ((2 * spectralEntropy(SabVector)) - spectralEntropy(sum1MatchedIntensities) - spectralEntropy(sum1PredIntensities)) / Math.log(4));
    }

    public double unweightedSpectralEntropy(float[] predicted, float[] matched) {
        if (predicted.length < 2) {
            return 0;
        }

        float[] SabVector = new float[predicted.length];
        int numFrags = 0;
        for (float j : matched) {
            if (j != 0) {
                numFrags += 1;
            }
        }

        if (numFrags < 2) {
            return 0;
        } else {
            for (int i = 0; i < SabVector.length; i++) {
                SabVector[i] = (predicted[i] + matched[i]) / 2;
            }
        }

        return 1 - ( ((2 * spectralEntropy(SabVector)) - spectralEntropy(matched) - spectralEntropy(predicted)) / Math.log(4));
    }

    public double weightedSpectralEntropy() {
        double unweighted = unweightedSpectralEntropy();
        if (unweighted == 0) {
            return 0;
        }
        return unweighted * Math.pow(spectralEntropy(sum1PredIntensities), 0.5);
    }

    public double heuristicSpectralEntropy() {
        if (sum1PredIntensities == null) {
            oneNormalizeIntensities(predIntensities, matchedIntensities);
        }

        double predEntropy = spectralEntropy(sum1PredIntensities);
        if (predEntropy < 1.75) {
            //reweighting
            double power = predEntropy / 2.75;

            float[] heuristicPred = new float[predIntensities.length];
            float[] heuristicMatched = new float[matchedIntensities.length];

            for (int i = 0; i < predIntensities.length; i++) {
                heuristicPred[i] = (float) Math.pow(predIntensities[i], power);
                heuristicMatched[i] = (float) Math.pow(matchedIntensities[i], power);
            }

            return unweightedSpectralEntropy(oneNormalize(heuristicPred), oneNormalize(heuristicMatched));
        } else {
            return unweightedSpectralEntropy();
        }
    }

    //top 24
    public double hypergeometricProbability() {
        if (predIntensities.length <= 1) {
            return 0;
        }
        this.getAllMatchedIntensities();

        if (allMatchedIntensities.length == 0) {
            HashSet<String> ftypes = new HashSet<>(Arrays.asList(predFragmentIonTypes));
            if (ftypes.size() == 1 && ftypes.contains("unknown")) {
                return 0; //predictions only include unannotated fragments
            } else {
                printInfo(pepObj.name + " contains only fragment ion types " + ftypes);
            }
        }

        int matchedIons = 0;
        for (float f : allMatchedIntensities) {
            if (f != 0) { //TODO: consider filtering for values above intensity threshold, same as predIntensity filtering
                matchedIons += 1;
            }
        }

        int successes = 0;
        int possible = 0; //number of predicted fragments that aren't unknown
        for (int i = 0; i < matchedIntensities.length; i++) {
            if (!Objects.equals(predFragmentIonTypes[i], "unknown")) {
                possible++;
                float f = matchedIntensities[i];
                if (f != 0) {
                    successes++;
                }
            }
        }
        if (successes > matchedIons) { //sig fig issue
            successes = matchedIons;
        }

        HypergeometricDistribution hgd = new HypergeometricDistribution(rng,
                allMatchedIntensities.length, matchedIons, possible);

        return -1 * Math.log10(hgd.upperCumulativeProbability(successes));
    }

    //top 12
    //biased towards longer peptides?
    public double intersection() {
        if (predIntensities.length <= 1) {
            return 0;
        }

        //calculate overlap of top predicted mzs and top matched mzs
        HashSet<Float> predSet = new HashSet<>();
        for (float f : matchedIntensities) {
            predSet.add(f);
        }

        this.getAllMatchedIntensities();
        float intersection = 0;
        float[] sortedIntensities = new float[allMatchedIntensities.length];
        for (int i = 0; i < allMatchedIntensities.length; i++) {
            sortedIntensities[i] = allMatchedIntensities[i];
        }
        Arrays.sort(sortedIntensities);
        float iters = 0;
        for (int i = sortedIntensities.length - 1; i > 0; i--) {
            if (sortedIntensities[i] == 0) { //using intensities directly, not m/z, since predicted m/z not saved
                break;
            }
            if (predSet.contains(sortedIntensities[i])) {
                intersection += 1;
            }
            iters += 1;
            if (iters >= Constants.topFragments) { //do we want to impose a ceiling?
                break;
            }
        }

        //return intersection / (matchedI.length + iters - intersection); //this would be jaccard
        return intersection;
    }

    //generic way of getting score
    public double getScore(String score) throws IOException, URISyntaxException {
        double returnScore = 0;
        switch (score) {
            case "brayCurtis":
                returnScore = brayCurtis();
                break;
            case "cosineSimilarity":
                returnScore = cosineSimilarity();
                break;
            case "spectralContrastAngle":
                returnScore = spectralContrastAngle();
                break;
            case "euclideanDistance":
                returnScore = euclideanDistance();
                break;
            case "pearsonCorr":
                returnScore = pearsonCorr();
                break;
            case "spearmanCorr":
                returnScore = spearmanCorr();
                break;
            case "hypergeometricProbability":
                returnScore = hypergeometricProbability();
                break;
            case "intersection":
                returnScore = intersection();
                break;
            case "dotProduct":
                returnScore = dotProduct();
                break;
            case "unweightedSpectralEntropy":
                returnScore = unweightedSpectralEntropy();
                break;
            case "weightedSpectralEntropy":
                returnScore = weightedSpectralEntropy();
                break;
            case "heuristicSpectralEntropy":
                returnScore = heuristicSpectralEntropy();
                break;
            default:
                printError("No score called " + score + ". Exiting");
                System.exit(1);
        }
        return returnScore;
    }

    public static void main(String[] args) {
        FragmentIonConstants.makeFragmentIonHierarchy();
        FragmentIonConstants.primaryFragmentIonTypes.add("y");
        Constants.matchWithDaltons = true;
        Constants.matchWithDaltonsAux = true;
        SpectrumComparison sc = new SpectrumComparison(new PeptideObj(),
                new float[]{10f, 20f, 30f}, new float[]{1f, 1f, 1f},
                new float[]{9.9997f, 19.9997f, 30.04f}, new float[]{1f, 1f, 1f}, new String[]{"y", "y", "y-NL"}, false);
        System.out.println(Arrays.toString(sc.matchedIntensities));
    }
}
