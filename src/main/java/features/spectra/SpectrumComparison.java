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
import mainsteps.PeptideObj;
import org.apache.commons.math3.distribution.HypergeometricDistribution;
import org.apache.commons.math3.random.Well19937c;
import org.apache.commons.math3.stat.correlation.PearsonsCorrelation;
import org.apache.commons.math3.stat.correlation.SpearmansCorrelation;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;
import java.util.LinkedHashSet;
import java.util.Random;

import static utils.NumericUtils.floatToDouble;
import static utils.Print.printError;

//TODO: also square root intensities? Squaring intensities may help for single cell data
public class SpectrumComparison {
    float[] predMZs;
    float[] predIntensities;
    float[] matchedIntensities;
    float[] unitNormMatchedIntensities;
    float[] unitNormPredIntensities;
    float[] sum1MatchedIntensities;
    float[] sum1PredIntensities;
    float[] allMatchedIntensities;
    public HashMap<String, Double> scores = new HashMap<>();
    int length;
    int matchedIons;
    public LinkedHashSet<Integer> matchedIdx = new LinkedHashSet<Integer>();
    private static final PearsonsCorrelation pc = new PearsonsCorrelation();
    private static final SpearmansCorrelation sc = new SpearmansCorrelation();
    public ArrayList<SpectrumComparison> spectrumComparisons = new ArrayList<>();
    public static Well19937c rng = new Well19937c(123);
    public PeptideObj pepObj;

    public SpectrumComparison(PeptideObj pepObj, float[] eMZs, float[] eIntensities,
                              float[] pMZs, float[] pIntensities, int length) {
        predMZs = pMZs;
        predIntensities = pIntensities;
        matchedIntensities = this.getMatchedIntensities(eMZs, eIntensities, predMZs, predIntensities);
        this.length = length;
        predMZs = null;

        this.pepObj = pepObj;
        if (Constants.features.contains("adjacent") || Constants.features.contains("bestScan")) {
            MassCalculator mc = new MassCalculator(pepObj.name.split("\\|")[0], pepObj.charge);
            pepObj.precursorMz = (mc.mass + pepObj.charge * mc.proton) / pepObj.charge;
        }
    }

    //TODO: this could be replaced with just also passing pred aux spectra arrays
    public SpectrumComparison(PeptideObj pepObj, float[] eMZs, float[] eIntensities,
                              float[] pMZs, float[] pIntensities, int length,
                              String[] fragmentIonTypes) {
        predMZs = pMZs;
        predIntensities = pIntensities;

        if (Constants.divideFragments.equals("0")) {
            matchedIntensities = this.getMatchedIntensities(eMZs, eIntensities, predMZs, predIntensities);
        } else {
            String[] fragmentsSplit = Constants.divideFragments.split(";");

            //get fragments that match the allowed
            for (String fragments : fragmentsSplit) {
                String[] allowedTypes = fragments.split("_");
                HashSet<String> allowedTypesSet = new HashSet<String>(Arrays.asList(allowedTypes));
                ArrayList<Integer> acceptedIdx = new ArrayList<>();

                for (int i = 0; i < fragmentIonTypes.length; i++) {
                    String ion = fragmentIonTypes[i];
                    if (allowedTypesSet.contains(ion)) {
                        acceptedIdx.add(i);
                    }
                }
                float[] mzs = new float[acceptedIdx.size()];
                float[] ints = new float[acceptedIdx.size()];
                for (int i = 0; i < acceptedIdx.size(); i++) {
                    mzs[i] = predMZs[acceptedIdx.get(i)];
                    ints[i] = predIntensities[acceptedIdx.get(i)];
                }

                //create new spectrumComparison obj and add to list
                spectrumComparisons.add(new SpectrumComparison(pepObj, eMZs, eIntensities,
                        mzs, ints, length));
            }
        }
        this.length = length;
        predMZs = null;

        this.pepObj = pepObj;
        if (Constants.features.contains("adjacent") || Constants.features.contains("bestScan")) {
            MassCalculator mc = new MassCalculator(pepObj.name.split("\\|")[0], pepObj.charge);
            pepObj.precursorMz = (mc.mass + pepObj.charge * mc.proton) / pepObj.charge;
        }
    }

    //TODO: if ever reimplement adjacent similarity, think of solution that does not require another constructor
    public SpectrumComparison(PeptideObj peptideObj, float[] eMZs, float[] eIntensities,
                              float[] pMZs, float[] pIntensities, int length, boolean willReload) {
        pepObj = peptideObj;
        predMZs = pMZs;
        predIntensities = pIntensities;

        matchedIntensities = this.getMatchedIntensities(eMZs, eIntensities, predMZs, predIntensities);
        this.length = length;
        if (! willReload) {
            predMZs = null;
        }
    }
    private SpectrumComparison() {}

    //get new scan read in
    public void reload(PeptideObj pobj, float[] eMZs, float[] eIntensities) {
        pepObj = pobj;

        unitNormMatchedIntensities = null;
        unitNormPredIntensities = null;
        sum1MatchedIntensities = null;
        sum1PredIntensities = null;

        allMatchedIntensities = null;
        matchedIntensities = getMatchedIntensities(eMZs, eIntensities, predMZs, predIntensities);
    }

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

    private float[] getMatchedIntensities(float[] expMZs, float[] expIntensities,
                                          float[] predMZs, float[] predIntensities) {
        if (predIntensities.length == 1) {
            return predIntensities;
        }
        int startPos = 0;
        int matchedNum = 0;
        float[] matchedInts = new float[predMZs.length];
        if (! Constants.matchWithDaltons) {
            double ppm = Constants.ppmTolerance / 1000000;

        /* Get best peaks from experimental spectrum that match to predicted peaks.
           Same experimental peak may match to the multiple predicted peaks,
              if they're close enough and experimental peak is strong.
           Unmatched peaks assigned 0
         */
            for (double mz : predMZs) {
                //see if any experimental peaks in vicinity
                //double fragmentError = ppm * mz;
                double fragmentMin = mz * (1 - ppm);
                double fragmentMax = mz * (1 + ppm);

                float predInt = 0;
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

                        if (potentialInt > predInt) { //new maximum intensity
                            predInt = potentialInt;
                        }
                        pastStart += 1;
                    } else { //outside of fragment tolerance range again
                        break;
                    }
                }

                matchedInts[matchedNum] = predInt;
                matchedNum += 1;
            }
        } else {
            for (double mz : predMZs) { //TODO: because fragments besides unknown have correct m/z, is this unnecessary?
                //see if any experimental peaks in vicinity
                //double fragmentError = ppm * mz;
                double fragmentMin = mz - Constants.DaTolerance;
                double fragmentMax = mz + Constants.DaTolerance;

                float predInt = 0;
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

                        if (potentialInt > predInt) { //new maximum intensity
                            predInt = potentialInt;
                        }
                        pastStart += 1;
                    } else { //outside of fragment tolerance range again
                        break;
                    }
                }

                matchedInts[matchedNum] = predInt;
                matchedNum += 1;
            }
        }
        return matchedInts;
    }

    private void getAllMatchedIntensities() {
        if (allMatchedIntensities == null) {
            MassCalculator mc = new MassCalculator(pepObj.name.split("\\|")[0], pepObj.charge);
            //calculate y and b m/zs
            float[] mzs = new float[4 * (mc.peptide.length() - 1)];
            String[] flags = {"y", "b"};
            int i = 0;
            for (int num = 1; num < mc.peptide.length(); num++) {
                for (int charge = 1; charge < 3; charge++) {
                    for (String flag : flags) {
                        mzs[i] = mc.calcMass(num, flag, charge, 0);
                        i += 1;
                    }
                }
            }
            //sort
            Arrays.sort(mzs);

            allMatchedIntensities = getMatchedIntensities(pepObj.scanNumObj.getExpMZs(), pepObj.scanNumObj.getExpIntensities(),
                    mzs, new float[mzs.length]);
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
        printError("Reimplement filterFragments!");
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

    //top 24
    public double hyperGeometricProbability() {
        this.getAllMatchedIntensities();
        matchedIons = 0;
        for (float f : allMatchedIntensities) {
            if (f != 0) {
                matchedIons += 1;
            }
        }

        HypergeometricDistribution hgd = new HypergeometricDistribution(rng,
                4 * (length - 1), matchedIons, predIntensities.length);
        int successes = 0;
        for (float f : matchedIntensities) {
            if (f != 0) {
                successes += 1;
            }
        }
        return -1 * Math.log10(hgd.upperCumulativeProbability(successes));
    }

    //top 12
    public double intersection() {
        //calculate overlap of top predicted mzs and top matched mzs
        HashSet<Float> predSet = new HashSet<>();
        for (float f : matchedIntensities) {
            predSet.add(f);
        }

        this.getAllMatchedIntensities();
        float intersection = 0;
        Arrays.sort(allMatchedIntensities);
        float iters = 0;
        for (int i = allMatchedIntensities.length - 1; i > 0; i--) {
            if (allMatchedIntensities[i] == 0) {
                break;
            }
            if (predSet.contains(allMatchedIntensities[i])) {
                intersection += 1;
            }
            iters += 1;
            if (iters >= Constants.topFragments) {
                break;
            }
        }

        //return intersection / (matchedI.length + iters - intersection); //this would be jaccard
        return intersection;
    }

    public void clearArrays() {
        predMZs = null;
        predIntensities = null;
        matchedIntensities = null;
        unitNormMatchedIntensities = null;
        unitNormPredIntensities = null;
        sum1MatchedIntensities = null;
        sum1PredIntensities = null;
        allMatchedIntensities = null;
    }
}
