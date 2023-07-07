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

import org.apache.commons.math3.distribution.HypergeometricDistribution;
import org.apache.commons.math3.random.Well19937c;
import org.apache.commons.math3.stat.correlation.PearsonsCorrelation;
import org.apache.commons.math3.stat.correlation.SpearmansCorrelation;

import java.util.*;
import java.util.stream.IntStream;

import static Features.FloatUtils.floatToDouble;
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
    HashMap<String, Double> scores = new HashMap<>();
    int length;
    int matchedIons;
    LinkedHashSet<Integer> matchedIdx = new LinkedHashSet<Integer>();
    private static ArrayList<Float> tmpMZs = new ArrayList<Float>();
    private static ArrayList<Float> tmpInts = new ArrayList<Float>();
    private static PearsonsCorrelation pc = new PearsonsCorrelation();
    private static SpearmansCorrelation sc = new SpearmansCorrelation();
    public ArrayList<SpectrumComparison> spectrumComparisons = new ArrayList<>();
    public static Well19937c rng = new Well19937c(123);

    public PeptideObj pepObj;
    ArrayList<Integer> sortedIndicesList = new ArrayList<>();

    public SpectrumComparison(PeptideObj pepObj, float[] eMZs, float[] eIntensities,
                              float[] pMZs, float[] pIntensities, int length,
                              boolean filterTop, boolean filterBase) {
        predMZs = pMZs;
        predIntensities = pIntensities;

//        if (Constants.sqrtPredIntensities) {
//            this.squareRootPredIntensities();
//        }
        if (filterBase) {
            this.filterIntensitiesByPercentage(Constants.percentBasePeak);
        }
        if (filterTop) {
            this.filterTopFragments();
        }

        int[] sortedIndices = IntStream.range(0, predMZs.length)
                .boxed().sorted((k, j) -> Float.compare(predMZs[k], predMZs[j]))
                .mapToInt(ele -> ele).toArray();
        for (int i : sortedIndices) {
            sortedIndicesList.add(i);
        }
        pMZs = predMZs;
        pIntensities = predIntensities;
        predMZs = new float[predMZs.length];
        predIntensities = new float[predIntensities.length];
        for (int i = 0; i < sortedIndices.length; i++) {
            predMZs[i] = pMZs[sortedIndices[i]];
            predIntensities[i] = pIntensities[sortedIndices[i]];
        }

        matchedIntensities = this.getMatchedIntensities(eMZs, eIntensities, predMZs, predIntensities);
        this.length = length;
        tmpMZs.clear();
        tmpInts.clear();
        predMZs = null;

        this.pepObj = pepObj;
    }

    public SpectrumComparison(PeptideObj pepObj, float[] eMZs, float[] eIntensities,
                              float[] pMZs, float[] pIntensities, int length,
                              boolean filterTop, boolean filterBase, String[] fragmentIonTypes) {
        predMZs = pMZs;
        predIntensities = pIntensities;

        if (Constants.divideFragments.equals("0")) {
//            if (Constants.sqrtPredIntensities) {
//                this.squareRootPredIntensities();
//            }
            if (filterBase) {
                this.filterIntensitiesByPercentage(Constants.percentBasePeak);
            }
            if (filterTop) {
                this.filterTopFragments(); //set higher in these cases, then use another top peaks when making fragment ion subsets
            }

            int[] sortedIndices = IntStream.range(0, predMZs.length)
                    .boxed().sorted((k, j) -> Float.compare(predMZs[k], predMZs[j]))
                    .mapToInt(ele -> ele).toArray();
            for (int i : sortedIndices) {
                sortedIndicesList.add(i);
            }
            pMZs = predMZs;
            pIntensities = predIntensities;
            //String[] fIonTypes = fragmentIonTypes;

            predMZs = new float[predMZs.length];
            predIntensities = new float[predIntensities.length];
            //fragmentIonTypes = new String[fragmentIonTypes.length];
            for (int i = 0; i < sortedIndices.length; i++) {
                predMZs[i] = pMZs[sortedIndices[i]];
                predIntensities[i] = pIntensities[sortedIndices[i]];
                //fragmentIonTypes[i] = fIonTypes[sortedIndices[i]];
            }

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

                //change between ppm and Da tolerance when matching.
                // Not the most accurate, as it assumes when using PredFull + another tool, PredFull is not used for it y/b predictions
//                if (! Constants.spectraRTPredModel.equals("PredFull")) {
//                    allowedTypesSet.remove("y");
//                    allowedTypesSet.remove("b");
//                    if (allowedTypesSet.size() != 0) { //other types only provided by PredFull are being used
//                        Constants.matchWithDaltons = true;
//                    } else {
//                        Constants.matchWithDaltons = false;
//
//                    }
//                }

                //create new spectrumComparison obj and add to list
                spectrumComparisons.add(new SpectrumComparison(pepObj, eMZs, eIntensities,
                        mzs, ints, length, Constants.useTopFragments, Constants.useBasePeak));
            }
        }
        this.length = length;
        tmpMZs.clear();
        tmpInts.clear();
        predMZs = null;
    }

    private void filterTopFragments() {
        //stick with arraylist because finding minimum will be faster than linkedlist due to indexing
        //skip if shorter
        if (predMZs.length > Constants.topFragments) {
            tmpInts.clear();
            tmpMZs.clear();

            for (float i : predIntensities) {
                tmpInts.add(i);
            }
            for (float i : predMZs) {
                tmpMZs.add(i);
            }

            predIntensities = new float[Constants.topFragments];
            predMZs = new float[Constants.topFragments];

            for (int i = 0; i < Constants.topFragments; i++) {
                int index = tmpInts.indexOf(Collections.max(tmpInts));
                predIntensities[i] = tmpInts.get(index);
                predMZs[i] = tmpMZs.get(index);
                tmpInts.set(index, -1f); //no longer highest intensity peak
            }
        }
    }

    public void filterIntensitiesByValue(float min) {
        tmpMZs.clear();
        tmpInts.clear();

        for (int i = 0; i < predIntensities.length; i++) {
            float intensity = predIntensities[i];

            if (intensity >= min) {
                tmpInts.add(intensity);
                tmpMZs.add(predMZs[i]);
            }
        }

        predMZs = new float[tmpMZs.size()];
        predIntensities = new float[tmpInts.size()];

        for (int i = 0; i < tmpInts.size(); i++) {
            predIntensities[i] = tmpInts.get(i);
            predMZs[i] = tmpMZs.get(i);
        }
    }

    public void filterIntensitiesByPercentage(float percentage) { //should be < 100
        if (percentage > 100) {
            System.out.println("percentBasePeak must be <= 100 but is set to " + Constants.percentBasePeak);
            System.exit(-1);
        } else {
            //get max intensity
            float maxIntensity = 0f;
            for (float f : predIntensities) {
                if (f > maxIntensity) {
                    maxIntensity = f;
                }
            }

            //make cutoff by percentage
            float cutoff = maxIntensity * percentage / 100f;
            filterIntensitiesByValue(cutoff);
        }
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
                        mzs[i] = mc.calcMass(num, flag, charge);
                        i += 1;
                    }
                }
            }
            //sort
            Arrays.sort(mzs);

            //modify getMatchedIntensities to be more general, get how many nonzero matched intensities
            //don't need length anymore?
            matchedIdx.clear();

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
        top = Math.min(top, sortedIndicesList.size());

        //calculate
        float[] newMatched = new float[top];
        float[] newPred = new float[top];

        if (Constants.adaptiveFragmentNum && predIntensities.length > top) {
            //function to filter vectors, and methods will work with these filtered vectors
            for (int i = 0; i < top; i++) {
                int index = sortedIndicesList.indexOf(i);
                newMatched[i] = matchedIntensities[index];
                newPred[i] = predIntensities[index];
            }
        } else {
            newPred = predIntensities;
            newMatched = matchedIntensities;
        }

        return new float[][]{newPred, newMatched};
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
            return 0;
        }
        if (Arrays.stream(floatToDouble(matchedIntensities)).sum() == 0 || matchedIntensities.length == 1) {
            return -1; //minimum pearson correlation
        } else {
            //uses Apache
            return pc.correlation(floatToDouble(matchedIntensities), floatToDouble(predIntensities));
        }
    }

    public double spearmanCorr() {
        if (predIntensities.length < 2) {
            return 0;
        }

        float[][] vectors = filterFragments(36);
        float[] predI = vectors[0];
        float[] matchedI = vectors[1];

        if (Arrays.stream(floatToDouble(matchedI)).sum() == 0 || matchedI.length == 1) {
            return -1;
        } else {
            //uses Apache
            return sc.correlation(floatToDouble(matchedI), floatToDouble(predI));
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

    public double unweightedSpectralEntropy() { //from https://www.nature.com/articles/s41592-021-01331-z
        if (predIntensities.length < 2) {
            return 0;
        }

        float[][] vectors = filterFragments(12);
        float[] predI = vectors[0];
        float[] matchedI = vectors[1];

        if (sum1PredIntensities == null) {
            oneNormalizeIntensities(predI, matchedI);
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

    public double hyperGeometricProbability() {
        this.getAllMatchedIntensities();
        matchedIons = 0;
        for (float f : allMatchedIntensities) {
            if (f != 0) {
                matchedIons += 1;
            }
        }

//        for (int idx = matchedIdx.size() - 1; idx > -1; idx--) {
//            pepObj.scanNumObj.expMZs = ArrayUtils.remove(pepObj.scanNumObj.expMZs, idx);
//            pepObj.scanNumObj.expIntensities = ArrayUtils.remove(pepObj.scanNumObj.expIntensities, idx);
//        }

        //calculate
        float[][] vectors = filterFragments(24);
        float[] predI = vectors[0];
        float[] matchedI = vectors[1];

        HypergeometricDistribution hgd = new HypergeometricDistribution(rng, 4 * (length - 1), matchedIons, predI.length);
        int successes = 0;
        for (float f : matchedI) {
            if (f != 0) {
                successes += 1;
            }
        }
        return -1 * Math.log10(hgd.upperCumulativeProbability(successes));
    }

    public double intersection() {
        //calculate overlap of top predicted mzs and top matched mzs
        int top = Constants.topFragments;
        if (Constants.adaptiveFragmentNum) {
            top = 12;
        }

        float[][] vectors = filterFragments(top);
        float[] matchedI = vectors[1];
        HashSet<Float> predSet = new HashSet<>();
        for (float f : matchedI) {
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
            if (iters >= top) {
                break;
            }
        }

        //return intersection / (matchedI.length + iters - intersection); //this would be jaccard
        return intersection;
    }
}
