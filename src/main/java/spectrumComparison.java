import org.apache.commons.lang.ArrayUtils;
import org.apache.commons.math.stat.correlation.PearsonsCorrelation;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;

public class spectrumComparison {
    double[] expMZs;
    double[] expIntensities;
    double[] predMZs;
    double[] predIntensities;
    double ppm;
    double[] matchedIntensities;
    double[] unitNormMatchedIntensities;
    double[] unitNormPredIntensities;

    public spectrumComparison(double[] eMZs, double[] eIntensities,
                              double[] pMZs, double[] pIntensities,
                              double ppmTolerance) {
        expMZs = eMZs;
        expIntensities = eIntensities;
        predMZs = pMZs;
        predIntensities = pIntensities;
        ppm =  ppmTolerance / 1000000;
        //could also consider just getting rid of low pred intensities before matching
        //this.filterIntensities(0.01);

        matchedIntensities = this.getMatchedIntensities();

        //predIntensities = this.rankIntensities(predIntensities);
        //matchedIntensities = this.rankIntensities(matchedIntensities);
    }

    public double[] getMatchedIntensities() {
        int startPos = 0;
        int matchedNum = 0;
        double[] matchedInts = new double[predMZs.length];

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

            double predInt = 0;
            int pastStart = 0;

            while (startPos + pastStart < expMZs.length) {
                double startMass = expMZs[startPos + pastStart];

                if (startMass < fragmentMin) { //yet to reach peak within fragment tolerance
                    startPos += 1;
                } else if (startMass <= fragmentMax) { //peak within fragment tolerance
                    double potentialInt = expIntensities[startPos + pastStart];

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
        return matchedInts;
    }

    public double[] getWeights(double[] freqs, double bin) {
        //will need to use predMZs
        int maxIndex = freqs.length;

        double[] weights = new double[predMZs.length];
        for (int i = 0; i < predMZs.length; i++) {
            int binIndex = (int) Math.floor(predMZs[i] / bin);
            if (binIndex < maxIndex) {
                weights[i] = freqs[binIndex];
            } else { //detected too big a fragment
                weights[i] = 0; //arbitrary, just ignore
            }
        }

        return weights;
    }

    private static double[] subNormalize(double[] vector) {
        //if we wish to normalize to unit vector
        double magnitude = 0;
        for (double i : vector) {
            magnitude += i * i;
        }
        magnitude =  Math.sqrt(magnitude);

        double[] newVector = new double[vector.length];

        if (magnitude != 0) { //fixes Bray curtis
            for (int i = 0; i < newVector.length; i++) {
                newVector[i] = vector[i] / magnitude;
            }
        }

        return newVector;
    }

    public void unitNormalize() { //only use after filtering
        unitNormPredIntensities = subNormalize(predIntensities);
        unitNormMatchedIntensities = subNormalize(matchedIntensities);
    }

    public void filterIntensities(double min) {

        ArrayList<Double> removableIntensities = new ArrayList<>();
        ArrayList<Double> removableMZs = new ArrayList<>();

        for (int i = 0; i < predIntensities.length; i++) {
            double intensity = predIntensities[i];
            double mz = predMZs[i];

            if (intensity < min) {
                removableIntensities.add(intensity);
                removableMZs.add(mz);
            }
        }

        for (int i = 0; i < removableIntensities.size(); i++) {
            predIntensities = ArrayUtils.removeElement(predIntensities, removableIntensities.get(i));
            predMZs = ArrayUtils.removeElement(predMZs, removableMZs.get(i));
        }
    }

    public double[] rankIntensities(double[] vector) { //need to adapt for unitNorm vectors too
        double[] ranks = new double[vector.length];

        double[] sortedVector = vector;
        Arrays.sort(sortedVector);

        for (int i = 0; i < sortedVector.length; i++) {
            ranks[i] = ArrayUtils.indexOf(sortedVector, vector[i]) + 1;
        }

        return ranks;
    }

    public double cosineSimilarity() {

        //numerator
        double num = 0;
        for (int i = 0; i < predMZs.length; i++) {
            num += predIntensities[i] * matchedIntensities[i];
        }

        //denominator
        double a = 0;
        double b = 0;
        for (int i = 0; i < predMZs.length; i++) {
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

        //numerator
        double num = 0;
        for (int i = 0; i < predMZs.length; i++) {
            num += predIntensities[i] * matchedIntensities[i] * weights[i];
        }

        //denominator
        double a = 0;
        double b = 0;
        for (int i = 0; i < predMZs.length; i++) {
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
        double cosSim = this.cosineSimilarity();
        return 1 - (2 * Math.acos(cosSim) / Math.PI);
    }

    public double weightedSpectralContrastAngle(double[] weights) {
        double cosSim = this.weightedCosineSimilarity(weights);
        return 1 - (2 * Math.acos(cosSim) / Math.PI);
    }

    public double euclideanDistance() {
        if (unitNormPredIntensities == null) {
            this.unitNormalize();
        }

        //max distance between two points in the positive quadrant with unit vectors is sqrt(2)
        if (Arrays.stream(unitNormMatchedIntensities).sum() == 0) {
            return 1 - Math.sqrt(2);
        } else {
            double numSum = 0;
            for (int i = 0; i < predMZs.length; i++) {
                double diff = unitNormPredIntensities[i] - unitNormMatchedIntensities[i];
                double square = diff * diff;
                numSum += square;
            }
            return 1 - Math.sqrt(numSum);
        }
    }

    public double weightedEuclideanDistance(double[] weights) {
        if (unitNormPredIntensities == null) {
            this.unitNormalize();
        }

        if (Arrays.stream(unitNormMatchedIntensities).sum() == 0) {
            return 1 - Math.sqrt(2);
        } else {
            //now we need to normalize again by weights
            double[] newNormPred = new double[unitNormPredIntensities.length];
            double[] newNormMatched = new double[unitNormMatchedIntensities.length];
            for (int i = 0; i < unitNormPredIntensities.length; i++) {
                newNormPred[i] = weights[i] * unitNormPredIntensities[i];
                newNormMatched[i] = weights[i] * unitNormMatchedIntensities[i];
            }

            newNormPred = subNormalize(newNormPred);
            newNormMatched = subNormalize(newNormMatched);

            //now just do euclidean distance
            double numSum = 0;
            for (int i = 0; i < predMZs.length; i++) {
                double diff = newNormPred[i] - newNormMatched[i];
                double square = diff * diff;
                numSum += square;
            }
            return 1 - Math.sqrt(numSum);
        }
    }

    public double brayCurtis() {
        if (unitNormPredIntensities == null) {
            this.unitNormalize();
        }

        //check if no matched peaks
        if (Arrays.stream(unitNormMatchedIntensities).sum() == 0) {
            return 0;
        } else {
            double num = 0;
            double den = 0;
            for (int i = 0; i < predMZs.length; i++) {
                double exp = unitNormMatchedIntensities[i];
                double pred = unitNormPredIntensities[i];

                num += Math.abs(exp - pred);
                den += (exp + pred);
            }
            return 1 - (num / den);
        }
    }

    public double weightedBrayCurtis(double[] weights) {
        if (unitNormPredIntensities == null) {
            this.unitNormalize();
        }

        //check if no matched peaks
        if (Arrays.stream(unitNormMatchedIntensities).sum() == 0) {
            return 0;
        } else {
            double num = 0;
            double den = 0;
            for (int i = 0; i < predMZs.length; i++) {
                double exp = unitNormMatchedIntensities[i];
                double pred = unitNormPredIntensities[i];

                num += weights[i] * Math.abs(exp - pred);
                den += weights[i] * (exp + pred);
            }
            return 1 - (num / den);
        }
    }

    public double pearsonCorr() {
        if (Arrays.stream(matchedIntensities).sum() == 0) {
            return -1; //minimum pearson correlation
        } else {
            //uses Apache
            PearsonsCorrelation pc = new PearsonsCorrelation();
            return pc.correlation(matchedIntensities, predIntensities);
        }
    }

    public double weightedPearsonCorr(double[] weights) {
        if (Arrays.stream(matchedIntensities).sum() == 0) {
            return -1; //minimum pearson correlation
        } else {

            double[] newPred = new double[predIntensities.length];
            double[] newMatched = new double[matchedIntensities.length];
            for (int i = 0; i < predIntensities.length; i++) {
                newPred[i] = weights[i] * predIntensities[i];
                newMatched[i] = weights[i] * matchedIntensities[i];
            }

            PearsonsCorrelation pc = new PearsonsCorrelation();
            return pc.correlation(newMatched, newPred);
        }
    }

    public double dotProduct() {
        if (Arrays.stream(matchedIntensities).sum() == 0) {
            return 0;
        } else {
            double predMax = 1 / Arrays.stream(predIntensities).max().getAsDouble();
            double matchedMax = 1 / Arrays.stream(matchedIntensities).max().getAsDouble();
            double multiplier = predMax * matchedMax;

            double num = 0;
            for (int i = 0; i < predMZs.length; i++) {
                num += predIntensities[i] * matchedIntensities[i] * multiplier;
            }

            return num;
        }
    }

    public double weightedDotProduct(double[] weights) {
        if (Arrays.stream(matchedIntensities).sum() == 0) {
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
            for (int i = 0; i < predMZs.length; i++) {
                num += newPred[i] * newMatched[i] * multiplier;
            }

            return num;
        }
    }

    public HashMap<String, Double> getAllSimilarities(double[] mzFreqs, int binwidth) {

        HashMap<String, Double> sims = new HashMap<>();

        sims.put("cosine", this.cosineSimilarity());
        sims.put("contrast", this.spectralContrastAngle());
        sims.put("euclidean", this.euclideanDistance());
        sims.put("bray-curtis", this.brayCurtis());
        sims.put("pearson", this.pearsonCorr());
        sims.put("dot", this.dotProduct());


        //weighted
        double[] weights = this.getWeights(mzFreqs, binwidth);
        sims.put("weightCosine", this.weightedCosineSimilarity(weights));
        sims.put("weightContrast", this.weightedSpectralContrastAngle(weights));
        sims.put("weightEuclidean", this.weightedEuclideanDistance(weights));
        sims.put("weightBray-curtis", this.weightedBrayCurtis(weights));
        sims.put("weightPearson", this.weightedPearsonCorr(weights));
        sims.put("weightDot", this.weightedDotProduct(weights));

        return sims;
    }

    public static void main(String[] args) {

    }
}
