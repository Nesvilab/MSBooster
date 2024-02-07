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

import com.github.sanity.pav.PairAdjacentViolators;
import com.github.sanity.pav.Point;
import kotlin.jvm.functions.Function1;
import org.apache.commons.math3.analysis.interpolation.LoessInterpolator;
import umontreal.ssj.gof.KernelDensity;
import umontreal.ssj.probdist.EmpiricalDist;
import umontreal.ssj.probdist.NormalDist;
import umontreal.ssj.randvar.KernelDensityGen;

import java.util.*;
import java.util.stream.IntStream;

public class StatMethods {

    long meanCounts = 0;
    long medianCounts = 0;
    double mean = 0;
    double M2 = 0;
    TreeMap<Double, Integer> medianHist = new TreeMap<>();
    final int medianBins = 10000;
    Double std = Double.NaN;
    Double median = Double.NaN;

    public StatMethods() {}

    public void updateVariance(double x) {
        meanCounts += 1;
        double delta = x - mean;
        mean += delta / meanCounts;
        M2 += delta * (x - mean);
    }

    public double getMean() {
        return mean;
    }

    public double getStd() {
        if (std.isNaN()) {
            std = Math.sqrt(M2 / (meanCounts - 1));
        }
        return std;
    }

    public void updateMedian(double x) {
        double key = Math.floor(x * medianBins);
        if (medianHist.containsKey(key)) {
            medianHist.put(key, medianHist.get(key) + 1);
        } else {
            medianHist.put(key, 1);
        }
        medianCounts += 1;
    }

    public double getMedian() {
        if (median.isNaN()) {
            int goal = (int) (medianCounts / 2);
            int counts = 0;
            for (Map.Entry<Double, Integer> entry : medianHist.entrySet()) {
                counts += entry.getValue();
                if (counts >= goal) {
                    median = entry.getKey() / medianBins;
                    medianHist.clear();
                    break;
                }
            }
        }
        return median;
    }

    public static float mean(double[] vector) {
        float meanX = 0;
        for (double x : vector) {
            meanX += x;
        }
        return meanX / (float) vector.length;
    }

    public static float mean(float[] vector) {
        float meanX = 0;
        for (float x : vector) {
            meanX += x;
        }
        return meanX / (float) vector.length;
    }

    public static float mean(ArrayList<Float> vector) {
        float meanX = 0;
        for (float x : vector) {
            meanX += x;
        }
        return meanX / vector.size();
    }

    public static float variance(double[] vector) {
        int vecLength = vector.length;
        float mean = mean(vector);
        float var = 0;

        for (double v : vector) {
            var += Math.pow(v - mean, 2);
        }

        return var / (vecLength - 1);
    }

    public static float variance(double[] vector, float mean) {
        int vecLength = vector.length;
        float var = 0;

        for (double v : vector) {
            var += Math.pow(v - mean, 2);
        }

        return var / (vecLength - 1);
    }

    public static float variance(ArrayList<Float> vector, float mean) {
        int vecLength = vector.size();
        float var = 0;

        for (float v : vector) {
            var += Math.pow(v - mean, 2);
        }

        return var / (vecLength - 1);
    }

    //covariance
    public static float variance(double[] vectorX, double[] vectorY) {
        int vecLength = vectorX.length;
        float meanX = mean(vectorX);
        float meanY = mean(vectorY);

        float var = 0;

        for (int i = 0; i < vecLength; i++) {
            var += (vectorX[i] - meanX) * (vectorY[i] - meanY);
        }

        return var / (vecLength - 1);
    }

    public static float variance(double[] vectorX, float meanX, double[] vectorY, float meanY) {
        int vecLength = vectorX.length;

        float var = 0;

        for (int i = 0; i < vecLength; i++) {
            var += (vectorX[i] - meanX) * (vectorY[i] - meanY);
        }

        return var / (vecLength - 1);
    }

    public static float[] linearRegression(double[] x, double[] y) {
        assert x.length == y.length : "vectors must be of same length";

        //returns beta 0 and beta 1
        float meanX = mean(x);
        float meanY = mean(y);

        float varX = variance(x, meanX);

        float covar = variance(x, meanX, y, meanY);

        float beta1 = covar / varX;
        float beta0 = meanY - (beta1 * meanX);

        return new float[] {beta0, beta1};
    }

    public static float zscore(float x, float mean, float sd) {
        return ((x - mean) / sd);
    }

    //TODO: automatic bandwidth selection
    public static Function1<Double, Double> LOESS(double[][] bins, double bandwidth, int robustIters) {
        if (bandwidth <= 0) {
            System.out.println("bandwidth is set to " + bandwidth + " but it must be greater than 0. Setting it to 0.05");
            bandwidth = 0.05;
        }
        if (bandwidth > 1) {
            System.out.println("bandwidth is set to " + bandwidth + " but maximum allowed is 1. Setting it to 1");
            bandwidth = 1;
        }

        //need to sort arrays (DIA-U mzml not in order)
        int[] sortedIndices = IntStream.range(0, bins[0].length)
                .boxed().sorted(Comparator.comparingDouble(k -> bins[0][k])).mapToInt(ele -> ele).toArray();
        double[] newX = new double[sortedIndices.length];
        double[] newY = new double[sortedIndices.length];
        for (int i = 0; i < sortedIndices.length; i++) {
            newX[i] = bins[0][sortedIndices[i]];
        }
        for (int i = 0; i < sortedIndices.length; i++) {
            newY[i] = bins[1][sortedIndices[i]];
        }

        //solve monotonicity issue
        double compare = -1;
        for (int i = 0; i < newX.length; i++) {
            double d = newX[i];
            if (d == compare) {
                newX[i] = newX[i - 1] + 0.00000001; //arbitrary increment to handle smooth method
            } else {
                compare = d;
            }
        }

        //fit loess
        if (newX.length < 50) {
            bandwidth = 1d;
        }

        LoessInterpolator loessInterpolator = new LoessInterpolator(bandwidth, robustIters);
        double[][] results = fittingRound(loessInterpolator, newX, newY, null, false);
        double[] fittedY = results[0];
        newX = results[1];
        newY = results[2];

        //remove Nan
        ArrayList<Integer> nanIdx = new ArrayList<>();
        int i = 0;
        for (double yval : fittedY) {
            if (Double.isNaN(yval)) {
                nanIdx.add(i);
            }
            i += 1;
        }

        fittedY = removeIdx(fittedY, nanIdx);
        newX = removeIdx(newX, nanIdx);
        newY = removeIdx(newY, nanIdx);

        Function1<Double, Double> ir = isotonicRegressor(newX, fittedY);
        for (i = 0; i < fittedY.length; i++) {
            fittedY[i] = ir.invoke(newX[i]);
        }
        results = fittingRound(loessInterpolator, newX, newY, fittedY, true);
        fittedY = results[0];
        newX = results[1];
        newY = results[2];

        return isotonicRegressor(newX, fittedY);
    }

    private static Function1<Double, Double> isotonicRegressor(double[] x, double[] y) {
        List<Point> points = new LinkedList<>();
        for (int i = 0; i < y.length; i++) {
            points.add(new Point(x[i], y[i]));
        }
        PairAdjacentViolators pav = new PairAdjacentViolators(points);

        return pav.interpolator(); //extrapolationStrategy can use flat or tangent (default tangent)
    }

    private static double[][] fittingRound(LoessInterpolator lr, double[] newX, double[] newY, double[] smoothedY,
                                           boolean extra) {
        double[] y;
        if (smoothedY != null) {
            y = smoothedY;
        } else {
            y = lr.smooth(newX, newY);
        }
        if (y.length > 100) {
            if (extra) {
                for (int i = 0; i < y.length; i++) {
                    y[i] = y[i] - newY[i];
                }

                //remove outliers
                float meanDiff = mean(y);
                float stdDiff = (float) Math.sqrt(variance(y));
                ArrayList<Integer> outliersIdx = new ArrayList<>();

                for (int i = 0; i < y.length; i++) {
                    if (Math.abs(zscore((float) y[i], meanDiff, stdDiff)) > 2) {
                        outliersIdx.add(i);
                    }
                }
                y = removeIdx(y, outliersIdx);
                newX = removeIdx(newX, outliersIdx);
                newY = removeIdx(newY, outliersIdx);
            }

            double[] weights = new double[y.length];
            for (int i = 0; i < y.length; i++) {
                weights[i] = Math.abs(median(Arrays.copyOfRange(y, Math.max(i - y.length/100, 0),
                        Math.min(i + y.length/100, y.length))));
            }

            //redo loess with weights
            y = lr.smooth(newX, newY, weights); //weighting so ends become better fit
        }
        double[][] results = new double[3][];
        results[0] = y;
        results[1] = newX;
        results[2] = newY;
        return results;
    }

    private static double[] removeIdx(double[] array, ArrayList<Integer> idx) {
        if (!idx.isEmpty()) {
            ArrayList<Double> newX = new ArrayList<>();
            for (int i = 0; i < array.length; i++) {
                if (! idx.contains(i)) {
                    newX.add(array[i]);
                }
            }
            array = new double[newX.size()];
            for (int i = 0; i < newX.size(); i++) {
                array[i] = newX.get(i);
            }
        }
        return array;
    }

    public static double residualSumOfSquares(double[] x, double[] y) {
        double sum = 0;
        for (int i = 0; i < x.length; i++) {
            sum += Math.pow(x[i] - y[i], 2);
        }
        return sum;
    }

    public static float[] movingAverage(float[] array, int windowSize) {
        int bl = array.length;
        float[] newStats = new float[bl];

        for (int i = 0; i < bl; i++) {
            float[] bin = Arrays.copyOfRange(array, Math.max(i - windowSize, 0), i + windowSize + 1);
            newStats[i] = mean(bin);
        }
        return newStats;
    }

    //get means and standard devs and interquartile ranges
    public static float[][] characterizebins(ArrayList<Float>[] bins, float IQR) {
        float[][] binStats = new float[bins.length][3]; //index by expRT, then mean or standard deviation
        for (int i = 0; i < bins.length; i++) {
            if (bins[i].size() == 0) {
                continue;
            } else if (bins[i].size() == 1) {
                binStats[i][0] = bins[i].get(0);
                binStats[i][1] = Float.MAX_VALUE;
                binStats[i][2] = Float.MAX_VALUE;
                continue;
            }

            //mean
            float m = StatMethods.mean(bins[i]);

            //standard dev
            float sd = (float) Math.sqrt(StatMethods.variance(bins[i], m));
            if (sd == 0f) {
                sd = Float.MAX_VALUE;
            }

            //interquartile range
            float plusMinus = IQR / 200f;
            float percentileIncrement = 1f / (float) (bins[i].size() - 1);
            int startIdx = Math.round(plusMinus / percentileIncrement);
            int endIdx = Math.round((plusMinus + .50f) / percentileIncrement);
            Collections.sort(bins[i]);
            float iqr = bins[i].get(endIdx) - bins[i].get(startIdx);
            if (iqr == 0f) {
                iqr = Float.MAX_VALUE;
            }

            binStats[i][0] = m;
            binStats[i][1] = sd;
            binStats[i][2] = iqr;
        }
        return binStats;
    }

    public static double probability(float exp, float pred, EmpiricalDist[] bins) {
        //get right bin to search
        EmpiricalDist ed = bins[Math.round(exp)];
        //check probability at point
        try {
            double p = KernelDensity.computeDensity(ed, new NormalDist(),
                    KernelDensityGen.getBaseBandwidth(ed), new double[]{pred})[0];
            if (Double.isNaN(p)) {
                p = 0;
            }
            return p;
        } catch (Exception e) { //nothing in bin
            return 0;
        }
    }

    public static float probabilityWithUniformPrior(int unifPriorSize, float unifProb, int binSize, float empiricalProb) {
        float w1 = (float) unifPriorSize / (float) (unifPriorSize + binSize);
        float w2 = (float) binSize / (float) (unifPriorSize + binSize);

        return w1 * unifProb + w2 * empiricalProb;
    }

    public static float median(ArrayList<Float> alist) {
        Collections.sort(alist);
        float median;
        if (alist.size() % 2 == 0)
            median = (alist.get(alist.size() / 2 - 1) + alist.get(alist.size() / 2)) / 2;
        else
            median = alist.get(alist.size() / 2);
        return median;
    }

    public static double median(double[] array) {
        ArrayList<Double> alist = new ArrayList<>();
        for (double f : array) {
            alist.add(f);
        }

        Collections.sort(alist);
        double median;
        if (alist.size() % 2 == 0)
            median = (alist.get(alist.size() / 2 - 1) + alist.get(alist.size() / 2)) / 2;
        else
            median = alist.get(alist.size() / 2);
        return median;
    }

    public static double medianDouble(ArrayList<Double> alist) {
        Collections.sort(alist);
        double median;
        if (alist.size() % 2 == 0)
            median = (alist.get(alist.size() / 2 - 1) + alist.get(alist.size() / 2)) / 2;
        else
            median = alist.get(alist.size() / 2);
        return median;
    }

    public static int consecutiveMatches(float[] array) {
        int maxConsecutive = 0;
        int currentScore = 0;

        for (float f : array) {
            if (f > 0f) {
                currentScore += 1;
            } else {
                if (currentScore > maxConsecutive) {
                    maxConsecutive = currentScore;
                }
                currentScore = 0;
            }
        }
        if (currentScore > maxConsecutive) {
            maxConsecutive = currentScore;
        }

        return maxConsecutive;
    }

    public static double meanSquaredError(double[] a, double[] b) {
        double mse = 0f;
        for(int i = 0; i < a.length; i++) {
            mse += Math.pow(a[i] - b[i], 2);
        }
        return mse / a.length;
    }

    public static ArrayList<double[][][]> trainTestSplit(double[][] expAndPred) {
        //arraylist of N splits
        //train-test,exp-pred,data
        double[] exp = expAndPred[0];
        double[] pred = expAndPred[1];

        ArrayList<double[][][]> splits = new ArrayList<>();

        for (int rep = 0; rep < Constants.regressionSplits; rep++) {
            //sort into train and test
            List<Double> trainExpList = new ArrayList<>();
            List<Double> trainPredList = new ArrayList<>();
            List<Double> testExpList = new ArrayList<>();
            List<Double> testPredList = new ArrayList<>();

            for (int i = 0; i < exp.length; i++) {
                if (i % Constants.regressionSplits != rep) {
                    trainExpList.add(exp[i]);
                    trainPredList.add(pred[i]);
                } else {
                    testExpList.add(exp[i]);
                    testPredList.add(pred[i]);
                }
            }

            double[][][] split = new double[2][2][];
            split[0][0] = trainExpList.stream().mapToDouble(Double::doubleValue).toArray();
            split[0][1] = trainPredList.stream().mapToDouble(Double::doubleValue).toArray();
            split[1][0] = testExpList.stream().mapToDouble(Double::doubleValue).toArray();
            split[1][1] = testPredList.stream().mapToDouble(Double::doubleValue).toArray();
            splits.add(split);
        }
        return splits;
    }
}