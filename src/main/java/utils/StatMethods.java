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

package utils;

import umontreal.ssj.gof.KernelDensity;
import umontreal.ssj.probdist.EmpiricalDist;
import umontreal.ssj.probdist.NormalDist;
import umontreal.ssj.randvar.KernelDensityGen;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.Map;
import java.util.TreeMap;


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

        return medianDouble(alist);
    }

    public static float median(float[] array) {
        ArrayList<Float> alist = new ArrayList<>();
        for (float f : array) {
            alist.add(f);
        }

        return median(alist);
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
        for (int i = 0; i < a.length; i++) {
            mse += Math.pow(a[i] - b[i], 2);
        }
        return mse / a.length;
    }

    public static float meanSquaredError(float[] a) {
        double mse = 0f;
        for (float value : a) {
            mse += Math.pow(value, 2);
        }
        return (float) (mse / a.length);
    }
}