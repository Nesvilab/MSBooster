package Features;

import com.github.sanity.pav.PairAdjacentViolators;
import com.github.sanity.pav.Point;
import kotlin.jvm.functions.Function1;
import org.apache.commons.math3.analysis.interpolation.LoessInterpolator;
import smile.stat.distribution.KernelDensity;

import java.util.*;

public class StatMethods {
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

    public static Function1<Double, Double> LOESS(double[][] bins, double bandwidth, int robustIters) {
        //solve monotonicity issue
        double compare = -1;
        for (int i = 0; i < bins[0].length; i++) {
            double d = bins[0][i];
            if (d == compare) {
                bins[0][i] = bins[0][i - 1] + 0.00000001; //arbitrary increment to handle smooth method
            } else {
                compare = d;
            }
        }

        //fit loess
        //may not need if sanity version uses spline
        LoessInterpolator loessInterpolator = new LoessInterpolator(bandwidth, robustIters);
        double[] y = loessInterpolator.smooth(bins[0], bins[1]);

        //isotonic regression
        List<Point> points = new LinkedList<>();
        for (int i = 0; i < y.length; i++) {
            points.add(new Point(bins[0][i], y[i]));
        }
        PairAdjacentViolators pav = new PairAdjacentViolators(points);

        return pav.interpolator(); //extrapolationStrategy can use flat or tangent (default tangent)
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

    public static double probability(float exp, float pred, KernelDensity[] bins) {
        //get right bin to search
        KernelDensity kd = bins[Math.round(exp)];

        //check probability at point
        try {
            return kd.p(pred);
        } catch (Exception e) { //nothing in bin
            return 0;
        }
    }

    public static float probabilityWithUniformPrior(int unifPriorSize, float unifProb,
                                                      int binSize, float empiricalProb) {
        float w1 = (float) unifPriorSize / (float) (unifPriorSize + binSize);
        float w2 = (float) binSize / (float) (unifPriorSize + binSize);

        return w1 * unifProb + w2 * empiricalProb;
    }

    public static void main(String[] args) {
        System.out.println(variance(new double[] {1.0, 1.0}));
        System.out.println( 1f / Float.POSITIVE_INFINITY);
    }
}
