package Features;

import java.util.ArrayList;

public class StatMethods {
    public static float mean(double[] vector) {
        float meanX = 0;
        for (double x : vector) {
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

    public static void main(String[] args) {
        System.out.println(variance(new double[] {1.0, 1.0}));
        System.out.println( 1f / Float.POSITIVE_INFINITY);
    }
}
