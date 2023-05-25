package Features;

public class SimilarityMethods {
    float[] vector1;
    float[] vector2;

    public SimilarityMethods(float[] vector1, float[] vector2) {
        this.vector1 = oneNormalize(vector1);
        this.vector2 = oneNormalize(vector2);
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
        float[] SabVector = new float[vector1.length];
        int numFrags = 0;
        for (float j : vector2) {
            if (j != 0) {
                numFrags += 1;
            }
        }

        if (numFrags < 2) {
            return 0;
        } else {
            for (int i = 0; i < SabVector.length; i++) {
                SabVector[i] = (vector1[i] + vector2[i]) / 2;
            }
        }

        return 1 - ( ((2 * spectralEntropy(SabVector)) - spectralEntropy(vector2) - spectralEntropy(vector1)) / Math.log(4));
    }
}
