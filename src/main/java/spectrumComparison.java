public class spectrumComparison {
    double[] expMZs;
    double[] expIntensities;
    double[] predMZs;
    double[] predIntensities;
    double ppm;
    double[] matchedIntensities;

    public spectrumComparison(double[] eMZs, double[] eIntensities,
                              double[] pMZs, double[] pIntensities,
                              double ppmTolerance) {
        expMZs = eMZs;
        expIntensities = eIntensities;
        predMZs = pMZs;
        predIntensities = pIntensities;
        ppm =  ppmTolerance / 1000000;
        matchedIntensities = this.getMatchedIntensities();
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

    public void normalize() {
        //if we wish to normalize
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
        return num / den;
    }
}
