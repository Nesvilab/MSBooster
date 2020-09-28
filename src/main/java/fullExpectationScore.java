import java.util.Arrays;

public class fullExpectationScore {
    double[] scoreList;
    double[] topScoreList;
    double[] sx;
    double beta0;
    double beta1;

    public double fullExpectationScore(double[] scoreList, double PSMscore) {
        //sorting score list
        this.scoreList = scoreList;
        Arrays.sort(this.scoreList);

        //if less than min
        double minScore = scoreList[0];
        if (PSMscore <= minScore) {
            return (double) scoreList.length;
        }

        //model with highest scoring part of score list, so extract those scores and log them
        this.topScoreList = this.topScoringPortion();

        //set log survival function
        this.sx = this.survivalFunction();

        //do linear regression and set y intercept and slope
        double[] betas = this.linearRegression();
        this.beta0 = betas[0];
        this.beta1 = betas[1];

        //normalize PSMscore
        PSMscore += (1.0 - scoreList[0]);

        //estimate survival function with y = mx + b
        //multiply by n to get expectation score
        return Math.min(Math.exp(beta0 + beta1 * PSMscore) * topScoreList.length, scoreList.length);
    }

    private double[] topScoringPortion() {
        int topN = (int) Math.ceil(scoreList.length * constants.highScoringProp);
        double[] myList = Arrays.copyOfRange(scoreList, topN, scoreList.length + 1);

        double minScore = scoreList[0];
        for (int i = 0; i < myList.length; i++) {
            //avoid 0s when taking log
            //normalize so minimum score is constants.addToZero
            myList[i] += (1.0 - minScore);
            myList[i] = Math.log(myList[i]);
        }

        return myList;
    }

    private double[] survivalFunction() {
        int tsl_length = topScoreList.length;
        double doubleListLength = (double) scoreList.length;

        double[] sx = new double[tsl_length];
        for (int i = 0; i < tsl_length; i++) {
            double num = (double) (tsl_length - i);
            sx[i] = Math.log(num / doubleListLength);
        }

        return sx;
    }

    private double[] linearRegression() {
        //returns beta 0 and beta 1
        double meanX = mean(topScoreList);
        double meanY = mean(sx);

        double varX = variance(topScoreList, meanX);

        double covar = variance(topScoreList, meanX, sx, meanY);

        double beta1 = covar / varX;
        double beta0 = meanY - (beta1 * meanX);

        return new double[] {beta0, beta1};
    }

    private static double mean(double[] vector) {
        double meanX = 0;
        for (double x : vector) {
            meanX += x;
        }
        return meanX / (double) vector.length;
    }

    private static double variance(double[] vector) {
        int vecLength = vector.length;
        double mean = mean(vector);
        double var = 0;

        for (double v : vector) {
            var += Math.pow(v - mean, 2);
        }

        return var / (vecLength - 1);
    }

    private static double variance(double[] vector, double mean) {
        int vecLength = vector.length;
        double var = 0;

        for (double v : vector) {
            var += Math.pow(v - mean, 2);
        }

        return var / (vecLength - 1);
    }

    //covariance
    private static double variance(double[] vectorX, double[] vectorY) {
        int vecLength = vectorX.length;
        double meanX = mean(vectorX);
        double meanY = mean(vectorY);

        double var = 0;

        for (int i = 0; i < vecLength; i++) {
            var += (vectorX[i] - meanX) * (vectorY[i] - meanY);
        }

        return var / (vecLength - 1);
    }

    private static double variance(double[] vectorX, double meanX, double[] vectorY, double meanY) {
        int vecLength = vectorX.length;

        double var = 0;

        for (int i = 0; i < vecLength; i++) {
            var += (vectorX[i] - meanX) * (vectorY[i] - meanY);
        }

        return var / (vecLength - 1);
    }

}
