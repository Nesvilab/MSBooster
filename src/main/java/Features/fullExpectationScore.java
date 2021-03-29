//package Features;
//
//import Features.constants;
//
//import java.util.Arrays;
//
//public class fullExpectationScore {
//    double[] scoreList;
//    double[] topScoreList;
//    double[] sx;
//    double beta0;
//    double beta1;
//
//    public double fullExpectationScore(double[] scoreList, double PSMscore) {
//        //sorting score list
//        this.scoreList = scoreList;
//        Arrays.sort(this.scoreList);
//
//        //if less than min
//        double minScore = scoreList[0];
//        if (PSMscore <= minScore) {
//            return (double) scoreList.length;
//        }
//
//        //model with highest scoring part of score list, so extract those scores and log them
//        this.topScoreList = this.topScoringPortion();
//
//        //set log survival function
//        this.sx = this.survivalFunction();
//
//        //do linear regression and set y intercept and slope
//        double[] betas = this.linearRegression();
//        this.beta0 = betas[0];
//        this.beta1 = betas[1];
//
//        //normalize PSMscore
//        PSMscore += (1.0 - scoreList[0]);
//
//        //estimate survival function with y = mx + b
//        //multiply by n to get expectation score
//        return Math.min(Math.exp(beta0 + beta1 * PSMscore) * topScoreList.length, scoreList.length);
//    }
//
//    private double[] topScoringPortion() {
//        int topN = (int) Math.ceil(scoreList.length * constants.highScoringProp);
//        double[] myList = Arrays.copyOfRange(scoreList, topN, scoreList.length + 1);
//
//        double minScore = scoreList[0];
//        for (int i = 0; i < myList.length; i++) {
//            //avoid 0s when taking log
//            //normalize so minimum score is constants.addToZero
//            myList[i] += (1.0 - minScore);
//            myList[i] = Math.log(myList[i]);
//        }
//
//        return myList;
//    }
//
//    private double[] survivalFunction() {
//        int tsl_length = topScoreList.length;
//        double doubleListLength = (double) scoreList.length;
//
//        double[] sx = new double[tsl_length];
//        for (int i = 0; i < tsl_length; i++) {
//            double num = (double) (tsl_length - i);
//            sx[i] = Math.log(num / doubleListLength);
//        }
//
//        return sx;
//    }
//
//    private double[] linearRegression() {
//        //returns beta 0 and beta 1
//        double meanX = StatMethods.mean(topScoreList);
//        double meanY = StatMethods.mean(sx);
//
//        double varX = StatMethods.variance(topScoreList, meanX);
//
//        double covar = StatMethods.variance(topScoreList, meanX, sx, meanY);
//
//        double beta1 = covar / varX;
//        double beta0 = meanY - (beta1 * meanX);
//
//        return new double[] {beta0, beta1};
//    }
//
//}
