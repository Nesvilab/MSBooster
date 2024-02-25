package Features;

import org.knowm.xchart.BitmapEncoder;
import org.knowm.xchart.CategoryChart;
import org.knowm.xchart.CategoryChartBuilder;
import org.knowm.xchart.Histogram;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.Set;

public class ScoreHistogram {

    Set<String> logScaleFeatures = Set.of("delta_RT_loess", "delta_RT_loess_normalized", "RT_probability_unif_prior",
            "hypergeometric_probability");

    public ScoreHistogram(PinReader pinReader, ArrayList<String> features) throws IOException {
        HashMap<String, ArrayList<Double>> targetScores = new HashMap<>();
        HashMap<String, ArrayList<Double>> decoyScores = new HashMap<>();
        HashMap<String, Double> scoreMax = new HashMap<>();
        HashMap<String, Double> scoreMin = new HashMap<>();

        for (int i = 0; i < features.size(); i++) {
            String feature = features.get(i);
            feature = Constants.camelToUnderscore.get(feature);
            features.set(i, feature);
            targetScores.put(feature, new ArrayList<>());
            decoyScores.put(feature, new ArrayList<>());
            scoreMax.put(feature, -1 * Double.MAX_VALUE);
            scoreMin.put(feature, Double.MAX_VALUE);
        }

        //get scores
        while (pinReader.next(true)) {
            for (String feature : features) {
                Double score = Double.valueOf(pinReader.getColumn(feature));
                if (logScaleFeatures.contains(feature)) {
                    score = Math.log10(score + 0.01);
                }
                if (pinReader.getTD() == 1) {
                    targetScores.get(feature).add(score);
                } else {
                    decoyScores.get(feature).add(score);
                }
                if (score > scoreMax.get(feature)) {
                    scoreMax.put(feature, score);
                }
                if (score < scoreMin.get(feature)) {
                    scoreMin.put(feature, score);
                }
            }
        }
        pinReader.reset();

        for (String feature : features) {
            String xAxisLabel = feature;
            if (logScaleFeatures.contains(feature)) {
                xAxisLabel = "log(" + xAxisLabel + " + 0.01)";
            }
            CategoryChart chart = new CategoryChartBuilder().width(1200).height(800).
                    xAxisTitle(xAxisLabel).yAxisTitle("PSMs").build();
            chart.getStyler().setPlotGridLinesVisible(false);
            //chart.getStyler().setXAxisTickMarkSpacingHint(10);
            chart.getStyler().setXAxisMaxLabelCount(10);
            chart.getStyler().setYAxisTickMarkSpacingHint(80);
            chart.getStyler().setXAxisLabelRotation(45);
            chart.getStyler().setXAxisDecimalPattern("#.##");

            //plot histogram
            try {
                Histogram histT = new Histogram(targetScores.get(feature), 100,
                        scoreMin.get(feature), scoreMax.get(feature));
                chart.addSeries("targets", histT.getxAxisData(), histT.getyAxisData());
            } catch (Exception e) {e.printStackTrace();}
            try {
                Histogram histD = new Histogram(decoyScores.get(feature), 100,
                        scoreMin.get(feature), scoreMax.get(feature));
                chart.addSeries("decoys", histD.getxAxisData(), histD.getyAxisData());
            } catch (Exception e) {e.printStackTrace();}

            String pinPath = new File(pinReader.name).getParent();
            String name = new File(pinReader.name).getName();
            String dir = pinPath + File.separator + "MSBooster_plots";
            if (!new File(dir).exists()) {
                new File(dir).mkdirs();
            }
            if (!new File(dir + File.separator + "score_histograms").exists()) {
                new File(dir + File.separator + "score_histograms").mkdirs();
            }
            BitmapEncoder.saveBitmap(chart, dir + File.separator + "score_histograms" + File.separator +
                            name.substring(0, name.length() - 4) + "_" + feature,
                    BitmapEncoder.BitmapFormat.PNG);
        }
    }
}
