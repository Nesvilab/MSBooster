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

import org.knowm.xchart.*;

import java.io.File;
import java.io.IOException;
import java.util.*;

public class ScoreHistogram {

    Set<String> logScaleFeatures = Set.of("delta_RT_loess", "delta_RT_loess_normalized", "RT_probability_unif_prior",
            "hypergeometric_probability");
    Set<String> logYScaleFeatures = Set.of("peptide_counts");
    Set<String> intScores = Set.of("peptide_counts"); //TODO: gets crowded for many bins, can consider binning by 10s

    public ScoreHistogram(PinReader pinReader, ArrayList<String> fs) throws IOException {
        HashMap<String, ArrayList<Double>> targetScores = new HashMap<>();
        HashMap<String, ArrayList<Double>> decoyScores = new HashMap<>();
        HashMap<String, Double> scoreMax = new HashMap<>();
        HashMap<String, Double> scoreMin = new HashMap<>();

        ArrayList<String> features = new ArrayList<>(fs.size());
        for (int i = 0; i < fs.size(); i++) {
            String feature = fs.get(i);
            feature = Constants.camelToUnderscore.get(feature);
            features.add(feature);
            targetScores.put(feature, new ArrayList<>());
            decoyScores.put(feature, new ArrayList<>());
            scoreMax.put(feature, -1 * Double.MAX_VALUE);
            scoreMin.put(feature, Double.MAX_VALUE);
        }

        //get scores
        while (pinReader.next(true)) {
            for (String feature : features) {
                double score = Double.parseDouble(pinReader.getColumn(feature));
                if (logScaleFeatures.contains(feature)) {
                    score = Math.log10(score + 0.01);
                }
                if (intScores.contains(feature)) {
                    if (pinReader.getTD() == 1) {
                        targetScores.get(feature).add(score + 0.5);
                    } else {
                        decoyScores.get(feature).add(score + 0.5);
                    }
                } else {
                    if (pinReader.getTD() == 1) {
                        targetScores.get(feature).add(score);
                    } else {
                        decoyScores.get(feature).add(score);
                    }
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
            String yAxisLabel = "PSMs";
            if (logYScaleFeatures.contains(feature)) {
                yAxisLabel = "log(PSMs + 1)";
                for (Double j = scoreMin.get(feature); j < scoreMax.get(feature) + 1; j++) {
                    targetScores.get(feature).add(j + 0.5);
                    decoyScores.get(feature).add(j + 0.5);
                }
            }
            CategoryChart chart = new CategoryChartBuilder().width(1200).height(800).
                    xAxisTitle(xAxisLabel).yAxisTitle(yAxisLabel).build();
            chart.getStyler().setPlotGridLinesVisible(false);
            //chart.getStyler().setXAxisTickMarkSpacingHint(10);
            chart.getStyler().setXAxisMaxLabelCount(10);
            chart.getStyler().setYAxisTickMarkSpacingHint(80);
            chart.getStyler().setXAxisLabelRotation(45);
            chart.getStyler().setXAxisDecimalPattern("#.##");

            //plot histogram
            int numBins = 100;
            if (intScores.contains(feature)) {
                numBins = (int) (scoreMax.get(feature) - scoreMin.get(feature) + 1);
            }
            try {
                Histogram histT;
                if (intScores.contains(feature)) {
                    histT = new Histogram(targetScores.get(feature), numBins,
                            scoreMin.get(feature), scoreMax.get(feature) + 1);
                } else {
                    histT = new Histogram(targetScores.get(feature), numBins,
                            scoreMin.get(feature), scoreMax.get(feature));
                }
                List<Double> xData = histT.getxAxisData();
                List<Double> yData = histT.getyAxisData();
                if (logYScaleFeatures.contains(feature)) {
                    xData.replaceAll(aDouble -> aDouble - 0.5);
                    yData.replaceAll(Math::log10);
                }
                chart.addSeries("targets", xData, yData);
            } catch (Exception e) {e.printStackTrace();}
            try {
                Histogram histD;
                if (intScores.contains(feature)) {
                    histD = new Histogram(decoyScores.get(feature), numBins,
                            scoreMin.get(feature), scoreMax.get(feature) + 1);
                } else {
                    histD = new Histogram(decoyScores.get(feature), numBins,
                            scoreMin.get(feature), scoreMax.get(feature));
                }
                List<Double> xData = histD.getxAxisData();
                List<Double> yData = histD.getyAxisData();
                if (logYScaleFeatures.contains(feature)) {
                    xData.replaceAll(aDouble -> aDouble - 0.5);
                    yData.replaceAll(Math::log10);
                }
                chart.addSeries("decoys", xData, yData);
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
