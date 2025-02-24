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

package figures;

import allconstants.Constants;
import org.knowm.xchart.CategoryChart;
import org.knowm.xchart.CategoryChartBuilder;
import org.knowm.xchart.Histogram;
import org.knowm.xchart.style.Styler;
import readers.datareaders.PinReader;

import java.awt.Color;
import java.awt.Font;
import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Set;

import static allconstants.Constants.figureDirectory;
import static figures.ExtensionPlotter.plot;

public class ScoreHistogram {
    HashSet<String> logScaleFeatures = new HashSet<>(Set.of(
            "delta_RT_loess", "delta_RT_loess_normalized", "RT_probability_unif_prior",
            "hypergeometric_probability", "delta_IM_loess"
    ));
    HashSet<String> intScores = new HashSet<>(Set.of("peptide_counts")); //TODO: gets crowded for many bins, can consider binning by 10s

    public ScoreHistogram(PinReader pinReader, ArrayList<String> fs) throws IOException {
        HashMap<String, ArrayList<Double>> targetScores = new HashMap<>();
        HashMap<String, ArrayList<Double>> decoyScores = new HashMap<>();
        HashMap<String, Double> scoreMax = new HashMap<>();
        HashMap<String, Double> scoreMin = new HashMap<>();

        ArrayList<String> features = new ArrayList<>(fs.size());
        ArrayList<Boolean> logScaleBooleans = new ArrayList<>(fs.size());

        for (String feature : fs) {
            String baseFeature = feature.split("\\^")[0];
            features.add(feature);
            if (logScaleFeatures.contains(baseFeature)) {
                logScaleBooleans.add(true);
            } else {
                logScaleBooleans.add(false);
            }
            targetScores.put(feature, new ArrayList<>());
            decoyScores.put(feature, new ArrayList<>());
            scoreMax.put(feature, -1 * Double.MAX_VALUE);
            scoreMin.put(feature, Double.MAX_VALUE);
        }

        //get scores
        while (pinReader.next(true)) {
            for (int i = 0; i < features.size(); i++) {
                String feature = features.get(i);
                boolean useLogScale = logScaleBooleans.get(i);

                double score = Double.parseDouble(pinReader.getColumn(feature));
                if (useLogScale) {
                    score = Math.log10(score + 0.01);
                }
                if (intScores.contains(feature)) {
                    if (pinReader.getTD() == 1) {
                        targetScores.get(feature).add(score - 0.5);
                    } else {
                        decoyScores.get(feature).add(score - 0.5);
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

        for (int i = 0; i < features.size(); i++) {
            String feature = features.get(i);
            String xAxisLabel = features.get(i);
            boolean useLogScale = logScaleBooleans.get(i);

            if (useLogScale) {
                xAxisLabel = "log(" + xAxisLabel + " + 0.01)";
            }
            String yAxisLabel = "PSMs";

            int width = 510;
            int height = 170;
            CategoryChart chart = null;
            if (Constants.plotExtension.equalsIgnoreCase("png")) {
                chart = new CategoryChartBuilder().width(width * 3).height(height * 3).
                        xAxisTitle(xAxisLabel).yAxisTitle(yAxisLabel).build();
            } else if (Constants.plotExtension.equalsIgnoreCase("pdf")) {
                chart = new CategoryChartBuilder().width(width).height(height).
                        xAxisTitle(xAxisLabel).yAxisTitle(yAxisLabel).build();
            }

            chart.getStyler().setPlotGridLinesVisible(false);
            //chart.getStyler().setXAxisTickMarkSpacingHint(10);
            chart.getStyler().setXAxisMaxLabelCount(10);
            //chart.getStyler().setYAxisTickMarkSpacingHint(80);
            chart.getStyler().setXAxisLabelRotation(45);
            chart.getStyler().setXAxisDecimalPattern("#.##");
            chart.getStyler().setLegendLayout(Styler.LegendLayout.Horizontal);
            chart.getStyler().setLegendPosition(Styler.LegendPosition.valueOf("OutsideS"));
            chart.getStyler().setLegendBorderColor(Color.white);
            chart.getStyler().setPlotBorderVisible(false);
            if (Constants.plotExtension.equalsIgnoreCase("png")) {
                chart.getStyler().setAxisTitleFont(new Font("Helvetica", Font.PLAIN, 21));
                chart.getStyler().setLegendFont(new Font("Helvetica", Font.PLAIN, 21));
                chart.getStyler().setAxisTickLabelsFont(new Font("Helvetica", Font.PLAIN, 15));
            } else if (Constants.plotExtension.equalsIgnoreCase("pdf")) {
                chart.getStyler().setAxisTitleFont(new Font("Helvetica", Font.PLAIN, 7));
                chart.getStyler().setLegendFont(new Font("Helvetica", Font.PLAIN, 7));
                chart.getStyler().setAxisTickLabelsFont(new Font("Helvetica", Font.PLAIN, 5));
            }
            chart.getStyler().setChartBackgroundColor(Color.white);

            //plot histogram
            int numBins = 100;
            if (intScores.contains(feature)) {
                numBins = (int) (scoreMax.get(feature) - scoreMin.get(feature) + 1);
            }
            try {
                Histogram histT;
                if (intScores.contains(feature)) {
                    histT = new Histogram(targetScores.get(feature), numBins,
                            scoreMin.get(feature) - 0.5, scoreMax.get(feature) + 1 - 0.5);
                } else {
                    histT = new Histogram(targetScores.get(feature), numBins,
                            scoreMin.get(feature), scoreMax.get(feature));
                }
                List<Double> xData = histT.getxAxisData();
                List<Double> yData = histT.getyAxisData();
                chart.addSeries("targets", xData, yData);
            } catch (Exception e) {e.printStackTrace();}
            try {
                Histogram histD;
                if (intScores.contains(feature)) {
                    histD = new Histogram(decoyScores.get(feature), numBins,
                            scoreMin.get(feature) - 0.5, scoreMax.get(feature) + 1 - 0.5);
                } else {
                    histD = new Histogram(decoyScores.get(feature), numBins,
                            scoreMin.get(feature), scoreMax.get(feature));
                }
                List<Double> xData = histD.getxAxisData();
                List<Double> yData = histD.getyAxisData();
                chart.addSeries("decoys", xData, yData);
            } catch (Exception e) {e.printStackTrace();}

            String name = new File(pinReader.name).getName();
            if (!new File(figureDirectory).exists()) {
                new File(figureDirectory).mkdirs();
            }
            if (!new File(figureDirectory + File.separator + "score_histograms").exists()) {
                new File(figureDirectory + File.separator + "score_histograms").mkdirs();
            }

            plot(chart, figureDirectory + File.separator + "score_histograms" + File.separator +
                    name.substring(0, name.length() - 4) + "_" + feature);
        }
        pinReader.close();
    }
}
