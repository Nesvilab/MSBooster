package figures;

import org.knowm.xchart.*;
import utils.StatMethods;

import java.io.File;
import java.io.IOException;
import java.util.*;

import static allconstants.Constants.figureDirectory;
import static figures.ExtensionPlotter.plot;

public class CumulativeLinePlots {
    public CumulativeLinePlots(HashMap<String, List<Double>> scoreDescriptors, boolean descending,
                               String scoreType, String scoreName)
            throws IOException {
        class MzmlScores implements Comparable<MzmlScores>{
            final int size;
            final double median;
            final double lowerQuartile;
            final double upperQuartile;
            final double lower5th;
            final double upper95th;
            MzmlScores(List<Double> scores) {
                size = scores.size();
                Collections.sort(scores);
                median = StatMethods.median(scores);

                float quartile = (float) size / 4f;
                lowerQuartile = scores.get((int) Math.floor(quartile));
                upperQuartile = scores.get((int) Math.floor(3f * quartile));

                float fifths = (float) size / 20f;
                lower5th = scores.get((int) Math.floor(fifths));
                upper95th = scores.get((int) Math.floor(19f * fifths));
            }

            @Override
            public int compareTo(MzmlScores other) {
                return Double.compare(this.median, other.median);
            }
        }

        //get scores
        List<Map.Entry<String, MzmlScores>> sortedEntries = new ArrayList<>();
        for (Map.Entry<String, List<Double>> entry : scoreDescriptors.entrySet()) {
            sortedEntries.add(new AbstractMap.SimpleEntry<>(
                    entry.getKey(), new MzmlScores(entry.getValue())));
        }

        // Sort by median
        if (descending) {
            sortedEntries.sort(Comparator.comparingDouble((Map.Entry<String, MzmlScores> e) ->
                    e.getValue().median).reversed());
        } else {
            sortedEntries.sort(Comparator.comparingDouble((Map.Entry<String, MzmlScores> e) ->
                    e.getValue().median));
        }

        ArrayList<Double> medians = new ArrayList<>();
        ArrayList<Double> lowerQuartiles = new ArrayList<>();
        ArrayList<Double> upperQuartiles = new ArrayList<>();
        ArrayList<Double> lowerFifths = new ArrayList<>();
        ArrayList<Double> upperNinetyfifths = new ArrayList<>();
        ArrayList<String> fileNames = new ArrayList<>();
        for (Map.Entry<String, MzmlScores> entry : sortedEntries) {
            medians.add(entry.getValue().median);
            lowerQuartiles.add(entry.getValue().lowerQuartile);
            upperQuartiles.add(entry.getValue().upperQuartile);
            lowerFifths.add(entry.getValue().lower5th);
            upperNinetyfifths.add(entry.getValue().upper95th);
            fileNames.add(entry.getKey() + "(n=" + entry.getValue().size + ")");
        }

        //plot
        CategoryChart chart = new CategoryChartBuilder()
                .width(1500)
                .height(1000)
                .title("Cumulative " + scoreType + " Similarity Score QC")
                .xAxisTitle("mzML files")
                .yAxisTitle(scoreName)
                .build();

        chart.getStyler().setDefaultSeriesRenderStyle(CategorySeries.CategorySeriesRenderStyle.Line);
        chart.getStyler().setMarkerSize(6);
        chart.getStyler().setLegendVisible(true);
        chart.getStyler().setXAxisLabelRotation(90);
        chart.getStyler().setOverlapped(true);
        chart.addSeries("5th percentile", fileNames, lowerFifths);
        chart.addSeries("25th percentile", fileNames, lowerQuartiles);
        chart.addSeries("median", fileNames, medians);
        chart.addSeries("75th percentile", fileNames, upperQuartiles);
        chart.addSeries("95th percentile", fileNames, upperNinetyfifths);

        plot(chart, figureDirectory + File.separator + "cumulativeQC" +
                File.separator + "cumulative_" + scoreType + "_lineplot");
    }
}
