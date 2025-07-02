package figures;

import org.knowm.xchart.*;
import org.knowm.xchart.style.lines.SeriesLines;
import org.knowm.xchart.style.markers.SeriesMarkers;
import utils.StatMethods;

import java.awt.*;
import java.io.File;
import java.io.IOException;
import java.util.*;

import static allconstants.Constants.figureDirectory;
import static figures.ExtensionPlotter.plot;

public class CumulativeAbsDeltaLinePlots {
    static int binSize = 20;
    static double windowDivider = 60d;
    public CumulativeAbsDeltaLinePlots(HashMap<String, HashMap<String, double[][]>> absDeltas) throws IOException {
        class deltaRT implements Comparable<deltaRT>{
            final double exp;
            final double delta;

            deltaRT(double exp, double delta) {
                this.exp = exp;
                this.delta = delta;
            }

            @Override
            public int compareTo(deltaRT other) {
                return Double.compare(this.exp, other.exp);
            }
        }

        for (Map.Entry<String, HashMap<String, double[][]>> entryMass : absDeltas.entrySet()) {
            String mass = entryMass.getKey();
            HashMap<String, double[][]> entries = entryMass.getValue();
            double minX = Double.MAX_VALUE, maxX = -1 * Double.MAX_VALUE;

            //for each mzml, save medians, 5th percentiles, and 95th percentiles for each bin
            TreeMap<String, double[][]> stats = new TreeMap<>(); //median, 5th, 95th
            for (Map.Entry<String, double[][]> entry : entries.entrySet()) {
                deltaRT[] deltaRTs = new deltaRT[entry.getValue()[0].length];
                for (int i = 0; i < deltaRTs.length; i++) {
                    deltaRTs[i] = new deltaRT(entry.getValue()[0][i], entry.getValue()[1][i]);
                }
                Arrays.sort(deltaRTs);
                if (deltaRTs[0].exp < minX) {
                    minX = deltaRTs[0].exp;
                }
                if (deltaRTs[deltaRTs.length - 1].exp > maxX) {
                    maxX = deltaRTs[deltaRTs.length - 1].exp;
                }

                //we want small bins
                int numBins = deltaRTs.length / binSize;
                if (numBins < 2) {
                    numBins = deltaRTs.length;
                }

                ArrayList<deltaRT>[] bins = new ArrayList[numBins];
                for (int i = 0; i < numBins; i++) {
                    bins[i] = new ArrayList<>();
                }

                for (int i = 0; i < deltaRTs.length; i++) {
                    bins[(int) Math.floor((double) i / (double) deltaRTs.length * (double) numBins)].add(deltaRTs[i]);
                }

                //get stats for each bin
                double[] medians = new double[numBins];
                double[] lower = new double[numBins];
                double[] upper = new double[numBins];
                double[] avgExp = new double[numBins];

                for (int i = 0; i < numBins; i++) {
                    ArrayList<deltaRT> bin = bins[i];
                    double[] values = new double[bin.size()];
                    double[] exps = new double[bin.size()];
                    for (int j = 0; j < values.length; j++) {
                        values[j] = bin.get(j).delta;
                        exps[j] = bin.get(j).exp;
                    }

                    Arrays.sort(values);
                    medians[i] = StatMethods.median(values);
                    int twentieth = values.length / 20;
                    lower[i] = values[twentieth];
                    upper[i] = values[19 * twentieth];
                    avgExp[i] = StatMethods.mean(exps);
                }

                //sliding window average
                medians = StatMethods.movingAverage(medians, (int) Math.ceil((double) medians.length / windowDivider));
                lower = StatMethods.movingAverage(lower, (int) Math.ceil((double) lower.length / windowDivider));
                upper = StatMethods.movingAverage(upper, (int) Math.ceil((double) upper.length / windowDivider));

                //save
                stats.put(entry.getKey(), new double[][]{avgExp, medians, lower, upper});
            }

            //plot
            String title = "Cumulative Absolute Delta RT Trends";
            String baseName = figureDirectory + File.separator + "cumulativeQC" + File.separator + "cumulative_absDeltaRT_trend";
            if (!mass.isEmpty()) {
                title = title + " -- " + mass.split("&")[0];
                baseName = baseName + "_" + mass.split("&")[0];
            }
            XYChart chart = new XYChartBuilder()
                    .width(1500)
                    .height(1000)
                    .title(title)
                    .xAxisTitle("Experimental RT")
                    .yAxisTitle("Calibrated RT - Experimental RT")
                    .build();

            chart.getStyler().setMarkerSize(6);
            chart.getStyler().setLegendVisible(true);

            //set line at 0
            double[] xRange = new double[] {minX, maxX};
            double[] yZero = new double[] {0.0, 0.0};

            XYSeries zeroLine = chart.addSeries("delta RT = 0", xRange, yZero);
            zeroLine.setLineColor(Color.BLACK);
            zeroLine.setLineWidth(3.0f);  // Make it bold
            zeroLine.setMarker(SeriesMarkers.NONE);

            int colorIdx = 1;
            for (Map.Entry<String, double[][]> entry : stats.entrySet()) {
                BasicStroke bs = lineStyles[colorIdx % lineStyles.length];

                XYSeries p5 = chart.addSeries(entry.getKey() + " 5th percentile",
                        entry.getValue()[0], entry.getValue()[2]);
                p5.setMarker(SeriesMarkers.DIAMOND);
                p5.setLineColor(FigureUtils.getColorWithAlpha(colorIdx, 0.5f));
                p5.setMarkerColor(FigureUtils.getColorWithAlpha(colorIdx, 0.5f));
                p5.setLineStyle(bs);
                if (colorIdx != 1) {
                    p5.setShowInLegend(false);
                }

                XYSeries median = chart.addSeries(entry.getKey() + " median",
                        entry.getValue()[0], entry.getValue()[1]);
                median.setMarker(SeriesMarkers.CIRCLE);
                median.setLineColor(FigureUtils.getColor(colorIdx));
                median.setMarkerColor(FigureUtils.getColor(colorIdx));
                median.setLineStyle(bs);

                XYSeries p95 = chart.addSeries(entry.getKey() + " 95th percentile",
                        entry.getValue()[0], entry.getValue()[3]);
                p95.setMarker(SeriesMarkers.SQUARE);
                p95.setLineColor(FigureUtils.getColorWithAlpha(colorIdx, 0.5f));
                p95.setMarkerColor(FigureUtils.getColorWithAlpha(colorIdx, 0.5f));
                p95.setLineStyle(bs);
                if (colorIdx != 1) {
                    p95.setShowInLegend(false);
                }

                colorIdx++;
            }

            plot(chart, baseName);
        }
    }

    BasicStroke[] lineStyles = new BasicStroke[] {
            SeriesLines.SOLID,
            SeriesLines.DASH_DASH,
            SeriesLines.DASH_DOT,
            SeriesLines.DOT_DOT,
    };
}
