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

import kotlin.jvm.functions.Function1;
import org.knowm.xchart.*;
import org.knowm.xchart.style.Styler;

import java.awt.*;
import java.io.File;
import java.io.IOException;
import java.util.List;
import java.util.*;

public class CalibrationFigure {
    String folderString;
    String mode;
    String charge;
    public CalibrationFigure() {}

    //TODO: repeat but with single entries at a time, with different outFile names
    public void plotFigure(MzmlReader mzml, String outFile, float opacity,
                           HashMap<String, double[][]> massToData,
                           HashMap<String, Function1<Double, Double>> loessFunctions) throws IOException {
        String pinPath = new File(outFile).getParent();
        String pinName = new File(outFile).getName();
        if (mode.equals("IM")) {
            pinName = "charge" + charge + "_" + pinName;
        }

        String dir = pinPath + File.separator + "MSBooster_plots";
        if (! new File(dir).exists()) {
            new File(dir).mkdirs();
        }
        if (! new File(dir + File.separator + folderString).exists()) {
            new File(dir + File.separator + folderString).mkdirs();
        }

        XYChart chart = new XYChartBuilder().width(1000).height(1000).build();
        chart.getStyler().setDefaultSeriesRenderStyle(XYSeries.XYSeriesRenderStyle.Scatter);
        chart.getStyler().setChartTitleVisible(true);
        chart.setTitle(mzml.pathStr);
        chart.setXAxisTitle("experimental " + mode);
        chart.setYAxisTitle("predicted " + mode);
        chart.getStyler().setLegendVisible(false);
        chart.getStyler().setMarkerSize(8);
        chart.getStyler().setYAxisGroupPosition(0, Styler.YAxisPosition.Right);
        chart.getStyler().setLegendVisible(true);
        chart.getStyler().setLegendPosition(Styler.LegendPosition.InsideNW);

        // Series
        List<Float> xData = new ArrayList<Float>();
        List<Float> yData = new ArrayList<Float>();

        //for PTMs besides oxM and C57
        List<List<Float>> xDataMod = new ArrayList<>();
        List<List<Float>> yDataMod = new ArrayList<>();

        ArrayList<String> massesList = new ArrayList<>();
        massesList.addAll(Arrays.asList(Constants.massesForLoessCalibration.split(",")));
        massesList.remove("");
        for (int i = 0; i < massesList.size(); i++) {
            xDataMod.add(new ArrayList<>());
            yDataMod.add(new ArrayList<>());
        }

        //set y lim
        int modIdx = 0;
        float minVal = Float.MAX_VALUE;
        float maxVal = Float.MIN_VALUE;
        for (Map.Entry<String, double[][]> entry : massToData.entrySet()) {
            if (entry.getKey().isEmpty() || entry.getKey().equals("others")) {
                for (double d : entry.getValue()[0]) {
                    xData.add((float) d);
                    if (d < minVal) {
                        minVal = (float) d;
                    }
                    if (d > maxVal) {
                        maxVal = (float) d;
                    }
                }
                for (double d : entry.getValue()[1]) {
                    yData.add((float) d);
                }
            } else {
                try { //might be empty if no regression
                    for (double d : entry.getValue()[0]) {
                        xDataMod.get(modIdx).add((float) d);
                        if (d < minVal) {
                            minVal = (float) d;
                        }
                        if (d > maxVal) {
                            maxVal = (float) d;
                        }
                    }
                    for (double d : entry.getValue()[1]) {
                        yDataMod.get(modIdx).add((float) d);
                    }
                } catch (Exception ignored) {}
                modIdx++;
            }
        }

        double ymax = 0;
        for (float f : yData) {
            if (f > ymax) {
                ymax = f;
            }
        }
        for (List<Float> listf : yDataMod) {
            for (float f : listf) {
                if (f > ymax) {
                    ymax = f;
                }
            }
        }
        chart.getStyler().setYAxisMax(ymax);

        if (!xData.isEmpty()) {
            XYSeries series = chart.addSeries("scatter", xData, yData);
            series.setMarkerColor(new Color(0, 0, 0, opacity));
        }

        for (int i = 0; i < xDataMod.size(); i++) {
            List<Float> ix = xDataMod.get(i);
            List<Float> iy = yDataMod.get(i);
            if (!ix.isEmpty()) {
                XYSeries seriesMod = chart.addSeries("scatterMods - " + massesList.get(i), ix, iy);
                seriesMod.setMarkerColor(new Color(65 * (i + 1) % 255, 105 * (i + 1) % 255, 225 * (i + 1) % 255));
            }
        }

        //loess regression
        // generates Log data
        int j = 1;
        for (String mass : loessFunctions.keySet()) {
            if (loessFunctions.get(mass) == null) {
                continue;
            }

            List<Float> x1Data = new ArrayList<Float>();
            List<Double> y1Data = new ArrayList<Double>();
            for (float i = minVal; i < maxVal; i = i + (maxVal - minVal) / 1000f) {
                x1Data.add(i);
                double y = loessFunctions.get(mass).invoke((double) i);
                y1Data.add(y);
            }
            XYSeries series1;
            if (mass.isEmpty()) {
                series1 = chart.addSeries("regression", x1Data, y1Data);
            } else {
                series1 = chart.addSeries("regression - " + mass, x1Data, y1Data);
            }
            series1.setMarkerColor(new Color(243 * j % 255, 9 * j % 255, 9 * j % 255));
            j += 2;
        }

        if (Constants.plotExtension.equalsIgnoreCase("png")) {
            BitmapEncoder.saveBitmap(chart, dir + File.separator + folderString +
                            File.separator + pinName.substring(0, pinName.length() - 4),
                    BitmapEncoder.BitmapFormat.PNG);
        } else if (Constants.plotExtension.equalsIgnoreCase("pdf")) {
            VectorGraphicsEncoder.saveVectorGraphic(chart,
                    dir + File.separator + folderString +
                            File.separator + pinName.substring(0, pinName.length() - 4),
                    VectorGraphicsEncoder.VectorGraphicsFormat.PDF);
        }
    }
}

