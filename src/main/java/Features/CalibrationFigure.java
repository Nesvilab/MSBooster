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
import org.knowm.xchart.style.markers.SeriesMarkers;

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
        chart.getStyler().setMarkerSize(8);
        chart.getStyler().setYAxisGroupPosition(0, Styler.YAxisPosition.Right);
        chart.getStyler().setLegendVisible(true);
        chart.getStyler().setLegendPosition(Styler.LegendPosition.InsideNW);

        // Series
        List<Float> xData = new ArrayList<Float>();
        List<Float> yData = new ArrayList<Float>();

        //for massOffsets or other masses specified
        HashMap<String, List<Float>> xDataMod = new HashMap<>();
        HashMap<String, List<Float>> yDataMod = new HashMap<>();

        ArrayList<String> massesList = new ArrayList<>(massToData.keySet());
        massesList.remove("");
        massesList.remove("others");
        for (String mass : massesList) {
            xDataMod.put(mass, new ArrayList<>());
            yDataMod.put(mass, new ArrayList<>());
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
                        xDataMod.get(entry.getKey()).add((float) d);
                        if (d < minVal) {
                            minVal = (float) d;
                        }
                        if (d > maxVal) {
                            maxVal = (float) d;
                        }
                    }
                    for (double d : entry.getValue()[1]) {
                        yDataMod.get(entry.getKey()).add((float) d);
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
        for (List<Float> listf : yDataMod.values()) {
            for (float f : listf) {
                if (f > ymax) {
                    ymax = f;
                }
            }
        }
        chart.getStyler().setYAxisMax(ymax);

        ArrayList<List<Float>> allXdata = new ArrayList<>();
        ArrayList<List<Float>> allYdata = new ArrayList<>();
        ArrayList<List<Color>> allColorData = new ArrayList<>();
        ArrayList<List<String>> allNameData = new ArrayList<>();
        if (!xData.isEmpty()) {
            allXdata.add(xData);
            allYdata.add(yData);

            List<Color> colorList = new ArrayList<>();
            List<String> nameList = new ArrayList<>();
            for (int i = 0; i < xData.size(); i++) {
                nameList.add("scatter");
                colorList.add(new Color(0, 0, 0, opacity));
            }
            allColorData.add(colorList);
            allNameData.add(nameList);
        }

        int colorI = 0;
        for (String mass : massesList) {
            List<Float> ix = xDataMod.get(mass);
            List<Float> iy = yDataMod.get(mass);
            if (!ix.isEmpty()) {
                allXdata.add(ix);
                allYdata.add(iy);

                List<Color> colorList = new ArrayList<>();
                List<String> nameList = new ArrayList<>();
                for (int j = 0; j < ix.size(); j++) {
                    nameList.add("scatterMods - " + mass);
                    colorList.add(new Color(
                            65 * (colorI + 1) % 255, 105 * (colorI + 1) % 255, 225 * (colorI + 1) % 255,
                            (int) (Math.min(opacity * 2, 1) * 255f)));
                }
                allColorData.add(colorList);
                allNameData.add(nameList);

                colorI++;
            }
        }

        ArrayList<Float> oneX = mixSeriesFloat(allXdata);
        ArrayList<Float> oneY = mixSeriesFloat(allYdata);
        ArrayList<Color> oneColor = mixSeriesColor(allColorData);
        ArrayList<String> oneName = mixSeriesString(allNameData);

        HashSet<String> includedSeries = new HashSet<>();
        for (int i = 0; i < oneX.size(); i++) {
            boolean showInLegend = false;
            String seriesName = oneName.get(i);
            if (!includedSeries.contains(seriesName)) {
                showInLegend = true;
                includedSeries.add(seriesName);
            } else {
                seriesName = String.valueOf(i);
            }

            XYSeries series = chart.addSeries(seriesName,
                    Collections.singletonList(oneX.get(i)), Collections.singletonList(oneY.get(i)));
            series.setMarkerColor(oneColor.get(i));
            series.setMarker(SeriesMarkers.CIRCLE);
            series.setShowInLegend(showInLegend);
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
            series1.setMarkerColor(new Color(255, 50 * j % 255, 0));
            j++;
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

    //mixes series so that they appear evenly dispersed in plot
    //do this for x and y separately
    private ArrayList<Float> mixSeriesFloat(ArrayList<List<Float>> series) {
        int maxLength = 0;
        for (List<Float> s : series) {
            if (s.size() > maxLength) {
                maxLength = s.size();
            }
        }

        ArrayList<Float> allPoints = new ArrayList<>();
        for (int i = 0; i < maxLength; i++) {
            for (List<Float> s : series) {
                if (i < s.size()) {
                    allPoints.add(s.get(i));
                }
            }
        }

        Collections.reverse(allPoints);
        return allPoints;
    }

    private ArrayList<Color> mixSeriesColor(ArrayList<List<Color>> series) {
        int maxLength = 0;
        for (List<Color> s : series) {
            if (s.size() > maxLength) {
                maxLength = s.size();
            }
        }

        ArrayList<Color> allPoints = new ArrayList<>();
        for (int i = 0; i < maxLength; i++) {
            for (List<Color> s : series) {
                if (i < s.size()) {
                    allPoints.add(s.get(i));
                }
            }
        }

        Collections.reverse(allPoints);
        return allPoints;
    }

    private ArrayList<String> mixSeriesString(ArrayList<List<String>> series) {
        int maxLength = 0;
        for (List<String> s : series) {
            if (s.size() > maxLength) {
                maxLength = s.size();
            }
        }

        ArrayList<String> allPoints = new ArrayList<>();
        for (int i = 0; i < maxLength; i++) {
            for (List<String> s : series) {
                if (i < s.size()) {
                    allPoints.add(s.get(i));
                }
            }
        }

        Collections.reverse(allPoints);
        return allPoints;
    }
}

