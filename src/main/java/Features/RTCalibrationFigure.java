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

import org.knowm.xchart.BitmapEncoder;
import org.knowm.xchart.XYChart;
import org.knowm.xchart.XYChartBuilder;
import org.knowm.xchart.XYSeries;
import org.knowm.xchart.style.Styler;

import java.awt.*;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.*;
import java.util.List;

public class RTCalibrationFigure {
    public RTCalibrationFigure(MzmlReader mzml, String outFile, float opacity) throws IOException {
        String pinPath = new File(outFile).getParent();
        String pinName = new File(outFile).getName();

        String dir = pinPath + File.separator + "MSBooster_plots";
        if (! new File(dir).exists()) {
            new File(dir).mkdirs();
        }
        if (! new File(dir + File.separator + "RT_calibration_curves").exists()) {
            new File(dir + File.separator + "RT_calibration_curves").mkdirs();
        }

        XYChart chart = new XYChartBuilder().width(1000).height(1000).build();
        chart.getStyler().setDefaultSeriesRenderStyle(XYSeries.XYSeriesRenderStyle.Scatter);
        chart.getStyler().setChartTitleVisible(true);
        chart.setTitle(mzml.pathStr);
        chart.setXAxisTitle("experimental RT");
        chart.setYAxisTitle("predicted RT");
        chart.getStyler().setLegendVisible(false);
        chart.getStyler().setMarkerSize(8);
        chart.getStyler().setYAxisGroupPosition(0, Styler.YAxisPosition.Right);
        chart.getStyler().setLegendVisible(true);
        chart.getStyler().setLegendPosition(Styler.LegendPosition.InsideNW);

        // Series
        List<Float> xData = new ArrayList<Float>();
        List<Float> yData = new ArrayList<Float>();
        List<Float> eScores = new ArrayList<Float>();

        //for PTMs besides oxM and C57
        List<List<Float>> xDataMod = new ArrayList<>();
        List<List<Float>> yDataMod = new ArrayList<>();

        ArrayList<String> massesList = new ArrayList<>();
        massesList.addAll(Arrays.asList(Constants.RTmassesForCalibration.split(",")));
        if (! Constants.RTfigure_masses.equals("")) {
            for (String mass : Constants.RTfigure_masses.split(",")) {
                if (! massesList.contains(mass)) {
                    massesList.add(mass);
                }
            }
        }
        massesList.remove("");
        for (int i = 0; i < massesList.size(); i++) {
            xDataMod.add(new ArrayList<>());
            yDataMod.add(new ArrayList<>());
        }

        HashMap<String, Float> bestEScore = new HashMap<>();
//        for (int scanNum : new TreeSet<Integer>(mzml.scanNumberObjects.keySet())) {
//            MzmlScanNumber scanNumObj = mzml.getScanNumObject(scanNum);
//            float rt = scanNumObj.RT; //experimental RT for this scan
//            if (rt < minRT) {
//                minRT = rt;
//            }
//            if (rt > maxRT) {
//                maxRT = rt;
//            }
//
//            for (int i = 1; i < scanNumObj.peptideObjects.size() + 1; i++) {
//                PeptideObj pep = scanNumObj.getPeptideObject(i);
//                //only get best ones
//                float currentScore = Float.parseFloat(pep.escore);
////                if (currentScore < Constants.RTescoreCutoff && pep.spectralSimObj.predIntensities[0] != 0f) {
//                if (currentScore < Constants.RTescoreCutoff) {
//                    if (Constants.plotBestPSMPerPeptide) {
//                        if (bestEScore.containsKey(pep.name)) {
//                            if (bestEScore.get(pep.name) > currentScore) {
//                                bestEScore.put(pep.name, currentScore);
//                            }
//                        } else {
//                            bestEScore.put(pep.name, currentScore);
//                        }
//                    } else {
//                        eScores.add(Float.parseFloat(pep.escore));
//                    }
//                } else {
//                    break;
//                }
//            }
//        }
//        if (Constants.plotBestPSMPerPeptide) {
//            eScores.addAll(bestEScore.values());
//        }
//
//        eScores.sort(Comparator.naturalOrder());
//        //this will plot more than numRTplot if PSMs are tied
//        float maxEScore = eScores.get(Math.min(Constants.numRTplot, eScores.size()) - 1);
//        for (int scanNum : new TreeSet<Integer>(mzml.scanNumberObjects.keySet())) {
//            MzmlScanNumber scanNumObj = mzml.getScanNumObject(scanNum);
//            float rt = scanNumObj.RT; //experimental RT for this scan
//            for (int i = 1; i < scanNumObj.peptideObjects.size() + 1; i++) {
//                PeptideObj pep = scanNumObj.getPeptideObject(i);
//                //only get best ones
//                float parsedScore = Float.parseFloat(pep.escore);
////                if (parsedScore <= maxEScore && pep.spectralSimObj.predIntensities[0] != 0f) {
//                if (parsedScore <= maxEScore) {
//                    if (Constants.plotBestPSMPerPeptide && parsedScore != bestEScore.get(pep.name)) {
//                        break;
//                    }
//                    boolean containsMod = false;
//                    int modIndex = 0;
//                    for (int j = 0; j < massesList.size(); j++) {
//                        String[] masses = massesList.get(j).split("/");
//                        for (String mass : masses) {
//                            if (pep.name.contains(mass)) {
//                                containsMod = true;
//                                modIndex = j;
//                                break;
//                            }
//                        }
//                    }
//                    if (! containsMod) {
//                        xData.add(rt);
//                        yData.add(pep.RT);
//                    } else {
//                        xDataMod.get(modIndex).add(rt);
//                        yDataMod.get(modIndex).add(pep.RT);
//                    }
//                } else {
//                    break;
//                }
//            }
//        }

        //set y lim
        int modIdx = 0;
        float minRT = Float.MAX_VALUE;
        float maxRT = Float.MIN_VALUE;
        for (Map.Entry<String, double[][]> entry : mzml.expAndPredRTs.entrySet()) {
            if (entry.getKey().isEmpty() || entry.getKey().equals("others")) {
                for (double d : entry.getValue()[0]) {
                    xData.add((float) d);
                    if (d < minRT) {
                        minRT = (float) d;
                    }
                    if (d > maxRT) {
                        maxRT = (float) d;
                    }
                }
                for (double d : entry.getValue()[1]) {
                    yData.add((float) d);
                }
            } else {
                for (double d : entry.getValue()[0]) {
                    xDataMod.get(modIdx).add((float) d);
                    if (d < minRT) {
                        minRT = (float) d;
                    }
                    if (d > maxRT) {
                        maxRT = (float) d;
                    }
                }
                for (double d : entry.getValue()[1]) {
                    yDataMod.get(modIdx).add((float) d);
                }
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

        if (xData.size() > 0) {
            XYSeries series = chart.addSeries("scatter", xData, yData);
            series.setMarkerColor(new Color(0, 0, 0, opacity));
        }

        for (int i = 0; i < xDataMod.size(); i++) {
            List<Float> ix = xDataMod.get(i);
            List<Float> iy = yDataMod.get(i);
            if (ix.size() > 0) {
                XYSeries seriesMod = chart.addSeries("scatterMods - " + massesList.get(i), ix, iy);
                seriesMod.setMarkerColor(new Color(65 * (i + 1) % 255, 105 * (i + 1) % 255, 225 * (i + 1) % 255));
            }
        }

        //loess regression
        // generates Log data
        int j = 1;
        for (String mass : mzml.RTLOESS.keySet()) {
            BufferedWriter calibrationPoints = null;
            if (Constants.writeCalibration) {
                calibrationPoints = new BufferedWriter(new FileWriter(
                        pinPath + File.separator + "MSBooster_RTplots" + File.separator +
                                pinName.substring(0, pinName.length() - 4) + "_calibration" + mass + ".csv"));
                calibrationPoints.write("experimental RT,predicted RT\n");
            }

            List<Float> x1Data = new ArrayList<Float>();
            List<Double> y1Data = new ArrayList<Double>();
            for (float i = minRT; i < maxRT; i = i + (maxRT - minRT) / 1000f) {
                x1Data.add(i);
                double y = mzml.RTLOESS.get(mass).invoke((double) i); //TODO adapt
                y1Data.add(y);
                if (Constants.writeCalibration) {
                    calibrationPoints.write(i + "," + y + "\n");
                }
            }
            if (Constants.writeCalibration) {
                calibrationPoints.close();
            }
            XYSeries series1;
            if (mass.equals("")) {
                series1 = chart.addSeries("regression", x1Data, y1Data);
            } else {
                series1 = chart.addSeries("regression - " + mass, x1Data, y1Data);
            }
            series1.setMarkerColor(new Color(243 * j % 255, 9 * j % 255, 9 * j % 255));
            j += 2;
        }

        BitmapEncoder.saveBitmap(chart, pinPath + File.separator + "MSBooster_plots" + File.separator +
                        "RT_calibration_curves" + File.separator + pinName.substring(0, pinName.length() - 4),
                BitmapEncoder.BitmapFormat.PNG);
    }
}
