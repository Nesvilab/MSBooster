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
import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;
import java.util.TreeSet;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;

public class RTCalibrationFigure {
    public RTCalibrationFigure(mzMLReader mzml, String outFile, float opacity) throws IOException {
        XYChart chart = new XYChartBuilder().width(1000).height(1000).build();
        //new SwingWrapper<>(chart).displayChart();


        chart.getStyler().setDefaultSeriesRenderStyle(XYSeries.XYSeriesRenderStyle.Scatter);
        chart.getStyler().setChartTitleVisible(true);
        chart.setTitle(mzml.pathStr);
        chart.setXAxisTitle("experimental RT");
        chart.setYAxisTitle("predicted RT");
        chart.getStyler().setLegendVisible(false);
        chart.getStyler().setMarkerSize(3);
        chart.getStyler().setYAxisGroupPosition(0, Styler.YAxisPosition.Right);
        // Series
        List<Float> xData = new ArrayList<Float>();
        List<Float> yData = new ArrayList<Float>();
        float minRT = Float.MAX_VALUE;
        float maxRT = Float.MIN_VALUE;
        for (int scanNum : new TreeSet<Integer>(mzml.scanNumberObjects.keySet())) {
            mzmlScanNumber scanNumObj = mzml.getScanNumObject(scanNum);
            float rt = scanNumObj.RT; //experimental RT for this scan
            if (rt < minRT) {
                minRT = rt;
            }
            if (rt > maxRT) {
                maxRT = rt;
            }

            for (int i = 1; i < scanNumObj.peptideObjects.size() + 1; i++) {
                peptideObj pep = scanNumObj.getPeptideObject(i);
                //only get best ones
                if (Float.parseFloat(pep.escore) < Constants.RTescoreCutoff && pep.spectralSimObj.predMZs[0] != 0f) {
                    xData.add(rt);
                    yData.add(pep.RT);
                } else {
                    break;
                }
            }
        }

        XYSeries series = chart.addSeries("scatter", xData, yData);
        series.setMarkerColor(new Color(0, 0, 0, opacity));

        //loess regression
        // generates Log data
        List<Float> x1Data = new ArrayList<Float>();
        List<Double> y1Data = new ArrayList<Double>();
        for (float i = minRT; i < maxRT; i = i + (maxRT - minRT) / 1000f) {
            x1Data.add(i);
            y1Data.add(mzml.RTLOESS.invoke((double) i));
        }
        XYSeries series1 = chart.addSeries("regression", x1Data, y1Data);
        series1.setMarkerColor(new Color(243, 9, 9));
        //series.setLineWidth(1);
        //series.setLineColor(XChartSeriesColors.RED);

//        BitmapEncoder.saveBitmap(chart, "C:/Users/yangkl/Downloads/proteomics/hla/Sample_Chart",
//                BitmapEncoder.BitmapFormat.PNG);
        String pinPath = new File(outFile).getParent();
        if (! new File(pinPath + File.separator + "MSBooster_RTplots").exists()) {
            new File(pinPath + File.separator + "MSBooster_RTplots").mkdirs();
        }
        String pinName = new File(outFile).getName();
        BitmapEncoder.saveBitmap(chart, pinPath + File.separator + "MSBooster_RTplots" + File.separator +
                        pinName.substring(0, pinName.length() - 4),
                BitmapEncoder.BitmapFormat.PNG);
    }
}
