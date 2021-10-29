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
        chart.setXAxisTitle("experimental RT (min)");
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
                if (Float.parseFloat(pep.escore) < Constants.RTescoreCutoff && pep.targetORdecoy == 1) {
                    xData.add(rt);
                    yData.add(pep.RT);
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
        File mzmlPath = new File(mzml.pathStr);
        if (! new File(outFile + File.separator + "RTplots").exists()) {
            new File(outFile + File.separator + "RTplots").mkdirs();
        }
        BitmapEncoder.saveBitmap(chart, outFile + File.separator + "RTplots" + File.separator + "RTplot_" + mzmlPath.getName(),
                BitmapEncoder.BitmapFormat.PNG);
    }

    public static void main(String[] args) throws Exception {
//        mzMLReader mzml = new mzMLReader("C:/Users/yangkl/Downloads/proteomics/hla/" +
//                "AG20210507_SR_HLA_A0201_W632_2IPs_Fxn01.mzml");
//        pinReader pin = new pinReader("C:/Users/yangkl/Downloads/proteomics/hla/" +
//                "AG20210507_SR_HLA_A0201_W632_2IPs_Fxn01.pin");
//        SpectralPredictionMapper predictedSpectra = SpectralPredictionMapper.createSpectralPredictionMapper(
//                "C:/Users/yangkl/Downloads/proteomics/hla/" +
//                "spectraRT.predicted.bin");
        ExecutorService executorService = Executors.newFixedThreadPool(11);
        mzMLReader mzml = new mzMLReader(new mgfFileReader("C:/Users/kevin/OneDriveUmich/proteomics/mzml/" +
                "20180819_TIMS2_12-2_AnBr_SA_200ng_HeLa_50cm_120min_100ms_11CT_3_A1_01_2769.mgf", true, executorService));
        pinReader pin = new pinReader("C:/Users/kevin/Downloads/proteomics/timsTOF/" +
                "20180819_TIMS2_12-2_AnBr_SA_200ng_HeLa_50cm_120min_100ms_11CT_3_A1_01_2769.pin");
        SpectralPredictionMapper predictedSpectra = SpectralPredictionMapper.createSpectralPredictionMapper(
                "C:/Users/kevin/Downloads/proteomics/timsTOF/" +
                        "spectraRT.predicted.bin");
        mzml.setPinEntries(pin, predictedSpectra);
        mzml.setLOESS(Constants.RTregressionSize, Constants.bandwidth, Constants.robustIters, "RT");
        RTCalibrationFigure fig = new RTCalibrationFigure(mzml, "C:/Users/kevin/Downloads/proteomics/narrow/", 0.2f);
        executorService.shutdown();
    }
}
