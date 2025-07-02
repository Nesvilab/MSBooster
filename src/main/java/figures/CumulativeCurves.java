package figures;

import org.knowm.xchart.XYChart;
import org.knowm.xchart.XYChartBuilder;
import org.knowm.xchart.XYSeries;
import org.knowm.xchart.style.Styler;

import java.io.File;
import java.io.IOException;
import java.util.List;
import java.util.Map;
import java.util.TreeMap;

import static allconstants.Constants.figureDirectory;
import static figures.ExtensionPlotter.plot;

public class CumulativeCurves {
    public CumulativeCurves(String name, TreeMap<String, List<Float>[]> curves, String mode) throws IOException {
        XYChart chart = new XYChartBuilder().width(1500).height(1000).build();
        chart.getStyler().setDefaultSeriesRenderStyle(XYSeries.XYSeriesRenderStyle.Scatter);
        chart.getStyler().setChartTitleVisible(true);
        chart.setTitle(name);
        chart.setXAxisTitle("experimental " + mode);
        chart.setYAxisTitle("predicted " + mode);
        chart.getStyler().setMarkerSize(6);
        chart.getStyler().setYAxisGroupPosition(0, Styler.YAxisPosition.Right);
        chart.getStyler().setLegendVisible(true);
        chart.getStyler().setLegendPosition(Styler.LegendPosition.InsideNW);

        String pathSuffix = "";
        if (name.contains("&")) {
            pathSuffix = "_" + name.split("&")[0];
        }

        int colorIdx = 1;
        for (Map.Entry<String, List<Float>[]> entry : curves.entrySet()) {
            XYSeries series = chart.addSeries(entry.getKey(), entry.getValue()[0], entry.getValue()[1]);
            series.setLineColor(FigureUtils.getColor(colorIdx));
            series.setMarkerColor(FigureUtils.getColor(colorIdx));
            colorIdx++;
        }

        plot(chart, figureDirectory + File.separator + "cumulativeQC" +
                File.separator + mode + pathSuffix + "_regression_curves");
    }
}
