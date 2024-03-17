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

package External;

import static utils.Print.printInfo;

import Features.Constants;
import Features.MainClass;
import Features.MyFileUtils;
import Features.PeptideObj;
import Features.PredictionEntry;
import Features.StatMethods;
import com.google.common.util.concurrent.AtomicDouble;
import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.TreeMap;
import java.util.concurrent.ConcurrentHashMap;
import java.util.concurrent.ExecutionException;
import java.util.concurrent.Future;
import java.util.concurrent.atomic.AtomicInteger;
import org.knowm.xchart.BitmapEncoder;
import org.knowm.xchart.BoxChart;
import org.knowm.xchart.BoxChartBuilder;
import org.knowm.xchart.style.BoxStyler;

public class NCEcalibrator {
    public static Object[] calibrateNCE(String currentModel,
                                        ArrayList<String> models, KoinaMethods km, String jsonOutFolder)
            throws IOException, ExecutionException, InterruptedException {
        //correct to CID model
        //TODO do in opposite direction. Or more general method to get right fragmentation model
        if (Constants.autoSwitchFragmentation) {
            if (Constants.FragmentationType.equals("CID") &&
                    Constants.spectraRTPredModel.contains("Prosit_2020_intensity_HCD")) {
                String oldModel = Constants.spectraModel;
                Constants.spectraRTPredModel =
                        Constants.spectraRTPredModel.replace("Prosit_2020_intensity_HCD",
                                "Prosit_2020_intensity_CID");
                Constants.spectraModel =
                        Constants.spectraModel.replace("Prosit_2020_intensity_HCD",
                                "Prosit_2020_intensity_CID");
                if (!oldModel.equals(Constants.spectraModel)) {
                    printInfo("Replacing Prosit_2020_intensity_HCD with " +
                            "Prosit_2020_intensity_CID");
                    currentModel = "Prosit_2020_intensity_CID";
                    List<String> modelsList = Arrays.asList(Constants.spectraRTPredModel.split(","));
                    models = new ArrayList<>(modelsList);
                    if (models.size() != 1) {
                        if (models.get(0).equals(models.get(1))) {
                            models.remove(1);
                            Constants.spectraRTPredModel = models.get(0);
                            if (Constants.spectraRTPredModel.equals("DIA-NN")) {
                                Constants.useKoina = false;
                            }
                        }
                    }
                }
            }
        }

        TreeMap<Integer, ArrayList<Double>> similarities = new TreeMap<>();
        double bestMedianDouble = 0d;
        int bestNCEint = 0;
        if (Constants.nceModels.contains(currentModel)) {
            printInfo("Calibrating NCE");

            //create nce calibration directory
            MyFileUtils.createWholeDirectory(jsonOutFolder);

            //write full.tsv
            //get peptides formatted for jsonwriter
            HashSet<String> allHits = km.writeFullPeptideFile(
                    jsonOutFolder + File.separator + "NCE_calibration_full.tsv", currentModel);

            //iterate through every NCE
            //TODO: make this into one method, taking a model and iterating through all NCEs
            ConcurrentHashMap<Integer, ArrayList<Double>> concurrentSimilarities = new ConcurrentHashMap<>();
            AtomicDouble bestMedian = new AtomicDouble(0);
            AtomicInteger bestNCE = new AtomicInteger();
            boolean old = Constants.removeRankPeaks;
            Constants.removeRankPeaks = false;

            List<Future> futureList = new ArrayList<>(Constants.numThreads);
            String finalCurrentModel = currentModel;
            for (int NCE = Constants.minNCE; NCE < Constants.maxNCE + 1; NCE++) {
                int finalNCE = NCE;
                futureList.add(MainClass.executorService.submit(() -> {
                    ConcurrentHashMap<String, PredictionEntry> allPreds =
                            km.getKoinaPredictions(allHits, finalCurrentModel, finalNCE,
                                    jsonOutFolder + File.separator + finalNCE,
                                    jsonOutFolder + File.separator + "NCE_calibration_full.tsv");

                    //compare pred and exp and set NCE
                    //make mzmlsscannumber objects and set peptide objects and calculate similarities
                    //can move this after reading in mzml files
                    ArrayList<Double> similarity = new ArrayList<>();
                    for (PeptideObj peptideObj : km.getPeptideObjects(allPreds)) {
                        similarity.add(peptideObj.spectralSimObj.unweightedSpectralEntropy());
                    }

                    //calculate median
                    concurrentSimilarities.put(finalNCE, similarity);
                    double median = StatMethods.medianDouble(similarity);
                    if (median > bestMedian.get()) {
                        bestMedian.set(median);
                        bestNCE.set(finalNCE);
                    }
                }));
            }
            for (Future future : futureList) {
                future.get();
            }
            similarities = new TreeMap<>(concurrentSimilarities);
            Constants.removeRankPeaks = old;

            printInfo("Best NCE for " + currentModel + " after calibration is " + bestNCE);
            bestMedianDouble = bestMedian.get();
            printInfo("Median similarity is " + String.format("%.4f", bestMedianDouble));
            if (bestNCE.get() == Constants.minNCE) {
                printInfo("Consider lowering minNCE below " + Constants.minNCE);
            } else if (bestNCE.get() == Constants.maxNCE) {
                printInfo("Consider increasing maxNCE above " + Constants.maxNCE);
            }
            bestNCEint = bestNCE.get();
        }
        return new Object[]{currentModel, models, similarities, bestMedianDouble, bestNCEint};
    }

    public static void plotNCEchart(TreeMap<Integer, ArrayList<Double>> similarities) {
        //create figure
        BoxChart chart = new BoxChartBuilder().title("NCE calibration").
                width(200 * (Constants.maxNCE - Constants.minNCE)).height(1000).
                yAxisTitle("unweightedSpectralEntropy").build(); //TODO: adapt to any MS2 feature?
        chart.getStyler().setBoxplotCalCulationMethod(BoxStyler.BoxplotCalCulationMethod.N_LESS_1);

        for (Map.Entry<Integer, ArrayList<Double>> entry : similarities.entrySet()) {
            chart.addSeries(String.valueOf(entry.getKey()), entry.getValue());
        }

        try {
            BitmapEncoder.saveBitmap(chart, Constants.outputDirectory + File.separator +
                            "MSBooster_plots" + File.separator + "NCE_calibration.png",
                    BitmapEncoder.BitmapFormat.PNG);
        } catch (IOException e) {
            e.printStackTrace();
        }
    }
}
