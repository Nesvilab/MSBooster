/*
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or
 * implied. See the License for the specific language governing
 * permissions and limitations under the License.
 */

package koinaclasses;

import allconstants.Constants;
import allconstants.NceConstants;
import com.google.common.util.concurrent.AtomicDouble;
import mainsteps.MainClass;
import mainsteps.PeptideObj;
import org.knowm.xchart.BoxChart;
import org.knowm.xchart.BoxChartBuilder;
import org.knowm.xchart.style.BoxStyler;
import peptideptmformatting.PeptideFormatter;
import predictions.PredictionEntryHashMap;
import utils.MyFileUtils;
import utils.ProgressReporter;
import utils.StatMethods;
import writers.PeptideFileCreator;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.LinkedList;
import java.util.List;
import java.util.TreeMap;
import java.util.concurrent.ConcurrentHashMap;
import java.util.concurrent.ExecutionException;
import java.util.concurrent.Future;
import java.util.concurrent.atomic.AtomicInteger;

import static allconstants.Constants.figureDirectory;
import static figures.ExtensionPlotter.plot;
import static utils.Print.printInfo;

public class NCEcalibrator {
    public static Object[] calibrateNCE(String currentModel,
                                        KoinaMethods km, String jsonOutFolder, boolean verbose)
            throws IOException, ExecutionException, InterruptedException {
        ArrayList<PeptideFormatter> peptideFormatterArrayList = km.peptideArraylist;
        HashMap<String, LinkedList<Integer>> scanNums = km.scanNums;
        HashMap<String, LinkedList<PeptideFormatter>> peptides = km.peptides;

        TreeMap<Integer, ArrayList<Double>> similarities = new TreeMap<>();
        double bestMedianDouble = 0d;
        int bestNCEint = 0;
        if (NceConstants.nceModels.contains(currentModel)) {
            printInfo("Calibrating NCE for " + currentModel);

            //create nce calibration directory
            MyFileUtils.createWholeDirectory(jsonOutFolder);

            //write full.tsv
            //get peptides formatted for jsonwriter
            PeptideFileCreator.createPartialFile(jsonOutFolder + File.separator + "NCE_calibration_full.tsv",
                    currentModel, peptideFormatterArrayList);
            HashSet<String> allHits = KoinaMethods.createPartialKoinaSet(currentModel, peptideFormatterArrayList);

            //iterate through every NCE
            //TODO: make this into one method, taking a model and iterating through all NCEs
            ConcurrentHashMap<Integer, ArrayList<Double>> concurrentSimilarities = new ConcurrentHashMap<>();
            AtomicDouble bestMedian = new AtomicDouble(0);
            AtomicInteger bestNCE = new AtomicInteger();
            ProgressReporter pr = new ProgressReporter(NceConstants.maxNCE - NceConstants.minNCE + 1);

            List<Future> futureList = new ArrayList<>(Constants.numThreads);
            for (int NCE = NceConstants.minNCE; NCE < NceConstants.maxNCE + 1; NCE++) {
                int finalNCE = NCE;
                futureList.add(MainClass.executorService.submit(() -> {
                    PredictionEntryHashMap allPreds =
                            km.getKoinaPredictions(allHits, currentModel, finalNCE,
                                    jsonOutFolder + File.separator + finalNCE,
                                    jsonOutFolder + File.separator + "NCE_calibration_full.tsv");

                    //compare pred and exp and set NCE
                    //make mzmlsscannumber objects and set peptide objects and calculate similarities
                    //can move this after reading in mzml files
                    ArrayList<Double> similarity = new ArrayList<>();
                    try {
                        for (PeptideObj peptideObj : km.getPeptideObjects(allPreds, scanNums, peptides)) {
                            if (peptideObj != null) {
                                similarity.add(peptideObj.spectralSimObj.unweightedSpectralEntropy());
                            }
                        }
                    } catch (Exception e) {
                        throw new RuntimeException(e);
                    }

                    //calculate median
                    concurrentSimilarities.put(finalNCE, similarity);
                    double median = StatMethods.medianDouble(similarity);
                    if (median > bestMedian.get()) {
                        bestMedian.set(median);
                        bestNCE.set(finalNCE);
                    }
                    pr.progress();
                }));
            }
            for (Future future : futureList) {
                future.get();
            }
            similarities = new TreeMap<>(concurrentSimilarities);

            printInfo("Best NCE for " + currentModel + " after calibration is " + bestNCE);
            bestMedianDouble = bestMedian.get();
            if (verbose) {
                printInfo("Median similarity for " + currentModel + " is " + String.format("%.4f", bestMedianDouble));
            }
            if (bestNCE.get() == NceConstants.minNCE) {
                printInfo("Consider lowering minNCE below " + NceConstants.minNCE);
            } else if (bestNCE.get() == NceConstants.maxNCE) {
                printInfo("Consider increasing maxNCE above " + NceConstants.maxNCE);
            }
            bestNCEint = bestNCE.get();
        }
        return new Object[]{similarities, bestMedianDouble, bestNCEint};
        //similarities is ordered map of NCE and list of similarities for PSMs
        //bestMedianDouble is best median similarity across NCEs
        //bestNCEint is NCE at which highest median similarity is achieved
    }

    public static void plotNCEchart(String currentModel, TreeMap<Integer, ArrayList<Double>> similarities) {
        try {
            if (!new File(figureDirectory).exists()) {
                new File(figureDirectory).mkdirs();
            }

            int currentNCE = NceConstants.minNCE;
            int added = 0;

            while (currentNCE <= NceConstants.maxNCE) {
                //create figure
                BoxChart chart = new BoxChartBuilder().title("NCE calibration").
                        width(510).height(340).
                        yAxisTitle("unweightedSpectralEntropy").build(); //TODO: adapt to any MS2 feature?
                chart.getStyler().setBoxplotCalCulationMethod(BoxStyler.BoxplotCalCulationMethod.N_LESS_1);

                int startNCE = currentNCE;
                while (added < 6 && currentNCE <= NceConstants.maxNCE) {
                    chart.addSeries(String.valueOf(currentNCE), similarities.get(currentNCE));
                    added++;
                    currentNCE++;
                }
                int endNCE = currentNCE - 1;
                added = 0;
                plot(chart, figureDirectory + File.separator + "NCE_calibration_" +
                        currentModel + startNCE + "to" + endNCE);
            }
        } catch (IOException e) {
            e.printStackTrace();
        }
    }
}
