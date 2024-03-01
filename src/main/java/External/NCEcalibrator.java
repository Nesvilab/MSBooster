package External;

import Features.*;
import com.google.common.util.concurrent.AtomicDouble;
import org.apache.commons.io.FileUtils;
import org.knowm.xchart.BitmapEncoder;
import org.knowm.xchart.BoxChart;
import org.knowm.xchart.BoxChartBuilder;
import org.knowm.xchart.style.BoxStyler;
import umich.ms.fileio.exceptions.FileParsingException;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.*;
import java.util.concurrent.*;
import java.util.concurrent.atomic.AtomicInteger;

public class NCEcalibrator {
    public static Object[] calibrateNCE(PinMzmlMatcher pmMatcher, String currentModel,
                                        ArrayList<String> models, KoinaMethods km, String jsonOutFolder)
            throws IOException, FileParsingException, ExecutionException, InterruptedException {
        //correct to CID model
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
                System.out.println("Replacing Prosit_2020_intensity_HCD with " +
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

        TreeMap<Integer, ArrayList<Double>> similarities = new TreeMap<>();
        double bestMedianDouble = 0d;
        int bestNCEint = 0;
        if (Constants.nceModels.contains(currentModel)) {
            System.out.println("Calibrating NCE");

            //create nce calibration directory
            if (Files.exists(Paths.get(jsonOutFolder))) {
                try {
                    FileUtils.cleanDirectory(new File(jsonOutFolder));
                } catch (IOException e) {}
            } else {
                Files.createDirectories(Paths.get(jsonOutFolder));
            }

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

            System.out.println("Best NCE for " + currentModel + " after calibration is " + bestNCE);
            bestMedianDouble = bestMedian.get();
            System.out.println("Median similarity is " + String.format("%.4f", bestMedianDouble));
            if (bestNCE.get() == Constants.minNCE) {
                System.out.println("Consider lowering minNCE below " + Constants.minNCE);
            } else if (bestNCE.get() == Constants.maxNCE) {
                System.out.println("Consider increasing maxNCE above " + Constants.maxNCE);
            }
            bestNCEint = bestNCE.get();

            //clean directories
            try {
                FileUtils.cleanDirectory(new File(jsonOutFolder));
                Files.deleteIfExists(Paths.get(jsonOutFolder));
            } catch (Exception e) {
                try {
                    FileUtils.cleanDirectory(new File(jsonOutFolder));
                    Files.deleteIfExists(Paths.get(jsonOutFolder));
                } catch (Exception e2) {}
            }
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
