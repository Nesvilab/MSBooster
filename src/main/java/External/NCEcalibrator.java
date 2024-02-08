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
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.*;
import java.util.concurrent.*;

public class NCEcalibrator {
    public static Object[] calibrateNCE(PinMzmlMatcher pmMatcher, String currentModel, ArrayList<String> models)
            throws IOException, FileParsingException, ExecutionException, InterruptedException {
        //first load mzml files
        for (int j = 0; j < pmMatcher.mzmlReaders.length; j++) {
            if (pmMatcher.mzmlReaders[j] == null) {
                MzmlReader mzml = new MzmlReader(pmMatcher.mzmlFiles[j].getCanonicalPath());
                pmMatcher.mzmlReaders[j] = mzml;
            }

            if (Constants.FragmentationType.isEmpty()) {
                try {
                    Set<String> fragTypes = pmMatcher.mzmlReaders[j]
                            .scanNumberObjects.firstEntry().getValue().NCEs.keySet();
                    if (fragTypes.contains("CID")) {
                        Constants.FragmentationType = "CID";
                    } else {
                        Constants.FragmentationType = "HCD";
                    }
                } catch (Exception e) {
                    System.out.println("Setting fragmentation type to HCD. " +
                            "You can specify this with '--FragmentationType' via the command line " +
                            "or 'FragmentationType=' in the param file.");
                    Constants.FragmentationType = "HCD";
                }
            }
        }

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

        if (Constants.nceModels.contains(Constants.spectraModel)) {
            System.out.println("Calibrating NCE");

            //need to collect top 1000 peptides for calibration
            //approximate by doing a subset per pin
            int numTopPSMs = (int) Math.ceil((float) Constants.numPSMsToCalibrateNCE /
                    (float) pmMatcher.pinFiles.length);
            HashSet<String> peptideSet = new HashSet<>();
            HashMap<String, LinkedList<Integer>> scanNums = new HashMap<>();
            HashMap<String, LinkedList<String>> peptides = new HashMap<>();
            for (int j = 0; j < pmMatcher.pinFiles.length; j++) {
                File pinFile = pmMatcher.pinFiles[j];
                PinReader pinReader = new PinReader(pinFile.getAbsolutePath());
                LinkedList[] topPSMs = pinReader.getTopPSMs(numTopPSMs);
                peptideSet.addAll(topPSMs[0]);
                scanNums.put(pmMatcher.mzmlFiles[j].getName(), topPSMs[1]);
                peptides.put(pmMatcher.mzmlFiles[j].getName(), topPSMs[0]);
                if (Constants.instrument.isEmpty()) {
                    pinReader.attachMzml(pmMatcher.mzmlReaders[j]);
                    Constants.instrument = pinReader.getInstrument();
                }
            }

            //write to file
            HashSet<String> allHits = new HashSet<>();
            String jsonOutFolder = Constants.outputDirectory + File.separator + "NCE_calibration";
            if (Files.exists(Paths.get(jsonOutFolder))) {
                try {
                    FileUtils.cleanDirectory(new File(jsonOutFolder));
                } catch (IOException e) {}
            } else {
                Files.createDirectories(Paths.get(jsonOutFolder));
            }
            FileWriter myWriter = new FileWriter(Constants.outputDirectory +
                    File.separator + "NCE_calibration" + File.separator + "NCE_calibration_full.tsv");
            for (String peptide : peptideSet) {
                String[] peptFormats = peptide.split(",");

                String stripped = peptFormats[2];

                if ((currentModel.contains("Prosit") || currentModel.contains("ms2pip") || currentModel.contains("Deeplc"))
                        && stripped.contains("U")) { // no peptides with U
                    continue;
                }
                if (currentModel.contains("ms2pip") && stripped.length() > 30) { //peptide has length limit
                    continue;
                }
                if (currentModel.contains("Prosit") && currentModel.contains("TMT") && stripped.length() > 30) {
                    continue;
                }

                String pep = peptFormats[1].replace("UniMod", "UNIMOD"); //need diann here
                if (pep.contains("[TMT]")) {
                    pep = pep.replace("[TMT]", "[UNIMOD:737]");
                }

                if (pep.startsWith("[")) { //this is the case for all n term mods //TODO deal with c term mods
                    int splitpoint = pep.indexOf("]");
                    if (currentModel.contains("Prosit")) {
                        if (currentModel.contains("TMT") && pep.startsWith("[UNIMOD:737]")) {
                            pep = pep.substring(0, splitpoint + 1) + "-" + pep.substring(splitpoint + 1);
                        } else {
                            pep = pep.substring(splitpoint + 1);
                        }
                    } else {
                        pep = pep.substring(0, splitpoint + 1) + "-" + pep.substring(splitpoint + 1);
                    }
                }
                if (currentModel.contains("Prosit") && currentModel.contains("TMT")) {
                    pep = pep.replace("S[UNIMOD:737]", "S");
                    if (!pep.startsWith("[")) {
                        pep = "[UNIMOD:737]-" + pep;
                    }
                }
                String[] baseCharge = peptFormats[0].split("\\|");
                allHits.add(pep + "," + baseCharge[1]);

                //need to generate full.tsv
                myWriter.write(baseCharge[0] + "\t" + baseCharge[1] + "\n");
            }
            myWriter.close();


            //iterate through every NCE
            ConcurrentHashMap<Integer, ArrayList<Double>> concurrentSimilarities = new ConcurrentHashMap<>();
            AtomicDouble bestMedian = new AtomicDouble(0);
            boolean old = Constants.removeRankPeaks;
            Constants.removeRankPeaks = false;

            List<Future> futureList = new ArrayList<>(Constants.numThreads);
            String finalCurrentModel = currentModel;
            for (int NCE = Constants.minNCE; NCE < Constants.maxNCE + 1; NCE++) {
                int finalNCE = NCE;
                futureList.add(MainClass.executorService.submit(() -> {
                    ScheduledThreadPoolExecutor executorService = new ScheduledThreadPoolExecutor(1);

                    HashSet<String> NCEhits = new HashSet<>();
                    for (String s : allHits) {
                        NCEhits.add(s + "," + finalNCE + "," + Constants.instrument + "," + Constants.FragmentationType);
                    }

                    JSONWriter jw = new JSONWriter(finalCurrentModel, NCEhits);
                    String jsonFolder = "";
                    try {
                        jsonFolder = jw.write(true, "NCE_calibration" + File.separator + finalNCE,
                                executorService);
                    } catch (Exception e) {
                        e.printStackTrace();
                        System.exit(1);
                    }

                    //send predictions to Koina
                    KoinaLibReader klr = new KoinaLibReader();
                    KoinaModelCaller kmc = new KoinaModelCaller();
                    kmc.callModel(finalCurrentModel, klr, jsonFolder, executorService, false, false);
                    executorService.shutdown();
                    try {
                        kmc.assignMissingPeptidePredictions(klr, jsonOutFolder +
                                File.separator + "NCE_calibration_full.tsv");
                    } catch (IOException e) {
                        throw new RuntimeException(e);
                    }
                    ConcurrentHashMap<String, PredictionEntry> allPreds = klr.allPreds;

                    //compare pred and exp and set NCE
                    //make mzmlsscannumber objects and set peptide objects and calculate similarities
                    //can move this after reading in mzml files
                    ArrayList<Double> similarity = new ArrayList<>();

                    for (int j = 0; j < pmMatcher.mzmlReaders.length; j++) {
                        MzmlReader mzmlReader = pmMatcher.mzmlReaders[j];
                        LinkedList<Integer> thisScanNums = scanNums.get(pmMatcher.mzmlFiles[j].getName());
                        LinkedList<String> thisPeptides = peptides.get(pmMatcher.mzmlFiles[j].getName());

                        for (int k = 0; k < thisScanNums.size(); k++) {
                            int scanNum = thisScanNums.get(k);
                            MzmlScanNumber msn = mzmlReader.getScanNumObject(scanNum);

                            String peptide = thisPeptides.get(k).split(",")[0];
                            String[] baseCharge = peptide.split("\\|");

                            PeptideObj pobj = msn.setPeptideObject(
                                    new PeptideFormatter(baseCharge[0], baseCharge[1], "base"),
                                    1, 1, "0", allPreds, false);
                            similarity.add(pobj.spectralSimObj.unweightedSpectralEntropy());
                        }
                    }

                    //calculate median
                    concurrentSimilarities.put(finalNCE, similarity);
                    double median = StatMethods.medianDouble(similarity);
                    if (median > bestMedian.get()) {
                        bestMedian.set(median);
                        Constants.NCE = String.valueOf(finalNCE);
                    }
                }));
            }
            for (Future future : futureList) {
                future.get();
            }
            TreeMap<Integer, ArrayList<Double>> similarities = new TreeMap<>(concurrentSimilarities);
            Constants.removeRankPeaks = old;

            System.out.println("Best NCE after calibration is " + Constants.NCE);
            if (Constants.NCE.equals(String.valueOf(Constants.minNCE))) {
                System.out.println("Consider lowering minNCE below " + Constants.minNCE);
            } else if (Constants.NCE.equals(String.valueOf(Constants.maxNCE))) {
                System.out.println("Consider increasing maxNCE above " + Constants.maxNCE);
            }

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
                                "MSBooster_plots" + File.separator + "calibration.png",
                        BitmapEncoder.BitmapFormat.PNG);
            } catch (IOException e) {
                e.printStackTrace();
            }
            try {
                FileUtils.cleanDirectory(new File(Constants.outputDirectory + File.separator + "NCE_calibration"));
                Files.deleteIfExists(Paths.get(Constants.outputDirectory + File.separator + "NCE_calibration"));
            } catch (Exception e) {
                try {
                    FileUtils.cleanDirectory(new File(Constants.outputDirectory + File.separator + "NCE_calibration"));
                    Files.deleteIfExists(Paths.get(Constants.outputDirectory + File.separator + "NCE_calibration"));
                } catch (Exception e2) {}
            }
        }
        return new Object[]{currentModel, models};
    }
}
