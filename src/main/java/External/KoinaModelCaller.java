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

import static utils.Print.printError;
import static utils.Print.printInfo;

import Features.*;
import com.google.gson.Gson;
import com.google.gson.stream.JsonReader;
import org.knowm.xchart.BitmapEncoder;
import org.knowm.xchart.VectorGraphicsEncoder;
import org.knowm.xchart.XYChart;
import org.knowm.xchart.XYChartBuilder;
import org.knowm.xchart.style.Styler;

import java.awt.*;
import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.*;
import java.util.concurrent.ConcurrentHashMap;
import java.util.concurrent.Future;
import java.util.concurrent.ScheduledThreadPoolExecutor;
import java.util.concurrent.TimeUnit;
import java.util.concurrent.atomic.AtomicBoolean;
import java.util.concurrent.atomic.AtomicInteger;
import java.util.concurrent.atomic.AtomicLong;

public class KoinaModelCaller {
    private static final int AlphaPeptDeepMzIdx = 1;
    private static final int AlphaPeptDeepIntIdx = 0;
    private static final int AlphaPeptDeepFragIdx = 2;
    private static final int PrositMzIdx = 1;
    private static final int PrositIntIdx = 2;
    private static final int PrositFragIdx = 0;
    private static final int ms2pipMzIdx = 1;
    private static final int ms2pipIntIdx = 2;
    private static final int ms2pipFragIdx = 0;
    private static String modelType;
    private static String finalModel;
    private static AtomicBoolean emptyUrl = new AtomicBoolean(false);

    public KoinaModelCaller(){}

    //TODO: add functions to clean this up
    public void callModel(String model, KoinaLibReader klr, String jsonFolder, ScheduledThreadPoolExecutor executorService,
                          boolean verbose, boolean makeFigure) {
        if (Constants.KoinaURL.isEmpty() && !emptyUrl.get()) {
            emptyUrl.set(true);
            printError("You must specify a URL to use the Koina service via the --KoinaURL parameter");
            printError("Exiting");
            System.exit(1);
        }
        modelType = model.toLowerCase().split("_")[0];
        finalModel = model;

        if (verbose) {
            printInfo("Calling " + model + " model");
        }
        long startTime = System.currentTimeMillis();

        String property = null;
        //decide if this is RT or MS2 model
        if (Constants.KoinaRTmodels.contains(model)) {
            property = "rt";
        } else if (Constants.KoinaMS2models.contains(model)) {
            property = "ms2";
        } else if (Constants.KoinaIMmodels.contains(model)) {
            property = "im";
        } else {
            printError(model + " not in Koina models");
            System.exit(1);
        }

        try {
            //pass json files to curl http request
            File[] fileArray = new File(jsonFolder).listFiles();
            ArrayList<String> filenameArraylist = new ArrayList<>();
            for (File f : fileArray) {
                String fname = f.toString();
                if (fname.endsWith(property + ".json")) {
                    filenameArraylist.add(fname);
                }
            }
            int numProcesses = filenameArraylist.size();

            ProgressReporter pr = new ProgressReporter(numProcesses);
            AtomicLong waitTime;
            if (property.equals("rt") || property.equals("im")) {
                waitTime = new AtomicLong(Constants.initialKoinaMillisecondsToWaitRtIm);
            } else {
                waitTime = new AtomicLong(Constants.initialKoinaMillisecondsToWaitMs2);
            }

            ConcurrentHashMap<Integer, Long> completionTimes = new ConcurrentHashMap<>();
            KoinaTask[] tasks = new KoinaTask[numProcesses];
            Future[] futureList = new Future[numProcesses];
            boolean[] needsToBeRun = new boolean[numProcesses]; //default is false
            Arrays.fill(needsToBeRun, true);
            AtomicLong nextAllowedSubmissionTime = new AtomicLong(System.currentTimeMillis());
            long jobStart = System.currentTimeMillis();

            for (int i = 0; i < numProcesses; i++) {
                KoinaTask task = new KoinaTask(filenameArraylist.get(i), property, model,
                        klr, waitTime);
                tasks[i] = task;
            }

            AtomicInteger finishedJobs = new AtomicInteger(1);
            AtomicInteger attemptedJobs = new AtomicInteger(0);

            int allowedFailedAttempts = 0; //finish all tasks that have not been attempted before going to ones that have failed once, etc

            while (finishedJobs.get() < tasks.length + 1) {
                for (int i = 0; i < tasks.length; i++) {
                    KoinaTask task = tasks[i];
                    if (needsToBeRun[i] && task.failedAttempts == allowedFailedAttempts) {
                        //0.01 seconds between submissions
                        long delay = Math.max(10L - (System.currentTimeMillis() - nextAllowedSubmissionTime.get()), 0);
                        futureList[i] = executorService.schedule(task, delay, TimeUnit.MILLISECONDS);
                        nextAllowedSubmissionTime.set(System.currentTimeMillis() + delay + 10L);

                        needsToBeRun[i] = false;
                        attemptedJobs.incrementAndGet();
                    }

                    Future future = futureList[i];
                    if (task.failedAttempts == allowedFailedAttempts) {
                        if (!task.completed && future.isDone()) {
                            boolean success = (boolean) future.get();
                            if (success) {
                                if (verbose) {
                                    pr.progress();
                                }
                                completionTimes.put(finishedJobs.getAndIncrement(), System.currentTimeMillis() - jobStart);
                                task.completed = true;
                            } else {
                                waitTime.addAndGet(5000);
                                needsToBeRun[i] = true;
                                futureList[i] = null;
                                task.failedAttempts++;
                            }
                            attemptedJobs.decrementAndGet();
                        }
                    }
                }
                if (attemptedJobs.get() == 0) {
                    allowedFailedAttempts++;
                }
                Thread.sleep(100);
            }

            long endTime = System.currentTimeMillis();
            long elapsedTime = endTime - startTime;
            if (verbose) {
                printInfo("cURL and parse time in milliseconds: " + elapsedTime);
            }

            //create plot
            if (makeFigure) {
                double[] xData = new double[numProcesses + 1];
                double[] yData = new double[numProcesses + 1];
                xData[0] = 0;
                yData[0] = 0;

                for (int i = 1; i < numProcesses + 1; i++) {
                    xData[i] = completionTimes.get(i) / 1000d;
                    yData[i] = i;
                }
                Arrays.sort(xData);
                XYChart chart = new XYChartBuilder().width(500).height(500).build();
                chart.setTitle("Koina timing");
                chart.setXAxisTitle("Time (sec)");
                chart.setYAxisTitle("Peptides predicted in " + JSONWriter.maxJsonLength + "'s");
                chart.getStyler().setLegendPosition(Styler.LegendPosition.InsideNW);
                chart.getStyler().setYAxisDecimalPattern("0");
                chart.addSeries(model, xData, yData);

                // Increase the font sizes
                chart.getStyler().setAxisTitleFont(new Font("Helvetica", Font.PLAIN, 21));
                chart.getStyler().setLegendFont(new Font("Helvetica", Font.PLAIN, 21));
                chart.getStyler().setAxisTickLabelsFont(new Font("Helvetica", Font.PLAIN, 15));

                String dir = Constants.outputDirectory + File.separator + "MSBooster_plots";
                if (!new File(dir).exists()) {
                    new File(dir).mkdirs();
                }

                if (Constants.plotExtension.equalsIgnoreCase("png")) {
                    BitmapEncoder.saveBitmap(chart,
                            dir + File.separator + "Koina_timing_" + model + ".png",
                            BitmapEncoder.BitmapFormat.PNG);
                } else if (Constants.plotExtension.equalsIgnoreCase("pdf")) {
                    VectorGraphicsEncoder.saveVectorGraphic(chart,
                            dir + File.separator + "Koina_timing_" + model + ".pdf",
                            VectorGraphicsEncoder.VectorGraphicsFormat.PDF);
                }
            }
        } catch (Exception e) {
            e.printStackTrace();
        }
    }

    public static void parseKoinaOutput(String fileName, String koinaString, String property, String model,
                                        KoinaLibReader klr) throws IOException {
        if (property.equalsIgnoreCase("rt")) {
            String rts = koinaString.split("data")[2];
            String[] results = rts.substring(3, rts.length() - 4).split(",");
            float[] parsedResults = new float[results.length];
            for (int i = 0; i < results.length; i++) {
                parsedResults[i] = Float.parseFloat(results[i]);
            }

            assignRTs(fileName, parsedResults, klr);
        } if (property.equalsIgnoreCase("im")) {
            String ims = koinaString.split("data")[2];
            String[] results = ims.substring(3, ims.length() - 4).split(",");
            float[] parsedResults = new float[results.length];
            for (int i = 0; i < results.length; i++) {
                parsedResults[i] = Float.parseFloat(results[i]);
            }

            assignIMs(fileName, parsedResults, klr);
        } else if (property.equalsIgnoreCase("ms2")) {
            //get indices for processing
            int mzIdx = 0;
            int intIdx = 0;
            int fragIdx = 0;
            if (model.contains("AlphaPept")) {
                mzIdx = AlphaPeptDeepMzIdx;
                intIdx = AlphaPeptDeepIntIdx;
                fragIdx = AlphaPeptDeepFragIdx;
            } else if (model.contains("Prosit")) {
                mzIdx = PrositMzIdx;
                intIdx = PrositIntIdx;
                fragIdx = PrositFragIdx;
            } else if (model.contains("ms2pip")) {
                mzIdx = ms2pipMzIdx;
                intIdx = ms2pipIntIdx;
                fragIdx = ms2pipFragIdx;
            }

            String msInfo = koinaString.split("outputs")[1];
            int numPeptides = Integer.parseInt(msInfo.split("shape")[1].split(",")[0].substring(3));

            msInfo = msInfo.substring(3, msInfo.length() - 3);
            String[] dataResults = msInfo.split("},");

            //intensities
            msInfo = dataResults[intIdx].split("data\":\\[")[1];
            msInfo = msInfo.substring(0, msInfo.length() - 1);
            String[] results = msInfo.split(",");
            int vectorLength = results.length / numPeptides;
            HashSet<Integer>[] acceptedIdx = new HashSet[numPeptides]; //dealing with -1 values
            float[][] allIntensities = new float[numPeptides][];
            for (int i = 0; i < numPeptides; i++) {
                ArrayList<Float> intensities = new ArrayList<>();
                HashSet<Integer> accepted = new HashSet<>();
                for (int j = i * vectorLength; j < (i + 1) * vectorLength; j++) {
                    String result = results[j];
                    float intensity = Float.parseFloat(result);
                    if (intensity > 0) {
                        intensities.add(intensity);
                        accepted.add(j);
                    }
                }
                float[] intensitiesArray = new float[intensities.size()];
                for (int j = 0; j < intensities.size(); j++) {
                    intensitiesArray[j] = intensities.get(j);
                }
                allIntensities[i] = intensitiesArray;
                acceptedIdx[i] = accepted;
            }

            //TODO: when mz is calculated correcly by Koina for PTMs, can add back in
            //mz
//            msInfo = dataResults[mzIdx].split("data\":\\[")[1];
//            msInfo = msInfo.substring(0, msInfo.length() - 1);
//            results = msInfo.split(",");
//            float[][] allMZs = new float[numPeptides][];
//            for (int i = 0; i < numPeptides; i++) {
//                ArrayList<Float> mz = new ArrayList<>();
//                for (int j = i * vectorLength; j < (i + 1) * vectorLength; j++) {
//                    String result = results[j];
//                    if (acceptedIdx[i].contains(j)) {
//                        mz.add(Float.parseFloat(result));
//                    }
//                }
//                float[] mzArray = new float[mz.size()];
//                for (int j = 0; j < mz.size(); j++) {
//                    mzArray[j] = mz.get(j);
//                }
//                allMZs[i] = mzArray;
//            }

            //fragment annotations
            msInfo = dataResults[fragIdx].split("data\":\\[")[1];
            msInfo = msInfo.substring(0, msInfo.length() - 1);
            results = msInfo.split(",");
            String[][] allFragmentIonTypes = new String[numPeptides][];
            int[][] allFragNums = new int[numPeptides][];
            int[][] allCharges = new int[numPeptides][];
            for (int i = 0; i < numPeptides; i++) {
                ArrayList<String> fragmentIonTypes = new ArrayList<>();
                ArrayList<Integer> fragNums = new ArrayList<>();
                ArrayList<Integer> charges = new ArrayList<>();
                for (int j = i * vectorLength; j < (i + 1) * vectorLength; j++) {
                    String result = results[j];
                    result = result.substring(1, result.length() - 1);
                    if (acceptedIdx[i].contains(j)) {
                        String[] info = result.split("\\+");
                        charges.add(Integer.parseInt(info[1]));
                        fragmentIonTypes.add(info[0].substring(0, 1));
                        fragNums.add(Integer.parseInt(info[0].substring(1)));
                    }
                }
                String[] fragmentIonTypesArray = new String[fragmentIonTypes.size()];
                int[] fragNumsArray = new int[fragNums.size()];
                int[] chargesArray = new int[charges.size()];
                for (int j = 0; j < fragmentIonTypes.size(); j++) {
                    fragmentIonTypesArray[j] = fragmentIonTypes.get(j);
                    fragNumsArray[j] = fragNums.get(j);
                    chargesArray[j] = charges.get(j);
                }
                allFragmentIonTypes[i] = fragmentIonTypesArray;
                allFragNums[i] = fragNumsArray;
                allCharges[i] = chargesArray;
            }

            assignMS2(fileName, allIntensities, allFragmentIonTypes, allFragNums, allCharges, klr);
        }
    }

    private static void assignRTs(String fileName, float[] RTs, KoinaLibReader klr) throws IOException {
        String[] peptides = readJSON(fileName, RTs.length);
        PredictionEntryHashMap preds = klr.getPreds();
        for (int i = 0; i < peptides.length; i++) {
            PeptideFormatter pf = new PeptideFormatter(peptides[i], 1, "diann");
            String peptide = pf.getBase() + "|";
            if (modelType.equals("prosit")) {
                peptide = peptide.replace("C[" + PTMhandler.carbamidomethylationMass + "]", "C");
                peptide = peptide.replace("C", "C[" + PTMhandler.carbamidomethylationMass + "]");
            }

            int entries = 0; //in case RT prediction is available but not MS2
            for (int charge = Constants.minPrecursorCharge; charge < Constants.maxPrecursorCharge + 1; charge++) {
                String peptideCharge = peptide + charge;
                if (preds.containsKey(peptideCharge)) {
                    PredictionEntry pe = preds.get(peptideCharge);
                    pe.setRT(RTs[i]);
                    preds.put(peptideCharge, pe);
                    entries++;
                }
            }
            if (entries == 0) { //RT was predicted but not MS2
                for (int charge = Constants.minPrecursorCharge; charge < Constants.maxPrecursorCharge + 1; charge++) {
                    String peptideCharge = peptide + charge;
                    PredictionEntry pe = new PredictionEntry();
                    pe.setRT(RTs[i]);
                    preds.put(peptideCharge, pe);
                }
            }
        }
    }

    private static void assignIMs(String fileName, float[] IMs, KoinaLibReader klr) throws IOException {
        double CCS_IM_COEF = 1059.62245;
        double IM_GAS_MASS = 28.0;

        String[] peptides = readJSON(fileName, IMs.length);
        PredictionEntryHashMap preds = klr.getPreds();
        for (int i = 0; i < peptides.length; i++) {
            String[] peptideSplit = peptides[i].split("\\|");
            PeptideFormatter pf = new PeptideFormatter(peptideSplit[0], peptideSplit[1], "diann");
            String peptide = pf.getBaseCharge();

            PredictionEntry pe;
            if (preds.containsKey(peptide)) {
                pe = preds.get(peptide);
            } else {
                pe = new PredictionEntry();
            }
            float im = IMs[i];
            if (Constants.KoinaCCSmodels.contains(finalModel)) {
                //convert from ccs to 1/K0
                //Mason Schampp according to https://github.com/MannLabs/alphabase/blob/bbaecc380ae157d0f4cc87fffec097ecb7a8ceca/alphabase/peptide/mobility.py#L18
                MassCalculator mc = new MassCalculator(pf.getBase(), peptideSplit[1]);
                double reducedMass = mc.mass + mc.charge * mc.proton;
                reducedMass = reducedMass * IM_GAS_MASS / (reducedMass + IM_GAS_MASS);
                im = (float) (im * Math.sqrt(reducedMass) / mc.charge / CCS_IM_COEF);
            }
            pe.setIM(im);
            preds.put(peptide, pe);
        }
    }

    private static void assignMS2(String fileName, float[][] intensities,
                                  String[][] fragmentIonTypes, int[][] fragNums, int[][] charges, KoinaLibReader klr)
            throws IOException {
        String[] peptides = readJSON(fileName, intensities.length);
        PredictionEntryHashMap preds = klr.getPreds();
        for (int i = 0; i < peptides.length; i++) {
            String[] peptideSplit = peptides[i].split("\\|");
            PeptideFormatter pf = new PeptideFormatter(peptideSplit[0], peptideSplit[1], "diann");
            String peptide = pf.getBaseCharge();
            if (modelType.equals("prosit")) {
                peptide = peptide.replace("C[" + PTMhandler.carbamidomethylationMass + "]", "C");
                peptide = peptide.replace("C", "C[" + PTMhandler.carbamidomethylationMass + "]");
            }

            //problem with Koina not including PTM mass in m/z. Need to calculate here
            String[] pepSplit = peptide.split("\\|");
            MassCalculator mc = new MassCalculator(pepSplit[0], pepSplit[1]);
            float[] mzs = new float[intensities[i].length];
            for (int j = 0; j < intensities[i].length; j++) {
                mzs[j] = mc.calcMass(fragNums[i][j], fragmentIonTypes[i][j], charges[i][j]);
            }

            PredictionEntry pe = new PredictionEntry(mzs, intensities[i], fragNums[i],
                    charges[i], fragmentIonTypes[i], new int[0]);

            if (preds.containsKey(peptide)) {
                PredictionEntry oldPe = preds.get(peptide);
                pe.setRT(oldPe.getRT());
                pe.setIM(oldPe.getIM());
            }
            preds.put(peptide, pe);
        }
    }

    public void assignMissingPeptidePredictions(KoinaLibReader klr, String fulltsv)
            throws IOException {
        BufferedReader TSVReader = new BufferedReader(new FileReader(fulltsv));
        String l;
        String[] line;
        PredictionEntryHashMap preds = klr.getPreds();

        while ((l = TSVReader.readLine()) != null) {
            line = l.split("\t");
            if (!preds.containsKey(line[0] + "|" + line[1])) {
                //get predictionEntry
                PeptideFormatter pf;
                PredictionEntry tmp = null;
                String stripped = "";
                String baseCharge = "";
                switch (modelType) {
                    case "alphapept":
                    case "ms2pip":
                    case "deeplc":
                        pf = new PeptideFormatter(
                                new PeptideFormatter(line[0], line[1], "base").getDiann(), line[1], "diann");
                        baseCharge = pf.getBaseCharge();
                        tmp = preds.get(baseCharge);
                        stripped = pf.getStripped();
                        break;
                    case "prosit":
                        if (finalModel.contains("TMT")) {
                            pf = new PeptideFormatter(
                                    new PeptideFormatter(line[0], line[1], "base").getPrositTMT(), line[1], "prosit");
                        } else {
                            pf = new PeptideFormatter(
                                    new PeptideFormatter(line[0], line[1], "base").getProsit(), line[1], "prosit");
                        }
                        baseCharge = pf.getBaseCharge();
                        tmp = preds.get(baseCharge);
                        stripped = pf.getStripped();
                        break;
                    default:
                        printError(modelType + " not supported by Koina");
                        System.exit(1);
                }
                MassCalculator mc = new MassCalculator(line[0], line[1]);
                try {
                    float[] newMZs = new float[tmp.getMzs().length];
                    for (int i = 0; i < newMZs.length; i++) {
                        newMZs[i] = mc.calcMass(tmp.getFragNums()[i],
                                Constants.flagTOion.get(tmp.getFlags()[i]), tmp.getCharges()[i]);
                    }

                    //add to hashmap
                    PredictionEntry newPred = new PredictionEntry(newMZs, tmp.getIntensities(),
                            tmp.getFragNums(), tmp.getCharges(), tmp.getFragmentIonTypes(), tmp.getFlags());
                    newPred.setRT(tmp.getRT());
                    newPred.setIM(tmp.getIM());
                    preds.put(mc.fullPeptide, newPred);
                } catch (Exception e) {
                    if (!Constants.foundBest && klr.failed) { //allow it to run without error
                        new PredictionEntry();
                        PredictionEntry newPred = new PredictionEntry(new float[]{0}, new float[]{0},
                                new int[]{0}, new int[]{0}, new String[]{"y"}, new int[]{1});
                        newPred.setRT(0);
                        newPred.setIM(0);
                        preds.put(mc.fullPeptide, newPred);
                    } else if (! PeptideSkipper.skipPeptide(stripped, baseCharge.split("\\|")[1])) {
                        e.printStackTrace();
                        printError("Missing peptide to transfer prediction onto " + l + ": " + baseCharge);
                        printError("Exiting now.");
                        System.exit(1);
                    }
                }
            }
        }
    }

    private static String[] readJSON(String fileName, int length) throws IOException {
        String[] peptides = new String[length];
        String[] charges = new String[length];
        ArrayList<String> names = new ArrayList<>();

        Gson gson = new Gson();
        //first pass to get order of names
        JsonReader jr = gson.newJsonReader(new FileReader(fileName));
        jr.beginObject();
        while (jr.hasNext()) {
            String name = jr.nextName();
            if (name.equals("id")) {
                jr.skipValue();
            } else if (name.equals("inputs")) {
                jr.beginArray();
                while (jr.hasNext()) {
                    jr.beginObject();
                    while (jr.hasNext()) {
                        name = jr.nextName();
                        if (name.equals("shape") || name.equals("datatype") || name.equals("data")) {
                            jr.skipValue();
                        } else if (name.equals("name")) {
                            names.add(jr.nextString());
                        }
                    }
                    jr.endObject();
                }
                jr.endArray();
            }
        }
        jr.endObject();

        //extract peptides (and charges) from JSON
        jr = gson.newJsonReader(new FileReader(fileName));
        jr.beginObject();
        int i = 0;
        int inputsIdx = 0;
        while (jr.hasNext()) {
            String name = jr.nextName();
            if (name.equals("id")) {
                jr.skipValue();
            } else if (name.equals("inputs")) {
                jr.beginArray();
                while (jr.hasNext()) {
                    jr.beginObject();
                    String currentName = names.get(inputsIdx);
                    while (jr.hasNext()) {
                        name = jr.nextName();
                        if (name.equals("shape") || name.equals("datatype") || name.equals("name")) {
                            jr.skipValue();
                        } else if (name.equals("data")) {
                            if (currentName.equals("peptide_sequences")) {
                                jr.beginArray();
                                while (jr.hasNext()) {
                                    jr.beginArray();
                                    peptides[i] = jr.nextString();
                                    i++;
                                    jr.endArray();
                                }
                                jr.endArray();
                                i = 0;
                                inputsIdx++;
                            } else if (currentName.equals("precursor_charges")) {
                                jr.beginArray();
                                while (jr.hasNext()) {
                                    jr.beginArray();
                                    charges[i] = jr.nextString();
                                    i++;
                                    jr.endArray();
                                }
                                jr.endArray();
                                i = 0;
                                inputsIdx++;
                            } else {
                                jr.skipValue();
                            }
                        }
                    }
                    jr.endObject();
                }
                jr.endArray();
            }
        }
        jr.endObject();

        if (names.contains("precursor_charges")) { //for ms2 info
            for (i = 0; i < peptides.length; i++) {
                peptides[i] = peptides[i] + "|" + charges[i];
            }
        }
        return peptides;
    }
}
