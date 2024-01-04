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

import Features.*;
import com.google.gson.Gson;
import com.google.gson.stream.JsonReader;

import java.io.*;
import java.util.*;
import java.util.concurrent.ConcurrentHashMap;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Future;
import java.util.concurrent.atomic.AtomicInteger;

public class KoinaModelCaller {
    private final int AlphaPeptDeepMzIdx = 1;
    private final int AlphaPeptDeepIntIdx = 0;
    private final int AlphaPeptDeepFragIdx = 2;
    private final int PrositMzIdx = 1;
    private final int PrositIntIdx = 2;
    private final int PrositFragIdx = 0;
    private final int ms2pipMzIdx = 1;
    private final int ms2pipIntIdx = 2;
    private final int ms2pipFragIdx = 0;
    private String modelType;
    private final int maxJobs = Math.min(150, Constants.numThreads * 2); //stable up to 700, according to Koina

    public KoinaModelCaller(){}

    //TODO: add functions to clean this up
    public void callModel(String model, KoinaLibReader klr, ExecutorService executorService) {
        this.modelType = model.toLowerCase().split("_")[0];

        System.out.println("Calling " + model + " model");
        long startTime = System.currentTimeMillis();

        String property = null;
        //decide if this is RT or MS2 model
        if (Constants.KoinaRTmodels.contains(model)) {
            property = "rt";
        } else if (Constants.KoinaMS2models.contains(model)) {
            property = "ms2";
        } else {
            System.out.println(model + " not in Koina models");
            System.exit(1);
        }

        try {
            //pass json files to curl http request
            File[] fileArray = new File(Constants.JsonDirectory).listFiles();
            ArrayList<String> filenameArraylist = new ArrayList<>();
            for (File f : fileArray) {
                String fname = f.toString();
                if (fname.endsWith(property + ".json")) {
                    filenameArraylist.add(fname);
                }
            }
            String[] filenameArray = new String[filenameArraylist.size()];
            for (int i = 0; i < filenameArray.length; i++) {
                filenameArray[i] = filenameArraylist.get(i);
            }

            Process[] processes = new Process[filenameArray.length];
            BufferedReader[] readers = new BufferedReader[filenameArray.length];
            String[] commands = new String[filenameArray.length];
            int numProcesses = 0;

            //https://www.baeldung.com/linux/curl-parallel-download#multiple-curl-instances-with-a-for-loop
            for (String fileString : filenameArray) {
                String command = "curl -s -H content-type:application/json -d @" + fileString +
                        " https://koina.proteomicsdb.org/v2/models/" + model + "/infer";
                commands[numProcesses] = command;

                //modified so that there are only 2 active jobs per thread at once
                if (numProcesses < maxJobs) {
                    ProcessBuilder builder = new ProcessBuilder(command.split(" "));
                    builder.redirectErrorStream(true);
                    processes[numProcesses] = builder.start();
                    readers[numProcesses] = new BufferedReader(new InputStreamReader(processes[numProcesses].getInputStream()));
                }
                numProcesses += 1;
            }

            ProgressReporter pr = new ProgressReporter(numProcesses);
            List<Future> futureList = new ArrayList<>();
            AtomicInteger jobIdx = new AtomicInteger(maxJobs);

            for (int i = 0; i < numProcesses; i++) {
                String finalProperty = property;
                int finalI = i;
                int finalNumProcesses = numProcesses;
                futureList.add(executorService.submit(() -> {
                    int attempts = 0;
                    while (attempts < Constants.numKoinaAttempts) {
                        StringBuilder koinaSb = new StringBuilder();
                        try {
                            Process p = processes[finalI];
                            BufferedReader reader = readers[finalI];

                            String line = "";
                            while ((line = reader.readLine()) != null) {
                                koinaSb.append(line);
                            }
                            reader.close();
                            p.waitFor();
                            p.destroy();

                            parseKoinaOutput(filenameArray[finalI], koinaSb.toString(),
                                    finalProperty, model, klr);
                            break;
                        } catch (Exception e) {
                            attempts++;
                            if (attempts == Constants.numKoinaAttempts) {
                                System.out.println(commands[finalI]);
                                System.out.println(filenameArray[finalI] + " had output that ended in: ");
                                System.out.println(koinaSb.toString().substring(Math.max(0, koinaSb.toString().length() - 1000)));
                                System.out.println("Retried calling " + filenameArray[finalI] + " " + attempts +
                                        " times. Exiting.");
                                e.printStackTrace();
                                System.exit(1);
                                break;
                            }

                            //set up again
                            String command = "curl -s -H content-type:application/json -d @" + filenameArray[finalI] +
                                    " https://koina.proteomicsdb.org/v2/models/" + model + "/infer";

                            ProcessBuilder builder = new ProcessBuilder(command.split(" "));
                            builder.redirectErrorStream(true);
                            try {
                                processes[finalI] = builder.start();
                            } catch (IOException ioException) {
                                ioException.printStackTrace();
                            }
                            readers[finalI] = new BufferedReader(new InputStreamReader(processes[finalI].getInputStream()));
                        }
                    }
                    pr.progress();

                    int stoppingPoint = jobIdx.get();
                    if (stoppingPoint < finalNumProcesses) {
                        ProcessBuilder builder = new ProcessBuilder(commands[stoppingPoint].split(" "));
                        builder.redirectErrorStream(true);
                        try {
                            processes[stoppingPoint] = builder.start();
                        } catch (IOException e) {
                            throw new RuntimeException(e);
                        }
                        readers[stoppingPoint] = new BufferedReader(new InputStreamReader(processes[stoppingPoint].getInputStream()));

                        jobIdx.incrementAndGet();
                    }
                }));
            }
            ConcurrentHashMap<Integer, Integer> completedTasks = new ConcurrentHashMap<>();
            while (completedTasks.size() < futureList.size()) {
                for (int i = 0; i < jobIdx.get(); i++) {
                    if (!completedTasks.containsKey(i) && futureList.get(i).isDone()) {
                        futureList.get(i).get();
                        completedTasks.put(i, 0);
                    }
                }
                Thread.sleep(100);
            }

            long endTime = System.currentTimeMillis();
            long elapsedTime = endTime - startTime;
            System.out.println("cURL and parse time in milliseconds: " + elapsedTime);

        } catch (Exception e) {
//            System.out.println(); //TODO print error message from bufferedreader
            e.printStackTrace();
        }
    }

    public void parseKoinaOutput(String fileName, String koinaString, String property, String model,
                                        KoinaLibReader klr) throws IOException {
        if (property.toLowerCase().equals("rt")) {
            String rts = koinaString.split("data")[2];
            String[] results = rts.substring(3, rts.length() - 4).split(",");
            float[] parsedResults = new float[results.length];
            for (int i = 0; i < results.length; i++) {
                parsedResults[i] = Float.parseFloat(results[i]);
            }

            assignRTs(fileName, parsedResults, klr);
        } else if (property.toLowerCase().equals("ms2")) {
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

            //mz
            msInfo = dataResults[mzIdx].split("data\":\\[")[1];
            msInfo = msInfo.substring(0, msInfo.length() - 1);
            results = msInfo.split(",");
            float[][] allMZs = new float[numPeptides][];
            for (int i = 0; i < numPeptides; i++) {
                ArrayList<Float> mz = new ArrayList<>();
                for (int j = i * vectorLength; j < (i + 1) * vectorLength; j++) {
                    String result = results[j];
                    if (acceptedIdx[i].contains(j)) {
                        mz.add(Float.parseFloat(result));
                    }
                }
                float[] mzArray = new float[mz.size()];
                for (int j = 0; j < mz.size(); j++) {
                    mzArray[j] = mz.get(j);
                }
                allMZs[i] = mzArray;
            }

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

            assignMS2(fileName, allMZs, allIntensities, allFragmentIonTypes, allFragNums, allCharges, klr);
        }
    }

    private void assignRTs(String fileName, float[] RTs, KoinaLibReader klr) throws IOException {
        String[] peptides = readJSON(fileName, RTs.length);
        ConcurrentHashMap<String, PredictionEntry> preds = klr.getPreds();
        for (int i = 0; i < peptides.length; i++) {
            PeptideFormatter pf = new PeptideFormatter(peptides[i], 1, "diann");
            int entries = 0; //in case RT prediction is available but not MS2
            for (int charge = Constants.minPrecursorCharge; charge < Constants.maxPrecursorCharge + 1; charge++) {
                String peptide = pf.getBase() + "|" + charge;
                if (modelType.equals("prosit")) {
                    peptide = peptide.replace("C[" + PTMhandler.carbamidomethylationMass + "]", "C");
                    peptide = peptide.replace("C", "C[" + PTMhandler.carbamidomethylationMass + "]");
                }
                if (preds.containsKey(peptide)) {
                    PredictionEntry pe = preds.get(peptide);
                    pe.setRT(RTs[i]);
                    preds.put(peptide, pe);
                    entries++;
                }
            }
            if (entries == 0) { //RT was predicted but not MS2
                for (int charge = Constants.minPrecursorCharge; charge < Constants.maxPrecursorCharge + 1; charge++) {
                    String peptide = pf.getBase() + "|" + charge;
                    if (modelType.equals("prosit")) {
                        peptide = peptide.replace("C[" + PTMhandler.carbamidomethylationMass + "]", "C");
                        peptide = peptide.replace("C", "C[" + PTMhandler.carbamidomethylationMass + "]");
                    }
                    PredictionEntry pe = new PredictionEntry();
                    pe.setRT(RTs[i]);
                    preds.put(peptide, pe);
                }
            }
        }
    }

    private void assignMS2(String fileName, float[][] mzs, float[][] intensities,
                           String[][] fragmentIonTypes, int[][] fragNums, int[][] charges, KoinaLibReader klr)
            throws IOException {
        String[] peptides = readJSON(fileName, mzs.length);
        ConcurrentHashMap<String, PredictionEntry> preds = klr.getPreds();
        for (int i = 0; i < peptides.length; i++) {
            PeptideFormatter pf = new PeptideFormatter(peptides[i].split("\\|")[0],
                    peptides[i].split("\\|")[1], "diann");
            String peptide = pf.getBaseCharge();
            if (modelType.equals("prosit")) {
                peptide = peptide.replace("C[" + PTMhandler.carbamidomethylationMass + "]", "C");
                peptide = peptide.replace("C", "C[" + PTMhandler.carbamidomethylationMass + "]");
            }
            PredictionEntry pe = new PredictionEntry(mzs[i], intensities[i], fragNums[i],
                    charges[i], fragmentIonTypes[i], true);

            if (preds.containsKey(peptide)) {
                pe.setRT(preds.get(peptide).getRT());
            }
            preds.put(peptide, pe);
        }
    }

    public void assignMissingPeptidePredictions(KoinaLibReader klr) throws IOException {
        BufferedReader TSVReader = new BufferedReader(new FileReader(
                Constants.spectraRTPredInput.substring(0, Constants.spectraRTPredInput.length() - 4) + "_full.tsv"));
        String l;
        String[] line;
        ConcurrentHashMap<String, PredictionEntry> preds = klr.getPreds();

        while ((l = TSVReader.readLine()) != null) {
            line = l.split("\t");
            if (!preds.containsKey(line[0] + "|" + line[1])) {
                //get predictionEntry
                PeptideFormatter pf;
                PredictionEntry tmp = new PredictionEntry();
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
                        pf = new PeptideFormatter(
                                new PeptideFormatter(line[0], line[1], "base").getProsit(), line[1], "prosit");
                        baseCharge = pf.getBaseCharge();
                        tmp = preds.get(baseCharge);
                        stripped = pf.getStripped();
                        break;
                    default:
                        System.out.println(modelType + " not supported by Koina");
                        System.exit(1);
                }
                try {
                    MassCalculator mc = new MassCalculator(line[0], line[1]);
                    float[] newMZs = new float[tmp.getMzs().length];
                    for (int i = 0; i < newMZs.length; i++) {
                        newMZs[i] = mc.calcMass(tmp.getFragNums()[i],
                                Constants.flagTOion.get(tmp.getFlags()[i]), tmp.getCharges()[i]);
                    }

                    //add to hashmap
                    PredictionEntry newPred = new PredictionEntry();
                    newPred.setMzs(newMZs);
                    newPred.setIntensities(tmp.getIntensities());
                    newPred.setRT(tmp.getRT());
                    newPred.setIM(tmp.getIM());
                    newPred.setFragmentIonTypes(tmp.getFragmentIonTypes());
                    newPred.setFragNums(tmp.getFragNums());
                    newPred.setFlags(tmp.getFlags());
                    preds.put(mc.fullPeptide, newPred);
                } catch (Exception e) {
                    if (! PeptideSkipper.skipPeptide(stripped, baseCharge.split("\\|")[1])) {
                        e.printStackTrace();
                        System.out.println("Missing peptide to transfer prediction onto " + l + ": " + baseCharge);
                        System.out.println("Exiting now.");
                        System.exit(1);
                    }
                }
            }
        }
    }

    private String[] readJSON(String fileName, int length) throws IOException {
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
