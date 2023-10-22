package External;

import Features.Constants;
import Features.KoinaLibReader;
import Features.PredictionEntry;
import Features.ProgressReporter;
import com.google.gson.Gson;
import com.google.gson.stream.JsonReader;

import java.io.*;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.concurrent.ConcurrentHashMap;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.Future;

public class KoinaModelCaller {
    private static final int AlphaPeptDeepMzIdx = 1;
    private static final int AlphaPeptDeepIntIdx = 0;
    private static final int PrositMzIdx = 1;
    private static final int PrositIntIdx = 2;
    private static final int ms2pipMzIdx = 1;
    private static final int ms2pipIntIdx = 2;

    public static void callModel(String model, KoinaLibReader klr) {
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
            int numProcesses = 0;

            for (String fileString : filenameArray) {
                String command = "curl -s --parallel --parallel-immediate --parallel-max 36 " +
                        "-H content-type:application/json -d @" + fileString +
                        " https://koina.proteomicsdb.org/v2/models/" + model + "/infer";

                ProcessBuilder builder = new ProcessBuilder(command.split(" "));
                builder.redirectErrorStream(true);
                processes[numProcesses] = builder.start();
                readers[numProcesses] = new BufferedReader(new InputStreamReader(processes[numProcesses].getInputStream()));
                numProcesses += 1;
            }

            List<Future> futureList = new ArrayList<>(Constants.numThreads);
            ExecutorService executorService = Executors.newFixedThreadPool(Constants.numThreads);
            for (int j = 0; j < Constants.numThreads; j++) {
                int start = (int) (numProcesses * (long) j) / Constants.numThreads;
                int end = (int) (numProcesses * (long) (j + 1)) / Constants.numThreads;
                String finalProperty = property;
                futureList.add(executorService.submit(() -> {
                    ProgressReporter pr = new ProgressReporter(end - start);
                    for (int i = start; i < end; i++) {
                        int attempts = 0;
                        while (attempts < 3) {
                            StringBuilder koinaSb = new StringBuilder();
                            try {
                                Process p = processes[i];
                                BufferedReader reader = readers[i];

                                String line = "";
                                while ((line = reader.readLine()) != null) {
                                    koinaSb.append(line);
                                }
                                reader.close();
                                p.waitFor();
                                p.destroy();
                            } catch (Exception e) {
                                e.printStackTrace();
                            }

                            try {
                                KoinaModelCaller.parseKoinaOutput(filenameArray[i], koinaSb.toString(),
                                        finalProperty, model, klr);
                                break;
                            } catch (Exception e) {
                                attempts++;
                                if (attempts == 3) {
                                    System.out.println(filenameArray[i] + " had output that ended in: ");
                                    System.out.println(koinaSb.toString().substring(Math.max(0, koinaSb.toString().length() - 1000)));
                                    System.out.println("Retried calling " + filenameArray[i] + " " + attempts +
                                            " times. Moving on.");
                                    e.printStackTrace();
                                    break;
                                }

                                //set up again
                                String command = "curl -s --parallel --parallel-immediate --parallel-max 36 " +
                                        "-H content-type:application/json -d @" + filenameArray[i] +
                                        " https://koina.proteomicsdb.org/v2/models/" + model + "/infer";

                                ProcessBuilder builder = new ProcessBuilder(command.split(" "));
                                builder.redirectErrorStream(true);
                                try {
                                    processes[i] = builder.start();
                                } catch (IOException ioException) {
                                    ioException.printStackTrace();
                                }
                                readers[i] = new BufferedReader(new InputStreamReader(processes[i].getInputStream()));
                            }
                        }
                        if (start == 0) {
                            pr.progress();
                        }
                    }
                }));
            }
            for (Future future : futureList) {
                future.get();
            }

            long endTime = System.currentTimeMillis();
            long elapsedTime = endTime - startTime;
            System.out.println("cURL and parse time in milliseconds: " + elapsedTime);

        } catch (Exception e) {
            System.out.println(); //TODO print error message from bufferedreader
            e.printStackTrace();
        }
    }

    public static void parseKoinaOutput(String fileName, String koinaString, String property, String model,
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
            if (model.contains("AlphaPept")) {
                mzIdx = AlphaPeptDeepMzIdx;
                intIdx = AlphaPeptDeepIntIdx;
            } else if (model.contains("Prosit")) {
                mzIdx = PrositMzIdx;
                intIdx = PrositIntIdx;
            } else if (model.contains("ms2pip")) {
                mzIdx = ms2pipMzIdx;
                intIdx = ms2pipIntIdx;
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
            float[][] allIntensities = new float[numPeptides][];
            for (int i = 0; i < numPeptides; i++) {
                ArrayList<Float> intensities = new ArrayList<>();
                for (int j = i * vectorLength; j < (i + 1) * vectorLength; j++) {
                    String result = results[j];
                    if (result.equals("-1.0")) {
                        break;
                    } else {
                        intensities.add(Float.parseFloat(result));
                    }
                }
                float[] intensitiesArray = new float[intensities.size()];
                for (int j = 0; j < intensities.size(); j++) {
                    intensitiesArray[j] = intensities.get(j);
                }
                allIntensities[i] = intensitiesArray;
            }

            //mz
            msInfo = dataResults[mzIdx].split("data\":\\[")[1];
            msInfo = msInfo.substring(0, msInfo.length() - 1);
            results = msInfo.split(",");
            vectorLength = results.length / numPeptides;
            float[][] allMZs = new float[numPeptides][];
            for (int i = 0; i < numPeptides; i++) {
                ArrayList<Float> mz = new ArrayList<>();
                for (int j = i * vectorLength; j < (i + 1) * vectorLength; j++) {
                    String result = results[j];
                    if (result.equals("-1.0")) {
                        break;
                    } else {
                        mz.add(Float.parseFloat(result));
                    }
                }
                float[] mzArray = new float[mz.size()];
                for (int j = 0; j < mz.size(); j++) {
                    mzArray[j] = mz.get(j);
                }
                allMZs[i] = mzArray;
            }

            assignMS2(fileName, allMZs, allIntensities, klr);
        }
    }

    private static void assignRTs(String fileName, float[] RTs, KoinaLibReader klr) throws IOException {
        String[] peptides = readJSON(fileName, RTs.length);
        ConcurrentHashMap<String, PredictionEntry> preds = klr.getPreds();
        for (int i = 0; i < peptides.length; i++) {
            String peptide = peptides[i];
            for (int charge = Constants.minPrecursorCharge; charge < Constants.maxPrecursorCharge + 1; charge++) {
                if (preds.containsKey(peptide + "|" + charge)) {
                    PredictionEntry pe = preds.get(peptide + "|" + charge);
                    pe.setRT(RTs[i]);
                    preds.put(peptide + "|" + charge, pe);
                }
            }
        }
    }

    private static void assignMS2(String fileName, float[][] mzs, float[][] intensities, KoinaLibReader klr)
            throws IOException {
        String[] peptides = readJSON(fileName, mzs.length);
        ConcurrentHashMap<String, PredictionEntry> preds = klr.getPreds();
        for (int i = 0; i < peptides.length; i++) {
            String peptide = peptides[i];
            PredictionEntry pe;
            if (preds.containsKey(peptide)) {
                pe = preds.get(peptide);
            } else {
                pe = new PredictionEntry();
            }
            pe.setMzs(mzs[i]);
            pe.setIntensities(intensities[i]);
            preds.put(peptide, pe);
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
