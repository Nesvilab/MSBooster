package External;

import Features.Constants;
import com.ctc.wstx.shaded.msv_core.verifier.jarv.Const;

import java.io.*;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
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

    public static void callModel(String model) throws IOException {
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
            int numProcesses = 0;
            for (File file : fileArray) {
                String fileString = file.toString();
                if (fileString.endsWith(".json")) {
                    numProcesses += 1;
                }
            }

            Process[] processes = new Process[numProcesses];
            BufferedReader[] readers = new BufferedReader[numProcesses];
            numProcesses = 0;

            for (File file : fileArray) {
                String fileString = file.toString();
                if (fileString.endsWith(".json")) {
                    String command = "curl -s --parallel --parallel-immediate --parallel-max 36 " +
                            "-H content-type:application/json -d @" + fileString +
                            " https://koina.proteomicsdb.org/v2/models/" + model + "/infer";

                    ProcessBuilder builder = new ProcessBuilder(command.split(" "));
                    builder.redirectErrorStream(true);
                    processes[numProcesses] = builder.start();
                    readers[numProcesses] = new BufferedReader(new InputStreamReader(processes[numProcesses].getInputStream()));
                    numProcesses += 1;
                }
            }

            List<Future> futureList = new ArrayList<>(Constants.numThreads);
            ExecutorService executorService = Executors.newFixedThreadPool(Constants.numThreads);
            for (int j = 0; j < Constants.numThreads; j++) {
                int start = (int) (numProcesses * (long) j) / Constants.numThreads;
                int end = (int) (numProcesses * (long) (j + 1)) / Constants.numThreads;
                String finalProperty = property;
                futureList.add(executorService.submit(() -> {
                    for (int i = start; i < end; i++) {
                        Process p = processes[i];
                        BufferedReader reader = readers[i];
                        StringBuilder koinaSb = new StringBuilder();

                        String line = "";
                        while (true) {
                            try {
                                if ((line = reader.readLine()) == null) break;
                            } catch (IOException e) {
                                e.printStackTrace();
                            }
                            koinaSb.append(line);
                        }
                        try {
                            reader.close();
                        } catch (IOException e) {
                            e.printStackTrace();
                        }
                        try {
                            p.waitFor();
                        } catch (InterruptedException e) {
                            e.printStackTrace();
                        }
                        p.destroy();

                        //TODO parse
                        KoinaModelCaller.parseKoinaOutput(koinaSb.toString(), finalProperty, model);
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

    public static void parseKoinaOutput(String koinaString, String property, String model) {
        if (property.toLowerCase().equals("rt")) {
            String rts = koinaString.split("data")[2];
            String[] results = rts.substring(3, rts.length() - 4).split(",");
            float[] parsedResults = new float[results.length];
            for (int i = 0; i < results.length; i++) {
                parsedResults[i] = Float.parseFloat(results[i]);
            }
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
            System.out.println(numPeptides);

//            msInfo = msInfo.substring(3, msInfo.length() - 3);
//            String[] dataResults = msInfo.split("},");
//
//            //intensities
//            msInfo = dataResults[intIdx].split("data\":\\[")[1];
//            msInfo = msInfo.substring(0, msInfo.length() - 1);
//            String[] results = msInfo.split(",");
//            int vectorLength = results.length / numPeptides;
//            for (int i = 0; i < numPeptides; i++) {
//                ArrayList<Float> intensities = new ArrayList<>();
//                for (int j = i * vectorLength; j < (i + 1) * vectorLength; j++) {
//                    String result = results[j];
//                    if (result.equals("-1.0")) {
//                        break;
//                    } else {
//                        intensities.add(Float.parseFloat(result));
//                    }
//                }
//            }
//
//            //mz
//            rawString = dataResults[mzIdx].split("data\":\\[")[1];
//            rawString = rawString.substring(0, rawString.length() - 1);
//            results = rawString.split(",");
//            vectorLength = results.length / numPeptides;
//            for (int i = 0; i < numPeptides; i++) {
//                ArrayList<Float> mz = new ArrayList<>();
//                for (int j = i * vectorLength; j < (i + 1) * vectorLength; j++) {
//                    String result = results[j];
//                    if (result.equals("-1.0")) {
//                        break;
//                    } else {
//                        mz.add(Float.parseFloat(result));
//                    }
//                }
//            }
//
//            endTime = System.currentTimeMillis();
//            elapsedTime = endTime - startTime;
//            if (rep == 0) {
//                System.out.println("Parsing time in milliseconds: " + elapsedTime);
//            }
        }
    }
}
