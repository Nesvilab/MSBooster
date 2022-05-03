package External;

import Features.Constants;

import java.io.*;
import java.nio.file.Files;

public class ExternalModelCaller {
    //TODO: can we make the repeated parts more concise?
    public static void callModel(Runtime run, String model) {
        long startTime = System.nanoTime();
        switch (model) {
            case "DIA-NN":
                try {
                    //DIA-NN command
                    System.out.println("Generating DIA-NN predictions");
                    ProcessBuilder builder = new ProcessBuilder(Constants.DiaNN,
                            "--lib",
                            Constants.spectraRTPredInput,
                            "--predict",
                            "--threads",
                            String.valueOf(Constants.numThreads),
                            "--strip-unknown-mods");
                    System.out.println(String.join(" ", builder.command()));
                    builder.redirectErrorStream(true);
                    Process process = builder.start();
                    InputStream is = process.getInputStream();
                    BufferedReader reader = new BufferedReader(new InputStreamReader(is));

                    //print DIA-NN output while running
                    String line = null;

                    //limit time
                    long start = System.currentTimeMillis();
                    long end = start + Constants.predictionTimeLimit * 60 * 1000; //stop after 2 hours

                    while ((! (line = reader.readLine()).contains("Finished")) &&
                            System.currentTimeMillis() < end) {
                        if (line.equals("")) {
                            System.out.println("null");
                        } else {
                            System.out.println(line);
                        } //wait for the last line to print instead of using waitfor?
                        //also check if readline gives null at some point
                    }
                    System.out.println(line);

                    Constants.spectraRTPredFile =
                            Constants.spectraRTPredInput.substring(0, Constants.spectraRTPredInput.length() - 4) +
                            ".predicted.bin";
                    int DIANNtermination = process.waitFor();

                    if (DIANNtermination != 0) {
                        System.out.println("Abnormal DIANN termination: " + DIANNtermination);
                    }
                    File predFile = new File(Constants.spectraRTPredFile);
                    if (Files.isReadable(predFile.toPath())) {
                        System.out.println("Done generating DIA-NN predictions");
                    } else {
                        System.out.println("Cannot find DIA-NN's output. Please rerun MSBooster");
                        System.exit(-1);
                    }

                } catch (IOException | InterruptedException e) {
                    e.printStackTrace();
                    System.exit(1);
                }
                break;
            case "DeepMSPeptide":
                try {
                    System.out.println("Generating DeepMSPeptide predictions");
                    ProcessBuilder builder = new ProcessBuilder("." + File.separator + "DeepMSPeptideRevised.exe",
                            Constants.detectPredInput);
                    System.out.println(String.join(" ", builder.command()));
                    builder.redirectErrorStream(true);
                    Process process = builder.start();
                    InputStream is = process.getInputStream();
                    BufferedReader reader = new BufferedReader(new InputStreamReader(is));

                    String line = null;
                    while ((line = reader.readLine()) != null) {
                        System.out.println(line);
                    }

                    Constants.detectPredFile =
                            Constants.detectPredInput.substring(0, Constants.detectPredInput.length() - 4) +
                                    "_Predictions.txt";
                    System.out.println("Done generating DeepMSPeptide predictions");
                } catch (IOException e) {
                    e.printStackTrace();
                    System.exit(1);
                }
                break;
        }
        long endTime = System.nanoTime();
        long duration = (endTime - startTime);
        System.out.println("Model running took " + duration / 1000000 +" milliseconds");
    }
}
