package External;

import Features.Constants;

import java.io.*;

public class ExternalModelCaller {
    //TODO: can we make the repeated parts more concise?
    public static void callModel(Runtime run, String model) {
        long startTime = System.nanoTime();
        switch (model) {
            case "DIA-NN":
                try {
                    System.out.println("Generating DIA-NN predictions");
                    ProcessBuilder builder = new ProcessBuilder(Constants.DiaNN,
                            "--lib",
                            Constants.spectraRTPredInput,
                            "--predict",
                            "--threads",
                            String.valueOf(Constants.numThreads));
                    System.out.println(String.join(" ", builder.command()));
                    builder.redirectErrorStream(true);
                    Process process = builder.start();
                    InputStream is = process.getInputStream();
                    BufferedReader reader = new BufferedReader(new InputStreamReader(is));

                    String line = null;
                    while ((line = reader.readLine()) != null) {
                        System.out.println(line);
                    }

                    Constants.spectraRTPredFile =
                            Constants.spectraRTPredInput.substring(0, Constants.spectraRTPredInput.length() - 4) +
                            ".predicted.bin";
                    System.out.println("Done generating DIA-NN predictions");
                } catch (IOException e) {
                    e.printStackTrace();
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
                }
                break;
        }
        long endTime = System.nanoTime();
        long duration = (endTime - startTime);
        System.out.println("Model running took " + duration / 1000000 +" milliseconds");
    }
}
