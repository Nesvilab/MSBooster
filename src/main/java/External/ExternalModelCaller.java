package External;

import Features.Constants;

import java.io.IOException;

public class ExternalModelCaller {
    public static void callModel(Runtime run, String model) {
        switch (model) {
            case "DIA-NN":
                try {
                    System.out.println("Generating DIA-NN predictions");
                    // Command to create an external process
                    String command = Constants.DiaNN + " --lib " + Constants.spectraRTPredInput + " --predict --threads ";

                    // Running the above command
                    Process proc = run.exec(command + Constants.numThreads);
                    proc.waitFor();
                    Constants.spectraRTPredFile =
                            Constants.spectraRTPredInput.substring(0, Constants.spectraRTPredInput.length() - 4) +
                            ".predicted.bin";
                    System.out.println("Done generating DIA-NN predictions");
                } catch (IOException | InterruptedException e) {
                    e.printStackTrace();
                }
                break;
        }
    }
}
