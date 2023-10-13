package External;

import Features.Constants;

import java.io.*;

public class DeepMSPeptideModelCaller {
    public static void callModel() {
        long startTime = System.nanoTime();
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
        long endTime = System.nanoTime();
        long duration = (endTime - startTime);
        System.out.println("Model running took " + duration / 1000000 +" milliseconds");
    }
}
