package External;

import Features.Constants;

import java.io.BufferedReader;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;

public class AlphaPeptDeepModelCaller {
    public static void callModel(String mode) {
        long startTime = System.nanoTime();
        try {
            if (mode.equals("transfer")) {
                System.out.println("Generating AlphaPeptDeep transfer learned model");
                ProcessBuilder builder = new ProcessBuilder(Constants.AlphaPeptDeep,
                        "transfer", Constants.paramsList);
                System.out.println(String.join(" ", builder.command()));
                builder.redirectErrorStream(true);
                Process process = builder.start();
                InputStream is = process.getInputStream();
                BufferedReader reader = new BufferedReader(new InputStreamReader(is));

                String line = null;
                while ((line = reader.readLine()) != null) {
                    System.out.println(line);
                }
            } else if (mode.equals("predict")) {
                System.out.println("Generating AlphaPeptDeep predictions");
                ProcessBuilder builder = new ProcessBuilder(Constants.AlphaPeptDeep,
                        "predict", Constants.paramsList);
                System.out.println(String.join(" ", builder.command()));
                builder.redirectErrorStream(true);
                Process process = builder.start();
                InputStream is = process.getInputStream();
                BufferedReader reader = new BufferedReader(new InputStreamReader(is));

                String line = null;
                while ((line = reader.readLine()) != null) {
                    System.out.println(line);
                }
            }
        } catch (IOException e) {
            e.printStackTrace();
        }
        long endTime = System.nanoTime();
        long duration = (endTime - startTime);
        System.out.println("Model running took " + duration / 1000000 +" milliseconds");
    }
}
