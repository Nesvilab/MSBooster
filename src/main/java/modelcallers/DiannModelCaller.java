/*
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or
 * implied. See the License for the specific language governing
 * permissions and limitations under the License.
 */

package modelcallers;

import allconstants.Constants;
import peptideptmformatting.PTMhandler;

import java.io.*;
import java.nio.file.Files;

import static utils.Print.printError;
import static utils.Print.printInfo;

public class DiannModelCaller {
    static boolean retry = false;
    public static String callModel(String inputFile, boolean verbose) {
        long startTime = System.nanoTime();
        String predFileString = null;
        try {
            boolean useTMT = false;
            //DIA-NN command
            if (verbose) {
                printInfo("Generating DIA-NN predictions");
            }
            predFileString = inputFile.substring(0, inputFile.length() - 4) + ".predicted.bin";
            String line;

            //check for TMT
            BufferedReader br = new BufferedReader(new FileReader(inputFile));
            while ((line = br.readLine()) != null) {
                if (line.contains("[TMT]")) {
                    useTMT = true;
                    break;
                }
            }
            br.close();

            for (int i = 0; i < Constants.splitPredInputFile; i++) {
                String inputString = inputFile;
                if (Constants.splitPredInputFile > 1) {
                    inputString += i;
                    printInfo("Predicting batch " + (i + 1));
                }

                //actual prediction
                ProcessBuilder builder;
                if (useTMT) {
                    builder = new ProcessBuilder(Constants.DiaNN,
                            "--lib",
                            inputString,
                            "--predict",
                            "--threads",
                            String.valueOf(Constants.numThreads),
                            "--strip-unknown-mods",
                            "--predict-n-frag",
                            "100",
                            "--mod",
                            //"TMT,229.1629",
                            "TMT," + PTMhandler.getTmtMass(),
                            "--original-mods");
                } else {
                    builder = new ProcessBuilder(Constants.DiaNN,
                            "--lib",
                            inputString,
                            "--predict",
                            "--threads",
                            String.valueOf(Constants.numThreads),
                            "--strip-unknown-mods",
                            "--predict-n-frag",
                            "100");
                }
                if (verbose) {
                    printInfo(String.join(" ", builder.command()));
                }
                builder.redirectErrorStream(true);
                Process process = builder.start();
                InputStream is = process.getInputStream();
                BufferedReader reader = new BufferedReader(new InputStreamReader(is));

                //print DIA-NN output while running
                while ((line = reader.readLine()) != null) {
                    if (verbose) {
                        printInfo(line);
                    }
                }

                int DIANNtermination = process.waitFor();

                if (DIANNtermination == -1073741515) {
                    printError("Microsoft Visual C++ Redistributable is missing. Please download at " +
                            "https://learn.microsoft.com/en-us/cpp/windows/latest-supported-vc-redist");
                    System.exit(1);
                }
                if (DIANNtermination == 137) {
                    printError("Out of memory during DIA-NN prediction. " +
                            "Please allocate more memory, or increase splitPredInputFile " +
                            "parameter until successfully predicted.");
                    System.exit(1);
                }
                if (DIANNtermination == -1073741819) {
                    if (retry) {
                        retry = false;
                        printError("Encountered segmentation fault/access violation.");
                        System.exit(1);
                    }
                    printError("Encountered segmentation fault/access violation. Retrying.");
                    return callModelRetry(inputFile, verbose);
                }
                if (DIANNtermination != 0) {
                    printError("Abnormal DIANN termination: " + DIANNtermination + ", please run the " +
                            "following command from the command line for more information\n" +
                            String.join(" ", builder.command()));
                    System.exit(1);
                }

                if (Constants.splitPredInputFile != 1) {
                    File inputf = new File(inputString);
                    inputf.delete();

                    // Concatenate binary files with buffering
                    try (FileOutputStream outfile = new FileOutputStream(predFileString + ".total", true);
                         BufferedOutputStream bos = new BufferedOutputStream(outfile)) {

                        File filename = new File(predFileString);

                        try (FileInputStream infile = new FileInputStream(filename);
                             BufferedInputStream bis = new BufferedInputStream(infile)) {

                            byte[] buffer = new byte[8192];  // 8KB buffer
                            int bytesRead;
                            while ((bytesRead = bis.read(buffer)) != -1) {
                                bos.write(buffer, 0, bytesRead);
                            }
                        }
                    } catch (Exception e) {
                        e.printStackTrace();
                    }
                }
            }

            File predFile = new File(predFileString);
            //move total file to typical name, if total file exists
            if (Constants.splitPredInputFile != 1) {
                predFile.delete();
                File oldFile = new File(predFileString + ".total");
                oldFile.renameTo(predFile);
            }

            if (Files.isReadable(predFile.toPath())) {
                if (verbose) {
                    printInfo("Done generating DIA-NN predictions");
                }
            } else {
                printError("Cannot find DIA-NN's output. Please rerun MSBooster");
                System.exit(1);
            }

        } catch (IOException | InterruptedException e) {
            e.printStackTrace();
            System.exit(1);
        }
        if (verbose) {
            long endTime = System.nanoTime();
            long duration = (endTime - startTime);
            printInfo("Model running took " + duration / 1000000 + " milliseconds");
        }
        
        return predFileString;
    }

    private static String callModelRetry(String inputFile, boolean verbose) {
        retry = true;
        String result = callModel(inputFile, verbose);
        retry = false;
        return result;
    }
}
