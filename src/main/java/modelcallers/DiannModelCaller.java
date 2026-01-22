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

package modelcallers;

import allconstants.Constants;
import peptideptmformatting.PTMhandler;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.io.RandomAccessFile;
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

                    //concatenate files together
                    //adapted from https://stackoverflow.com/questions/2243073/java-multiple-connection-downloading/2243731#2243731
                    int data = 0;
                    try {
                        File filename = new File(predFileString + ".total");
                        FileWriter outfile = new FileWriter(filename, true);

                        filename = new File(predFileString);
                        RandomAccessFile infile = new RandomAccessFile(filename, "r");
                        data = infile.read();
                        while (data != -1) {
                            outfile.write(data);
                            data = infile.read();
                        }
                        infile.close();
                        outfile.close();
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
