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

import Features.Constants;

import java.io.*;
import java.nio.file.Files;
import java.util.concurrent.TimeUnit;

public class ExternalModelCaller {
    //TODO: can we make the repeated parts more concise?
    public static void callModel(String model, String mode) {
        long startTime = System.nanoTime();
        switch (model) {
            case "DIA-NN":
                try {
                    //DIA-NN command
                    System.out.println("Generating DIA-NN predictions");
                    Constants.spectraRTPredFile =
                            Constants.spectraRTPredInput.substring(0, Constants.spectraRTPredInput.length() - 4) +
                                    ".predicted.bin";
                    String line;

                    //get num PSMs
                    int subFileSize = 0;
                    int linenumTotal = 0;
                    if (Constants.splitPredInputFile != 1) {
                        BufferedReader br = new BufferedReader(new FileReader(Constants.spectraRTPredInput));
                        int linenum = -1;
                        while ((line = br.readLine()) != null) {
                            linenum += 1;
                        }
                        br.close();

                        subFileSize = linenum / Constants.splitPredInputFile;
                        linenumTotal = linenum;
                    }

                    for (int i = 1; i < Constants.splitPredInputFile + 1; i++) {
                        String inputString = Constants.spectraRTPredInput;

                        //splitting in case large input file
                        if (Constants.splitPredInputFile != 1) {
                            //get new input string name
                            inputString += i;

                            //get row splits to go in each file
                            int startRow = (i - 1) * subFileSize;
                            int endRow = i * subFileSize;
                            if (i == Constants.splitPredInputFile) {
                                endRow = linenumTotal;
                            }

                            BufferedReader br = new BufferedReader(new FileReader(Constants.spectraRTPredInput));
                            line = br.readLine();
                            FileWriter myWriter = new FileWriter(inputString);
                            myWriter.write(line + "\n");
                            int linenum = 0;
                            while ((line = br.readLine()) != null) {
                                if (linenum >= startRow && linenum < endRow) {
                                    myWriter.write(line + "\n");
                                }
                                linenum += 1;
                            }
                            myWriter.close();
                            br.close();
                        }

                        //actual prediction
                        ProcessBuilder builder = new ProcessBuilder(Constants.DiaNN,
                                "--lib",
                                inputString,
                                "--predict",
                                "--threads",
                                String.valueOf(Constants.numThreads),
                                "--strip-unknown-mods",
                                "--full-unimod",
                                "--predict-n-frag",
                                "100");
                        System.out.println(String.join(" ", builder.command()));
                        builder.redirectErrorStream(true);
                        Process process = builder.start();
                        InputStream is = process.getInputStream();
                        BufferedReader reader = new BufferedReader(new InputStreamReader(is));

                        //print DIA-NN output while running
                        while ((line = reader.readLine()) != null) {
                            System.out.println(line);
                        }

                        int DIANNtermination = process.waitFor();

                        if (DIANNtermination == -1073741515) {
                            System.out.println("Microsoft Visual C++ Redistributable is missing. Please download at " +
                                    "https://learn.microsoft.com/en-us/cpp/windows/latest-supported-vc-redist");
                            System.exit(-1);
                        }
                        if (DIANNtermination == 137) {
                            System.out.println("Out of memory during DIA-NN prediction. " +
                                    "Please allocate more memory, or increase splitPredInputFile " +
                                    "parameter until successfully predicted.");
                            System.exit(-1);
                        }
                        if (DIANNtermination != 0) {
                            System.out.println("Abnormal DIANN termination: " + DIANNtermination + ", please run the " +
                                    "following command from the command line for more information\n" +
                                    String.join(" ", builder.command()));
                            System.exit(-1);
                        }

                        if (Constants.splitPredInputFile != 1) {
                            File inputFile = new File(inputString);
                            inputFile.delete();

                            //concatenate files together
                            //adapted from https://stackoverflow.com/questions/2243073/java-multiple-connection-downloading/2243731#2243731
                            int data = 0;
                            try {
                                File filename = new File(Constants.spectraRTPredFile + ".total");
                                FileWriter outfile = new FileWriter(filename, true);

                                filename = new File(Constants.spectraRTPredFile);
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

                    File predFile = new File(Constants.spectraRTPredFile);
                    //move total file to typical name, if total file exists
                    if (Constants.splitPredInputFile != 1) {
                        predFile.delete();
                        File oldFile = new File(Constants.spectraRTPredFile + ".total");
                        oldFile.renameTo(predFile);
                    }

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
            case "alphapeptdeep":
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
        }
        long endTime = System.nanoTime();
        long duration = (endTime - startTime);
        System.out.println("Model running took " + duration / 1000000 +" milliseconds");
    }
}
