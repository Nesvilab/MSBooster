package External;

import Features.Constants;

import java.io.*;
import java.nio.file.Files;

public class DiannModelCaller {
    public static void callModel() {
        long startTime = System.nanoTime();
        try {
            boolean useTMT = false;
            //DIA-NN command
            System.out.println("Generating DIA-NN predictions");
            Constants.spectraRTPredFile =
                    Constants.spectraRTPredInput.substring(0, Constants.spectraRTPredInput.length() - 4) +
                            ".predicted.bin";
            String line;

            //check for TMT
            BufferedReader br = new BufferedReader(new FileReader(Constants.spectraRTPredInput));
            while ((line = br.readLine()) != null) {
                if (line.contains("[TMT]")) {
                    useTMT = true;
                    break;
                }
            }
            br.close();

            //get num PSMs
            int subFileSize = 0;
            int linenumTotal = 0;
            if (Constants.splitPredInputFile != 1) {
                br = new BufferedReader(new FileReader(Constants.spectraRTPredInput));
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

                    br = new BufferedReader(new FileReader(Constants.spectraRTPredInput));
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
                ProcessBuilder builder;
                if (useTMT) {
                    builder = new ProcessBuilder(Constants.DiaNN,
                            "--lib",
                            inputString,
                            "--predict",
                            "--threads",
                            String.valueOf(Constants.numThreads),
                            "--strip-unknown-mods",
                            "--full-unimod",
                            "--predict-n-frag",
                            "100",
                            "--mod",
                            "TMT,229.1629");
                } else {
                    builder = new ProcessBuilder(Constants.DiaNN,
                            "--lib",
                            inputString,
                            "--predict",
                            "--threads",
                            String.valueOf(Constants.numThreads),
                            "--strip-unknown-mods",
                            "--full-unimod",
                            "--predict-n-frag",
                            "100");
                }
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
        long endTime = System.nanoTime();
        long duration = (endTime - startTime);
        System.out.println("Model running took " + duration / 1000000 +" milliseconds");
    }
}
