package transferlearn;

import allconstants.Constants;
import peptideptmformatting.PTMhandler;
import readers.datareaders.FastaReader;
import speclib.io.ParquetToSpecLib;
import utils.Print;

import java.io.*;
import java.net.HttpURLConnection;
import java.net.URL;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.sql.Connection;
import java.sql.DriverManager;
import java.sql.Statement;
import java.util.Arrays;
import java.util.HashMap;
import java.util.Locale;
import java.util.stream.Collectors;

import static allconstants.Constants.versionNumber;
import static peptideptmformatting.PTMhandler.*;
import static transferlearn.Helpers.*;
import static utils.Print.printInfo;

public class Predictor {
    static Long waitTime = 5000L;

    private static void errorMessage() {
        Print.printError("Usage: java -cp MSBooster.jar src.main.java.transferlearn.Predictor " +
                "--paramsList <msbooster parameters> --api-key <key> " +
                "--url <server url> --model <path/to/model/weights> " +
                "optional: --peptide-list-to-dict <path/to/csv> --basename <output base name>  --output-dir <output directory> " +
                "--ms2 <predict ms2> --rt <predict rt> --im <predict ccs> " +
                "--min-charge <int> max-charge <int> --output-format <'parquet' or 'mgf' or 'librarytsv'>");
        Print.printError("Example: java -cp MSBooster.jar src.main.java.transferlearn.Predictor " +
                "--paramsList msbooster_params.txt --url http://localhost:8001 --model model.zip " +
                "--ms2 true --rt true --im false --output_format mgf " +
                "(rather than the default library.tsv in parquet format, write to mgf)");
        Print.printError("Example: java -cp MSBooster.jar src.main.java.transferlearn.Predictor " +
                "--paramsList msbooster_params.txt --url http://localhost:8001 --model model.zip " +
                "(all properties set to true)");
        Print.printError("Example: java -cp MSBooster.jar src.main.java.transferlearn.Predictor " +
                "--paramsList msbooster_params.txt --url http://localhost:8001 --model model.zip " +
                "--basename predictions " +
                "(returns predictions.parquet instead of default model.parquet)");
        Print.printError("Example: java -cp MSBooster.jar src.main.java.transferlearn.Predictor " +
                "--paramsList msbooster_params.txt --url http://localhost:8001 --peptide-list-to-dict peptides.csv " +
                "--min-charge 1 --max-charge 4 " +
                "(For all peptides in peptides.csv, predicts precursors with charge 1 to 4, rather than default 2 to 3." +
                "In peptides.csv, columns must be peptides,proteins,is_decoy. " +
                "Peptides must be formatted as PEPTIDE[PTM mass]R." +
                "'is_decoy' is either 0 or 1)");
        System.exit(1);
    }

    public static void main(String[] args) throws Exception {
        Locale.setDefault(Locale.US);
        Thread.setDefaultUncaughtExceptionHandler((t, e) -> {
			e.printStackTrace();
			System.exit(1);
		});

        printInfo(versionNumber + " Transfer Learned Prediction");

        //parse arguments
        if (args.length % 2 != 0) {
            Print.printError("Malformed arguments, args of length " + args.length);
            errorMessage();
        }

        String params = "";
        String peptideList = "";
        String customMods = "";
        String url = "";
        String model = "";
        String apiKey = "";
        String ms2 = "true";
        String rt = "true";
        String im = "true";
        String basename = "";
        int minCharge = 2;
        int maxCharge = 3;
        String fasta = "";
        String outputFormat = "parquet";
        String outputDir = "";
        String keepDecoys = "1";

        for (int i = 0; i < args.length; i++) {
            switch (args[i]) {
                case "--paramsList":
                    params =  args[i + 1];
                    break;
                case "--peptide-list-to-predict":
                    peptideList = args[i + 1];
                    File file = new File(args[i + 1]);
                    Constants.peptideListDirectory = file.getParent();
                    break;
                case "--custom-mods":
                    customModsStringToTsv(args[i + 1]);
                    customMods = Constants.additionalModsFile;
                    PTMhandler.unimodToModMassAll = makeUnimodToModMassAll(false);
                    AAunimodToModMassAll = makeUnimodToModMassAll(true);
                    PTMhandler.AAunimodToModMassAllKeys = AAunimodToModMassAll.keySet();
                    PTMhandler.aamassToAlphapeptdeep = PTMhandler.makeModAAmassToAlphapeptdeep();
                    break;
                case "--url":
                    url =  args[i + 1];
                    break;
                case "--model":
                    model =  args[i + 1];
                    break;
                case "--api-key":
                    apiKey = args[i + 1];
                    break;
                case "--ms2":
                    ms2 = args[i + 1];
                    break;
                case "--rt":
                    rt = args[i + 1];
                    break;
                case "--im":
                    im = args[i + 1];
                    break;
                case "--basename":
                    basename =  args[i + 1];
                    break;
                case "--min-charge":
                    minCharge = Integer.parseInt(args[i + 1]);
                    break;
                case "--max-charge":
                    maxCharge = Integer.parseInt(args[i + 1]);
                    break;
                case "--fasta":
                    fasta = args[i + 1];
                    break;
                case "--keep-decoys":
                    keepDecoys = args[i + 1];
                    break;
                case "--decoy-tag":
                    Constants.decoyPrefix = args[i + 1];
                    break;
                case "--output-format":
                    outputFormat = args[i + 1];
                    if (outputFormat.equals("parquet") ||
                            outputFormat.equals("mgf") ||
                            outputFormat.equals("librarytsv") ||
                            outputFormat.equals("speclib")) {} else {
                        Print.printError("Unknown output format: " + outputFormat + ". " +
                                "Must be one of mgf, parquet, speclib, or librarytsv");
                        errorMessage();
                    }
                    break;
                case "--output-dir":
                    outputDir =  args[i + 1];
                    File directory = new File(outputDir);
                    if (! directory.exists()) {
                        directory.mkdirs();
                    }
                    break;
            }
        }

        if (params.isEmpty() || url.isEmpty() || model.isEmpty() || apiKey.isEmpty() ||
                (fasta.isEmpty() && (outputFormat.equals("librarytsv") || outputFormat.equals("speclib")))) {
            Print.printError("At least one of params, url, model, fasta, or apiKey is missing.");
            errorMessage();
        }

        //use precursor charge states from MSFragger results
        if (peptideList.isEmpty()) {
            minCharge = 0;
            maxCharge = 0;
        }

        File[] inputFiles = PredictUtils.createPredictInputFiles(peptideList, params, keepDecoys);
        if (outputDir.isEmpty()) {
            outputDir = inputFiles[0].getParent();
        }

        //prediction
        Print.printInfo("Generating predictions");
        URL uploadUrl = new URL(url + "/predict/upload");
        File modelZip = new File(model);
        if (modelZip.getName().contains("_")) {
            Print.printError(modelZip.getName() + " cannot contain the underscore character. " +
                    "Please replace them with dashes and try again.");
            System.exit(1);
        }

        if (basename.isEmpty()) {
            String zipName = modelZip.getName();
            int lastIndex = zipName.lastIndexOf(".zip");
            basename = (lastIndex != -1)
                    ? zipName.substring(0, lastIndex)
                    : zipName;
        }

        //verify custom mods exists
        FileInputStream customFis = null;
        if (!customMods.isEmpty()) {
            if (!Files.exists(Paths.get(customMods))) {
                Print.printError("Custom mods file does not exist");
                System.exit(1);
            }
            customFis = new FileInputStream(customMods);
        }

        EndJob[] endJobs = new EndJob[inputFiles.length];
        String[] jobIds = new String[inputFiles.length];
        for (int i = 0; i < inputFiles.length; i++) {
            String jobId = PredictUtils.sendPredRequest(url, uploadUrl, apiKey, inputFiles[i], modelZip,
                    outputFormat, ms2, rt, im, minCharge, maxCharge, customFis, customMods);
            jobIds[i] = jobId;

            //set up shut down hook
            EndJob endJob = new Helpers.EndJob(url + "/predict/cancel/" + jobId);
            Runtime.getRuntime().addShutdownHook(endJob);
            endJobs[i] = endJob;
        }

        //read fasta if needed
        FastaReader fastaReader = new FastaReader(fasta);
        HashMap<String, String> protMap = new HashMap<>();
        if (outputFormat.equals("librarytsv") ||  outputFormat.equals("speclib")) {
            Print.printInfo("Reading FASTA file");
            fastaReader.mapProtToGene();

            if (! peptideList.isEmpty()) {
                protMap = mapProteinsListToGenes(peptideList, fastaReader.protToGene);
            } else {
                protMap = mapProteinsListToGenes(inputFiles[0].getAbsolutePath(), fastaReader.protToGene);
            }
        }

        String downloadExtension = "";
        switch (outputFormat) {
            case "librarytsv":
                downloadExtension = ".tsv";
                break;
            case "speclib":
                downloadExtension = ".speclib";
                break;
            case "parquet":
                downloadExtension = ".parquet";
                break;
            case  "mgf":
                downloadExtension = ".mgf";
                break;
            default:
                break;
        }

        //download files and merge
        File totalPredFileInitiator = new File(outputDir + File.separator + basename + downloadExtension);
        if (totalPredFileInitiator.exists()) {
            totalPredFileInitiator.delete();
        }
        String[] subFilePaths = new String[inputFiles.length];

        //initiate total file
        String totalFilePath = outputDir + File.separator + basename + downloadExtension;
        BufferedWriter totalPredFile = null;
        if (outputFormat.equals("mgf")) {
            totalPredFile = new BufferedWriter(new FileWriter(totalFilePath, true));
        }

        try (Connection conn = DriverManager.getConnection("jdbc:duckdb:")) {
            for (int i = 0; i < inputFiles.length; i++) {
                String subFilePath = PredictUtils.downloadAndProcess(jobIds[i], url, outputFormat, outputDir,
                        basename + i, protMap, endJobs[i], totalFilePath, conn);
                subFilePaths[i] = subFilePath;
                Print.printInfo("Downloaded file " + (i + 1) + " of " + inputFiles.length);

                //start merging files together
                if (outputFormat.equals("mgf")) {
                    try (BufferedReader in = new BufferedReader(new FileReader(subFilePath))) {
                        String str;
                        while ((str = in.readLine()) != null) {
                            totalPredFile.write(str + "\n");
                        }
                    } catch (Exception e) {
                        e.printStackTrace();
                        totalPredFile.close();
                        System.exit(1);
                    }
                }

                //delete old file
                if (outputFormat.equals("mgf") || outputFormat.equals("librarytsv")) {
                    File miniFile = new File(subFilePath);
                    miniFile.delete();
                }
            }
            if (outputFormat.equals("mgf")) {
                totalPredFile.close();
            }
        } catch (Exception e) {
            e.printStackTrace();
            if (outputFormat.equals("mgf")) {
                totalPredFile.close();
            }
            System.exit(1);
        }

        if (outputFormat.equals("speclib") || outputFormat.equals("parquet")) {
            if (subFilePaths.length > 1) {
                Print.printInfo("Merging parquet files");
                try (Connection conn = DriverManager.getConnection("jdbc:duckdb:");
                     Statement stmt = conn.createStatement()) {

                    String fileList = Arrays.stream(subFilePaths)
                            .map(p -> "'" + p.replace("\\", "/") + "'")
                            .collect(Collectors.joining(", "));

                    stmt.execute("COPY (SELECT * FROM read_parquet([" + fileList + "])) " +
                            "TO '" + totalPredFileInitiator.getAbsolutePath().replace("\\", "/") +
                            "' (FORMAT PARQUET)");
                }

                //delete old files
                for (String subFilePath : subFilePaths) {
                    File subFile = new File(subFilePath);
                    if (subFile.exists()) {
                        subFile.delete();
                    }
                }
            } else {
                Files.move(Paths.get(subFilePaths[0]), Paths.get(totalFilePath));
            }
        }
        if (outputFormat.equals("speclib")) {
            Print.printInfo("Converting parquet to speclib format");
            ParquetToSpecLib ptsl = new ParquetToSpecLib(totalPredFileInitiator.getAbsolutePath(), protMap, -3,
                    true, true, true);
            String speclibPath = totalPredFileInitiator.getAbsolutePath().replace(".parquet", ".speclib");
            ptsl.convertAndWrite(speclibPath);

            Print.printInfo("Total speclib file at " + speclibPath);
        } else {
            Print.printInfo("Total file at " + totalFilePath);
        }

        System.exit(0);
    }
}
