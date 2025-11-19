package transferlearn;

import allconstants.Constants;
import allconstants.NceConstants;
import mainsteps.MainClass;
import mainsteps.ParameterUtils;
import readers.datareaders.FastaReader;
import speclib.io.ParquetToSpecLib;
import utils.Print;

import java.io.*;
import java.net.HttpURLConnection;
import java.net.URL;
import java.nio.charset.StandardCharsets;
import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Locale;

import static allconstants.Constants.versionNumber;
import static transferlearn.Helpers.*;
import static utils.Print.printInfo;

public class Predictor {
    static Long waitTime = 15000L;

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

        printInfo(versionNumber);

        //parse arguments
        if (args.length % 2 != 0) {
            Print.printError("Malformed arguments, args of length " + args.length);
            errorMessage();
        }

        String params = "";
        String peptideList = "";
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

        File inputFile;
        if (peptideList.isEmpty()) { //predict all peptides in pin files
            Constants.spectraModel = "alphapeptdeep";
            Constants.rtModel = "alphapeptdeep";
            Constants.imModel = "alphapeptdeep";
            Constants.createPredFileOnly = true;
            MainClass.main(new String[]{"--paramsList", params, "--keepDecoys", keepDecoys});
            inputFile = new File(Constants.spectraRTPrefix + ".csv");
            minCharge = 0;
            maxCharge = 0;

            //convert input to parquet
            printInfo("Converting AlphaPeptDeep input to parquet");
            String parquetPath = Constants.spectraRTPrefix + ".parquet";
            Helpers.convertCsvToParquet(inputFile.getPath(), parquetPath, true);
            inputFile = new File(parquetPath);
        } else { //predict everything in peptide list
            HashMap<String, String> paramsMap = new HashMap<>();
            ParameterUtils.processCommandLineInputs(new String[]{"--paramsList", params, "--requirePinMzml", "false"},
                    paramsMap);
            ParameterUtils.updateConstants(paramsMap);
            inputFile = new File(Constants.spectraRTPrefix + ".parquet");

            //if using peptide list, nce and instrument should be provided
            NceConstants.mzmlNCEs.put("HCD", String.valueOf(NceConstants.NCE));

            //adapt peptide list to alphapeptdeep format
            convertPeptideListToApdInput(peptideList, inputFile);
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

        HttpURLConnection connection = setUpConnection(url, uploadUrl);
        connection.setRequestMethod("POST");
        connection.setDoOutput(true);

        String boundary = "----JavaFormBoundary" + System.currentTimeMillis();
        connection.setRequestProperty("Content-Type", "multipart/form-data; boundary=" + boundary);
        connection.setRequestProperty("X-Api-Key", apiKey);

        //send request
        try (OutputStream os = connection.getOutputStream();
             PrintWriter writer = new PrintWriter(new OutputStreamWriter(os, StandardCharsets.UTF_8), true);
             FileInputStream fisInput = new FileInputStream(inputFile);
             FileInputStream fisModel = new FileInputStream(modelZip);
        ) {
            // input file
            writer.append("--").append(boundary).append("\r\n");
            writer.append("Content-Disposition: form-data; name=\"input_file\"; filename=\"")
                    .append(inputFile.getName()).append("\"\r\n");
            writer.append("Content-Type: application/vnd.apache.parquet\r\n\r\n");
            writer.flush();

            byte[] buffer = new byte[4096];
            int bytesRead;
            while ((bytesRead = fisInput.read(buffer)) != -1) {
                os.write(buffer, 0, bytesRead);
            }
            os.flush();
            writer.append("\r\n");

            //model zip
            writer.append("--").append(boundary).append("\r\n");
            writer.append("Content-Disposition: form-data; name=\"model_zip\"; filename=\"")
                    .append(modelZip.getName()).append("\"\r\n");
            writer.append("Content-Type: application/zip\r\n\r\n");
            writer.flush();

            buffer = new byte[4096];
            while ((bytesRead = fisModel.read(buffer)) != -1) {
                os.write(buffer, 0, bytesRead);
            }
            os.flush();
            writer.append("\r\n");

            // output format
            writer.append("--").append(boundary).append("\r\n");
            writer.append("Content-Disposition: form-data; name=\"output_format\"\r\n\r\n");
            writer.append(outputFormat).append("\r\n");

            // ms2
            writer.append("--").append(boundary).append("\r\n");
            writer.append("Content-Disposition: form-data; name=\"ms2\"\r\n\r\n");
            writer.append(ms2).append("\r\n");

            // rt
            writer.append("--").append(boundary).append("\r\n");
            writer.append("Content-Disposition: form-data; name=\"rt\"\r\n\r\n");
            writer.append(rt).append("\r\n");

            // ccs
            writer.append("--").append(boundary).append("\r\n");
            writer.append("Content-Disposition: form-data; name=\"ccs\"\r\n\r\n");
            writer.append(im).append("\r\n");

            // mincharge
            writer.append("--").append(boundary).append("\r\n");
            writer.append("Content-Disposition: form-data; name=\"mincharge\"\r\n\r\n");
            writer.append(Integer.toString(minCharge)).append("\r\n");

            // maxcharge
            writer.append("--").append(boundary).append("\r\n");
            writer.append("Content-Disposition: form-data; name=\"maxcharge\"\r\n\r\n");
            writer.append(Integer.toString(maxCharge)).append("\r\n");

            // end of multipart
            writer.append("--").append(boundary).append("--").append("\r\n");
            writer.flush();
        }

        //get response
        int responseCode = connection.getResponseCode();
        InputStream responseStream;
        String jobId = "";
        if (responseCode >= 200 && responseCode < 300) {
            responseStream = connection.getInputStream();

            try {
                HashMap<String, Object> map = readJsonResponse(responseStream);
                jobId = map.get("job_id").toString();

                //set up shut down hook
                Runtime.getRuntime().addShutdownHook(new Helpers.EndJob(url + "/predict/cancel/" + jobId));
            } catch (Exception e) { //success that we don't handle yet
                Print.printError(String.valueOf(e));
                BufferedReader in = new BufferedReader(new InputStreamReader(responseStream));
                String line;
                while ((line = in.readLine()) != null) {
                    Print.printError(line);
                }
                System.exit(1);
            }
        } else { //error
            responseStream = connection.getErrorStream();

            BufferedReader in = new BufferedReader(new InputStreamReader(responseStream));
            String line;
            while ((line = in.readLine()) != null) {
                Print.printError(line);
            }
            System.exit(1);
        }

        //check status of job id
        Print.printInfo("Job ID: " + jobId);
        URL statusUrl = new URL(url + "/predict/status/" + jobId);

        HashSet<String> nonResults = new HashSet<>(Arrays.asList("PENDING", "RECEIVED", "STARTED"));
        String status = "PENDING";
        HashMap<String, Object> map = null;
        int oldStdoutLen = 0;
        while (nonResults.contains(status.toUpperCase())) {
            Thread.sleep(waitTime);

            connection = setUpConnection(url, statusUrl);
            responseStream = connection.getInputStream();
            map = readJsonResponse(responseStream);
            status = map.get("status").toString();
            String stdout = map.get("stdout").toString();
            if (!stdout.isEmpty()) {
                System.out.print(stdout.substring(oldStdoutLen));
                oldStdoutLen = stdout.length();
            }
        }

        //download
        String downloadExtension = ".mgf";
        if (!outputFormat.equals("mgf")) {
            downloadExtension = ".parquet";
        }
        if (outputDir.isEmpty()) {
            outputDir = inputFile.getParent();
        }
        File downloadFile = new File(outputDir, basename + downloadExtension);

        if (status.equals("SUCCESS")) {
            URL downloadUrl = new URL(url + "/predict/download/" + jobId);
            connection = setUpConnection(url, downloadUrl);
            connection.setRequestMethod("GET");

            responseCode = connection.getResponseCode();
            if (responseCode != HttpURLConnection.HTTP_OK) {
                Print.printError("Server returned HTTP " + responseCode + " " + connection.getResponseMessage());
                Print.printError((String) map.get("stdout"));
                System.exit(1);
            }

            // Read input stream (file content) from server
            try (InputStream inputStream = connection.getInputStream();
                 FileOutputStream outputStream = new FileOutputStream(downloadFile)) {

                byte[] buffer = new byte[4096];
                int bytesRead;
                while ((bytesRead = inputStream.read(buffer)) != -1) {
                    outputStream.write(buffer, 0, bytesRead);
                }
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
                    protMap = mapProteinsListToGenes(inputFile.getAbsolutePath(), fastaReader.protToGene);
                }
            }

            String downloadPath = downloadFile.getAbsolutePath();
            if (outputFormat.equals("librarytsv")) {
                String tsvPath = downloadPath.replace(".parquet", ".tsv");
                convertParquetToLibraryTsv(downloadPath, tsvPath, protMap);

                Print.printInfo("File downloaded to: " + tsvPath);
            } else if (outputFormat.equals("speclib")) {
                ParquetToSpecLib ptsl = new ParquetToSpecLib(downloadPath, protMap);
                String speclibPath = downloadPath.replace(".parquet", ".speclib");
                ptsl.convertAndWrite(speclibPath);

                Print.printInfo("File downloaded to: " + speclibPath);
            } else {
                Print.printInfo("File downloaded to: " + downloadPath);
            }
            System.exit(0);
        } else {
            Print.printError(String.valueOf(map));
            Print.printError(connection.getResponseMessage());
            System.exit(1);
        }
    }
}
