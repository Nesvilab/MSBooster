package transferlearn;

import allconstants.Constants;
import allconstants.NceConstants;
import mainsteps.MainClass;
import mainsteps.ParameterUtils;
import mainsteps.PinMzmlMatcher;
import peptideptmformatting.PeptideFormatter;
import utils.Print;

import java.io.*;
import java.net.HttpURLConnection;
import java.net.URL;
import java.nio.charset.StandardCharsets;
import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Locale;

import static transferlearn.Helpers.*;

public class Predictor {
    static Long waitTime = 15000L;

    private static void errorMessage() {
        Print.printError("Usage: java -cp MSBooster.jar src.main.java.transferlearn.Predictor " +
                "--paramsList <msbooster parameters> --api-key <key> " +
                "--url <server url> --model <path/to/model/weights> " +
                "optional: --peptide-list-to-dict <path/to/csv> --basename <output base name> " +
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
        String outputFormat = "parquet";

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
                case "--output-format":
                    outputFormat = args[i + 1];
                    if (outputFormat.equals("parquet") ||
                    outputFormat.equals("mgf") ||
                    outputFormat.equals("librarytsv")) {} else {
                        Print.printError("Unknown output format: " + outputFormat + ". " +
                                "Must be one of mgf, parquet, or librarytsv");
                        errorMessage();
                    }
                    break;
            }
        }

        if (params.isEmpty() || url.isEmpty() || model.isEmpty() || apiKey.isEmpty()) {
            Print.printError("At least one of params, url, model, or apiKey is missing.");
            errorMessage();
        }

        File inputFile;
        if (peptideList.isEmpty()) { //predict all peptides in pin files
            Constants.spectraModel = "alphapeptdeep";
            Constants.rtModel = "alphapeptdeep";
            Constants.imModel = "alphapeptdeep";
            Constants.createPredFileOnly = true;
            MainClass.main(new String[]{"--paramsList", params});
            inputFile = new File(Constants.spectraRTPrefix + ".csv");
        } else { //predict everything in peptide list
            HashMap<String, String> paramsMap = new HashMap<>();
            ParameterUtils.processCommandLineInputs(new String[]{"--paramsList", params, "--requirePinMzml", "false"},
                    paramsMap);
            ParameterUtils.updateConstants(paramsMap);
            inputFile = new File(Constants.spectraRTPrefix + ".csv");

            //if using peptide list, nce and instrument should be provided
            NceConstants.mzmlNCEs.put("HCD", String.valueOf(NceConstants.NCE));

            //adapt peptide list to alphapeptdeep format
            try (BufferedReader reader = new BufferedReader(new FileReader(peptideList));
                 BufferedWriter writer = new BufferedWriter(new FileWriter(inputFile))
            ) {
                //skip old header, write new header
                String line = reader.readLine(); //peptide,proteins,is_decoy
                writer.write("sequence,mods,mod_sites,charge,nce,instrument,modified,proteins,is_decoy\n");
                while ((line = reader.readLine()) != null) {
                    String[] lineSplit = line.split(",");

                    String charge = "";
                    while (Character.isDigit(lineSplit[0].charAt(lineSplit[0].length() - 1))) {
                        charge += lineSplit[0].charAt(lineSplit[0].length() - 1);
                        lineSplit[0] = lineSplit[0].substring(0, lineSplit[0].length() - 1);
                    }

                    PeptideFormatter pf = new PeptideFormatter(lineSplit[0], charge, "apdpred");
                    if (charge.isEmpty()) {
                        for (int c = minCharge; c < maxCharge + 1; c++) {
                            pf.charge = String.valueOf(c);
                            writer.write(pf.getStripped() + "," + pf.getAlphapeptdeepMods() + "," +
                                    pf.getModPositions() + "," + pf.getCharge() + "," + NceConstants.getNCE() + "," +
                                    Constants.instrument + "," + pf.getLibrarytsv() + "," +
                                    lineSplit[1] + "," + lineSplit[2] + "\n");
                        }
                    } else {
                        writer.write(pf.getStripped() + "," + pf.getAlphapeptdeepMods() + "," +
                                pf.getModPositions() + "," + charge + "," + NceConstants.getNCE() + "," +
                                Constants.instrument + "," + pf.getLibrarytsv() + "," +
                                lineSplit[1] + "," + lineSplit[2] + "\n");
                    }

                }
            } catch (IOException e) {
                Print.printError("Error reading and writing input file for AlphaPeptDeep prediction: "
                        + e.getMessage());
                System.exit(1);
            }
        }

        //convert input to parquet
        String parquetPath = Constants.spectraRTPrefix + ".parquet";
        Helpers.convertCsvToParquet(inputFile.getPath(), parquetPath);
        inputFile = new File(parquetPath);

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
        while (nonResults.contains(status.toUpperCase())) {
            Thread.sleep(waitTime);

            connection = setUpConnection(url, statusUrl);
            responseStream = connection.getInputStream();
            map = readJsonResponse(responseStream);
            status = map.get("status").toString();
        }

        //download
        String downloadExtension = ".mgf";
        if (!outputFormat.equals("mgf")) {
            downloadExtension = ".parquet";
        }
        File downloadPath = new File(inputFile.getParent(), basename + downloadExtension);

        if (status.equals("SUCCESS")) {
            URL downloadUrl = new URL(url + "/predict/download/" + jobId);
            connection = setUpConnection(url, downloadUrl);
            connection.setRequestMethod("GET");

            responseCode = connection.getResponseCode();
            if (responseCode != HttpURLConnection.HTTP_OK) {
                Print.printError("Server returned HTTP " + responseCode + " " + connection.getResponseMessage());
                Print.printError((String) map.get("result"));
                System.exit(1);
            }

            // Read input stream (file content) from server
            try (InputStream inputStream = connection.getInputStream();
                 FileOutputStream outputStream = new FileOutputStream(downloadPath)) {

                byte[] buffer = new byte[4096];
                int bytesRead;
                while ((bytesRead = inputStream.read(buffer)) != -1) {
                    outputStream.write(buffer, 0, bytesRead);
                }
            }

            if (outputFormat.equals("librarytsv")) {
                String tsvPath = downloadPath.getAbsolutePath().replace(".parquet", ".tsv");
                convertParquetToCsv(String.valueOf(downloadPath), tsvPath);
                downloadPath.delete();

                Print.printInfo("File downloaded to: " + tsvPath);
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
