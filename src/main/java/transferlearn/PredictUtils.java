package transferlearn;

import allconstants.Constants;
import allconstants.NceConstants;
import mainsteps.MainClass;
import mainsteps.ParameterUtils;
import speclib.io.ParquetToSpecLib;
import utils.Print;

import java.io.*;
import java.net.HttpURLConnection;
import java.net.URL;
import java.nio.charset.StandardCharsets;
import java.sql.Connection;
import java.sql.SQLException;
import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;

import static transferlearn.Helpers.*;
import static transferlearn.Helpers.convertParquetToLibraryTsv;
import static transferlearn.Helpers.readJsonResponse;
import static transferlearn.Predictor.waitTime;
import static utils.Print.printInfo;

public class PredictUtils {
    public static File[] createPredictInputFiles(String peptideList, String params, String keepDecoys) throws Exception {
        File[] inputFiles;
        if (peptideList.isEmpty()) { //predict all peptides in pin files
            Constants.spectraModel = "alphapeptdeep";
            Constants.rtModel = "alphapeptdeep";
            Constants.imModel = "alphapeptdeep";
            Constants.createPredFileOnly = true;
            MainClass.main(new String[]{"--paramsList", params, "--keepDecoys", keepDecoys});
            File inputFile = new File(Constants.spectraRTPrefix + ".csv");

            //convert input to parquet
            printInfo("Converting AlphaPeptDeep input to parquet");
            String parquetPath = Constants.spectraRTPrefix + ".parquet";
            Helpers.convertCsvToParquet(inputFile.getPath(), parquetPath, true);
            inputFile = new File(parquetPath);
            inputFiles = new File[]{inputFile};
        } else { //predict everything in peptide list
            HashMap<String, String> paramsMap = new HashMap<>();
            ParameterUtils.processCommandLineInputs(new String[]{"--paramsList", params, "--requirePinMzml", "false"},
                    paramsMap);
            ParameterUtils.updateConstants(paramsMap);

            //if using peptide list, nce and instrument should be provided
            NceConstants.mzmlNCEs.put("HCD", String.valueOf(NceConstants.NCE));

            //adapt peptide list to alphapeptdeep format
            //max peptides is 500k per file
            int maxPeptides = 500000;
            inputFiles = convertPeptideListToApdInput(peptideList, Constants.spectraRTPrefix, maxPeptides);
        }
        return inputFiles;
    }

    public static String sendPredRequest(String url, URL uploadUrl, String apiKey, File inputFile, File modelZip,
                                       String outputFormat, String ms2, String rt, String im,
                                       int minCharge, int maxCharge, FileInputStream customFis, String customMods
                                       ) throws Exception {
        HttpURLConnection connection = setUpConnection(url, uploadUrl);
        connection.setRequestMethod("POST");
        connection.setDoOutput(true);

        String boundary = "----JavaFormBoundary" + System.currentTimeMillis();
        connection.setRequestProperty("Content-Type", "multipart/form-data; boundary=" + boundary);
        connection.setRequestProperty("X-Api-Key", apiKey);

        int bufferSize = 65536;
        connection.setChunkedStreamingMode(bufferSize); // 64KB chunks

        //send request
        try (OutputStream os = connection.getOutputStream();
             BufferedOutputStream bos = new BufferedOutputStream(os);
             PrintWriter writer = new PrintWriter(new OutputStreamWriter(bos, StandardCharsets.UTF_8), true);
             FileInputStream fisInput = new FileInputStream(inputFile);
             FileInputStream fisModel = new FileInputStream(modelZip);
        ) {
            // input file
            writer.append("--").append(boundary).append("\r\n");
            writer.append("Content-Disposition: form-data; name=\"input_file\"; filename=\"")
                    .append(inputFile.getName()).append("\"\r\n");
            writer.append("Content-Type: application/vnd.apache.parquet\r\n\r\n");
            writer.flush();

            byte[] buffer = new byte[bufferSize];
            int bytesRead;
            while ((bytesRead = fisInput.read(buffer)) != -1) {
                bos.write(buffer, 0, bytesRead);
            }
            bos.flush();
            writer.append("\r\n");

            //model zip
            writer.append("--").append(boundary).append("\r\n");
            writer.append("Content-Disposition: form-data; name=\"model_zip\"; filename=\"")
                    .append(modelZip.getName()).append("\"\r\n");
            writer.append("Content-Type: application/zip\r\n\r\n");
            writer.flush();

            buffer = new byte[4096];
            while ((bytesRead = fisModel.read(buffer)) != -1) {
                bos.write(buffer, 0, bytesRead);
            }
            bos.flush();
            writer.append("\r\n");

            // output format
            writer.append("--").append(boundary).append("\r\n");
            writer.append("Content-Disposition: form-data; name=\"output_format\"\r\n\r\n");
            writer.append(outputFormat).append("\r\n");

            // decoy tag
            writer.append("--").append(boundary).append("\r\n");
            writer.append("Content-Disposition: form-data; name=\"decoy_tag\"\r\n\r\n");
            writer.append(Constants.decoyPrefix).append("\r\n");

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

            // custom mods
            if (customFis != null) {
                writer.append("--").append(boundary).append("\r\n");
                writer.append("Content-Disposition: form-data; name=\"custom_mods\"; filename=\"")
                        .append(new File(customMods).getName()).append("\"\r\n");
                writer.append("Content-Type: text/tab-separated-values\r\n\r\n");
                writer.flush();

                buffer = new byte[4096];
                while ((bytesRead = customFis.read(buffer)) != -1) {
                    os.write(buffer, 0, bytesRead);
                }
                os.flush();
                writer.append("\r\n");
                customFis.close();
            }

            // end of multipart
            writer.append("--").append(boundary).append("--").append("\r\n");
            writer.flush();
        }

        //register job cancel
        String jobId = getJobId(connection);
        return jobId;
    }

    private static String getJobId(HttpURLConnection connection) throws IOException {
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
                e.printStackTrace();
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
        return jobId;
    }

    public static String downloadAndProcess(String jobId, String url, String outputFormat, String outputDir, String basename,
                                            HashMap<String, String> protMap, EndJob endJob,
                                            String totalPredFilePath, Connection conn, boolean verbose)
            throws IOException, InterruptedException, SQLException {

        //check status of job id
        Print.printInfo("Job ID: " + jobId);
        URL statusUrl = new URL(url + "/predict/status/" + jobId);

        HashSet<String> nonResults = new HashSet<>(Arrays.asList("PENDING", "RECEIVED", "STARTED"));
        String status = "PENDING";
        HashMap<String, Object> map = null;
        int oldStdoutLen = 0;
        int oldJobsAhead = 0;
        InputStream responseStream;
        HttpURLConnection connection = null;
        boolean firstCheck = true;
        while (nonResults.contains(status.toUpperCase())) {
            if (firstCheck) {
                firstCheck = false;
            } else {
                Thread.sleep(waitTime);
            }

            connection = setUpConnection(url, statusUrl);
            responseStream = connection.getInputStream();
            map = readJsonResponse(responseStream);
            status = map.get("status").toString();
            String stdout = map.get("stdout").toString();
            if (!stdout.isEmpty()) {
                if (verbose) {
                    System.out.print(stdout.substring(oldStdoutLen));
                }
                oldStdoutLen = stdout.length();
            } else { //print queue position
                int jobsAhead = ((Number) map.get("queue_position")).intValue();
                if (jobsAhead != oldJobsAhead) {
                    oldJobsAhead = jobsAhead;
                    if (jobsAhead == 1) {
                        System.out.println("Your job is queued. There is currently 1 job ahead of yours.");
                    } else {
                        System.out.println("Your job is queued. There are currently " +
                                jobsAhead + " jobs ahead of yours.");
                    }
                }
            }
        }

        //download
        String downloadExtension = ".mgf";
        if (!outputFormat.equals("mgf")) {
            downloadExtension = ".parquet";
        }

        File downloadFile = new File(outputDir, basename + downloadExtension);

        if (status.equals("SUCCESS")) {
            endJob.ended = true;

            URL downloadUrl = new URL(url + "/predict/download/" + jobId);
            connection = setUpConnection(url, downloadUrl);
            connection.setRequestMethod("GET");

            int responseCode = connection.getResponseCode();
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

            String downloadPath = downloadFile.getAbsolutePath();
            if (outputFormat.equals("librarytsv")) {
                String tsvPath = downloadPath.replace(".parquet", ".tsv");
                convertParquetToLibraryTsv(downloadPath, totalPredFilePath, protMap, conn);
                // Delete input parquet now that we're done reading it
                new File(downloadPath).delete();

                Print.printInfo("File downloaded to: " + tsvPath);
                return tsvPath;
            } else {
                Print.printInfo("File downloaded to: " + downloadPath);
                return downloadPath;
            }
        } else {
            endJob.ended = true;
            Print.printError(String.valueOf(map));
            Print.printError(connection.getResponseMessage());
            System.exit(1); //does this and all tasks?
            return "";
        }
    }
}
