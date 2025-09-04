package transferlearn;

import utils.Print;

import java.io.*;
import java.net.HttpURLConnection;
import java.net.URL;
import java.nio.charset.StandardCharsets;
import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;

import static transferlearn.Helpers.readJsonResponse;
import static transferlearn.Helpers.setUpConnection;

public class Predictor {
    static Long waitTime = 15000L;

    public static void main(String[] args) throws IOException, InterruptedException {
        //parse arguments
        if (args.length != 3 && args.length != 6) {
            Print.printError("Usage: java -cp MSBooster.jar src.main.java.transferlearn.Predictor " +
                    "<server url> <path/to/peptide/input/file> <path/to/model/weights> " +
                    "<predict ms2> <predict rt> <predict ccs>");
            Print.printError("Example: java -cp MSBooster.jar src.main.java.transferlearn.Predictor " +
                    "http://localhost:8001 input.csv model.zip true true false");
            Print.printError("Example: java -cp MSBooster.jar src.main.java.transferlearn.Predictor " +
                    "http://localhost:8001 input.csv model.zip (all properties set to true)");
            System.exit(1);
        }
        URL uploadUrl = new URL(args[0] + "/upload");
        File inputFile = new File(args[1]);
        File modelZip = new File(args[2]);
        if (modelZip.getName().contains("_")) {
            Print.printError(modelZip.getName() + " cannot contain the underscore character. " +
                    "Please replace them with dashes and try again.");
            System.exit(1);
        }
        String ms2 = "true";
        String rt = "true";
        String ccs = "true";
        if (args.length == 6) {
            ms2 = args[3];
            rt = args[4];
            ccs = args[5];
        }

        HttpURLConnection connection = setUpConnection(args[0], uploadUrl);
        connection.setRequestMethod("POST");
        connection.setDoOutput(true);

        String boundary = "----JavaFormBoundary" + System.currentTimeMillis();
        connection.setRequestProperty("Content-Type", "multipart/form-data; boundary=" + boundary);

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
            writer.append("Content-Type: text/tab-separated-values\r\n\r\n");
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
            writer.append(ccs).append("\r\n");

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
        URL statusUrl = new URL(args[0] + "/status/" + jobId);

        HashSet<String> nonResults = new HashSet<>(Arrays.asList("PENDING", "RECEIVED", "STARTED"));
        String status = "PENDING";
        HashMap<String, Object> map = null;
        while (nonResults.contains(status.toUpperCase())) {
            Thread.sleep(waitTime);

            connection = setUpConnection(args[0], statusUrl);
            responseStream = connection.getInputStream();
            map = readJsonResponse(responseStream);
            status = map.get("status").toString();
        }

        //download
        File downloadPath = new File(inputFile.getParent(), jobId + ".mgf");

        if (status.equals("SUCCESS")) {
            URL downloadUrl = new URL(args[0] + "/download/" + jobId);
            connection = setUpConnection(args[0], downloadUrl);
            connection.setRequestMethod("GET");

            responseCode = connection.getResponseCode();
            if (responseCode != HttpURLConnection.HTTP_OK) {
                throw new IOException("Server returned HTTP " + responseCode + " " + connection.getResponseMessage());
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

            Print.printInfo("File downloaded to: " + downloadPath);
        } else {
            Print.printError(String.valueOf(map));
            Print.printError(connection.getResponseMessage());
        }
    }
}
