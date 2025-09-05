package transferlearn;

import allconstants.Constants;
import mainsteps.MainClass;
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

    private static void errorMessage() {
        Print.printError("Usage: java -cp MSBooster.jar src.main.java.transferlearn.Predictor " +
                "--paramsList <msbooster parameters> " +
                "--url <server url> --model <path/to/model/weights> " +
                "optional: --ms2 <predict ms2> --rt <predict rt> --im <predict ccs> --basename <output base name>");
        Print.printError("Example: java -cp MSBooster.jar src.main.java.transferlearn.Predictor " +
                "--paramsList msbooster_params.txt --url http://localhost:8001 --model model.zip " +
                "--ms2 true --rt true --im false");
        Print.printError("Example: java -cp MSBooster.jar src.main.java.transferlearn.Predictor " +
                "--paramsList msbooster_params.txt --url http://localhost:8001 --model model.zip " +
                "(all properties set to true)");
        Print.printError("Example: java -cp MSBooster.jar src.main.java.transferlearn.Predictor " +
                "--paramsList msbooster_params.txt --url http://localhost:8001 --model model.zip " +
                "--basename predictions " +
                "(returns predictions.mgf instead of default mymodel.mgf)");
        System.exit(1);
    }

    public static void main(String[] args) throws Exception {
        //parse arguments
        if (args.length % 2 != 0) {
            errorMessage();
        }

        String params = "";
        String url = "";
        String model = "";
        String ms2 = "true";
        String rt = "true";
        String im = "true";
        String basename = "";

        for (int i = 0; i < args.length; i++) {
            switch (args[i]) {
                case "--paramsList":
                    params =  args[i + 1];
                    break;
                case "--url":
                    url =  args[i + 1];
                    break;
                case "--model":
                    model =  args[i + 1];
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
            }
        }

        if (params.isEmpty() || url.isEmpty() || model.isEmpty()) {
            errorMessage();
        }

        //create pred file only
        MainClass.main(new String[]{"--paramsList", params});
        File inputFile = new File(Constants.spectraRTPrefix + ".csv");

        //prediction
        Print.printInfo("Generating predictions");
        URL uploadUrl = new URL(url + "/upload");
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
        URL statusUrl = new URL(url + "/status/" + jobId);

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
        File downloadPath = new File(inputFile.getParent(), basename + ".mgf");

        if (status.equals("SUCCESS")) {
            URL downloadUrl = new URL(url + "/download/" + jobId);
            connection = setUpConnection(url, downloadUrl);
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
            System.exit(0);
        } else {
            Print.printError(String.valueOf(map));
            Print.printError(connection.getResponseMessage());
            System.exit(1);
        }
    }
}
