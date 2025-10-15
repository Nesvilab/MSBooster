package transferlearn;

import utils.Print;

import java.io.*;
import java.net.HttpURLConnection;
import java.net.URL;
import java.nio.charset.StandardCharsets;
import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;
import java.time.LocalDateTime;

import static transferlearn.Helpers.readJsonResponse;
import static transferlearn.Helpers.setUpConnection;

public class Trainer {
    static Long waitTime = 15000L;

    private static void errorMessage() {
        Print.printError("Usage: java -cp MSBooster.jar src.main.java.transferlearn.Trainer " +
                "--url <server url> --library <path/to/librarytsv> --api-key <key> " +
                "optional: --basename <output base name>");
        Print.printError("Example: java -cp MSBooster.jar src.main.java.transferlearn.Trainer " +
                "--url http://localhost:8000 --library library.tsv --basename weights (returns weights.zip)");
        Print.printError("Example: java -cp MSBooster.jar src.main.java.transferlearn.Trainer " +
                "--url http://localhost:8000 --library library.tsv (returns weights-<datetime>.zip)");
        System.exit(1);
    }

    public static void main(String[] args) throws IOException, InterruptedException {
        train(args);
    }
    public static String train(String[] args) throws IOException, InterruptedException {
        //parse arguments
        if (args.length % 2 != 0) {
            Print.printError("Malformed arguments, args of length " + args.length);
            errorMessage();
        }

        String url = "";
        String library = "";
        String apiKey = "";
        String basename = "";

        for (int i = 0; i < args.length; i++) {
            switch (args[i]) {
                case "--url":
                    url =  args[i + 1];
                    break;
                case "--library":
                    library =  args[i + 1];
                    break;
                case "--api-key":
                    apiKey = args[i + 1];
                    break;
                case "--basename":
                    basename =  args[i + 1];
                    break;
            }
        }

        if (url.isEmpty() || library.isEmpty() || apiKey.isEmpty()) {
            Print.printError("At least one of url, library, or apiKey is missing.");
            errorMessage();
        }

        //convert library to parquet
        String parquetPath = library.substring(0, library.length() - 3) + "parquet";
        Helpers.convertCsvToParquet(library, parquetPath);
        library = parquetPath;

        Print.printInfo("Transfer learning started");
        URL uploadUrl = new URL(url + "/train/upload");
        File libraryparquet = new File(library);
        String outputBaseName;
        if (! basename.isEmpty()) {
            outputBaseName = basename;
        } else {
            String datetime = String.valueOf(LocalDateTime.now());
            outputBaseName = "weights-" + datetime.split("\\.")[0].replace(":", "-");
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
             FileInputStream fis = new FileInputStream(libraryparquet)) {

            // Start multipart
            writer.append("--").append(boundary).append("\r\n");
            writer.append("Content-Disposition: form-data; name=\"file\"; filename=\"")
                    .append(libraryparquet.getName()).append("\"\r\n");
            writer.append("Content-Type: application/vnd.apache.parquet\r\n\r\n");
            writer.flush();

            byte[] buffer = new byte[4096];
            int bytesRead;
            while ((bytesRead = fis.read(buffer)) != -1) {
                os.write(buffer, 0, bytesRead);
            }
            os.flush();

            // End of multipart
            writer.append("\r\n").flush();
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
        URL statusUrl = new URL(url + "/train/status/" + jobId);

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
        File downloadPath = new File(libraryparquet.getParent(), outputBaseName + ".zip");

        if (status.equals("SUCCESS")) {
            URL downloadUrl = new URL(url + "/train/download/" + jobId);
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
            return downloadPath.getAbsolutePath();
        } else {
            Print.printError(String.valueOf(map));
            Print.printError(connection.getResponseMessage());
            return "";
        }
    }
}
