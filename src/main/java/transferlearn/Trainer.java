package transferlearn;

import allconstants.Constants;
import peptideptmformatting.PTMhandler;
import utils.Print;

import java.io.*;
import java.net.HttpURLConnection;
import java.net.URL;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Locale;
import java.time.LocalDateTime;

import static allconstants.Constants.versionNumber;
import static transferlearn.Helpers.*;
import static utils.Print.printError;
import static utils.Print.printInfo;

public class Trainer {
    static Long waitTime = 15000L;

    private static void errorMessage() {
        Print.printError("Usage: java -cp MSBooster.jar src.main.java.transferlearn.Trainer " +
                "--url <server url> --library <path/to/librarytsv> --api-key <key> " +
                "optional: --basename <output base name> --output-dir <output directory>");
        Print.printError("Example: java -cp MSBooster.jar src.main.java.transferlearn.Trainer " +
                "--url http://localhost:8000 --library library.tsv --basename weights (returns weights.zip)");
        Print.printError("Example: java -cp MSBooster.jar src.main.java.transferlearn.Trainer " +
                "--url http://localhost:8000 --library library.tsv (returns weights-<datetime>.zip)");
        System.exit(1);
    }

    public static void main(String[] args) throws IOException, InterruptedException {
        Locale.setDefault(Locale.US);
        Thread.setDefaultUncaughtExceptionHandler((t, e) -> {
			e.printStackTrace();
			System.exit(1);
		});

        train(args);
    }
    public static String train(String[] args) throws IOException, InterruptedException {
        printInfo(versionNumber + " Transfer Learning");

        //parse arguments
        if (args.length % 2 != 0) {
            Print.printError("Malformed arguments, args of length " + args.length);
            errorMessage();
        }

        String url = "";
        String library = "";
        String instrument = "";
        float nce = 30f;
        String customMods = "";
        String apiKey = "";
        String basename = "";
        String outputDir = "";

        for (int i = 0; i < args.length; i++) {
            switch (args[i]) {
                case "--url":
                    url =  args[i + 1];
                    break;
                case "--library":
                    library =  args[i + 1];
                    break;
                case "--instrument":
                    instrument =  args[i + 1];
                    break;
                case "--nce":
                    nce = Float.parseFloat(args[i + 1]);
                    break;
                case "--custom-mods":
                    customModsStringToTsv(args[i + 1]);
                    customMods = Constants.additionalModsFile;
                    PTMhandler.unimodToModMassAll = PTMhandler.makeUnimodToModMassAll(false);
                    PTMhandler.AAunimodToModMassAll = PTMhandler.makeUnimodToModMassAll(true);
                    PTMhandler.AAunimodToModMassAllKeys = PTMhandler.AAunimodToModMassAll.keySet();
                    PTMhandler.aamassToAlphapeptdeep = PTMhandler.makeModAAmassToAlphapeptdeep();
                    break;
                case "--api-key":
                    apiKey = args[i + 1];
                    break;
                case "--basename":
                    basename =  args[i + 1];
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

        if (url.isEmpty() || library.isEmpty() || apiKey.isEmpty()) {
            Print.printError("At least one of url, library, or apiKey is missing.");
            errorMessage();
        }

        //the model is downloaded beside the library unless --output-dir says where
        final File libraryFile = new File(library).getAbsoluteFile();
        if (outputDir.isEmpty()) {
            outputDir = libraryFile.getParent();
        }

        printInfo("Converting input to parquet");
        File libraryparquet = parquetForUpload(libraryFile, new File(outputDir));

        printInfo("Transfer learning started");
        URL uploadUrl = new URL(url + "/train/upload");
        String outputBaseName;
        if (! basename.isEmpty()) {
            outputBaseName = basename;
        } else {
            String datetime = String.valueOf(LocalDateTime.now());
            outputBaseName = "weights-" + datetime.split("\\.")[0].replace(":", "-");
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
            writer.append("\r\n");

            // instrument
            writer.append("--").append(boundary).append("\r\n");
            writer.append("Content-Disposition: form-data; name=\"instrument\"\r\n\r\n");
            writer.append(instrument).append("\r\n");

            // nce
            writer.append("--").append(boundary).append("\r\n");
            writer.append("Content-Disposition: form-data; name=\"nce\"\r\n\r\n");
            writer.append(Float.toString(nce)).append("\r\n");

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

            writer.append("--").append(boundary).append("--").append("\r\n");
            writer.flush();
        } finally {
            deleteUpload(libraryparquet);
        }

        //get response
        int responseCode = connection.getResponseCode();
        InputStream responseStream;
        String jobId = "";
        EndJob endJob = null;
        if (responseCode >= 200 && responseCode < 300) {
            responseStream = connection.getInputStream();

            try {
                HashMap<String, Object> map = readJsonResponse(responseStream);
                jobId = map.get("job_id").toString();

                //set up shut down hook
                endJob = new Helpers.EndJob(url + "/train/cancel/" + jobId);
                Runtime.getRuntime().addShutdownHook(endJob);
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
        int oldStdoutLen = 0;
        int oldJobsAhead = 0;
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
        File downloadPath = new File(outputDir, outputBaseName + ".zip");

        if (status.equals("SUCCESS")) {
            endJob.ended = true;

            URL downloadUrl = new URL(url + "/train/download/" + jobId);
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
            endJob.ended = true;
            Print.printError(String.valueOf(map));
            Print.printError(connection.getResponseMessage());
            System.exit(1);
            return "";
        }
    }

    /**
     * The TSV library {@code library} converted to Parquet for the upload, under the library's name
     * with .parquet in place of its extension, in a new folder of its own in {@code dir}. Never beside
     * the library: its folder can hold a Parquet library of that name already (FragSpecLib writes
     * library.parquet beside library.tsv, and FragDIA reads it, footer and all), which a converted
     * copy, with no footer, would replace. The file and its folder are deleted when the JVM exits,
     * unless {@link #deleteUpload} deleted them before.
     */
    static File parquetForUpload(File library, File dir) throws IOException {
        Files.createDirectories(dir.toPath());
        final Path uploadDir = Files.createTempDirectory(dir.toPath(), "library-upload-");
        final File parquet = uploadDir.resolve(parquetName(library.getName())).toFile();
        uploadDir.toFile().deleteOnExit(); // registered first, so deleted after the file in it
        parquet.deleteOnExit();
        Helpers.convertCsvToParquet(library.getPath(), parquet.getPath(), false);
        return parquet;
    }

    /** Deletes the Parquet library {@link #parquetForUpload} made, and its folder. */
    static void deleteUpload(File parquet) {
        try {
            Files.deleteIfExists(parquet.toPath());
            Files.deleteIfExists(parquet.toPath().getParent());
        } catch (IOException e) {
            printError("Could not delete " + parquet.getParent() + ", made for the upload: " + e);
        }
    }

    /** {@code name} with its extension, if any, replaced by .parquet: library.tsv becomes library.parquet. */
    static String parquetName(String name) {
        final int dot = name.lastIndexOf('.');
        return (dot > 0 ? name.substring(0, dot) : name) + ".parquet";
    }
}
