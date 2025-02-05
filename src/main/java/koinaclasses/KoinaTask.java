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

package koinaclasses;

import allconstants.Constants;
import modelcallers.KoinaModelCaller;
import readers.predictionreaders.KoinaLibReader;

import javax.net.ssl.HttpsURLConnection;
import java.io.BufferedReader;
import java.io.FileInputStream;
import java.io.IOException;
import java.io.InputStreamReader;
import java.io.OutputStream;
import java.net.ConnectException;
import java.net.HttpURLConnection;
import java.net.URL;
import java.nio.charset.StandardCharsets;
import java.util.Timer;
import java.util.TimerTask;
import java.util.concurrent.Callable;
import java.util.concurrent.atomic.AtomicLong;

import static allconstants.Constants.numKoinaAttempts;
import static utils.Print.printError;
import static utils.Print.printInfo;

public class KoinaTask implements Callable<Boolean> {
    private final String jsonFilePath;
    private final String property;
    private final String model;
    KoinaLibReader klr;
    AtomicLong waitTime;
    public int failedAttempts = 0;
    public boolean completed = false;

    public KoinaTask(String jsonFilePath, String model, KoinaLibReader klr, AtomicLong waitTime) {
        this.jsonFilePath = jsonFilePath;
        this.property = klr.property;
        this.model = model;
        this.klr = klr;
        this.waitTime = waitTime;
    }

    // Utility method to read the content of the JSON file into a byte array
    private static byte[] readFile(String file) throws Exception {
        FileInputStream fis = new FileInputStream(file);
        byte[] data = new byte[fis.available()];
        fis.read(data);
        fis.close();

        return data;
    }

    @Override
    public Boolean call() throws Exception {
        URL url = new URL(Constants.KoinaURL + model + "/infer");
        HttpURLConnection connection;

        // Read the JSON file's contents into a string
        String jsonInputString = new String(readFile(jsonFilePath), StandardCharsets.UTF_8);

        int numAttempts = 0;
        while (true) {
            if (Constants.KoinaURL.startsWith("http:")) {
                connection = (HttpURLConnection) url.openConnection();
            } else { //https:
                connection = (HttpsURLConnection) url.openConnection();
            }

            connection.setRequestMethod("POST");
            connection.setRequestProperty("Content-Type", "application/json; utf-8");
            connection.setDoOutput(true);

            try (OutputStream os = connection.getOutputStream()) { //this step may occasionally time out
//                Random random = new Random();
//                int randomNumber = random.nextInt(5);
//                if (randomNumber == 0) {
//                    throw new ConnectException();
//                }
                byte[] input = jsonInputString.getBytes(StandardCharsets.UTF_8);
                os.write(input, 0, input.length);
                break;
            } catch (ConnectException e) {
                if (numAttempts >= numKoinaAttempts) {
                    printError("Koina server is busy. Please retry later, or use DIA-NN for prediction. Exiting");
                    System.exit(1);
                }
                printInfo("Koina server may be running slow. Resending " + jsonFilePath);
                numAttempts++;
            }
        }

        StringBuilder response = new StringBuilder();
        try {
            long start = System.currentTimeMillis();
            BufferedReader in = new BufferedReader(new InputStreamReader(connection.getInputStream()));

            //Goal is to cut off 504 error early
            Timer timer = new Timer();
            TimerTask timerTask = new TimerTask() {
                @Override
                public void run() {
                    if (response.toString().isEmpty()) {
                        try {
                            in.close();
                        } catch (IOException ignored) {}
                    }
                }
            };
            timer.schedule(timerTask, waitTime.get());

            String inputLine;
            while ((inputLine = in.readLine()) != null) {
                response.append(inputLine);
            }
            in.close();

            KoinaModelCaller.parseKoinaOutput(jsonFilePath, response.toString(), property, model, klr);
            long timeDiff = System.currentTimeMillis() - start;
            long currentWaitTime = waitTime.get();
            waitTime.set(currentWaitTime + ((3 * timeDiff - currentWaitTime) / 30));
            return true;
        } catch (Exception e) {
            try {
                if (failedAttempts == numKoinaAttempts) {
                    BufferedReader inErrors = new BufferedReader(new InputStreamReader(connection.getErrorStream()));
                    String inputLine;
                    while ((inputLine = inErrors.readLine()) != null) {
                        response.append(inputLine);
                    }
                    inErrors.close();

                    printError(jsonFilePath + " had the following output: ");
                    printError(response.toString());

                    if (Constants.foundBest) {
                        printError("Retried calling " + jsonFilePath + " " + failedAttempts + " times.");
                        printError("Exiting");
                        System.exit(1);
                    } else {
                        klr.failed = true;
                        return true;
                    }
                }

                return false;
            } catch (Exception e2) {
                e.printStackTrace();
                System.exit(1);
            }
            return false;
        }
    }
}
