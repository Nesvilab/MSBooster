/*
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or
 * implied. See the License for the specific language governing
 * permissions and limitations under the License.
 */

package koinaclasses;

import allconstants.Constants;
import modelcallers.KoinaModelCaller;
import readers.predictionreaders.KoinaLibReader;

import javax.net.ssl.HttpsURLConnection;
import java.io.BufferedReader;
import java.io.FileInputStream;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.io.OutputStream;
import java.net.HttpURLConnection;
import java.net.SocketException;
import java.net.SocketTimeoutException;
import java.net.URL;
import java.nio.charset.StandardCharsets;
import java.util.List;
import java.util.Locale;
import java.util.Timer;
import java.util.TimerTask;
import java.util.concurrent.Callable;
import java.util.concurrent.atomic.AtomicLong;
import java.util.function.Supplier;

import static allconstants.Constants.numKoinaAttempts;
import static utils.Print.printError;
import static utils.Print.printInfo;

/**
 * Sends one request to a Koina model and reads the predictions it answers with into the KoinaLibReader. The request
 * is a request file JSONWriter wrote, or some of its precursors ({@link #slice}).
 * <ul>
 *   <li>{@link Result#PREDICTED}: the predictions are read.</li>
 *   <li>{@link Result#NOT_RECEIVED}: the response could not be received, so KoinaModelCaller sends the request
 *       again, until the last attempt; then {@link Result#GIVEN_UP} while the best model is looked for, and a
 *       {@link KoinaException} otherwise.</li>
 *   <li>{@link Result#REFUSED}: Koina answered HTTP 400, as it does for a precursor its model cannot take. A
 *       request of several precursors would be refused again, so it is not resent: {@link KoinaBisection} sends its
 *       halves. A request of one precursor is sent once more ({@link #confirmation}) before its precursor counts as
 *       refused, as Koina answers HTTP 400 also for a failure of its own.</li>
 * </ul>
 * Koina answers HTTP 400 also when its model instance fails, whatever the request holds. An answer that names such a
 * failure ({@link #isServerFailure}) is no refusal: it counts as a response not received. A CUDA failure leaves the
 * instance failing every request it gets until it is restarted, so cutting its requests would take precursors for
 * refused that Koina can predict.
 * A response that arrives but cannot be read throws a {@link KoinaException}, since Koina would answer the same
 * request the same way. A server that takes the request and does not answer is not waited for forever
 * ({@link #timeoutMs}): its response counts as not received.
 */
public class KoinaTask implements Callable<KoinaTask.Result> {
    /** What came of sending the request. */
    public enum Result {
        /** the predictions were read into the KoinaLibReader */
        PREDICTED,
        /** the response could not be received; the request may be sent again */
        NOT_RECEIVED,
        /** Koina refused the request with HTTP 400; {@link #refusal()} says what it answered */
        REFUSED,
        /** the response could not be received in the last attempt; the KoinaLibReader is marked failed */
        GIVEN_UP
    }

    private final String jsonFilePath;
    //what messages call the request: the file, and which of its precursors a slice holds
    private final String label;
    private final String model;
    KoinaLibReader klr;
    AtomicLong waitTime;
    public int failedAttempts = 0;
    public boolean completed = false;
    //the request, read from the file once, or given to a slice
    private String request;
    //the precursors of the request
    private KoinaBatch batch;
    //the request file's precursors a slice starts at
    private final int firstPrecursor;
    //a precursor sent alone to see whether Koina refuses every request of the model (KoinaBisection)
    private final boolean probe;
    //a request of one precursor Koina refused, sent once more to confirm it refuses it (KoinaBisection)
    private final boolean confirming;
    //how long the request waits before it is sent
    private final long pauseMs;
    //what Koina answered when it refused the request
    private String refusal;
    //why the last attempt did not predict the request, when its response counts as not received
    private String failure;

    /**
     * What Koina's answer says when its model instance failed rather than the request: a CUDA failure (seen as
     * "TritonModelException: Failed to record the event.: device-side assert triggered" from
     * Prosit_2025_intensity_MultiFrag), or memory running out. Lower case.
     */
    private static final List<String> SERVER_FAILURES = List.of("device-side assert", "cuda error",
            "failed to record the event", "out of memory");

    public KoinaTask(String jsonFilePath, String model, KoinaLibReader klr, AtomicLong waitTime) {
        this.jsonFilePath = jsonFilePath;
        this.label = jsonFilePath;
        this.model = model;
        this.klr = klr;
        this.waitTime = waitTime;
        this.firstPrecursor = 0;
        this.probe = false;
        this.confirming = false;
        this.pauseMs = 0;
    }

    private KoinaTask(KoinaTask parent, int from, int to, boolean probe) {
        this.jsonFilePath = parent.jsonFilePath;
        this.model = parent.model;
        this.klr = parent.klr;
        this.waitTime = parent.waitTime;
        this.request = KoinaBatch.sliceRequest(parent.request, from, to);
        this.batch = parent.batch.slice(from, to);
        this.firstPrecursor = parent.firstPrecursor + from;
        this.probe = probe;
        this.confirming = false;
        this.pauseMs = 0;
        this.label = jsonFilePath + (to - from == 1 ? " precursor " + (firstPrecursor + 1)
                : " precursors " + (firstPrecursor + 1) + "-" + (firstPrecursor + to - from));
    }

    //the same request as refused, sent again after a pause
    private KoinaTask(KoinaTask refused, long pauseMs) {
        this.jsonFilePath = refused.jsonFilePath;
        this.label = refused.label;
        this.model = refused.model;
        this.klr = refused.klr;
        this.waitTime = refused.waitTime;
        this.request = refused.request;
        this.batch = refused.batch;
        this.firstPrecursor = refused.firstPrecursor;
        this.probe = refused.probe;
        this.confirming = true;
        this.pauseMs = pauseMs;
    }

    /**
     * A request for precursors [from, to) of this one, which has been sent.
     *
     * @param probe whether it is sent to see whether Koina refuses every request of the model
     */
    public KoinaTask slice(int from, int to, boolean probe) {
        if (batch == null) {
            throw new IllegalStateException(label + " has not been read yet");
        }
        return new KoinaTask(this, from, to, probe);
    }

    /**
     * This request of one precursor, which Koina refused, to be sent once more after a pause. Koina answers HTTP 400
     * also for a failure of its own that the next request may not meet, so a precursor counts as refused only when
     * Koina refuses it twice.
     */
    public KoinaTask confirmation(long pauseMs) {
        if (batch == null || batch.size() != 1) {
            throw new IllegalStateException(label + " is not a request of one precursor that has been sent");
        }
        return new KoinaTask(this, pauseMs);
    }

    /** Whether this request is sent again to confirm that Koina refuses it ({@link #confirmation}). */
    public boolean isConfirmation() {
        return confirming;
    }

    /** How long the request waits before it is sent, in milliseconds. */
    public long pauseMs() {
        return pauseMs;
    }

    /** The precursors of the request; null before it is first sent. */
    public KoinaBatch batch() {
        return batch;
    }

    /** The number of precursors of the request, once it has been sent. */
    public int size() {
        return batch.size();
    }

    /** What Koina answered when it refused the request; null unless it did. */
    public String refusal() {
        return refusal;
    }

    public boolean isProbe() {
        return probe;
    }

    /** Why the last attempt did not predict the request, when it was {@link Result#NOT_RECEIVED}; else null. */
    public String failure() {
        return failure;
    }

    /**
     * Whether what Koina answered with HTTP 400 says its model instance failed, whatever the request holds
     * ({@link #SERVER_FAILURES}), rather than that it cannot take something in the request.
     */
    static boolean isServerFailure(String answer) {
        if (answer == null) {
            return false;
        }
        String lower = answer.toLowerCase(Locale.ROOT);
        return SERVER_FAILURES.stream().anyMatch(lower::contains);
    }

    /** The request file, and which of its precursors this request holds if not all. */
    public String label() {
        return label;
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
    public Result call() throws Exception {
        URL url = new URL(KoinaModelCheck.modelsUrl() + model + "/infer");
        HttpURLConnection connection;

        // Read the JSON file's contents into a string
        if (request == null) {
            request = new String(readFile(jsonFilePath), StandardCharsets.UTF_8);
        }
        if (batch == null) {
            try {
                batch = KoinaBatch.fromRequest(request);
            } catch (IllegalArgumentException e) {
                throw new KoinaException("MSBooster could not read its request file " + jsonFilePath + " for Koina "
                        + "model " + model + ": " + e.getMessage(), e);
            }
        }

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
            //a server that does not answer is not waited for forever: a connection or response that times out
            //throws, and the request is sent again as for any response not received
            int timeout = timeoutMs();
            connection.setConnectTimeout(timeout);
            connection.setReadTimeout(timeout);

            try (OutputStream os = connection.getOutputStream()) { //this step may occasionally time out
//                Random random = new Random();
//                int randomNumber = random.nextInt(5);
//                if (randomNumber == 0) {
//                    throw new ConnectException();
//                }
                byte[] input = request.getBytes(StandardCharsets.UTF_8);
                os.write(input, 0, input.length);
                break;
            } catch (SocketException | SocketTimeoutException e) {
                if (numAttempts >= numKoinaAttempts) {
                    throw new KoinaException("Koina server is busy: MSBooster could not send " + label
                            + " to Koina model " + model + " (" + e.getMessage() + "). Please retry later, or use "
                            + "DIA-NN for prediction");
                }
                printInfo("Koina server may be running slow. Resending " + label);
                numAttempts++;
            }
        }

        //receive the response. Only a failure here is sent back to be retried: the server may answer next time
        long start = System.currentTimeMillis();
        HttpURLConnection sent = connection;
        int code;
        try {
            code = sent.getResponseCode();
        } catch (IOException e) {
            return failedAttempt(describe(e), () -> errorBody(sent));
        }
        if (code == HttpURLConnection.HTTP_BAD_REQUEST) {
            String answer = errorBody(sent);
            if (isServerFailure(answer)) {
                //Koina's model instance failed, not the request: it is sent again whole, as for a response not
                //received
                return failedAttempt("HTTP 400 " + KoinaBisection.quote(answer), () -> answer);
            }
            //Koina cannot take something in the request, or failed on it: KoinaBisection decides what is sent next
            refusal = answer;
            return Result.REFUSED;
        }
        String response;
        try {
            response = readResponse(sent);
        } catch (Exception e) {
            return failedAttempt(describe(e), () -> errorBody(sent));
        }

        //read the predictions. A response that cannot be read is not retried, since Koina would answer the same
        //request the same way: parseKoinaOutput throws a KoinaException naming the model, which ends the run
        KoinaModelCaller.parseKoinaOutput(response, batch.peptides(), batch.charges(), model, klr);
        long timeDiff = System.currentTimeMillis() - start;
        long currentWaitTime = waitTime.get();
        waitTime.set(Math.min(60000, currentWaitTime + ((3 * timeDiff - currentWaitTime) / 30))); //don't wait longer than minute
        return Result.PREDICTED;
    }

    /**
     * How long a connection, and then the response, is waited for: waitTime, which follows how long the model's
     * requests take, but never less than it is at the start for the model's kind, so a server a little slower than
     * usual is still waited for.
     */
    int timeoutMs() {
        long initial = klr.spec.kind() == KoinaModelSpec.Kind.MS2 ? Constants.initialKoinaMillisecondsToWaitMs2
                : Constants.initialKoinaMillisecondsToWaitRtIm;
        return (int) Math.min(Integer.MAX_VALUE, Math.max(1, Math.max(waitTime.get(), initial)));
    }

    //the body of the response. The read is cut off when nothing has arrived after waitTime, to end a 504 early
    private String readResponse(HttpURLConnection connection) throws IOException {
        StringBuilder response = new StringBuilder();
        BufferedReader in = new BufferedReader(new InputStreamReader(connection.getInputStream(),
                StandardCharsets.UTF_8));
        Timer timer = new Timer(true);
        try {
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
        } finally {
            timer.cancel();
            in.close();
        }
        return response.toString();
    }

    //after a response that could not be received: NOT_RECEIVED to have the request sent again, until the last
    //attempt. Then the request is given up: while the best model is looked for, the model is marked failed;
    //otherwise the run cannot go on
    //why: what went wrong, for messages; answer: what the server answered, read only when the request is given up
    private Result failedAttempt(String why, Supplier<String> answer) {
        failure = why;
        if (failedAttempts < numKoinaAttempts) {
            return Result.NOT_RECEIVED;
        }
        printError(label + " had the following output: ");
        printError(answer.get());

        if (Constants.foundBest) {
            throw new KoinaException("Koina model " + model + " did not answer " + label + " (" + why + "). "
                    + "Retried it " + failedAttempts + " times. Please retry later, or use DIA-NN for prediction");
        }
        klr.failed = true;
        return Result.GIVEN_UP;
    }

    //an exception that ended an attempt, for messages
    private static String describe(Exception e) {
        return e.getMessage() == null ? e.getClass().getSimpleName()
                : e.getClass().getSimpleName() + ": " + e.getMessage();
    }

    //what the server said with its error status, if it said anything
    private static String errorBody(HttpURLConnection connection) {
        InputStream errorStream = connection.getErrorStream();
        if (errorStream == null) { //no error status: the connection failed, or the read was cut off
            return "(no response)";
        }
        StringBuilder body = new StringBuilder();
        try (BufferedReader inErrors = new BufferedReader(new InputStreamReader(errorStream,
                StandardCharsets.UTF_8))) {
            String inputLine;
            while ((inputLine = inErrors.readLine()) != null) {
                body.append(inputLine);
            }
        } catch (IOException e) {
            body.append(" (the rest could not be read: ").append(e.getMessage()).append(")");
        }
        return body.toString();
    }
}
