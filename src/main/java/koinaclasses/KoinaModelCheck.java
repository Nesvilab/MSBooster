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
import com.google.gson.JsonArray;
import com.google.gson.JsonElement;
import com.google.gson.JsonObject;
import com.google.gson.JsonParseException;
import com.google.gson.JsonParser;

import java.io.IOException;
import java.io.InputStream;
import java.net.HttpURLConnection;
import java.net.MalformedURLException;
import java.net.URI;
import java.net.URL;
import java.nio.charset.StandardCharsets;
import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import java.util.Set;
import java.util.concurrent.ConcurrentHashMap;
import java.util.concurrent.ExecutionException;
import java.util.concurrent.FutureTask;

import static utils.Print.printInfo;

/**
 * Asks the Koina server about a model before MSBooster sends it a request: {@code GET <KoinaURL><model>} gives the
 * model's metadata (KServe v2), with the names of its outputs, and {@code GET <KoinaURL><model>/config} its Triton
 * configuration, with the largest batch it takes ({@code max_batch_size}). Each model is asked about once per run
 * and Koina URL, the first time a request is written or sent for it.
 * <ul>
 *   <li>A model the server does not serve (HTTP 400, as Koina answers, or 404, or an answer that is not model
 *       metadata) ends the run, before any of its requests is sent. So does a candidate of the best-model search:
 *       the search cannot leave a model out, and a candidate with no predictions would end it later anyway.</li>
 *   <li>A model that does not answer with every output its spec reads ends the run, as reading its first response
 *       would.</li>
 *   <li>Requests carry at most the spec's batch size, or the server's max_batch_size when that is smaller and
 *       above 0.</li>
 *   <li>A server that does not answer is asked again, as a request is ({@link Constants#numKoinaAttempts} more
 *       times). If it still does not answer, the check is no stricter than the requests: they carry the spec's
 *       batch size, and they are sent, retried and given up on as they always were.</li>
 * </ul>
 * With no Koina URL, nothing is asked: KoinaModelCaller says that a URL is needed before any request is sent.
 */
public final class KoinaModelCheck {
    //pause between two attempts to ask a server that did not answer
    private static final long RETRY_PAUSE_MS = 1000;
    //characters of an answer quoted in a message
    private static final int QUOTE_LENGTH = 300;

    //the answers, by Koina URL and model; a FutureTask so each model is asked once while other threads wait for it
    private static final ConcurrentHashMap<String, FutureTask<Answer>> ANSWERS = new ConcurrentHashMap<>();
    //the answers already reported, so a report is printed once
    private static final Set<String> REPORTED = ConcurrentHashMap.newKeySet();

    private KoinaModelCheck() {}

    /** What the server said about a model. */
    enum Status {
        /** the server answered with the model's metadata */
        SERVED,
        /** the server answered that it has no such model, or answered with something that is not model metadata */
        NOT_SERVED,
        /** the server did not answer */
        NO_ANSWER
    }

    /**
     * @param outputs      the output names in the model's metadata; empty unless SERVED
     * @param maxBatchSize the max_batch_size of the model's config; -1 if it is not known
     * @param problem      NOT_SERVED and NO_ANSWER: why. SERVED: why the max_batch_size is not known, or null
     */
    record Answer(Status status, List<String> outputs, int maxBatchSize, String problem) {}

    /**
     * Asks Koina about the model, once per run, and ends the run if the run cannot use it (see the class). A server
     * that could not be asked passes.
     *
     * @throws KoinaException if the run cannot use the model
     */
    public static void check(String model) {
        batchSize(model);
    }

    /**
     * The most precursors one request to the model carries: the spec's batch size, or the server's max_batch_size
     * when that is smaller and above 0. The server is asked about the model once per run (see the class).
     *
     * @throws KoinaException if the run cannot use the model
     */
    public static int batchSize(String model) {
        KoinaModelSpec spec = KoinaModelRegistry.require(model);
        if (Constants.KoinaURL == null || Constants.KoinaURL.isEmpty()) {
            return spec.batchSize();
        }
        String modelsUrl = modelsUrl();
        String key = modelsUrl + spec.name();
        Answer answer = answer(spec, modelsUrl, key);

        switch (answer.status()) {
            case NOT_SERVED:
                throw new KoinaException("Koina at " + modelsUrl + " does not serve the model " + spec.name() + ": "
                        + answer.problem() + ". Check the model name, and that the Koina URL ends with /v2/models/");
            case NO_ANSWER:
                if (REPORTED.add(key)) {
                    printInfo("Could not ask Koina at " + modelsUrl + " about the model " + spec.name() + " ("
                            + answer.problem() + "). Its requests carry up to " + spec.batchSize() + " precursors");
                }
                return spec.batchSize();
            default: //SERVED
                List<String> missing = new ArrayList<>(spec.outputs());
                missing.removeAll(answer.outputs());
                if (!missing.isEmpty()) {
                    throw new KoinaException("The Koina model " + spec.name() + " at " + modelsUrl
                            + " does not answer with the output" + (missing.size() == 1 ? " " : "s ")
                            + String.join(", ", missing) + " that MSBooster reads. It answers with "
                            + (answer.outputs().isEmpty() ? "no outputs" : String.join(", ", answer.outputs())));
                }
                int batchSize = spec.batchSize();
                if (answer.maxBatchSize() > 0 && answer.maxBatchSize() < batchSize) {
                    batchSize = answer.maxBatchSize();
                    if (REPORTED.add(key)) {
                        printInfo("Koina takes batches of up to " + batchSize + " precursors for " + spec.name()
                                + ", so its requests carry " + batchSize + " instead of " + spec.batchSize());
                    }
                } else if (answer.maxBatchSize() < 0 && REPORTED.add(key)) {
                    printInfo("Koina did not say how large a batch " + spec.name() + " takes ("
                            + answer.problem() + "). Its requests carry up to " + batchSize + " precursors");
                }
                return batchSize;
        }
    }

    /**
     * Forgets what Koina answered, so the next check of a model asks again. A run asks once per model; this is for
     * tests, whose stub servers may come back on a port an earlier one had.
     */
    public static void forgetAnswers() {
        ANSWERS.clear();
        REPORTED.clear();
    }

    /** The Koina URL, ending with a slash, so a model name can be appended to it. */
    public static String modelsUrl() {
        String url = Constants.KoinaURL;
        return url.endsWith("/") ? url : url + "/";
    }

    //the answer about the model, asking the server if this is the first time
    private static Answer answer(KoinaModelSpec spec, String modelsUrl, String key) {
        FutureTask<Answer> task = ANSWERS.get(key);
        if (task == null) {
            FutureTask<Answer> asking = new FutureTask<>(() -> ask(spec.name(), modelsUrl));
            task = ANSWERS.putIfAbsent(key, asking);
            if (task == null) {
                task = asking;
                asking.run();
            }
        }
        try {
            return task.get();
        } catch (InterruptedException e) {
            Thread.currentThread().interrupt();
            throw new KoinaException("Asking Koina about the model " + spec.name() + " was interrupted", e);
        } catch (ExecutionException e) {
            Throwable cause = e.getCause();
            if (cause instanceof RuntimeException) {
                throw (RuntimeException) cause;
            }
            if (cause instanceof Error) {
                throw (Error) cause;
            }
            throw new KoinaException("Asking Koina about the model " + spec.name() + " failed: " + cause, cause);
        }
    }

    //asks the server for the model's metadata, then its config
    static Answer ask(String model, String modelsUrl) {
        String metadataUrl = modelsUrl + model;
        Reply metadata = get(metadataUrl);
        if (metadata.code() == 400 || metadata.code() == 404) {
            return new Answer(Status.NOT_SERVED, List.of(), -1,
                    "GET " + metadataUrl + " answered HTTP " + metadata.code() + " " + quote(metadata.body()));
        }
        if (metadata.code() != 200) {
            return new Answer(Status.NO_ANSWER, List.of(), -1, metadata.problem(metadataUrl));
        }
        List<String> outputs = outputNames(metadata.body());
        if (outputs == null) {
            return new Answer(Status.NOT_SERVED, List.of(), -1,
                    "GET " + metadataUrl + " answered with something that is not model metadata: "
                            + quote(metadata.body()));
        }

        String configUrl = metadataUrl + "/config";
        Reply config = get(configUrl);
        if (config.code() != 200) {
            return new Answer(Status.SERVED, outputs, -1, config.problem(configUrl));
        }
        Integer maxBatchSize = maxBatchSize(config.body());
        if (maxBatchSize == null) {
            return new Answer(Status.SERVED, outputs, -1,
                    "GET " + configUrl + " answered with no max_batch_size: " + quote(config.body()));
        }
        return new Answer(Status.SERVED, outputs, maxBatchSize, null);
    }

    //the names in the "outputs" of KServe v2 model metadata, or null if the body is not model metadata
    private static List<String> outputNames(String body) {
        try {
            JsonElement root = JsonParser.parseString(body);
            if (!root.isJsonObject() || !root.getAsJsonObject().has("outputs")
                    || !root.getAsJsonObject().get("outputs").isJsonArray()) {
                return null;
            }
            JsonArray outputs = root.getAsJsonObject().getAsJsonArray("outputs");
            List<String> names = new ArrayList<>();
            for (JsonElement output : outputs) {
                if (output.isJsonObject() && output.getAsJsonObject().has("name")) {
                    names.add(output.getAsJsonObject().get("name").getAsString());
                }
            }
            return Collections.unmodifiableList(names);
        } catch (JsonParseException | IllegalStateException | UnsupportedOperationException e) {
            return null;
        }
    }

    //the max_batch_size of a Triton model config (0: the model takes no batches), or null if it has none
    private static Integer maxBatchSize(String body) {
        try {
            JsonElement root = JsonParser.parseString(body);
            if (!root.isJsonObject()) {
                return null;
            }
            JsonObject config = root.getAsJsonObject();
            if (!config.has("max_batch_size") || !config.get("max_batch_size").isJsonPrimitive()) {
                return null;
            }
            int maxBatchSize = config.get("max_batch_size").getAsInt();
            return maxBatchSize < 0 ? null : maxBatchSize;
        } catch (JsonParseException | IllegalStateException | NumberFormatException
                 | UnsupportedOperationException e) {
            return null;
        }
    }

    /**
     * The status and body of a GET. code is -1 when there was no answer.
     *
     * @param failure why there was no answer, or null
     */
    private record Reply(int code, String body, String failure) {
        String problem(String url) {
            if (code < 0) {
                return "GET " + url + " got no answer: " + failure;
            }
            return "GET " + url + " answered HTTP " + code + " " + quote(body);
        }
    }

    //a GET, asked again while there is no answer or the server answers with an error of its own (HTTP 5xx and
    //others), as a request is. HTTP 200, 400 and 404 are answers
    private static Reply get(String url) {
        long timeout = Constants.initialKoinaMillisecondsToWaitRtIm;
        Reply reply = null;
        for (int attempt = 0; attempt <= Constants.numKoinaAttempts; attempt++) {
            if (attempt > 0) {
                try {
                    Thread.sleep(RETRY_PAUSE_MS);
                } catch (InterruptedException e) {
                    Thread.currentThread().interrupt();
                    return reply;
                }
            }
            reply = getOnce(url, (int) Math.min(Integer.MAX_VALUE, timeout));
            if (reply.code() == 200 || reply.code() == 400 || reply.code() == 404) {
                return reply;
            }
            timeout += 5000; //as KoinaModelCaller waits longer after a failed request
        }
        return reply;
    }

    private static Reply getOnce(String url, int timeoutMs) {
        URL target;
        try {
            target = URI.create(url).toURL();
        } catch (IllegalArgumentException | MalformedURLException e) {
            throw new KoinaException("The Koina URL " + Constants.KoinaURL + " is not a URL: " + e.getMessage(), e);
        }
        HttpURLConnection connection = null;
        try {
            connection = (HttpURLConnection) target.openConnection();
            connection.setRequestMethod("GET");
            connection.setConnectTimeout(timeoutMs);
            connection.setReadTimeout(timeoutMs);
            int code = connection.getResponseCode();
            InputStream in = code < 400 ? connection.getInputStream() : connection.getErrorStream();
            String body = "";
            if (in != null) {
                try (in) {
                    body = new String(in.readAllBytes(), StandardCharsets.UTF_8);
                }
            }
            return new Reply(code, body, null);
        } catch (IOException e) {
            return new Reply(-1, "", e.toString());
        } finally {
            if (connection != null) {
                connection.disconnect();
            }
        }
    }

    //a body as it can be quoted in one line of a message
    private static String quote(String body) {
        String line = body.replaceAll("\\s+", " ").trim();
        if (line.isEmpty()) {
            return "(no body)";
        }
        return line.length() > QUOTE_LENGTH ? line.substring(0, QUOTE_LENGTH) + "..." : line;
    }
}
