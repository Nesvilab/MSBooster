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

package testutils;

import com.sun.net.httpserver.HttpExchange;
import com.sun.net.httpserver.HttpServer;
import koinaclasses.KoinaBatch;
import koinaclasses.KoinaModelSpec;

import java.io.IOException;
import java.io.InputStream;
import java.io.OutputStream;
import java.net.InetAddress;
import java.net.InetSocketAddress;
import java.nio.charset.StandardCharsets;
import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.concurrent.ConcurrentHashMap;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.atomic.AtomicInteger;
import java.util.function.Function;
import java.util.stream.Collectors;

/**
 * A Koina server on 127.0.0.1 for tests, so no test calls the real one. It speaks the three routes
 * MSBooster uses under {@link #url()} ({@code .../v2/models/}): {@code POST <model>/infer} answers
 * with the response a test gave it, and {@code GET <model>} and {@code GET <model>/config} answer
 * with the model's metadata and config. A model it was given nothing for gets HTTP 400 and
 * Triton's "unknown model" error, as Koina answers; a path outside {@code /v2/models/} gets 404.
 * Every infer body is kept, and every GET counted, so a test can check what MSBooster sent.
 *
 * <p>An infer response is answered at most {@code maxAnswers} times per model, and after that the
 * stub answers HTTP 500, which KoinaTask takes for a response it could not receive and resends the
 * request for. A response it receives but cannot read is not resent: KoinaModelCaller throws a
 * KoinaException, so one answer is enough to see it. A model can instead be given a function that
 * answers each request ({@link #serveInfer(String, Function)}), so requests of any part of a batch
 * are answered.
 *
 * <p>Koina answers HTTP 400 to a request with a precursor its model cannot take. The stub does so
 * for any request with a peptide it was told to refuse ({@link #refusePeptide}), or for every
 * request ({@link #refuseAll}), whatever else it would answer. Koina also answers HTTP 400 for a
 * failure of its own that the next request may not meet, such as a model instance out of memory;
 * the stub can answer so the first time each request arrives ({@link #failNewRequests}). A CUDA
 * failure leaves the model instance failing every request it gets, whatever the request holds,
 * until it is restarted; the stub can answer so to the first requests of a model, whatever they
 * hold ({@link #breakInstance}).
 */
public final class KoinaStubServer implements AutoCloseable {
    private static final String ROOT = "/v2/models/";

    private final HttpServer server;
    private final ExecutorService executor;
    private final Map<String, String> inferResponses = new ConcurrentHashMap<>();
    private final Map<String, Function<String, String>> inferResponders = new ConcurrentHashMap<>();
    private final Map<String, Set<String>> refusedPeptides = new ConcurrentHashMap<>();
    private final Set<String> refusingAll = ConcurrentHashMap.newKeySet();
    private final Map<String, Set<String>> failingNew = new ConcurrentHashMap<>();
    private final Map<String, String> failingNewBodies = new ConcurrentHashMap<>();
    private final Map<String, AtomicInteger> brokenFor = new ConcurrentHashMap<>();
    private final Map<String, Integer> maxAnswers = new ConcurrentHashMap<>();
    private final Map<String, String> metadata = new ConcurrentHashMap<>();
    private final Map<String, String> configs = new ConcurrentHashMap<>();
    private final Map<String, List<String>> requests = new ConcurrentHashMap<>();
    private final Map<String, AtomicInteger> gets = new ConcurrentHashMap<>();
    private final AtomicInteger failingGets = new AtomicInteger();

    private KoinaStubServer() throws IOException {
        server = HttpServer.create(new InetSocketAddress(InetAddress.getLoopbackAddress(), 0), 0);
        executor = Executors.newCachedThreadPool(r -> {
            Thread t = new Thread(r, "koina-stub");
            t.setDaemon(true);
            return t;
        });
        server.setExecutor(executor);
        server.createContext(ROOT, this::handle);
        server.start();
    }

    public static KoinaStubServer start() throws IOException {
        return new KoinaStubServer();
    }

    /** The Koina URL to set as {@code Constants.KoinaURL}; it ends with {@code /v2/models/}. */
    public String url() {
        return "http://" + server.getAddress().getHostString() + ":" + server.getAddress().getPort() + ROOT;
    }

    /** Answers {@code POST <model>/infer} with this body, at most {@code maxAnswers} times. */
    public void serveInfer(String model, String responseBody, int maxAnswers) {
        inferResponses.put(model, responseBody);
        this.maxAnswers.put(model, maxAnswers);
    }

    /** Answers each {@code POST <model>/infer} with HTTP 200 and the body the function gives for the request. */
    public void serveInfer(String model, Function<String, String> responder) {
        inferResponders.put(model, responder);
    }

    /**
     * Answers HTTP 400, as Koina does for a precursor its model cannot take, to any infer request of
     * the model whose peptide_sequences has this peptide.
     */
    public void refusePeptide(String model, String peptide) {
        refusedPeptides.computeIfAbsent(model, k -> ConcurrentHashMap.newKeySet()).add(peptide);
    }

    /** Answers HTTP 400 to every infer request of the model, as Koina does for an input it cannot take. */
    public void refuseAll(String model) {
        refusingAll.add(model);
    }

    /**
     * Answers HTTP 400, as Koina does for a failure of its own ({@link #failureOf}), to each infer
     * request of the model the first time it arrives, whatever else it would answer. The same request
     * sent again is answered as it would be otherwise.
     */
    public void failNewRequests(String model) {
        failNewRequests(model, failureOf(model));
    }

    /** {@link #failNewRequests(String)}, answering HTTP 400 with this body. */
    public void failNewRequests(String model, String body) {
        failingNewBodies.put(model, body);
        failingNew.put(model, ConcurrentHashMap.newKeySet());
    }

    /**
     * Answers the next {@code requests} infer requests of the model with HTTP 400 and a CUDA failure
     * ({@link #brokenInstanceOf}), whatever they hold, as a Koina model instance does once a CUDA
     * failure has broken it, until it is restarted. {@code Integer.MAX_VALUE} never restarts it.
     */
    public void breakInstance(String model, int requests) {
        brokenFor.put(model, new AtomicInteger(requests));
    }

    /** What the stub answers with HTTP 400 for a failure of Koina's own ({@link #failNewRequests}). */
    public static String failureOf(String model) {
        return "{\"error\":\"Failed to process the request(s) for model instance '" + model
                + "_0', message: RuntimeException: CUDA out of memory (stub)\"}";
    }

    /** What the stub answers with HTTP 400 for a broken model instance ({@link #breakInstance}), as Koina did. */
    public static String brokenInstanceOf(String model) {
        return "{\"error\":\"Failed to process the request(s) for model instance '" + model + "', message: "
                + "TritonModelException: Failed to record the event.: device-side assert triggered\\n\\nAt:\\n "
                + "/models/repo/" + model + "/1/model.py(164): predict_batch\\n\"}";
    }

    /** A failure of Koina's own whose answer does not say what failed, for {@link #failNewRequests(String, String)}. */
    public static String unnamedFailureOf(String model) {
        return "{\"error\":\"Failed to process the request(s) for model instance '" + model
                + "_0', message: stub instance is not ready\"}";
    }

    /** What the stub answers with HTTP 400 for a request with a peptide it refuses. */
    public static String refusalOf(String model, String peptide) {
        return "{\"error\":\"Failed to process the request(s) for model instance '" + model
                + "_0', message: ValueError: stub cannot take " + peptide + "\"}";
    }

    /** What the stub answers with HTTP 400 for a model that refuses every request. */
    public static String refusalOfAll(String model) {
        return "{\"error\":\"Failed to process the request(s) for model instance '" + model
                + "_0', message: ValueError: stub takes no request\"}";
    }

    public void serveMetadata(String model, String json) {
        metadata.put(model, json);
    }

    public void serveConfig(String model, String json) {
        configs.put(model, json);
    }

    /**
     * Answers {@code GET <model>} with metadata that has the spec's inputs and outputs, and
     * {@code GET <model>/config} with this max_batch_size, as Koina does for a model it serves.
     */
    public void serveModel(KoinaModelSpec spec, int maxBatchSize) {
        List<String> inputs = new ArrayList<>(List.of("peptide_sequences"));
        for (KoinaModelSpec.Input input : spec.inputs()) {
            inputs.add(input.koinaName());
        }
        serveModel(spec.name(), inputs, spec.outputs(), maxBatchSize);
    }

    /** {@link #serveModel(KoinaModelSpec, int)} with the inputs and outputs given. */
    public void serveModel(String model, List<String> inputs, List<String> outputs, int maxBatchSize) {
        serveMetadata(model, "{\"name\":\"" + model + "\",\"versions\":[\"1\"],\"platform\":\"stub\","
                + "\"inputs\":" + tensors(inputs, "[-1,1]") + ",\"outputs\":" + tensors(outputs, "[-1,-1]") + "}");
        serveConfig(model, "{\"name\":\"" + model + "\",\"platform\":\"\",\"backend\":\"python\","
                + "\"max_batch_size\":" + maxBatchSize + "}");
    }

    private static String tensors(List<String> names, String shape) {
        return names.stream()
                .map(name -> "{\"name\":\"" + name + "\",\"datatype\":\"FP32\",\"shape\":" + shape + "}")
                .collect(Collectors.joining(",", "[", "]"));
    }

    /** Answers the next {@code times} GETs with HTTP 503, as a server that is busy. */
    public void failNextGets(int times) {
        failingGets.set(times);
    }

    /** How many GETs of this path under {@code /v2/models/} arrived, e.g. {@code "UniSpec/config"}. */
    public int gets(String path) {
        AtomicInteger count = gets.get(path);
        return count == null ? 0 : count.get();
    }

    /** The infer bodies received for this model, in arrival order. */
    public List<String> requests(String model) {
        List<String> received = requests.get(model);
        if (received == null) {
            return List.of();
        }
        synchronized (received) {
            return new ArrayList<>(received);
        }
    }

    private void handle(HttpExchange exchange) throws IOException {
        try (exchange) {
            String path = exchange.getRequestURI().getPath().substring(ROOT.length());
            String method = exchange.getRequestMethod();
            if (method.equals("POST") && path.endsWith("/infer")) {
                String model = path.substring(0, path.length() - "/infer".length());
                String body;
                try (InputStream in = exchange.getRequestBody()) {
                    body = new String(in.readAllBytes(), StandardCharsets.UTF_8);
                }
                List<String> received = requests.computeIfAbsent(model,
                        k -> Collections.synchronizedList(new ArrayList<>()));
                int answered;
                synchronized (received) {
                    received.add(body);
                    answered = received.size();
                }
                String response = inferResponses.get(model);
                Function<String, String> responder = inferResponders.get(model);
                String refused = refusedPeptide(model, body);
                Set<String> seen = failingNew.get(model);
                AtomicInteger broken = brokenFor.get(model);
                if (seen != null && seen.add(body)) {
                    send(exchange, 400, failingNewBodies.get(model));
                } else if (broken != null
                        && broken.getAndUpdate(n -> n == Integer.MAX_VALUE ? n : Math.max(0, n - 1)) > 0) {
                    send(exchange, 400, brokenInstanceOf(model));
                } else if (refusingAll.contains(model)) {
                    send(exchange, 400, refusalOfAll(model));
                } else if (refused != null) {
                    send(exchange, 400, refusalOf(model, refused));
                } else if (responder != null) {
                    send(exchange, 200, responder.apply(body));
                } else if (response == null) {
                    send(exchange, 400, unknownModel(model));
                } else if (answered > maxAnswers.get(model)) {
                    send(exchange, 500, "{\"error\":\"stub: " + model + " was already answered "
                            + maxAnswers.get(model) + " times\"}");
                } else {
                    send(exchange, 200, response);
                }
            } else if (method.equals("GET") && failingGets.getAndUpdate(n -> Math.max(0, n - 1)) > 0) {
                gets.computeIfAbsent(path, k -> new AtomicInteger()).incrementAndGet();
                send(exchange, 503, "stub: busy");
            } else if (method.equals("GET") && path.endsWith("/config")) {
                gets.computeIfAbsent(path, k -> new AtomicInteger()).incrementAndGet();
                String model = path.substring(0, path.length() - "/config".length());
                String config = configs.get(model);
                send(exchange, config == null ? 400 : 200, config == null ? unknownModel(model) : config);
            } else if (method.equals("GET")) {
                gets.computeIfAbsent(path, k -> new AtomicInteger()).incrementAndGet();
                String json = metadata.get(path);
                send(exchange, json == null ? 400 : 200, json == null ? unknownModel(path) : json);
            } else {
                send(exchange, 404, unknownModel(path));
            }
        }
    }

    //the first peptide of the request the model is told to refuse, or null
    private String refusedPeptide(String model, String body) {
        Set<String> refused = refusedPeptides.get(model);
        if (refused == null) {
            return null;
        }
        for (String peptide : KoinaBatch.fromRequest(body).peptides()) {
            if (refused.contains(peptide)) {
                return peptide;
            }
        }
        return null;
    }

    private static String unknownModel(String model) {
        return "{\"error\":\"Request for unknown model: '" + model + "' is not found\"}";
    }

    private static void send(HttpExchange exchange, int status, String body) throws IOException {
        byte[] bytes = body.getBytes(StandardCharsets.UTF_8);
        exchange.getResponseHeaders().set("Content-Type", "application/json");
        exchange.sendResponseHeaders(status, bytes.length == 0 ? -1 : bytes.length);
        if (bytes.length > 0) {
            try (OutputStream out = exchange.getResponseBody()) {
                out.write(bytes);
            }
        }
    }

    @Override
    public void close() {
        server.stop(0);
        executor.shutdownNow();
    }
}
