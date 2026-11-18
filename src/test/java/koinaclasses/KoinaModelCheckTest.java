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

import static org.junit.jupiter.api.Assertions.assertDoesNotThrow;
import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertThrows;
import static org.junit.jupiter.api.Assertions.assertTrue;

import allconstants.Constants;
import java.io.ByteArrayOutputStream;
import java.io.PrintStream;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.Collections;
import java.util.LinkedHashSet;
import java.util.List;
import java.util.Set;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.Future;
import java.util.concurrent.ScheduledThreadPoolExecutor;
import java.util.stream.Collectors;
import java.util.stream.Stream;
import modelcallers.KoinaModelCaller;
import org.json.JSONObject;
import org.junit.jupiter.api.AfterEach;
import org.junit.jupiter.api.BeforeEach;
import org.junit.jupiter.api.Test;
import org.junit.jupiter.api.io.TempDir;
import readers.predictionreaders.KoinaLibReader;
import testutils.ExitTrap;
import testutils.KoinaStubServer;
import writers.JSONWriter;

/**
 * Before MSBooster writes or sends a request for a Koina model, it asks the server about the model once
 * ({@code GET <model>} and {@code GET <model>/config}). A model the server does not serve, or that does not answer
 * with an output its spec reads, ends the run with a message naming the model and the URL; a request carries at
 * most the server's max_batch_size precursors. A candidate of the best-model search that the server does not serve
 * ends the run too, since the search cannot leave a model out, and a server that does not answer is left to the
 * requests, which retry and give up as they always have.
 */
class KoinaModelCheckTest {
    private static final String RT_MODEL = "Prosit_2019_irt";
    private static final String MS2_MODEL = "UniSpec";

    @TempDir
    Path dir;

    private String koinaUrl;
    private Boolean foundBest;
    private Integer numKoinaAttempts;
    private String auxSpectraModel;

    @BeforeEach
    void setUp() {
        koinaUrl = Constants.KoinaURL;
        foundBest = Constants.foundBest;
        numKoinaAttempts = Constants.numKoinaAttempts;
        auxSpectraModel = Constants.auxSpectraModel;
        Constants.auxSpectraModel = "";
        KoinaModelCheck.forgetAnswers();
    }

    @AfterEach
    void tearDown() {
        Constants.KoinaURL = koinaUrl;
        Constants.foundBest = foundBest;
        Constants.numKoinaAttempts = numKoinaAttempts;
        Constants.auxSpectraModel = auxSpectraModel;
        KoinaModelCheck.forgetAnswers();
    }

    @Test
    void aServedModelIsAskedAboutOncePerRun() throws Exception {
        try (KoinaStubServer stub = KoinaStubServer.start()) {
            stub.serveModel(KoinaModelRegistry.require(MS2_MODEL), 1000);
            Constants.KoinaURL = stub.url();
            Constants.foundBest = true;

            ExecutorService threads = Executors.newFixedThreadPool(8);
            try {
                List<Future<Integer>> sizes = new ArrayList<>();
                for (int i = 0; i < 16; i++) {
                    sizes.add(threads.submit(() -> KoinaModelCheck.batchSize(MS2_MODEL)));
                }
                for (Future<Integer> size : sizes) {
                    assertEquals(100, size.get(), "the spec's batch size, which is below the server's");
                }
            } finally {
                threads.shutdownNow();
            }
            KoinaModelCheck.check(MS2_MODEL.toLowerCase());

            assertEquals(1, stub.gets(MS2_MODEL), "the metadata is asked for once");
            assertEquals(1, stub.gets(MS2_MODEL + "/config"), "the config is asked for once");
        }
    }

    @Test
    void theServersSmallerBatchSizeCutsTheRequestFiles() throws Exception {
        try (KoinaStubServer stub = KoinaStubServer.start()) {
            stub.serveModel(KoinaModelRegistry.require("Prosit_2020_intensity_HCD"), 40);
            stub.serveModel(KoinaModelRegistry.require("Prosit_2019_intensity"), 0); //no batching: not a limit
            stub.serveModel(KoinaModelRegistry.require("Prosit_2020_intensity_CID"), 5000);
            stub.serveMetadata(RT_MODEL, "{\"name\":\"" + RT_MODEL + "\",\"outputs\":[{\"name\":\"irt\"}]}"); //no config
            Constants.KoinaURL = stub.url();
            Constants.foundBest = true;

            int batchSize = KoinaModelCheck.batchSize("Prosit_2020_intensity_HCD");
            assertEquals(40, batchSize);
            assertEquals(1000, KoinaModelCheck.batchSize("Prosit_2019_intensity"));
            assertEquals(1000, KoinaModelCheck.batchSize("Prosit_2020_intensity_CID"));
            assertEquals(1000, KoinaModelCheck.batchSize(RT_MODEL));

            Set<String> entries = new LinkedHashSet<>();
            for (int i = 0; i < 100; i++) {
                entries.add("PEPTIDE" + "ACDEFGHIK".substring(0, 1 + i % 9) + "R" + i + ",2,30,QE,HCD");
            }
            ExecutorService executor = Executors.newFixedThreadPool(2);
            try {
                new JSONWriter("Prosit_2020_intensity_HCD", entries, false, batchSize)
                        .write(true, dir.toString(), executor);
            } finally {
                executor.shutdownNow();
            }
            assertEquals(List.of(20, 40, 40), requestSizes(dir));
        }
    }

    @Test
    void theRequestsSentAreCutAtTheServersBatchSize() throws Exception {
        try (ExitTrap trap = ExitTrap.install(); KoinaStubServer stub = KoinaStubServer.start()) {
            stub.serveModel(KoinaModelRegistry.require(RT_MODEL), 2);
            stub.serveInfer(RT_MODEL, "{\"model_name\":\"" + RT_MODEL + "\",\"outputs\":[{\"name\":\"irt\","
                    + "\"datatype\":\"FP32\",\"shape\":[2,1],\"data\":[10.5,20.5]}]}", 10);
            Constants.KoinaURL = stub.url();
            Constants.foundBest = true;

            Set<String> entries = new LinkedHashSet<>(List.of("PEPTIDEK,2,30,QE,HCD", "LGEPTIDEK,2,30,QE,HCD",
                    "AGEPTIDEK,2,30,QE,HCD", "SGEPTIDEK,2,30,QE,HCD"));
            ScheduledThreadPoolExecutor executor = new ScheduledThreadPoolExecutor(2);
            KoinaLibReader klr = new KoinaLibReader(RT_MODEL);
            try {
                new JSONWriter(RT_MODEL, entries, false, KoinaModelCheck.batchSize(RT_MODEL))
                        .write(true, dir.toString(), executor);
                new KoinaModelCaller().callModel(RT_MODEL, klr, dir.toString(), executor, false, false);
            } finally {
                executor.shutdownNow();
            }

            List<String> sent = stub.requests(RT_MODEL);
            assertEquals(2, sent.size(), "four precursors in requests of two");
            for (String request : sent) {
                assertEquals(2, new JSONObject(request).getJSONArray("inputs").getJSONObject(0)
                        .getJSONArray("data").length());
            }
            for (String peptide : List.of("PEPTIDEK", "LGEPTIDEK", "AGEPTIDEK", "SGEPTIDEK")) {
                assertTrue(klr.getPreds().containsKey(peptide + "|2"), peptide + " " + klr.getPreds().keySet());
            }
            assertEquals(1, stub.gets(RT_MODEL), "callModel does not ask again");
        }
    }

    @Test
    void aModelKoinaDoesNotServeEndsTheRunBeforeAnyRequest() throws Exception {
        try (ExitTrap trap = ExitTrap.install(); KoinaStubServer stub = KoinaStubServer.start()) {
            //an infer response, so that a request would be answered if one were sent
            stub.serveInfer(RT_MODEL, "{\"outputs\":[]}", 10);
            Constants.KoinaURL = stub.url();
            Constants.foundBest = true;
            Files.writeString(dir.resolve("spectraRT_" + RT_MODEL + "rt.json"), "{\"id\":\"0\",\"inputs\":[{\"name\":"
                    + "\"peptide_sequences\",\"datatype\":\"BYTES\",\"shape\":[1,1],\"data\":[[\"PEPTIDEK\"]]}]}");

            KoinaException e = assertThrows(KoinaException.class, () -> KoinaModelCheck.batchSize(RT_MODEL));
            assertTrue(e.getMessage().contains("does not serve the model " + RT_MODEL), e.getMessage());
            assertTrue(e.getMessage().contains(stub.url()), e.getMessage());
            assertTrue(e.getMessage().contains("HTTP 400"), e.getMessage());
            assertTrue(e.getMessage().contains("unknown model"), e.getMessage());

            ScheduledThreadPoolExecutor executor = new ScheduledThreadPoolExecutor(1);
            try {
                KoinaLibReader klr = new KoinaLibReader(RT_MODEL);
                KoinaException fromCall = assertThrows(KoinaException.class, () -> new KoinaModelCaller()
                        .callModel(RT_MODEL, klr, dir.toString(), executor, false, false));
                assertEquals(e.getMessage(), fromCall.getMessage());
            } finally {
                executor.shutdownNow();
            }
            assertTrue(stub.requests(RT_MODEL).isEmpty(), "no request goes to a model Koina does not serve");
            assertEquals(1, stub.gets(RT_MODEL));
        }
    }

    //the search cannot leave a model out: a candidate without predictions would end it later, with a message that
    //does not say why
    @Test
    void aCandidateKoinaDoesNotServeEndsTheSearchBeforeAnyRequest() throws Exception {
        try (ExitTrap trap = ExitTrap.install(); KoinaStubServer stub = KoinaStubServer.start()) {
            stub.serveInfer(RT_MODEL, "{\"outputs\":[]}", 10); //would be answered, were a request sent
            Constants.KoinaURL = stub.url();
            Constants.foundBest = false; //the best model is searched for
            Files.writeString(dir.resolve("spectraRT_" + RT_MODEL + "rt.json"), "{\"id\":\"0\",\"inputs\":[{\"name\":"
                    + "\"peptide_sequences\",\"datatype\":\"BYTES\",\"shape\":[1,1],\"data\":[[\"PEPTIDEK\"]]}]}");

            KoinaException e = assertThrows(KoinaException.class, () -> KoinaModelCheck.batchSize(RT_MODEL));
            assertTrue(e.getMessage().contains("Koina at " + stub.url() + " does not serve the model " + RT_MODEL),
                    e.getMessage());
            ScheduledThreadPoolExecutor executor = new ScheduledThreadPoolExecutor(1);
            try {
                KoinaLibReader klr = new KoinaLibReader(RT_MODEL);
                KoinaException fromCall = assertThrows(KoinaException.class, () -> new KoinaModelCaller()
                        .callModel(RT_MODEL, klr, dir.toString(), executor, false, false));
                assertEquals(e.getMessage(), fromCall.getMessage());
            } finally {
                executor.shutdownNow();
            }
            assertTrue(stub.requests(RT_MODEL).isEmpty(), "its requests would all fail, so none is sent");
            assertEquals(1, stub.gets(RT_MODEL), "asked about once");
        }
    }

    @Test
    void aServedModelPassesTheCheckInBothModes() throws Exception {
        try (KoinaStubServer stub = KoinaStubServer.start()) {
            stub.serveModel(KoinaModelRegistry.require(RT_MODEL), 1000);
            Constants.KoinaURL = stub.url();
            for (boolean searching : new boolean[] {true, false}) {
                Constants.foundBest = !searching;
                assertDoesNotThrow(() -> KoinaModelCheck.check(RT_MODEL));
            }
        }
    }

    @Test
    void aKoinaUrlWithoutV2ModelsIsNamedInTheMessage() throws Exception {
        try (KoinaStubServer stub = KoinaStubServer.start()) {
            String root = stub.url().substring(0, stub.url().length() - "v2/models/".length());
            Constants.KoinaURL = root;
            Constants.foundBest = true;

            KoinaException e = assertThrows(KoinaException.class, () -> KoinaModelCheck.check(MS2_MODEL));
            assertTrue(e.getMessage().contains("GET " + root + MS2_MODEL + " answered HTTP 404"), e.getMessage());
            assertTrue(e.getMessage().contains("ends with /v2/models/"), e.getMessage());
        }
    }

    @Test
    void aModelWithoutAnOutputItsSpecReadsEndsTheRun() throws Exception {
        try (KoinaStubServer stub = KoinaStubServer.start()) {
            stub.serveModel(MS2_MODEL, List.of("peptide_sequences"), List.of("intensities", "mz", "annotations"), 1000);
            Constants.KoinaURL = stub.url();

            for (boolean searching : new boolean[] {true, false}) {
                Constants.foundBest = !searching; //its responses could not be read, so the search cannot go on either
                KoinaException e = assertThrows(KoinaException.class, () -> KoinaModelCheck.batchSize(MS2_MODEL));
                assertTrue(e.getMessage().contains("The Koina model " + MS2_MODEL + " at " + stub.url()),
                        e.getMessage());
                assertTrue(e.getMessage().contains("the output annotation that MSBooster reads"), e.getMessage());
                assertTrue(e.getMessage().contains("It answers with intensities, mz, annotations"), e.getMessage());
            }
        }
    }

    @Test
    void anAnswerThatIsNotModelMetadataIsNoModel() throws Exception {
        try (KoinaStubServer stub = KoinaStubServer.start()) {
            stub.serveMetadata(RT_MODEL, "<html><body>Welcome</body></html>");
            Constants.KoinaURL = stub.url();
            Constants.foundBest = true;

            KoinaException e = assertThrows(KoinaException.class, () -> KoinaModelCheck.check(RT_MODEL));
            assertTrue(e.getMessage().contains("not model metadata: <html><body>Welcome</body></html>"),
                    e.getMessage());
        }
    }

    @Test
    void aBusyServerIsAskedAgain() throws Exception {
        try (KoinaStubServer stub = KoinaStubServer.start()) {
            stub.serveModel(KoinaModelRegistry.require(MS2_MODEL), 30);
            stub.failNextGets(1);
            Constants.KoinaURL = stub.url();
            Constants.foundBest = true;
            Constants.numKoinaAttempts = 1;

            assertEquals(30, KoinaModelCheck.batchSize(MS2_MODEL));
            assertEquals(2, stub.gets(MS2_MODEL), "one HTTP 503, then the metadata");
        }
    }

    @Test
    void aServerThatDoesNotAnswerIsLeftToTheRequests() throws Exception {
        try (KoinaStubServer stub = KoinaStubServer.start()) {
            stub.serveModel(KoinaModelRegistry.require(MS2_MODEL), 30);
            stub.failNextGets(100);
            Constants.KoinaURL = stub.url();
            Constants.foundBest = true;
            Constants.numKoinaAttempts = 1;

            String printed = stdout(() -> assertEquals(100, KoinaModelCheck.batchSize(MS2_MODEL),
                    "the spec's batch size; the requests retry and give up as they always have"));
            assertTrue(printed.contains("[INFO] - Could not ask Koina at " + stub.url() + " about the model "
                    + MS2_MODEL), printed);
            assertTrue(printed.contains("HTTP 503"), printed);
            assertEquals(2, stub.gets(MS2_MODEL), "asked numKoinaAttempts + 1 times");
        }

        //nothing listens at the port of a closed server
        String closedUrl;
        try (KoinaStubServer stub = KoinaStubServer.start()) {
            closedUrl = stub.url();
        }
        Constants.KoinaURL = closedUrl;
        Constants.numKoinaAttempts = 0;
        String printed = stdout(() -> assertEquals(1000, KoinaModelCheck.batchSize(RT_MODEL)));
        assertTrue(printed.contains("got no answer"), printed);
    }

    @Test
    void theKoinaUrlMayLackItsLastSlash() throws Exception {
        try (KoinaStubServer stub = KoinaStubServer.start()) {
            stub.serveModel(KoinaModelRegistry.require(MS2_MODEL), 50);
            Constants.KoinaURL = stub.url().substring(0, stub.url().length() - 1);
            Constants.foundBest = true;

            assertEquals(stub.url(), KoinaModelCheck.modelsUrl());
            assertEquals(50, KoinaModelCheck.batchSize(MS2_MODEL));
            assertEquals(1, stub.gets(MS2_MODEL));
        }
    }

    @Test
    void withoutAKoinaUrlNothingIsAsked() {
        Constants.KoinaURL = "";
        Constants.foundBest = true;
        assertEquals(100, KoinaModelCheck.batchSize(MS2_MODEL));
    }

    @Test
    void aKoinaUrlThatIsNotAUrlEndsTheRun() {
        Constants.KoinaURL = "koina.wilhelmlab.org/v2/models/";
        for (boolean searching : new boolean[] {true, false}) {
            Constants.foundBest = !searching;
            KoinaException e = assertThrows(KoinaException.class, () -> KoinaModelCheck.check(MS2_MODEL));
            assertTrue(e.getMessage().contains("The Koina URL koina.wilhelmlab.org/v2/models/ is not a URL"),
                    e.getMessage());
        }
    }

    @Test
    void theStubAnswersAsKoinaDoes() throws Exception {
        try (KoinaStubServer stub = KoinaStubServer.start()) {
            KoinaModelCheck.Answer answer = KoinaModelCheck.ask("NoSuchModel", stub.url());
            assertEquals(KoinaModelCheck.Status.NOT_SERVED, answer.status());
            //Koina answers HTTP 400 and this body for a model it does not serve
            assertTrue(answer.problem().endsWith("answered HTTP 400 "
                    + "{\"error\":\"Request for unknown model: 'NoSuchModel' is not found\"}"), answer.problem());
        }
    }

    //the precursors in each request file of the folder, smallest first
    private static List<Integer> requestSizes(Path dir) throws Exception {
        List<Integer> sizes = new ArrayList<>();
        try (Stream<Path> files = Files.list(dir)) {
            for (Path file : files.filter(f -> f.toString().endsWith(".json")).collect(Collectors.toList())) {
                JSONObject request = new JSONObject(Files.readString(file));
                sizes.add(request.getJSONArray("inputs").getJSONObject(0).getJSONArray("data").length());
            }
        }
        Collections.sort(sizes);
        return sizes;
    }

    //what the code prints to System.out
    private static String stdout(Runnable code) {
        PrintStream out = System.out;
        ByteArrayOutputStream bytes = new ByteArrayOutputStream();
        try {
            System.setOut(new PrintStream(bytes, true, StandardCharsets.UTF_8));
            code.run();
        } finally {
            System.setOut(out);
        }
        return bytes.toString(StandardCharsets.UTF_8);
    }
}
