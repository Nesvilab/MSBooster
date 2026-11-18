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

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertFalse;
import static org.junit.jupiter.api.Assertions.assertThrows;
import static org.junit.jupiter.api.Assertions.assertTrue;

import allconstants.Constants;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.concurrent.ScheduledThreadPoolExecutor;
import java.util.concurrent.atomic.AtomicBoolean;
import java.util.concurrent.atomic.AtomicLong;
import modelcallers.KoinaModelCaller;
import org.junit.jupiter.api.AfterEach;
import org.junit.jupiter.api.BeforeEach;
import org.junit.jupiter.api.Test;
import org.junit.jupiter.api.io.TempDir;
import readers.predictionreaders.KoinaLibReader;
import testutils.ExitTrap;
import testutils.KoinaStubServer;

/**
 * KoinaTask sends a request again only when its response could not be received, also when the server did not
 * answer in time. A response that arrives but cannot be read ends the run at once with a KoinaException naming the
 * model, which reaches the thread that called KoinaModelCaller.callModel; Koina would answer the same request the
 * same way, so it is not resent.
 */
class KoinaTaskTest {
    private static final String MODEL = "Prosit_2019_intensity";
    private static final String REQUEST = "{\"id\":\"0\",\"inputs\":["
            + "{\"shape\":[2,1],\"data\":[[\"PEPTIDEK\"],[\"LGEPTIDEK\"]],\"datatype\":\"BYTES\","
            + "\"name\":\"peptide_sequences\"},"
            + "{\"shape\":[2,1],\"data\":[[2],[3]],\"datatype\":\"INT32\",\"name\":\"precursor_charges\"},"
            + "{\"shape\":[2,1],\"data\":[[30],[30]],\"datatype\":\"FP32\",\"name\":\"collision_energies\"}]}";

    @TempDir
    Path dir;

    private String koinaUrl;
    private Boolean foundBest;
    private ScheduledThreadPoolExecutor executor;

    @BeforeEach
    void setUp() throws Exception {
        koinaUrl = Constants.KoinaURL;
        foundBest = Constants.foundBest;
        executor = new ScheduledThreadPoolExecutor(2);
        KoinaModelCheck.forgetAnswers();
        Files.writeString(dir.resolve("spectraRT_" + MODEL + "_0ms2.json"), REQUEST, StandardCharsets.UTF_8);
    }

    @AfterEach
    void tearDown() {
        executor.shutdownNow();
        Constants.KoinaURL = koinaUrl;
        Constants.foundBest = foundBest;
        KoinaModelCheck.forgetAnswers();
    }

    @Test
    void aResponseThatCannotBeReadIsNotResent() throws Exception {
        try (ExitTrap trap = ExitTrap.install(); KoinaStubServer stub = KoinaStubServer.start()) {
            stub.serveInfer(MODEL, response("q2+1"), Constants.numKoinaAttempts);
            stub.serveModel(KoinaModelRegistry.require(MODEL), 1000);
            Constants.KoinaURL = stub.url();
            Constants.foundBest = false;
            KoinaLibReader klr = new KoinaLibReader(MODEL);

            KoinaException e = assertThrows(KoinaException.class, () -> new KoinaModelCaller()
                    .callModel(MODEL, klr, dir.toString(), executor, false, false));
            assertTrue(e.getMessage().contains(MODEL), e.getMessage());
            assertTrue(e.getMessage().contains("q2+1"), e.getMessage());
            assertEquals(1, stub.requests(MODEL).size(), "an unreadable response must not be resent");
        }
    }

    @Test
    void aReadableResponseIsReadAgainstTheRequestsPrecursors() throws Exception {
        try (ExitTrap trap = ExitTrap.install(); KoinaStubServer stub = KoinaStubServer.start()) {
            stub.serveInfer(MODEL, response("y2+1"), Constants.numKoinaAttempts);
            stub.serveModel(KoinaModelRegistry.require(MODEL), 1000);
            Constants.KoinaURL = stub.url();
            KoinaLibReader klr = new KoinaLibReader(MODEL);

            new KoinaModelCaller().callModel(MODEL, klr, dir.toString(), executor, false, false);
            assertEquals(1, stub.requests(MODEL).size());
            assertFalse(klr.failed);
            assertEquals(2, klr.getPreds().size(), "" + klr.getPreds().keySet());
            assertEquals(3, klr.getPreds().get("PEPTIDEK|2").numFragments());
            assertEquals(2, klr.getPreds().get("LGEPTIDEK|3").numFragments());
        }
    }

    @Test
    void aResponseThatCannotBeReceivedIsResentThenGivenUp() throws Exception {
        try (ExitTrap trap = ExitTrap.install(); KoinaStubServer stub = KoinaStubServer.start()) {
            stub.serveInfer(MODEL, response("y2+1"), 0); //HTTP 500 every time
            stub.serveModel(KoinaModelRegistry.require(MODEL), 1000);
            Constants.KoinaURL = stub.url();
            Constants.foundBest = false; //while the best model is looked for, the model is marked failed
            KoinaLibReader klr = new KoinaLibReader(MODEL);

            new KoinaModelCaller().callModel(MODEL, klr, dir.toString(), executor, false, false);
            assertEquals(Constants.numKoinaAttempts + 1, stub.requests(MODEL).size());
            assertTrue(klr.failed);
            assertTrue(klr.getPreds().isEmpty());
        }
    }

    @Test
    void aResponseThatIsNeverReceivedEndsTheRun() throws Exception {
        try (ExitTrap trap = ExitTrap.install(); KoinaStubServer stub = KoinaStubServer.start()) {
            stub.serveInfer(MODEL, response("y2+1"), 0); //HTTP 500 every time
            stub.serveModel(KoinaModelRegistry.require(MODEL), 1000);
            Constants.KoinaURL = stub.url();
            Constants.foundBest = true;
            KoinaLibReader klr = new KoinaLibReader(MODEL);

            KoinaException e = assertThrows(KoinaException.class, () -> new KoinaModelCaller()
                    .callModel(MODEL, klr, dir.toString(), executor, false, false));
            assertTrue(e.getMessage().contains(MODEL), e.getMessage());
            assertTrue(e.getMessage().contains("did not answer"), e.getMessage());
            assertEquals(Constants.numKoinaAttempts + 1, stub.requests(MODEL).size());
        }
    }

    //a server that takes the request and never answers would keep the run waiting forever: the response is waited
    //for as long as the requests' wait time, then counts as not received, and the request is sent again
    @Test
    void aServerThatDoesNotAnswerIsNotWaitedForForever() throws Exception {
        Integer initialWait = Constants.initialKoinaMillisecondsToWaitMs2;
        AtomicBoolean first = new AtomicBoolean(true);
        try (ExitTrap trap = ExitTrap.install(); KoinaStubServer stub = KoinaStubServer.start()) {
            stub.serveInfer(MODEL, request -> {
                if (first.getAndSet(false)) {
                    try {
                        Thread.sleep(30_000); //longer than the test waits
                    } catch (InterruptedException e) {
                        Thread.currentThread().interrupt();
                    }
                }
                return response("y2+1");
            });
            stub.serveModel(KoinaModelRegistry.require(MODEL), 1000);
            Constants.KoinaURL = stub.url();
            Constants.foundBest = true;
            Constants.initialKoinaMillisecondsToWaitMs2 = 500;
            KoinaLibReader klr = new KoinaLibReader(MODEL);

            long start = System.currentTimeMillis();
            new KoinaModelCaller().callModel(MODEL, klr, dir.toString(), executor, false, false);
            long elapsed = System.currentTimeMillis() - start;

            assertEquals(2, stub.requests(MODEL).size(), "the request that was not answered, then again");
            assertEquals(2, klr.getPreds().size(), "" + klr.getPreds().keySet());
            assertTrue(elapsed < 20_000, "waited " + elapsed + " ms");
        } finally {
            Constants.initialKoinaMillisecondsToWaitMs2 = initialWait;
        }
    }

    @Test
    void theWaitIsNeverShorterThanAtTheStart() {
        Integer initialWait = Constants.initialKoinaMillisecondsToWaitRtIm;
        try {
            Constants.initialKoinaMillisecondsToWaitRtIm = 30_000;
            KoinaLibReader klr = new KoinaLibReader("Prosit_2019_irt");
            AtomicLong waitTime = new AtomicLong(2_000); //as it is after many requests that took less than a second
            KoinaTask task = new KoinaTask(dir.resolve("a.json").toString(), "Prosit_2019_irt", klr, waitTime);
            assertEquals(30_000, task.timeoutMs());
            waitTime.set(45_000);
            assertEquals(45_000, task.timeoutMs());
        } finally {
            Constants.initialKoinaMillisecondsToWaitRtIm = initialWait;
        }
    }

    //a Prosit response for the request: PEPTIDEK y1+1, b2+1 and the given annotation; LGEPTIDEK y1+1 and b1+2
    private static String response(String thirdAnnotation) {
        return "{\"id\":\"0\",\"model_name\":\"" + MODEL + "\",\"model_version\":\"1\",\"outputs\":["
                + "{\"name\":\"annotation\",\"datatype\":\"BYTES\",\"shape\":[2,3],\"data\":"
                + "[\"y1+1\",\"b2+1\",\"" + thirdAnnotation + "\",\"y1+1\",\"b1+2\",\"y2+1\"]},"
                + "{\"name\":\"mz\",\"datatype\":\"FP32\",\"shape\":[2,3],\"data\":"
                + "[147.11281,227.10263,276.15509,147.11281,57.02602,-1.0]},"
                + "{\"name\":\"intensities\",\"datatype\":\"FP32\",\"shape\":[2,3],\"data\":"
                + "[1.0,0.25,0.5,0.5,0.75,-1.0]}]}";
    }
}
