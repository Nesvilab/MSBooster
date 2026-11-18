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
import static org.junit.jupiter.api.Assertions.assertNotNull;
import static org.junit.jupiter.api.Assertions.assertNull;
import static org.junit.jupiter.api.Assertions.assertTrue;

import allconstants.Constants;
import com.google.gson.JsonArray;
import com.google.gson.JsonObject;
import java.io.ByteArrayOutputStream;
import java.io.PrintStream;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.LinkedHashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.TreeMap;
import java.util.concurrent.ConcurrentHashMap;
import java.util.concurrent.ScheduledThreadPoolExecutor;
import java.util.stream.Collectors;
import java.util.stream.Stream;
import modelcallers.KoinaModelCaller;
import org.junit.jupiter.api.AfterEach;
import org.junit.jupiter.api.BeforeEach;
import org.junit.jupiter.api.Test;
import predictions.PredictionEntry;
import predictions.PredictionEntryHashMap;
import readers.datareaders.PinReader;
import readers.predictionreaders.KoinaLibReader;
import testutils.ExitTrap;
import testutils.KoinaGolden;
import testutils.KoinaGolden.Case;
import testutils.KoinaGolden.Psm;
import testutils.KoinaStubServer;
import writers.JSONWriter;

/**
 * Koina answers HTTP 400 to a request with a precursor its model cannot take. Such a request is never sent again as
 * it is: its parts are sent until the precursors Koina refuses are alone, and those are left without predictions and
 * scored as precursors the model cannot take, while the others are predicted as if the refused ones had never been
 * in the run. A precursor counts as refused only when Koina refuses it twice, since Koina also answers HTTP 400 for
 * a failure of its own. A model Koina refuses every request of ends the run, also while the best model is looked
 * for, but only once no request it may still predict is out.
 */
class KoinaBisectionTest {
    private static final String MS2_MODEL = "Prosit_2020_intensity_HCD";
    private static final String RT_MODEL = "Prosit_2019_irt";
    //a peptide the stub refuses, as the models are sent it
    private static final String MARKER = "GSWEDLAPK";

    private static final List<Psm> PSMS = List.of(
            new Psm("K.LGEPTIDEK.A", 2),
            new Psm("R.AVLDQGSTWPEM[15.9949]AGK.S", 2),
            new Psm("K.TEC[57.0215]LDHFAVNPR.G", 3),
            new Psm("K." + MARKER + ".T", 2),
            new Psm("R.GS[79.9663]PEDLHVTAR.L", 2),
            new Psm("K.n[42.0106]ADVTGLPHK.E", 2),
            new Psm("K.VFLEN[0.9840]GTLR.A", 2),
            new Psm("R.HEAPGLDVSTWKPEEVLR.S", 4),
            new Psm("K.GDUVLPAEK.A", 2), //residue U: never sent
            new Psm("K.ADGLR.T", 1),
            new Psm("K." + MARKER + ".T", 3));

    private long confirmPauseMs;

    @BeforeEach
    void setUp() {
        confirmPauseMs = KoinaBisection.confirmPauseMs;
        KoinaBisection.confirmPauseMs = 50; //not the 5 s a run waits before it sends a refused precursor again
    }

    @AfterEach
    void tearDown() {
        KoinaBisection.confirmPauseMs = confirmPauseMs;
    }

    @Test
    void aRefusedPrecursorIsFoundAndTheOthersArePredictedAsWithoutIt() throws Exception {
        //requests of up to 4 precursors, so the refused precursors share requests with others
        Run withMarker = run(MS2_MODEL, PSMS, 4, stub -> stub.refusePeptide(MS2_MODEL, MARKER), true);
        Run withoutMarker = run(MS2_MODEL, withoutMarker(), 4, stub -> {}, true);

        assertFalse(withMarker.klr.failed);
        assertEquals(2, withMarker.klr.refusedCount(), "the marker at charges 2 and 3");
        assertTrue(withMarker.klr.isRefused(MARKER, "2"));
        assertTrue(withMarker.klr.isRefused(MARKER, "3"));
        assertFalse(withMarker.klr.isRefused(MARKER, "4"));
        assertTrue(withMarker.printed.contains("Koina model " + MS2_MODEL + " refused 2 precursors (HTTP 400)"),
                withMarker.printed);
        assertTrue(withMarker.printed.contains(MARKER + " charge 2: "
                + KoinaStubServer.refusalOf(MS2_MODEL, MARKER)), withMarker.printed);

        //the others are predicted as in a run without the marker
        assertEquals(KoinaGolden.dump(withoutMarker.transferred), KoinaGolden.dump(without(withMarker.transferred)));
        //the marker is scored as a precursor the model cannot take
        for (String key : List.of(MARKER + "|2", MARKER + "|3")) {
            PredictionEntry pe = withMarker.transferred.get(key);
            assertNotNull(pe, key + " " + withMarker.transferred.keySet());
            assertEquals(1, pe.numFragments(), key);
            assertEquals(0f, pe.getIntensity(0), key);
        }

        //a refused request is never sent again as a whole
        for (String request : withMarker.written) {
            if (request.contains(MARKER)) {
                assertEquals(1, withMarker.sent.stream().filter(request::equals).count(), request);
            }
        }
        //every request without the marker is predicted
        for (String request : withMarker.sent) {
            assertTrue(request.contains(MARKER) || !withMarker.refusedBodies.contains(request), request);
        }
        //a precursor alone is refused twice before it counts as refused
        for (String charge : List.of("2", "3")) {
            String alone = sentAlone(withMarker.sent, MARKER, charge);
            assertEquals(2, withMarker.sent.stream().filter(alone::equals).count(), MARKER + " " + charge);
        }
    }

    @Test
    void aRefusedPeptideOfAnRtModelHasNoRtAtAnyCharge() throws Exception {
        //one request of every precursor: the first refusal comes before anything is predicted, so some of its
        //precursors are sent alone first, then it is cut
        Run withMarker = run(RT_MODEL, PSMS, 1000, stub -> stub.refusePeptide(RT_MODEL, MARKER), true);
        Run withoutMarker = run(RT_MODEL, withoutMarker(), 1000, stub -> {}, true);

        assertFalse(withMarker.klr.failed);
        assertEquals(1, withMarker.klr.refusedCount(), "an RT model is sent no charge, so the peptide is refused");
        assertTrue(withMarker.klr.isRefused(MARKER, "2"));
        assertTrue(withMarker.klr.isRefused(MARKER, "5"));
        assertTrue(withMarker.printed.contains("Koina model " + RT_MODEL + " refused 1 precursor (HTTP 400)"),
                withMarker.printed);
        assertTrue(withMarker.printed.contains("  " + MARKER + ": "), withMarker.printed);

        assertEquals(KoinaGolden.dump(withoutMarker.transferred), KoinaGolden.dump(without(withMarker.transferred)));
        assertEquals(0f, withMarker.transferred.get(MARKER + "|2").getRT());
        assertEquals(0f, withMarker.transferred.get(MARKER + "|3").getRT());
        assertEquals(1, withMarker.sent.stream().filter(withMarker.written.get(0)::equals).count(),
                "the refused request is sent once");
    }

    //the search cannot leave a model out, so a candidate Koina refuses every request of ends it too, after the
    //request, then 8 of its precursors alone, each twice
    @Test
    void aModelKoinaRefusesEveryRequestOfEndsTheRunAfterAFewRequests() throws Exception {
        for (boolean foundBest : new boolean[] {true, false}) {
            Run run = attempt(RT_MODEL, PSMS, 1000, stub -> stub.refuseAll(RT_MODEL), foundBest);
            assertNotNull(run.failure, "foundBest " + foundBest);
            String message = run.failure.getMessage();
            assertTrue(message.contains("Koina model " + RT_MODEL + " refused every request MSBooster sent it, "
                    + "with HTTP 400 " + KoinaStubServer.refusalOfAll(RT_MODEL)), message);
            assertTrue(message.contains("and each of the 8 of its precursors sent alone, twice"), message);
            assertEquals(1 + 2 * KoinaBisection.PROBES, run.sent.size(), "foundBest " + foundBest);
            assertFalse(run.printed.contains("goes on"), run.printed);
        }
    }

    @Test
    void aModelKoinaRefusesEveryPrecursorOfPredictsNothingAndEndsTheRun() throws Exception {
        //requests of one precursor each, all refused twice: each is a precursor Koina refuses, and none is left
        for (boolean foundBest : new boolean[] {true, false}) {
            Run run = attempt(RT_MODEL, PSMS.subList(0, 3), 1, stub -> stub.refuseAll(RT_MODEL), foundBest);
            assertNotNull(run.failure, "foundBest " + foundBest);
            assertTrue(run.failure.getMessage().contains("Koina model " + RT_MODEL + " refused all 3 precursors it "
                    + "was sent, with HTTP 400 " + KoinaStubServer.refusalOfAll(RT_MODEL) + ", so it predicts "
                    + "nothing"), run.failure.getMessage());
            assertEquals(6, run.sent.size(), "each precursor twice");
        }
    }

    //Koina refuses a request at once, before inference, while a request of hundreds of precursors takes seconds. A
    //request file all of whose precursors are refused, and its precursors sent alone, can so be back before another
    //file is: the model does not refuse every request while that file may still be predicted
    @Test
    void aModelIsNotSaidToRefuseEveryRequestWhileARequestIsStillOut() throws Exception {
        List<String> refused = List.of("LGEPTIDEK", "AGEPTIDEK", "SGEPTIDEK");
        List<String> slow = List.of("TGEPTIDEK", "VGEPTIDEK", "MGEPTIDEK");
        Path jsonDir = KoinaGolden.workDir("bisection-pending");
        Files.writeString(jsonDir.resolve("spectraRT_" + RT_MODEL + "_0rt.json"), rtRequest(refused));
        Files.writeString(jsonDir.resolve("spectraRT_" + RT_MODEL + "_1rt.json"), rtRequest(slow));

        KoinaGolden.Settings settings = KoinaGolden.Settings.apply(Case.of(RT_MODEL, "LUMOS", "HCD", "27", PSMS));
        ScheduledThreadPoolExecutor executor = new ScheduledThreadPoolExecutor(2);
        try (ExitTrap trap = ExitTrap.install(); KoinaStubServer stub = KoinaStubServer.start()) {
            stub.serveModel(KoinaModelRegistry.require(RT_MODEL), 1000);
            stub.serveInfer(RT_MODEL, request -> {
                if (request.contains(slow.get(0))) {
                    try {
                        Thread.sleep(1500);
                    } catch (InterruptedException e) {
                        Thread.currentThread().interrupt();
                    }
                }
                return answer(RT_MODEL, request);
            });
            for (String peptide : refused) {
                stub.refusePeptide(RT_MODEL, peptide);
            }
            Constants.KoinaURL = stub.url();
            Constants.foundBest = true;

            KoinaLibReader klr = new KoinaLibReader(RT_MODEL);
            new KoinaModelCaller().callModel(RT_MODEL, klr, jsonDir.toString(), executor, false, false);

            assertEquals(3, klr.refusedCount());
            for (String peptide : refused) {
                assertTrue(klr.isRefused(peptide, "2"), peptide);
            }
            for (String peptide : slow) {
                assertTrue(klr.getPreds().containsKey(peptide + "|2"), peptide + " " + klr.getPreds().keySet());
            }
            assertEquals(1, stub.requests(RT_MODEL).stream().filter(r -> r.contains(slow.get(0))).count());
        } finally {
            executor.shutdownNow();
            settings.restore();
        }
    }

    //Koina also answers HTTP 400 for a failure of its own, such as a model instance out of memory: a precursor that
    //is predicted when it is sent again is no refusal, and the run goes on as if Koina had predicted it at once. A
    //failure its answer names (CUDA) is sent again whole; one it does not name is cut, and confirmed when alone
    @Test
    void aRefusalKoinaDoesNotRepeatIsNone() throws Exception {
        for (String model : List.of(RT_MODEL, MS2_MODEL)) {
            for (String failure : List.of(KoinaStubServer.failureOf(model), KoinaStubServer.unnamedFailureOf(model))) {
                //each request is refused the first time it arrives, the probes and the halves too: had one refusal
                //counted, the RT model would refuse every request and the MS2 model every precursor
                for (int maxBatchSize : new int[] {1000, 4, 1}) {
                    Run failing = run(model, withoutMarker(), maxBatchSize,
                            stub -> stub.failNewRequests(model, failure), true);
                    Run answering = run(model, withoutMarker(), maxBatchSize, stub -> {}, true);
                    String where = model + " in requests of " + maxBatchSize + ", " + failure;

                    assertEquals(0, failing.klr.refusedCount(), where);
                    assertFalse(failing.printed.contains("refused"), where + "\n" + failing.printed);
                    assertEquals(KoinaGolden.dump(answering.transferred), KoinaGolden.dump(failing.transferred),
                            where);
                }
            }
        }
    }

    //a CUDA failure leaves a Koina model instance failing every request it gets, whatever the request holds, until
    //the instance is restarted (seen on Prosit_2025_intensity_MultiFrag: "device-side assert triggered"). Its answer
    //is no refusal: the request is sent again whole, as a response not received, and nothing is taken for refused
    @Test
    void aRequestKoinaFailsOnWithACudaErrorIsSentAgainWhole() throws Exception {
        for (String model : List.of(RT_MODEL, MS2_MODEL)) {
            Run failing = run(model, withoutMarker(), 4, stub -> stub.breakInstance(model, 3), true);
            Run answering = run(model, withoutMarker(), 4, stub -> {}, true);

            assertEquals(0, failing.klr.refusedCount(), model);
            assertFalse(failing.printed.contains("refused"), model + "\n" + failing.printed);
            assertEquals(KoinaGolden.dump(answering.transferred), KoinaGolden.dump(failing.transferred), model);
            //no request is cut: each one sent is a request file, and only the 3 failed ones are sent twice
            for (String request : failing.sent) {
                assertTrue(failing.written.contains(request), model + " sent a part of a request: " + request);
            }
            assertEquals(answering.sent.size() + 3, failing.sent.size(), model);
            //the log says why a request is sent again
            assertTrue(failing.printed.contains("device-side assert triggered"), model + "\n" + failing.printed);
        }
    }

    //a model instance that stays broken ends the run as a server that does not answer, not as refused precursors;
    //while the best model is looked for, the model is marked failed
    @Test
    void aModelKoinaKeepsFailingOnWithACudaErrorEndsTheRunWithoutRefusals() throws Exception {
        Run chosen = attempt(RT_MODEL, withoutMarker(), 4, stub -> stub.breakInstance(RT_MODEL, Integer.MAX_VALUE),
                true);
        assertNotNull(chosen.failure);
        assertTrue(chosen.failure.getMessage().contains("did not answer"), chosen.failure.getMessage());
        assertTrue(chosen.failure.getMessage().contains("device-side assert triggered"), chosen.failure.getMessage());
        assertEquals(0, chosen.klr.refusedCount());
        for (String request : chosen.sent) {
            assertTrue(chosen.written.contains(request), "sent a part of a request: " + request);
        }

        Run looking = attempt(RT_MODEL, withoutMarker(), 4, stub -> stub.breakInstance(RT_MODEL, Integer.MAX_VALUE),
                false);
        assertNull(looking.failure);
        assertTrue(looking.klr.failed);
        assertEquals(0, looking.klr.refusedCount());
    }

    //a model that remembers the UniMod codes it finds (PeptideFormatter.foundUnimods) may write a mass with two codes
    //one way in its request and the other way when transferKoinaPreds writes the peptide again
    @Test
    void aRefusalIsFoundWhicheverUniModCodeWritesItsMass() {
        for (String model : List.of("AlphaPept_ms2_generic", "Deeplc_hela_hf")) {
            KoinaModelSpec spec = KoinaModelRegistry.require(model);
            int[] codes = twoCodesOfOneMass(spec);
            KoinaLibReader klr = new KoinaLibReader(model);
            boolean sendsCharge = spec.sends(KoinaModelSpec.Input.PRECURSOR_CHARGES);

            assertTrue(klr.refuse("PEPTIDEK[UNIMOD:" + codes[0] + "]", sendsCharge ? "2" : null, "refused"));
            assertTrue(klr.isRefused("PEPTIDEK[UNIMOD:" + codes[1] + "]", "2"), model + " " + codes[1]);
            assertTrue(klr.isRefused("PEPTIDEK[UNIMOD:" + codes[0] + "]", "2"), model + " " + codes[0]);
            assertEquals(!sendsCharge, klr.isRefused("PEPTIDEK[UNIMOD:" + codes[1] + "]", "3"),
                    model + ": an RT model is sent no charge");
            assertFalse(klr.isRefused("PEPTIDEK", "2"), model);
            assertFalse(klr.refuse("PEPTIDEK[UNIMOD:" + codes[1] + "]", sendsCharge ? "2" : null, "again"),
                    model + ": the same precursor");
            assertEquals(1, klr.refusedCount(), model);
        }
    }

    @Test
    void aRequestIsCutIntoTheRequestsOfItsPrecursors() {
        String request = "{\"id\":\"0\",\"inputs\":["
                + "{\"name\":\"peptide_sequences\",\"shape\":[3,1],\"datatype\":\"BYTES\","
                + "\"data\":[[\"PEPTIDEK\"],[\"LGEPTIDEK\"],[\"AGEPTIDEK\"]]},"
                + "{\"name\":\"precursor_charges\",\"shape\":[3,1],\"datatype\":\"INT32\",\"data\":[[2],[3],[4]]},"
                + "{\"name\":\"collision_energies\",\"shape\":[3,1],\"datatype\":\"FP32\","
                + "\"data\":[[27.0],[27.5],[28]]}]}";
        assertEquals("{\"id\":\"0\",\"inputs\":["
                + "{\"name\":\"peptide_sequences\",\"shape\":[2,1],\"datatype\":\"BYTES\","
                + "\"data\":[[\"LGEPTIDEK\"],[\"AGEPTIDEK\"]]},"
                + "{\"name\":\"precursor_charges\",\"shape\":[2,1],\"datatype\":\"INT32\",\"data\":[[3],[4]]},"
                + "{\"name\":\"collision_energies\",\"shape\":[2,1],\"datatype\":\"FP32\","
                + "\"data\":[[27.5],[28]]}]}", KoinaBatch.sliceRequest(request, 1, 3));
        KoinaBatch batch = KoinaBatch.fromRequest(request).slice(1, 3);
        assertEquals(List.of("LGEPTIDEK", "AGEPTIDEK"), List.of(batch.peptides()));
        assertEquals(List.of("3", "4"), List.of(batch.charges()));
    }

    //what a run of a model on PSMs made
    private static final class Run {
        KoinaLibReader klr;
        PredictionEntryHashMap transferred;
        //the request files, as written
        List<String> written;
        //the infer bodies the stub received
        List<String> sent;
        //the bodies the stub refused
        Set<String> refusedBodies = new LinkedHashSet<>();
        String printed;
        //what ended the run, if anything did
        KoinaException failure;
    }

    @FunctionalInterface
    private interface StubSetup {
        void apply(KoinaStubServer stub);
    }

    /** {@link #attempt}, for a run that is to predict the PSMs: what ended it is thrown. */
    private static Run run(String model, List<Psm> psms, int maxBatchSize, StubSetup setup, boolean foundBest)
            throws Exception {
        Run run = attempt(model, psms, maxBatchSize, setup, foundBest);
        if (run.failure != null) {
            throw run.failure;
        }
        return run;
    }

    /**
     * Runs PSMs through a model as the main run does, on a stub that serves the model with this max_batch_size and
     * answers each request with predictions made from its precursors: PinReader.createJSON, JSONWriter,
     * KoinaModelCaller.callModel (verbose) and transferKoinaPreds onto the spectraRT_full.tsv lines of the PSMs. A
     * KoinaException that ends the run is kept in {@link Run#failure}, and nothing is transferred then.
     *
     * @param foundBest whether the run's models are chosen, or the best model is looked for
     */
    private static Run attempt(String model, List<Psm> psms, int maxBatchSize, StubSetup setup, boolean foundBest)
            throws Exception {
        Case c = Case.of(model, "LUMOS", "HCD", "27", psms);
        Path workDir = KoinaGolden.workDir("bisection-" + model);
        Path pin = KoinaGolden.writePin(psms, workDir.resolve("test.pin"));
        Path jsonDir = workDir.resolve("jsonFiles");
        Path fullTsv = workDir.resolve("spectraRT_full.tsv");

        Run run = new Run();
        KoinaGolden.Settings settings = KoinaGolden.Settings.apply(c);
        ScheduledThreadPoolExecutor executor = new ScheduledThreadPoolExecutor(2);
        PrintStream out = System.out;
        ByteArrayOutputStream printed = new ByteArrayOutputStream();
        try (ExitTrap trap = ExitTrap.install(); KoinaStubServer stub = KoinaStubServer.start()) {
            stub.serveModel(KoinaModelRegistry.require(model), maxBatchSize);
            stub.serveInfer(model, request -> answer(model, request));
            setup.apply(stub);
            Constants.KoinaURL = stub.url();
            Constants.foundBest = foundBest;
            System.setOut(new PrintStream(printed, true, StandardCharsets.UTF_8));

            Set<String> pinHits = new LinkedHashSet<>();
            PinReader reader = new PinReader(pin.toString());
            try {
                reader.createJSON(pinHits, new KoinaPreCheck(KoinaModelRegistry.require(model)));
            } finally {
                reader.close();
            }
            Set<String> hits = ConcurrentHashMap.newKeySet();
            hits.addAll(pinHits);
            new JSONWriter(model, hits, true, KoinaModelCheck.batchSize(model)).write(true, jsonDir.toString(), executor);
            try (Stream<Path> files = Files.list(jsonDir)) {
                run.written = new ArrayList<>();
                for (Path file : files.sorted().collect(Collectors.toList())) {
                    run.written.add(Files.readString(file, StandardCharsets.UTF_8));
                }
            }

            Set<String> fullHits = new LinkedHashSet<>();
            reader = new PinReader(pin.toString());
            try {
                reader.createFull(fullHits);
            } finally {
                reader.close();
            }
            Files.writeString(fullTsv, String.join("\n", fullHits) + "\n", StandardCharsets.UTF_8);

            run.klr = new KoinaLibReader(model);
            try {
                new KoinaModelCaller().callModel(model, run.klr, jsonDir.toString(), executor, true, false);
            } catch (KoinaException e) {
                run.failure = e;
            } finally {
                run.sent = stub.requests(model);
                for (String request : run.sent) {
                    if (answer(model, request) == null) {
                        run.refusedBodies.add(request);
                    }
                }
            }
            //as the main run does: what callModel leaves is transferred, whatever it marked, unless the model was given
            //up while the best model is looked for, which leaves the model out
            if (run.failure == null && !run.klr.failed) {
                run.transferred = new PredictionEntryHashMap();
                run.transferred.transferKoinaPreds(new ArrayList<>(List.of(run.klr)), fullTsv.toString());
            }
        } finally {
            System.setOut(out);
            executor.shutdownNow();
            settings.restore();
            run.printed = printed.toString(StandardCharsets.UTF_8);
            System.out.print(run.printed);
        }
        return run;
    }

    private static List<Psm> withoutMarker() {
        return PSMS.stream().filter(psm -> !psm.peptide().contains(MARKER)).collect(Collectors.toList());
    }

    //the predictions of the precursors other than the marker
    private static Map<String, PredictionEntry> without(PredictionEntryHashMap preds) {
        Map<String, PredictionEntry> others = new TreeMap<>(preds);
        others.keySet().removeIf(key -> key.startsWith(MARKER + "|"));
        return others;
    }

    //the request body sent with this precursor alone
    private static String sentAlone(List<String> sent, String peptide, String charge) {
        for (String request : sent) {
            KoinaBatch batch = KoinaBatch.fromRequest(request);
            if (batch.size() == 1 && batch.peptides()[0].equals(peptide) && batch.charges()[0].equals(charge)) {
                return request;
            }
        }
        throw new AssertionError(peptide + " charge " + charge + " was not sent alone");
    }

    //a request of these peptides to an RT model that is sent peptides only
    private static String rtRequest(List<String> peptides) {
        return "{\"id\":\"0\",\"inputs\":[{\"name\":\"peptide_sequences\",\"datatype\":\"BYTES\",\"shape\":["
                + peptides.size() + ",1],\"data\":["
                + peptides.stream().map(p -> "[\"" + p + "\"]").collect(Collectors.joining(",")) + "]}]}";
    }

    //two UniMod codes the spec reads as one mass, as the mass notation writes it
    private static int[] twoCodesOfOneMass(KoinaModelSpec spec) {
        Map<String, Integer> codeOfMass = new HashMap<>();
        for (Map.Entry<Integer, Double> e : new TreeMap<>(spec.modMassMap().unimodToMass()).entrySet()) {
            String mass = String.format("%.4f", e.getValue());
            Integer other = codeOfMass.putIfAbsent(mass, e.getKey());
            if (other != null) {
                return new int[] {other, e.getKey()};
            }
        }
        throw new AssertionError(spec.name() + " reads no two UniMod codes as one mass");
    }

    /**
     * The stub's answer to a request: predictions made from each precursor alone, so a precursor is predicted the
     * same whatever request it is in. Null for a request with the marker, which the stub refuses when told to.
     */
    private static String answer(String model, String request) {
        KoinaBatch batch = KoinaBatch.fromRequest(request);
        if (List.of(batch.peptides()).contains(MARKER)) {
            return null;
        }
        int n = batch.size();
        JsonArray outputs = new JsonArray();
        if (KoinaModelRegistry.require(model).kind() == KoinaModelSpec.Kind.RT) {
            JsonArray irt = new JsonArray();
            for (String peptide : batch.peptides()) {
                irt.add(peptide.length() * 2.5f + (peptide.hashCode() & 0xff) / 100f);
            }
            outputs.add(output("irt", "FP32", n, 1, irt));
        } else {
            JsonArray annotations = new JsonArray();
            JsonArray mzs = new JsonArray();
            JsonArray intensities = new JsonArray();
            String[] fragments = {"y1+1", "b2+1", "y2+1"};
            for (int i = 0; i < n; i++) {
                float h = ((batch.peptides()[i] + "|" + batch.charges()[i]).hashCode() & 0xffff) / 65535f;
                float[] values = {1f, 0.1f + 0.5f * h, 0.05f + 0.3f * h};
                for (int j = 0; j < fragments.length; j++) {
                    annotations.add(fragments[j]);
                    mzs.add(100f + 50f * j);
                    intensities.add(values[j]);
                }
            }
            outputs.add(output("annotation", "BYTES", n, 3, annotations));
            outputs.add(output("mz", "FP32", n, 3, mzs));
            outputs.add(output("intensities", "FP32", n, 3, intensities));
        }
        JsonObject response = new JsonObject();
        response.addProperty("id", "0");
        response.addProperty("model_name", model);
        response.add("outputs", outputs);
        return response.toString();
    }

    private static JsonObject output(String name, String datatype, int rows, int columns, JsonArray data) {
        JsonObject output = new JsonObject();
        output.addProperty("name", name);
        output.addProperty("datatype", datatype);
        JsonArray shape = new JsonArray();
        shape.add(rows);
        shape.add(columns);
        output.add("shape", shape);
        output.add("data", data);
        return output;
    }
}
