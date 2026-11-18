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
import static org.junit.jupiter.api.Assertions.fail;

import java.nio.file.Path;
import java.util.ArrayList;
import java.util.List;
import java.util.TreeSet;
import java.util.stream.Collectors;
import java.util.stream.Stream;
import org.junit.jupiter.api.DynamicTest;
import org.junit.jupiter.api.Test;
import org.junit.jupiter.api.TestFactory;
import testutils.KoinaGolden;
import testutils.KoinaGolden.Case;
import testutils.KoinaGolden.Outcome;

/**
 * Every Koina model keeps sending the request it sends today and keeps making the same predictions
 * of the response Koina gave to it. The goldens were made by today's code from a response recorded
 * once from the public Koina server; this test serves that response from a local stub, so it never
 * calls Koina. See {@link KoinaGolden} for the folder layout and for how to update or re-record a
 * model's goldens on purpose ({@code -Dkoina.golden.update=<Model>} or
 * {@code -Dkoina.golden.record=<Model>}).
 */
class KoinaGoldenTest {
    @TestFactory
    Stream<DynamicTest> everyModelMatchesItsGolden() throws Exception {
        return KoinaGolden.fixtureModels().stream()
                .map(model -> DynamicTest.dynamicTest(model, () -> check(model)));
    }

    private void check(String model) throws Exception {
        Case c = Case.load(model);
        Path workDir = KoinaGolden.workDir(model);
        boolean record = KoinaGolden.requested(KoinaGolden.RECORD_PROPERTY, model);
        boolean update = record || KoinaGolden.requested(KoinaGolden.UPDATE_PROPERTY, model);

        KoinaGolden.LiveKoina live = record ? KoinaGolden.LiveKoina.fromProperty() : null;
        KoinaGolden.Responder responder = record ? live : (m, request) -> c.read(KoinaGolden.RESPONSE);
        Outcome outcome = KoinaGolden.run(c, workDir, responder);

        if (update) {
            if (record) {
                c.write(KoinaGolden.RESPONSE, live.responses.get(model));
                c.write(KoinaGolden.METADATA, live.metadata.get(model));
                c.write(KoinaGolden.CONFIG, live.configs.get(model));
            }
            c.write(KoinaGolden.REQUEST, outcome.request);
            c.write(KoinaGolden.PARSED, outcome.parsed);
            c.write(KoinaGolden.TRANSFERRED, outcome.transferred);
            System.out.println("[KoinaGoldenTest] " + (record ? "recorded" : "updated") + " the goldens of "
                    + model + " in " + c.dir);
        }

        String requestDiff = KoinaGolden.compareRequests(c.read(KoinaGolden.REQUEST), outcome.request);
        if (requestDiff != null) {
            fail(model + ": the request differs from " + KoinaGolden.REQUEST + " (written to "
                    + KoinaGolden.writeActual(model, outcome) + "): " + requestDiff);
        }
        assertEquals(List.of(outcome.request), outcome.sent,
                model + ": KoinaModelCaller must send the request file once, as written");
        assertFalse(outcome.failed, model + ": KoinaModelCaller could not receive the recorded response "
                + "from the stub server");

        List<String> diffs = new ArrayList<>();
        String parsedDiff = KoinaGolden.compareDumps(c.read(KoinaGolden.PARSED), outcome.parsed);
        if (parsedDiff != null) {
            diffs.add(KoinaGolden.PARSED + " " + parsedDiff);
        }
        String transferredDiff = KoinaGolden.compareDumps(c.read(KoinaGolden.TRANSFERRED), outcome.transferred);
        if (transferredDiff != null) {
            diffs.add(KoinaGolden.TRANSFERRED + " " + transferredDiff);
        }
        if (!diffs.isEmpty()) {
            fail(model + ": the predictions differ from the golden (written to "
                    + KoinaGolden.writeActual(model, outcome) + "):\n" + String.join("\n", diffs));
        }
    }

    @Test
    void theComparisonsCatchAChangedPrediction() throws Exception {
        Case c = Case.load("Prosit_2019_intensity");
        String golden = c.read(KoinaGolden.PARSED).replace("\r\n", "\n");
        assertNull(KoinaGolden.compareDumps(golden, golden));

        //one intensity off by a part in a million
        int start = golden.indexOf("\nF\t") + 1;
        int end = golden.indexOf('\n', start);
        String[] fields = golden.substring(start, end).split("\t", -1);
        fields[2] = Float.toString(Float.parseFloat(fields[2]) * 1.000001f);
        String changed = golden.substring(0, start) + String.join("\t", fields) + golden.substring(end);
        assertNotNull(KoinaGolden.compareDumps(golden, changed), "a changed intensity must not pass");
        //a lost fragment
        String shorter = golden.substring(0, start) + golden.substring(end + 1);
        assertNotNull(KoinaGolden.compareDumps(golden, shorter), "a lost fragment must not pass");

        String request = c.read(KoinaGolden.REQUEST);
        assertNull(KoinaGolden.compareRequests(request, request.replace("    ", "  ")),
                "whitespace is not part of a request");
        assertNotNull(KoinaGolden.compareRequests(request, request.replace("LGEPTIDEK", "LGEPTIDEQ")),
                "a changed peptide must not pass");
    }

    //every model of the registry, the one list of the Koina models MSBooster can call
    @Test
    void everyKoinaModelHasACompleteGolden() throws Exception {
        TreeSet<String> fixtures = new TreeSet<>(KoinaGolden.fixtureModels());
        List<String> registry = KoinaModelRegistry.all().stream()
                .map(KoinaModelSpec::name)
                .collect(Collectors.toList());
        List<String> missing = registry.stream()
                .filter(model -> !fixtures.contains(model))
                .sorted()
                .collect(Collectors.toList());
        assertTrue(missing.isEmpty(), "Koina models without a golden: " + missing + ". Add "
                + KoinaGolden.fixtureRoot() + "/<Model>/" + KoinaGolden.CASE + " and " + KoinaGolden.PSMS
                + ", then record it with -D" + KoinaGolden.RECORD_PROPERTY + "=<Model>");

        List<String> stale = fixtures.stream()
                .filter(model -> !registry.contains(model))
                .collect(Collectors.toList());
        assertTrue(stale.isEmpty(), "golden folders that name no Koina model: " + stale);

        List<String> incomplete = new ArrayList<>();
        for (String model : fixtures) {
            List<String> files = KoinaGolden.missingFiles(model);
            if (!files.isEmpty()) {
                incomplete.add(model + " " + files);
            }
        }
        assertTrue(incomplete.isEmpty(), "golden folders missing files: " + incomplete
                + ". Record them with -D" + KoinaGolden.RECORD_PROPERTY + "=<Model>");
    }
}
