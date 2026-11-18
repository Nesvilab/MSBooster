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

import allconstants.Constants;
import allconstants.NceConstants;
import com.google.gson.JsonArray;
import com.google.gson.JsonElement;
import com.google.gson.JsonObject;
import com.google.gson.JsonParser;
import koinaclasses.KoinaModelCheck;
import koinaclasses.KoinaModelRegistry;
import koinaclasses.KoinaPreCheck;
import modelcallers.KoinaModelCaller;
import peptideptmformatting.PeptideFormatter;
import predictions.PredictionEntry;
import predictions.PredictionEntryHashMap;
import readers.datareaders.PinReader;
import readers.predictionreaders.KoinaLibReader;
import writers.JSONWriter;

import java.io.IOException;
import java.net.URI;
import java.net.http.HttpClient;
import java.net.http.HttpRequest;
import java.net.http.HttpResponse;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.time.Duration;
import java.util.ArrayList;
import java.util.Comparator;
import java.util.HashMap;
import java.util.HashSet;
import java.util.LinkedHashSet;
import java.util.List;
import java.util.Map;
import java.util.Properties;
import java.util.Set;
import java.util.TreeMap;
import java.util.TreeSet;
import java.util.concurrent.ConcurrentHashMap;
import java.util.concurrent.ScheduledThreadPoolExecutor;
import java.util.stream.Collectors;
import java.util.stream.Stream;

/**
 * Golden runs of the Koina models: what today's code sends to a model for a fixed set of synthetic
 * PSMs, and what it makes of the response Koina gave once for them. {@code KoinaGoldenTest} runs
 * every model folder under {@code src/test/resources/koina/} through {@link #run} and compares the
 * result with the files in that folder, so a refactor of the Koina code has to keep the requests
 * and the predictions as they are.
 *
 * <p>A model folder holds:
 * <ul>
 *   <li>{@value #CASE}: the run settings (instrument, fragmentation, nce);</li>
 *   <li>{@value #PSMS}: the PSMs, one pin-style peptide (MSFragger mass notation, flanking
 *       residues included) and charge per line, with a note on what each one is there for;</li>
 *   <li>{@value #REQUEST}: the request JSONWriter writes for them;</li>
 *   <li>{@value #RESPONSE}: the body Koina answered that request with, byte for byte;</li>
 *   <li>{@value #METADATA} and {@value #CONFIG}: what Koina answered {@code GET <model>} and
 *       {@code GET <model>/config} with when the response was recorded;</li>
 *   <li>{@value #PARSED}: the predictions KoinaModelCaller parses from the response, keyed as the
 *       KoinaLibReader stores them (the model's own peptide, back in mass notation);</li>
 *   <li>{@value #TRANSFERRED}: those predictions after transferKoinaPreds moved them onto the PSMs'
 *       own peptides, which is what the rest of MSBooster reads.</li>
 * </ul>
 *
 * <p>Nothing here calls Koina unless asked to. The goldens change only on purpose:
 * <ul>
 *   <li>{@code -D}{@value #UPDATE_PROPERTY}{@code =<Model>[,<Model>...]|all} rewrites
 *       {@value #REQUEST}, {@value #PARSED} and {@value #TRANSFERRED} from the current code and the
 *       stored response, offline. Use it when a change to the output is deliberate.</li>
 *   <li>{@code -D}{@value #RECORD_PROPERTY}{@code =<Model>[,<Model>...]|all} also sends the request
 *       to the live Koina server ({@code -D}{@value #URL_PROPERTY}, default
 *       {@value #DEFAULT_KOINA_URL}) and stores its answer, metadata and config. Use it for a new
 *       model folder, or when the request changes, since a stored response only belongs to the
 *       request it answered. Only synthetic peptides go into these folders.</li>
 * </ul>
 * For example: {@code mvn -B test -Dtest=KoinaGoldenTest -Dkoina.golden.update=UniSpec}. Review the
 * diff of the folder before committing it. A new model gets a folder with {@value #CASE} and
 * {@value #PSMS}, and is then recorded; KoinaGoldenTest fails for a Koina model without one.
 *
 * <p>Each run starts from a cleared {@code PeptideFormatter.foundUnimods}, the case's instrument,
 * fragmentation and NCE (as the model's calibrated NCE too), and puts every setting it touches back
 * afterwards. The PSMs are read from one pin file, the way PeptideFileCreator reads each pin.
 */
public final class KoinaGolden {
    public static final String CASE = "case.properties";
    public static final String PSMS = "psms.tsv";
    public static final String REQUEST = "request.json";
    public static final String RESPONSE = "response.json";
    public static final String METADATA = "metadata.json";
    public static final String CONFIG = "config.json";
    public static final String PARSED = "parsed.tsv";
    public static final String TRANSFERRED = "transferred.tsv";
    public static final List<String> RECORDED_FILES = List.of(REQUEST, RESPONSE, METADATA, CONFIG, PARSED,
            TRANSFERRED);

    public static final String UPDATE_PROPERTY = "koina.golden.update";
    public static final String RECORD_PROPERTY = "koina.golden.record";
    public static final String URL_PROPERTY = "koina.golden.url";
    public static final String DEFAULT_KOINA_URL = "https://koina.wilhelmlab.org:443/v2/models/";

    private static final String DUMP_HEADER =
            "# P\tkey\tRT\tIM\tfragments\tdaltonMatching\n" +
            "# F\tmz\tintensity\tionType\tfragNum\tcharge\tfullAnnotation\tisotope\n";

    private KoinaGolden() {}

    /** The folder the model folders live in, in the source tree so that updates land there. */
    public static Path fixtureRoot() {
        return Paths.get(System.getProperty("basedir", System.getProperty("user.dir")))
                .resolve("src").resolve("test").resolve("resources").resolve("koina");
    }

    /**
     * A scratch folder for one model's run, under target/. Not a JUnit temp folder: the Koina code
     * leaves spectraRT_full.tsv open (PredictionEntryHashMap.transferKoinaPreds), and Windows refuses
     * to delete an open file.
     */
    public static Path workDir(String model) throws IOException {
        Path dir = Paths.get(System.getProperty("basedir", System.getProperty("user.dir")))
                .resolve("target").resolve("koina-golden-work").resolve(model);
        if (Files.exists(dir)) {
            try (Stream<Path> files = Files.walk(dir)) {
                files.sorted(Comparator.reverseOrder()).forEach(p -> p.toFile().delete());
            }
        }
        Files.createDirectories(dir);
        return dir;
    }

    /** Every model folder, by name. */
    public static List<String> fixtureModels() throws IOException {
        try (Stream<Path> dirs = Files.list(fixtureRoot())) {
            return dirs.filter(d -> Files.isRegularFile(d.resolve(CASE)))
                    .map(d -> d.getFileName().toString())
                    .sorted()
                    .collect(Collectors.toList());
        }
    }

    /** Whether -D{@code property} names this model (or is "all"). */
    public static boolean requested(String property, String model) {
        String value = System.getProperty(property, "").trim();
        if (value.isEmpty()) {
            return false;
        }
        for (String name : value.split(",")) {
            if (name.trim().equalsIgnoreCase("all") || name.trim().equals(model)) {
                return true;
            }
        }
        return false;
    }

    public record Psm(String peptide, int charge) {}

    /** One model's PSMs and run settings. */
    public static final class Case {
        public final String model;
        public final Path dir;
        public final String instrument;
        public final String fragmentation;
        public final String nce;
        public final List<Psm> psms;

        private Case(String model, Path dir, String instrument, String fragmentation, String nce,
                     List<Psm> psms) {
            this.model = model;
            this.dir = dir;
            this.instrument = instrument;
            this.fragmentation = fragmentation;
            this.nce = nce;
            this.psms = psms;
        }

        /** A case made in a test, with no model folder: it cannot be read from or written to. */
        public static Case of(String model, String instrument, String fragmentation, String nce, List<Psm> psms) {
            return new Case(model, null, instrument, fragmentation, nce, List.copyOf(psms));
        }

        public static Case load(String model) throws IOException {
            Path dir = fixtureRoot().resolve(model);
            Properties settings = new Properties();
            try (var reader = Files.newBufferedReader(dir.resolve(CASE), StandardCharsets.UTF_8)) {
                settings.load(reader);
            }
            List<Psm> psms = new ArrayList<>();
            for (String line : Files.readAllLines(dir.resolve(PSMS), StandardCharsets.UTF_8)) {
                if (line.isBlank() || line.startsWith("#") || line.startsWith("peptide\t")) {
                    continue;
                }
                String[] fields = line.split("\t");
                psms.add(new Psm(fields[0].trim(), Integer.parseInt(fields[1].trim())));
            }
            return new Case(model, dir, required(settings, "instrument", dir),
                    required(settings, "fragmentation", dir), required(settings, "nce", dir), psms);
        }

        private static String required(Properties settings, String key, Path dir) {
            String value = settings.getProperty(key);
            if (value == null || value.isBlank()) {
                throw new IllegalArgumentException(dir.resolve(CASE) + " has no " + key);
            }
            return value.trim();
        }

        public String read(String file) throws IOException {
            return Files.readString(dir.resolve(file), StandardCharsets.UTF_8);
        }

        public void write(String file, String content) throws IOException {
            Files.writeString(dir.resolve(file), content, StandardCharsets.UTF_8);
        }
    }

    /** Gives the infer response for a request: the stored one, or the live server's. */
    @FunctionalInterface
    public interface Responder {
        String respond(String model, String request) throws Exception;
    }

    /** What the current code made of a case. */
    public static final class Outcome {
        /** The request file JSONWriter wrote, as text. */
        public String request;
        /** The bodies KoinaModelCaller sent to the stub server. */
        public List<String> sent;
        /** Whether the KoinaLibReader was marked failed. */
        public boolean failed;
        /** What the request writer left out of the PSMs, and why. */
        public KoinaPreCheck preCheck;
        public String parsed;
        public String transferred;
    }

    /**
     * Runs a case through today's code: the pin's PSMs through PinReader.createJSON and JSONWriter
     * (as PeptideFileCreator does for one pin, with the batch size KoinaModelCheck finds), the
     * request to {@code responder}, its answer back through KoinaModelCaller.callModel, and the
     * parsed predictions through PredictionEntryHashMap.transferKoinaPreds onto the
     * spectraRT_full.tsv lines PinReader.createFull writes for the same pin. Koina is a stub server
     * from the start: it answers {@code GET <model>} and {@code GET <model>/config} with the recorded
     * metadata and config (the live ones when recording), and the infer request with the
     * responder's answer. Where the Koina code calls System.exit, the run fails instead (see
     * {@link ExitTrap}); the code prints why before it exits. A response the code cannot read, or a
     * model whose metadata lacks an output the spec reads, ends the run with a KoinaException.
     */
    public static Outcome run(Case c, Path workDir, Responder responder) throws Exception {
        Path pin = writePin(c, workDir.resolve("golden.pin"));
        Path jsonDir = workDir.resolve("jsonFiles");
        Path fullTsv = workDir.resolve("spectraRT_full.tsv");
        Outcome outcome = new Outcome();

        Settings settings = Settings.apply(c);
        ScheduledThreadPoolExecutor executor = new ScheduledThreadPoolExecutor(Math.max(1, Constants.numThreads));
        try (ExitTrap trap = ExitTrap.install(); KoinaStubServer stub = KoinaStubServer.start()) {
            serveRecordedMetadata(stub, c, responder);
            Constants.KoinaURL = stub.url();

            //request: PeptideFileCreator.createPeptideFile for a Koina model, with one pin file
            Set<String> pinHits = new LinkedHashSet<>();
            outcome.preCheck = new KoinaPreCheck(KoinaModelRegistry.require(c.model));
            PinReader reader = new PinReader(pin.toString());
            try {
                reader.createJSON(pinHits, outcome.preCheck);
            } finally {
                reader.close();
            }
            Set<String> hSetHits = ConcurrentHashMap.newKeySet();
            hSetHits.addAll(pinHits);
            new JSONWriter(c.model, hSetHits, true, KoinaModelCheck.batchSize(c.model))
                    .write(true, jsonDir.toString(), executor);
            List<Path> written;
            try (Stream<Path> files = Files.list(jsonDir)) {
                written = files.filter(p -> p.toString().endsWith(".json")).sorted().collect(Collectors.toList());
            }
            if (written.size() != 1) {
                throw new AssertionError(c.model + ": a golden case must fit one request file, but JSONWriter wrote "
                        + written.size() + " (" + written + ")");
            }
            outcome.request = Files.readString(written.get(0), StandardCharsets.UTF_8);

            //spectraRT_full.tsv: PeptideFileCreator.createPeptideFile with "createFull"
            Set<String> fullHits = new LinkedHashSet<>();
            reader = new PinReader(pin.toString());
            try {
                reader.createFull(fullHits);
            } finally {
                reader.close();
            }
            StringBuilder full = new StringBuilder();
            for (String hit : fullHits) {
                full.append(hit).append('\n');
            }
            Files.writeString(fullTsv, full.toString(), StandardCharsets.UTF_8);

            //predictions, from the stub server
            stub.serveInfer(c.model, responder.respond(c.model, outcome.request), Constants.numKoinaAttempts);
            KoinaLibReader klr = new KoinaLibReader(c.model);
            new KoinaModelCaller().callModel(c.model, klr, jsonDir.toString(), executor, false, false);
            outcome.sent = stub.requests(c.model);
            outcome.failed = klr.failed;
            outcome.parsed = dump(klr.getPreds());

            PredictionEntryHashMap transferred = new PredictionEntryHashMap();
            transferred.transferKoinaPreds(new ArrayList<>(List.of(klr)), fullTsv.toString());
            outcome.transferred = dump(transferred);
        } catch (ExitTrap.ExitCalled e) {
            throw new AssertionError(c.model + ": the Koina code called System.exit(" + e.status
                    + "); the message it printed before is in the test output", e);
        } finally {
            executor.shutdownNow();
            settings.restore();
        }
        return outcome;
    }

    //lets the stub answer GET <model> and GET <model>/config as Koina did when the response was recorded, or as
    //it does now when it is recorded
    private static void serveRecordedMetadata(KoinaStubServer stub, Case c, Responder responder)
            throws IOException, InterruptedException {
        if (responder instanceof LiveKoina) {
            LiveKoina live = (LiveKoina) responder;
            live.fetchModel(c.model);
            stub.serveMetadata(c.model, live.metadata.get(c.model));
            stub.serveConfig(c.model, live.configs.get(c.model));
            return;
        }
        if (Files.isRegularFile(c.dir.resolve(METADATA))) {
            stub.serveMetadata(c.model, c.read(METADATA));
        }
        if (Files.isRegularFile(c.dir.resolve(CONFIG))) {
            stub.serveConfig(c.model, c.read(CONFIG));
        }
    }

    private static Path writePin(Case c, Path pin) throws IOException {
        return writePin(c.psms, pin);
    }

    /** A pin file of these PSMs, one target PSM of rank 1 per scan, as the golden runs read them. */
    public static Path writePin(List<Psm> psms, Path pin) throws IOException {
        StringBuilder sb = new StringBuilder("SpecId\tLabel\tScanNr\tExpMass\tretentiontime\trank\thyperscore\t"
                + "log10_evalue\tPeptide\tProteins\n");
        int scan = 0;
        for (Psm psm : psms) {
            scan++;
            sb.append("golden.").append(scan).append('.').append(scan).append('.').append(psm.charge).append("_1\t1\t")
                    .append(scan).append("\t1000.0\t").append(10.0 + scan).append("\t1\t30.0\t-8.0\t")
                    .append(psm.peptide).append("\tsp|GOLDEN|GOLDEN_SYNTH\n");
        }
        Files.createDirectories(pin.getParent());
        Files.writeString(pin, sb.toString(), StandardCharsets.UTF_8);
        return pin;
    }

    /** The settings a run changes, and what they were before. */
    public static final class Settings {
        private final String instrument = Constants.instrument;
        private final String fragmentationType = Constants.FragmentationType;
        private final String koinaUrl = Constants.KoinaURL;
        private final Boolean foundBest = Constants.foundBest;
        private final String auxSpectraModel = Constants.auxSpectraModel;
        private final String features = Constants.features;
        private final Map<String, String> mzmlNCEs = new HashMap<>(NceConstants.mzmlNCEs);
        private final Map<String, String> calibratedModels = new HashMap<>(NceConstants.calibratedModels);
        private final Set<String> foundUnimods = new HashSet<>(PeptideFormatter.foundUnimods);

        /**
         * Sets the case's instrument, fragmentation and NCE (as the model's calibrated NCE too), the best model
         * as not found yet, no aux model and no features, clears PeptideFormatter.foundUnimods, and forgets
         * what Koina answered about models.
         *
         * @return what to {@link #restore} afterwards
         */
        public static Settings apply(Case c) {
            Settings saved = new Settings();
            Constants.instrument = c.instrument;
            Constants.FragmentationType = c.fragmentation;
            Constants.foundBest = false; //a failed batch marks the reader failed instead of exiting
            Constants.auxSpectraModel = "";
            Constants.features = ""; //createFull counts peptides only for the peptideCounts feature
            NceConstants.mzmlNCEs.clear();
            NceConstants.mzmlNCEs.put(c.fragmentation, c.nce);
            NceConstants.calibratedModels.clear();
            NceConstants.calibratedModels.put(c.model, c.nce);
            PeptideFormatter.foundUnimods.clear();
            KoinaModelCheck.forgetAnswers(); //the stub server is new, and may have an earlier one's port
            return saved;
        }

        public void restore() {
            Constants.instrument = instrument;
            Constants.FragmentationType = fragmentationType;
            Constants.KoinaURL = koinaUrl;
            Constants.foundBest = foundBest;
            Constants.auxSpectraModel = auxSpectraModel;
            Constants.features = features;
            NceConstants.mzmlNCEs.clear();
            NceConstants.mzmlNCEs.putAll(mzmlNCEs);
            NceConstants.calibratedModels.clear();
            NceConstants.calibratedModels.putAll(calibratedModels);
            PeptideFormatter.foundUnimods.clear();
            PeptideFormatter.foundUnimods.addAll(foundUnimods);
            KoinaModelCheck.forgetAnswers();
        }
    }

    /** Sends a request to the live Koina server, and keeps its metadata and config. */
    public static final class LiveKoina implements Responder {
        private final String url;
        private final HttpClient client = HttpClient.newBuilder().connectTimeout(Duration.ofSeconds(30)).build();
        public final Map<String, String> responses = new HashMap<>();
        public final Map<String, String> metadata = new HashMap<>();
        public final Map<String, String> configs = new HashMap<>();

        public LiveKoina(String url) {
            this.url = url.endsWith("/") ? url : url + "/";
        }

        public static LiveKoina fromProperty() {
            return new LiveKoina(System.getProperty(URL_PROPERTY, DEFAULT_KOINA_URL));
        }

        /** Gets the model's metadata and config from the server and keeps them. */
        public void fetchModel(String model) throws IOException, InterruptedException {
            metadata.put(model, get(model));
            configs.put(model, get(model + "/config"));
        }

        @Override
        public String respond(String model, String request) throws Exception {
            HttpRequest post = HttpRequest.newBuilder(URI.create(url + model + "/infer"))
                    .timeout(Duration.ofMinutes(5))
                    .header("Content-Type", "application/json; utf-8")
                    .POST(HttpRequest.BodyPublishers.ofString(request, StandardCharsets.UTF_8))
                    .build();
            HttpResponse<String> response = client.send(post,
                    HttpResponse.BodyHandlers.ofString(StandardCharsets.UTF_8));
            if (response.statusCode() != 200) {
                throw new AssertionError(model + ": Koina answered HTTP " + response.statusCode() + ": "
                        + response.body());
            }
            responses.put(model, response.body());
            return response.body();
        }

        private String get(String path) throws IOException, InterruptedException {
            HttpRequest get = HttpRequest.newBuilder(URI.create(url + path))
                    .timeout(Duration.ofMinutes(1))
                    .GET()
                    .build();
            HttpResponse<String> response = client.send(get,
                    HttpResponse.BodyHandlers.ofString(StandardCharsets.UTF_8));
            if (response.statusCode() != 200) {
                throw new AssertionError("GET " + url + path + ": HTTP " + response.statusCode() + ": "
                        + response.body());
            }
            return response.body();
        }
    }

    /**
     * Predictions as text: one P line per precursor in key order, then one F line per fragment in
     * the order the entry holds them. Floats are written with Float.toString, so they read back
     * exactly.
     */
    public static String dump(Map<String, PredictionEntry> preds) {
        StringBuilder sb = new StringBuilder(DUMP_HEADER);
        for (Map.Entry<String, PredictionEntry> e : new TreeMap<>(preds).entrySet()) {
            dumpEntry(sb, "P", "F", e.getKey(), e.getValue());
            if (e.getValue().auxSpectra != null) {
                dumpEntry(sb, "A", "AF", e.getKey(), e.getValue().auxSpectra);
            }
        }
        return sb.toString();
    }

    private static void dumpEntry(StringBuilder sb, String precursorTag, String fragmentTag, String key,
                                  PredictionEntry pe) {
        int n = pe.numFragments();
        sb.append(precursorTag).append('\t').append(key).append('\t').append(Float.toString(pe.getRT()))
                .append('\t').append(Float.toString(pe.getIM())).append('\t').append(n)
                .append('\t').append(pe.daltonMatching).append('\n');
        for (int i = 0; i < n; i++) {
            sb.append(fragmentTag).append('\t').append(Float.toString(pe.getMz(i)))
                    .append('\t').append(Float.toString(pe.getIntensity(i)))
                    .append('\t').append(pe.getIonTypeString(i))
                    .append('\t').append(pe.getFragNum(i))
                    .append('\t').append(pe.getCharge(i))
                    .append('\t').append(i < pe.fullAnnotations.length ? pe.fullAnnotations[i] : "")
                    .append('\t').append(i < pe.isotopes.length ? String.valueOf(pe.isotopes[i]) : "")
                    .append('\n');
        }
    }

    /**
     * Compares two dumps line by line. Text must match exactly; a number may differ from the golden
     * by at most two ulps, which allows the same arithmetic in another order and nothing more.
     *
     * @return null when they match, else what differs first
     */
    public static String compareDumps(String expected, String actual) {
        String[] want = lines(expected);
        String[] got = lines(actual);
        for (int i = 0; i < Math.min(want.length, got.length); i++) {
            if (!sameLine(want[i], got[i])) {
                return "line " + (i + 1) + "\n  expected: " + want[i] + "\n  actual:   " + got[i];
            }
        }
        if (want.length != got.length) {
            int i = Math.min(want.length, got.length);
            return "expected " + want.length + " lines, got " + got.length + "; first extra line " + (i + 1) + ": "
                    + (want.length > got.length ? "expected " + want[i] : "actual " + got[i]);
        }
        return null;
    }

    private static String[] lines(String text) {
        String normalized = text.replace("\r\n", "\n");
        if (normalized.endsWith("\n")) {
            normalized = normalized.substring(0, normalized.length() - 1);
        }
        return normalized.split("\n", -1);
    }

    private static boolean sameLine(String expected, String actual) {
        if (expected.equals(actual)) {
            return true;
        }
        String[] want = expected.split("\t", -1);
        String[] got = actual.split("\t", -1);
        if (want.length != got.length) {
            return false;
        }
        for (int i = 0; i < want.length; i++) {
            if (!want[i].equals(got[i]) && !closeFloats(want[i], got[i])) {
                return false;
            }
        }
        return true;
    }

    private static boolean closeFloats(String expected, String actual) {
        float a;
        float b;
        try {
            a = Float.parseFloat(expected);
            b = Float.parseFloat(actual);
        } catch (NumberFormatException e) {
            return false;
        }
        return Math.abs(a - b) <= 2 * Math.ulp(Math.max(Math.abs(a), Math.abs(b)));
    }

    /**
     * Compares two requests as JSON: the same inputs in the same order with the same values, while
     * whitespace and the order of an object's keys may differ.
     *
     * @return null when they match, else what differs first
     */
    public static String compareRequests(String expected, String actual) {
        JsonElement want = JsonParser.parseString(expected);
        JsonElement got = JsonParser.parseString(actual);
        if (want.equals(got)) {
            return null;
        }
        if (want.isJsonObject() && got.isJsonObject()) {
            String difference = inputDifference(want.getAsJsonObject(), got.getAsJsonObject());
            if (difference != null) {
                return difference;
            }
        }
        return "expected\n" + want + "\nactual\n" + got;
    }

    //names the input, field and row where two requests first differ, or null if they differ elsewhere
    private static String inputDifference(JsonObject want, JsonObject got) {
        JsonArray wantInputs = want.getAsJsonArray("inputs");
        JsonArray gotInputs = got.getAsJsonArray("inputs");
        if (wantInputs == null || gotInputs == null) {
            return null;
        }
        List<String> wantNames = inputNames(wantInputs);
        List<String> gotNames = inputNames(gotInputs);
        if (!wantNames.equals(gotNames)) {
            return "inputs " + gotNames + ", expected " + wantNames;
        }
        for (int i = 0; i < wantInputs.size(); i++) {
            JsonObject a = wantInputs.get(i).getAsJsonObject();
            JsonObject b = gotInputs.get(i).getAsJsonObject();
            Set<String> keys = new TreeSet<>(a.keySet());
            keys.addAll(b.keySet());
            for (String key : keys) {
                JsonElement x = a.get(key);
                JsonElement y = b.get(key);
                if (x != null && x.equals(y)) {
                    continue;
                }
                String where = "input " + wantNames.get(i) + ", " + key;
                if (x != null && y != null && x.isJsonArray() && y.isJsonArray()) {
                    JsonArray xs = x.getAsJsonArray();
                    JsonArray ys = y.getAsJsonArray();
                    for (int row = 0; row < Math.min(xs.size(), ys.size()); row++) {
                        if (!xs.get(row).equals(ys.get(row))) {
                            return where + " row " + row + ": expected " + xs.get(row) + ", actual " + ys.get(row);
                        }
                    }
                    return where + ": expected " + xs.size() + " rows, actual " + ys.size();
                }
                return where + ": expected " + x + ", actual " + y;
            }
        }
        return null;
    }

    private static List<String> inputNames(JsonArray inputs) {
        List<String> names = new ArrayList<>();
        for (JsonElement input : inputs) {
            JsonElement name = input.isJsonObject() ? input.getAsJsonObject().get("name") : null;
            names.add(name == null ? "?" : name.getAsString());
        }
        return names;
    }

    /** Writes what a run produced next to the build output, for a look after a failure. */
    public static Path writeActual(String model, Outcome outcome) throws IOException {
        Path dir = Paths.get(System.getProperty("basedir", System.getProperty("user.dir")))
                .resolve("target").resolve("koina-golden-actual").resolve(model);
        Files.createDirectories(dir);
        Files.writeString(dir.resolve(REQUEST), outcome.request, StandardCharsets.UTF_8);
        Files.writeString(dir.resolve(PARSED), outcome.parsed, StandardCharsets.UTF_8);
        Files.writeString(dir.resolve(TRANSFERRED), outcome.transferred, StandardCharsets.UTF_8);
        return dir;
    }

    /** The recorded files a model folder is missing. */
    public static List<String> missingFiles(String model) {
        Path dir = fixtureRoot().resolve(model);
        List<String> missing = new ArrayList<>();
        for (String file : Stream.concat(Stream.of(CASE, PSMS), RECORDED_FILES.stream()).collect(Collectors.toList())) {
            if (!Files.isRegularFile(dir.resolve(file))) {
                missing.add(file);
            }
        }
        return missing;
    }
}
