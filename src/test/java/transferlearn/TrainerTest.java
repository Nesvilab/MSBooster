package transferlearn;

import static org.junit.jupiter.api.Assertions.assertArrayEquals;
import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertTrue;

import com.sun.net.httpserver.HttpExchange;
import com.sun.net.httpserver.HttpServer;
import java.io.ByteArrayOutputStream;
import java.io.IOException;
import java.io.InputStream;
import java.io.OutputStream;
import java.net.InetSocketAddress;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.Arrays;
import java.util.HashSet;
import java.util.Set;
import java.util.stream.Collectors;
import java.util.stream.Stream;
import org.junit.jupiter.api.AfterEach;
import org.junit.jupiter.api.BeforeEach;
import org.junit.jupiter.api.Test;
import org.junit.jupiter.api.io.TempDir;

/**
 * Training on the server, against a stand-in for it on this computer: the library goes up as Parquet,
 * converted in a folder of its own, so a Parquet library beside the TSV (FragSpecLib writes
 * library.parquet beside library.tsv, and FragDIA reads its footer) is left as it was.
 */
public class TrainerTest {

    private static final String LIBRARY_TSV = String.join("\n",
            "PrecursorMz\tProductMz\tModifiedPeptideSequence\tPrecursorCharge\tLibraryIntensity\tNormalizedRetentionTime",
            "500.25\t300.1\tPEPTIDEK\t2\t1.0\t10.5",
            "500.25\t401.2\tPEPTIDEK\t2\t0.5\t10.5",
            "");
    private static final byte[] LIBRARY_PARQUET = "FragSpecLib's library.parquet, footer and all"
            .getBytes(StandardCharsets.UTF_8);
    private static final byte[] MODEL = "the trained model".getBytes(StandardCharsets.UTF_8);

    private HttpServer server;
    private byte[] uploaded;
    private Long waitTime;

    @BeforeEach
    public void startTheServer() throws IOException {
        server = HttpServer.create(new InetSocketAddress("127.0.0.1", 0), 0);
        server.createContext("/train/upload", exchange -> {
            uploaded = readAll(exchange.getRequestBody());
            respond(exchange, "{\"job_id\": \"7\"}".getBytes(StandardCharsets.UTF_8));
        });
        server.createContext("/train/status/7", exchange -> respond(exchange,
                "{\"status\": \"SUCCESS\", \"stdout\": \"trained\\n\", \"queue_position\": 0}"
                        .getBytes(StandardCharsets.UTF_8)));
        server.createContext("/train/download/7", exchange -> respond(exchange, MODEL));
        server.start();
        waitTime = Trainer.waitTime;
        Trainer.waitTime = 1L;
    }

    @AfterEach
    public void stopTheServer() {
        server.stop(0);
        Trainer.waitTime = waitTime;
    }

    @Test
    public void theParquetLibraryBesideTheTsvIsLeftAsItWas(@TempDir Path dir) throws Exception {
        final Path tsv = library(dir);

        final String model = Trainer.train(args(tsv, "--basename", "model"));

        assertArrayEquals(LIBRARY_PARQUET, Files.readAllBytes(dir.resolve("library.parquet")));
        assertEquals(dir.resolve("model.zip").toString(), model, "the model is downloaded beside the library");
        assertArrayEquals(MODEL, Files.readAllBytes(dir.resolve("model.zip")));
        assertEquals(names("library.tsv", "library.parquet", "model.zip"), names(dir),
                "nothing made for the upload is left");
    }

    @Test
    public void theLibraryGoesUpAsParquetUnderItsOwnName(@TempDir Path dir) throws Exception {
        Trainer.train(args(library(dir), "--basename", "model"));

        final String body = new String(uploaded, StandardCharsets.ISO_8859_1);
        assertTrue(body.contains("filename=\"library.parquet\""), body);
        assertTrue(body.contains("PAR1"), "the upload is a Parquet file");
    }

    @Test
    public void withAnOutputDirTheLibraryFolderIsNotTouched(@TempDir Path dir) throws Exception {
        final Path tsv = library(dir.resolve("lib"));
        final Path out = dir.resolve("out");

        final String model = Trainer.train(args(tsv, "--basename", "model", "--output-dir", out.toString()));

        assertEquals(out.resolve("model.zip").toString(), model);
        assertEquals(names("library.tsv", "library.parquet"), names(tsv.getParent()));
        assertArrayEquals(LIBRARY_PARQUET, Files.readAllBytes(tsv.resolveSibling("library.parquet")));
        assertEquals(names("model.zip"), names(out), "nothing made for the upload is left");
    }

    @Test
    public void theUploadIsNamedAfterTheLibrary() {
        assertEquals("library.parquet", Trainer.parquetName("library.tsv"));
        assertEquals("lib.v2.parquet", Trainer.parquetName("lib.v2.tsv"));
        assertEquals("library.parquet", Trainer.parquetName("library"));
    }

    /** library.tsv in {@code dir}, with FragSpecLib's library.parquet beside it. */
    private static Path library(Path dir) throws IOException {
        Files.createDirectories(dir);
        final Path tsv = dir.resolve("library.tsv");
        Files.write(tsv, LIBRARY_TSV.getBytes(StandardCharsets.UTF_8));
        Files.write(dir.resolve("library.parquet"), LIBRARY_PARQUET);
        return tsv;
    }

    private String[] args(Path library, String... more) {
        final String[] base = {
                "--url", "http://127.0.0.1:" + server.getAddress().getPort(),
                "--library", library.toString(),
                "--api-key", "key"};
        final String[] all = Arrays.copyOf(base, base.length + more.length);
        System.arraycopy(more, 0, all, base.length, more.length);
        return all;
    }

    private static Set<String> names(String... names) {
        return new HashSet<>(Arrays.asList(names));
    }

    private static Set<String> names(Path dir) throws IOException {
        try (Stream<Path> files = Files.list(dir)) {
            return files.map(p -> p.getFileName().toString()).collect(Collectors.toSet());
        }
    }

    private static void respond(HttpExchange exchange, byte[] body) throws IOException {
        exchange.sendResponseHeaders(200, body.length);
        try (OutputStream os = exchange.getResponseBody()) {
            os.write(body);
        }
    }

    private static byte[] readAll(InputStream in) throws IOException {
        final ByteArrayOutputStream bytes = new ByteArrayOutputStream();
        final byte[] buffer = new byte[8192];
        int n;
        while ((n = in.read(buffer)) != -1) {
            bytes.write(buffer, 0, n);
        }
        return bytes.toByteArray();
    }
}
