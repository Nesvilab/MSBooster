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

package transferlearn;

import static org.junit.jupiter.api.Assertions.assertArrayEquals;
import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertFalse;
import static org.junit.jupiter.api.Assertions.assertTrue;

import java.io.File;
import java.io.FileInputStream;
import java.io.IOException;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.sql.Connection;
import java.sql.DriverManager;
import java.sql.ResultSet;
import java.sql.Statement;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.stream.Collectors;
import java.util.stream.Stream;

import org.junit.jupiter.api.Test;
import org.junit.jupiter.api.io.TempDir;
import testutils.ParquetFooters;

// The server predicts the peptides it is given and knows nothing of the enzyme, so the library it
// returns records none. FragDIA reads the enzyme from the Parquet footer and nowhere else - so
// whichever way the server path leaves its chunks, every Parquet file it leaves has to carry it.
public class PredictorTest {

    /** A chunk as the server returns it: a few library rows, and a footer of its own. */
    private static File chunk(Path dir, String name, int rows) throws Exception {
        File parquet = dir.resolve(name).toFile();
        try (Connection conn = DriverManager.getConnection("jdbc:duckdb:");
             Statement stmt = conn.createStatement()) {
            stmt.execute("COPY (SELECT CAST(400 + i AS FLOAT) AS PrecursorMz, 'PEPTIDEK' AS PeptideSequence " +
                    "FROM range(" + rows + ") t(i)) TO '" + path(parquet) + "' (FORMAT PARQUET, " +
                    "KV_METADATA {'written_by': 'the server'})");
        }
        return parquet;
    }

    private static Map<String, String> trypsin(String fasta) {
        EnzymeOptions enzyme = new EnzymeOptions();
        String[] args = {"--cut", "KR", "--no-cut", "P", "--sense", "C", "--num-enzyme-termini", "1",
                "--clip-n-term-m", "0"};
        for (int i = 0; i < args.length; i += 2) {
            assertTrue(enzyme.accept(args[i], args[i + 1]));
        }
        enzyme.validate();
        return enzyme.footer(fasta);
    }

    private static long rows(File parquet) throws Exception {
        try (Connection conn = DriverManager.getConnection("jdbc:duckdb:");
             Statement stmt = conn.createStatement();
             ResultSet rs = stmt.executeQuery("SELECT COUNT(*) FROM read_parquet('" + path(parquet) + "')")) {
            rs.next();
            return rs.getLong(1);
        }
    }

    private static Set<String> names(Path dir) throws IOException {
        try (Stream<Path> files = Files.list(dir)) {
            return files.map(p -> p.getFileName().toString()).collect(Collectors.toSet());
        }
    }

    private static String path(File file) {
        return file.getAbsolutePath().replace("\\", "/");
    }

    // One chunk is the library, under the run's name, and it records the enzyme and the FASTA - and
    // only them: what the server's writer kept there is replaced, not added to.
    @Test
    public void oneChunkBecomesTheLibraryWithTheFooter(@TempDir Path dir) throws Exception {
        File only = chunk(dir, "model0.parquet", 3);
        File total = dir.resolve("model.parquet").toFile();
        Map<String, String> footer = trypsin("C:\\fasta\\human.fasta");

        Predictor.finishParquet(new String[] {only.getAbsolutePath()}, total, false, footer);

        assertEquals(footer, ParquetFooters.read(total));
        assertEquals(3, rows(total));
        assertEquals(new HashSet<>(Collections.singletonList("model.parquet")), names(dir),
                "nothing staged is left beside the library");
    }

    // With nothing to record the chunk is only moved, exactly as it was downloaded: no pass over the
    // library is spent writing a footer that says nothing.
    @Test
    public void oneChunkWithNothingToRecordIsOnlyMoved(@TempDir Path dir) throws Exception {
        File only = chunk(dir, "model0.parquet", 3);
        byte[] downloaded = Files.readAllBytes(only.toPath());
        File total = dir.resolve("model.parquet").toFile();

        Predictor.finishParquet(new String[] {only.getAbsolutePath()}, total, false, Collections.emptyMap());

        assertArrayEquals(downloaded, Files.readAllBytes(total.toPath()));
        assertEquals(new HashSet<>(Collections.singletonList("model.parquet")), names(dir));
    }

    // Merged, the one file holds every chunk's rows and the footer, and the chunks are gone.
    @Test
    public void mergedChunksCarryTheFooter(@TempDir Path dir) throws Exception {
        String[] chunks = {chunk(dir, "model0.parquet", 2).getAbsolutePath(),
                chunk(dir, "model1.parquet", 5).getAbsolutePath()};
        File total = dir.resolve("model.parquet").toFile();
        Map<String, String> footer = trypsin("");

        Predictor.finishParquet(chunks, total, true, footer);

        assertEquals(footer, ParquetFooters.read(total));
        assertFalse(footer.containsKey("database_name"), "no --fasta, no database_name");
        assertEquals(7, rows(total));
        assertEquals(new HashSet<>(Collections.singletonList("model.parquet")), names(dir));
    }

    // Not merged, the chunks are the library's shards, and FragDIA compares the enzyme across a
    // directory's shards - so each one carries the same footer, and keeps every row.
    @Test
    public void unmergedChunksEachCarryTheFooter(@TempDir Path dir) throws Exception {
        List<File> chunks = Arrays.asList(chunk(dir, "model0.parquet", 2), chunk(dir, "model1.parquet", 5));
        File total = dir.resolve("model.parquet").toFile();
        Map<String, String> footer = trypsin("/data/human.fasta");

        Predictor.finishParquet(chunks.stream().map(File::getAbsolutePath).toArray(String[]::new), total,
                false, footer);

        assertEquals(footer, ParquetFooters.read(chunks.get(0)));
        assertEquals(footer, ParquetFooters.read(chunks.get(1)));
        assertEquals(2, rows(chunks.get(0)));
        assertEquals(5, rows(chunks.get(1)));
        assertEquals(new HashSet<>(Arrays.asList("model0.parquet", "model1.parquet")), names(dir));
    }

    // A rewrite that fails leaves the library as it was and nothing staged beside it.
    @Test
    public void aFooterThatCannotBeWrittenLeavesTheLibraryAsItWas(@TempDir Path dir) throws Exception {
        File notParquet = dir.resolve("model.parquet").toFile();
        Files.write(notParquet.toPath(), "not a parquet file".getBytes());

        org.junit.jupiter.api.Assertions.assertThrows(java.sql.SQLException.class,
                () -> Helpers.writeFooter(notParquet, trypsin("")));

        assertEquals("not a parquet file", new String(Files.readAllBytes(notParquet.toPath())));
        assertEquals(new HashSet<>(Collections.singletonList("model.parquet")), names(dir));
    }

    // The rewritten library is moved over the original in one step, never by deleting the original
    // first, so a move that fails leaves the library where it was. Windows refuses to move a file
    // another handle holds open, which is how the failure is brought about here; elsewhere that move
    // succeeds, and the library is then the replacement - one whole file either way, never neither.
    @Test
    public void aReplacementThatCannotBeMovedLeavesTheLibraryAsItWas(@TempDir Path dir) throws Exception {
        File library = dir.resolve("model.parquet").toFile();
        File staged = dir.resolve("model.parquet.footer.parquet").toFile();
        Files.write(library.toPath(), "the library".getBytes(StandardCharsets.UTF_8));
        Files.write(staged.toPath(), "its replacement".getBytes(StandardCharsets.UTF_8));

        boolean moved;
        try (FileInputStream held = new FileInputStream(staged)) {
            try {
                Helpers.replaceFile(staged, library);
                moved = true;
            } catch (IOException e) {
                moved = false;
            }
        }

        if (System.getProperty("os.name").startsWith("Windows")) {
            assertFalse(moved, "Windows moved a file another handle holds open");
            assertEquals("the library", text(library));
            assertEquals("its replacement", text(staged));
        } else {
            assertTrue(moved);
            assertEquals("its replacement", text(library));
            assertFalse(staged.exists());
        }
    }

    private static String text(File file) throws IOException {
        return new String(Files.readAllBytes(file.toPath()), StandardCharsets.UTF_8);
    }
}
