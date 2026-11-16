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
import java.util.ArrayList;
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

// The server predicts the peptides it is given and knows nothing of the --meta file, so the library
// it returns records none of it. FragDIA reads the enzyme from the Parquet footer and nowhere else,
// and FragPipe requires a library's shards to carry the same footer - so whichever way the server
// path leaves its chunks, every Parquet file it leaves has to carry the file's pairs, in its order.
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

    /**
     * The footer of a --meta file FragPipe writes for trypsin with no enzyme name, as the predictor
     * reads it: {@code database_name} only when {@code fasta} is not empty, as FragPipe writes none
     * without a FASTA.
     */
    private static Map<String, String> trypsin(Path dir, String fasta) throws IOException {
        List<String> lines = new ArrayList<>(Arrays.asList(
                "# Library footer, written by FragPipe",
                "search_enzyme_cut_1 = KR",
                "search_enzyme_nocut_1 = P",
                "search_enzyme_sense_1 = C",
                "num_enzyme_termini = 1",
                "clip_nTerm_M = 0"));
        if (!fasta.isEmpty()) {
            lines.add("database_name = " + fasta);
        }
        //read and then removed, so the directory holds only what finishParquet leaves there
        Path file = dir.resolve("library.meta");
        Files.write(file, lines, StandardCharsets.UTF_8);
        MetaFile meta = new MetaFile();
        meta.accept(file.toString());
        meta.read();
        Files.delete(file);
        return meta.footer();
    }

    /** A footer as an ordered list, so that a comparison sees the order of its pairs as well. */
    private static List<Map.Entry<String, String>> inOrder(Map<String, String> footer) {
        return new ArrayList<>(footer.entrySet());
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

    // One chunk is the library, under the run's name, and it records the meta file's pairs - and
    // only them: what the server's writer kept there is replaced, not added to.
    @Test
    public void oneChunkBecomesTheLibraryWithTheFooter(@TempDir Path dir) throws Exception {
        File only = chunk(dir, "model0.parquet", 3);
        File total = dir.resolve("model.parquet").toFile();
        Map<String, String> footer = trypsin(dir, "C:\\fasta\\human.fasta");

        Predictor.finishParquet(new String[] {only.getAbsolutePath()}, total, false, footer);

        assertEquals(inOrder(footer), inOrder(ParquetFooters.read(total)));
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
        Map<String, String> footer = trypsin(dir, "");

        Predictor.finishParquet(chunks, total, true, footer);

        assertEquals(inOrder(footer), inOrder(ParquetFooters.read(total)));
        assertFalse(ParquetFooters.read(total).containsKey("database_name"), "the file names no database");
        assertEquals(7, rows(total));
        assertEquals(new HashSet<>(Collections.singletonList("model.parquet")), names(dir));
    }

    // Not merged, the chunks are the library's shards, and FragPipe requires every shard's whole
    // footer to be the same - so each one carries the same pairs in the same order, and keeps every
    // row.
    @Test
    public void unmergedChunksEachCarryTheFooter(@TempDir Path dir) throws Exception {
        List<File> chunks = Arrays.asList(chunk(dir, "model0.parquet", 2), chunk(dir, "model1.parquet", 5));
        File total = dir.resolve("model.parquet").toFile();
        Map<String, String> footer = trypsin(dir, "/data/human.fasta");

        Predictor.finishParquet(chunks.stream().map(File::getAbsolutePath).toArray(String[]::new), total,
                false, footer);

        assertEquals(inOrder(footer), inOrder(ParquetFooters.read(chunks.get(0))));
        assertEquals(inOrder(footer), inOrder(ParquetFooters.read(chunks.get(1))));
        assertEquals(2, rows(chunks.get(0)));
        assertEquals(5, rows(chunks.get(1)));
        assertEquals(new HashSet<>(Arrays.asList("model0.parquet", "model1.parquet")), names(dir));
    }

    // A rewrite that fails leaves the library as it was and nothing staged beside it.
    @Test
    public void aFooterThatCannotBeWrittenLeavesTheLibraryAsItWas(@TempDir Path dir) throws Exception {
        Map<String, String> footer = trypsin(dir, "");
        File notParquet = dir.resolve("model.parquet").toFile();
        Files.write(notParquet.toPath(), "not a parquet file".getBytes());

        org.junit.jupiter.api.Assertions.assertThrows(java.sql.SQLException.class,
                () -> Helpers.writeFooter(notParquet, footer));

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
