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

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertFalse;
import static org.junit.jupiter.api.Assertions.assertThrows;
import static org.junit.jupiter.api.Assertions.assertTrue;

import java.io.ByteArrayOutputStream;
import java.io.File;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.sql.Connection;
import java.sql.DriverManager;
import java.sql.SQLException;
import java.sql.Statement;
import java.util.AbstractMap.SimpleEntry;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;

import org.junit.jupiter.api.Test;
import org.junit.jupiter.api.io.TempDir;
import testutils.ParquetFooters;

// FragPipe writes the meta file and every library builder - FragCast, FragSpecLib and the two
// predictors here - reads it, so each has to read the same file into the same pairs, and refuse the
// same files. The pairs go into the library's Parquet footer verbatim and in file order: FragDIA
// digests with the enzyme it finds there, and FragPipe requires a library's shards to carry the same
// footer, so a builder that read a line differently would record something nobody wrote.
public class MetaFileTest {

    private static final String PATH = "test.meta";

    private static final char BOM = 0xFEFF;

    /** {@code text}'s pairs, read as a meta file named {@link #PATH}. */
    private static LinkedHashMap<String, String> parse(String text) {
        return MetaFile.parse(PATH, text.getBytes(StandardCharsets.UTF_8));
    }

    /** The pairs as an ordered list, so that a comparison sees their order as well as their values. */
    private static List<Map.Entry<String, String>> inOrder(Map<String, String> pairs) {
        return new ArrayList<>(pairs.entrySet());
    }

    @SafeVarargs
    private static List<Map.Entry<String, String>> pairs(Map.Entry<String, String>... pairs) {
        return Arrays.asList(pairs);
    }

    private static Map.Entry<String, String> pair(String key, String value) {
        return new SimpleEntry<>(key, value);
    }

    private static String rejection(String text) {
        return assertThrows(IllegalArgumentException.class, () -> parse(text),
                "accepted " + text).getMessage();
    }

    /** What FragPipe writes for strict trypsin, searched against a FASTA. */
    private static final String TRYPSIN = String.join("\n",
            "# Library footer, written by FragPipe",
            "search_enzyme_name_1 = stricttrypsin",
            "search_enzyme_cut_1 = KR",
            "search_enzyme_nocut_1 = @",
            "search_enzyme_sense_1 = C",
            "num_enzyme_termini = 2",
            "clip_nTerm_M = 1",
            "database_name = C:\\fasta\\human.fasta",
            "");

    // --- what a line means ----------------------------------------------------------------------

    // Every pair is kept, in the order of the file's lines - which is the order of the footer, and
    // FragPipe compares whole footers. A key sorted or hashed into another place would not be the
    // footer it wrote.
    @Test
    public void thePairsAreKeptInFileOrder() {
        assertEquals(pairs(pair("search_enzyme_name_1", "stricttrypsin"), pair("search_enzyme_cut_1", "KR"),
                        pair("search_enzyme_nocut_1", "@"), pair("search_enzyme_sense_1", "C"),
                        pair("num_enzyme_termini", "2"), pair("clip_nTerm_M", "1"),
                        pair("database_name", "C:\\fasta\\human.fasta")),
                inOrder(parse(TRYPSIN)));
        assertEquals(pairs(pair("z", "1"), pair("a", "2"), pair("m", "3")),
                inOrder(parse("z = 1\na = 2\nm = 3\n")), "not sorted");
    }

    // The tools know no key names: whatever the file says is recorded, whatever it is called.
    @Test
    public void anyKeyIsTaken() {
        assertEquals(pairs(pair("anything at all", "x"), pair("search_enzyme_cut_9", "Z")),
                inOrder(parse("anything at all = x\nsearch_enzyme_cut_9 = Z")));
    }

    // Key and value are trimmed, the line before them; what is between the words is kept.
    @Test
    public void keyAndValueAreTrimmed() {
        assertEquals(pairs(pair("a key", "a  value")), inOrder(parse(" \t a key \t=\t a  value \t ")));
        assertEquals(pairs(pair("a", "1")), inOrder(parse("a=1")), "no spaces are needed around '='");
    }

    // Only the first '=' splits a line, so a value may hold its own.
    @Test
    public void aValueMayContainEquals() {
        assertEquals(pairs(pair("a", "b=c = d=")), inOrder(parse("a = b=c = d=")));
        assertEquals(pairs(pair("a", "=")), inOrder(parse("a==")));
    }

    // A comment is a whole line. A FASTA path can hold a '#', and a value cut short at one would be
    // recorded as a path nobody wrote - so a '#' anywhere but first is part of the key or the value.
    @Test
    public void aHashInsideAKeyOrValueIsNotAComment() {
        assertEquals(pairs(pair("database_name", "D:\\data #2\\human#1.fasta # really")),
                inOrder(parse("database_name = D:\\data #2\\human#1.fasta # really")));
        assertEquals(pairs(pair("key#", "v")), inOrder(parse("key# = v")));
    }

    @Test
    public void commentAndBlankLinesAreSkipped() {
        assertEquals(pairs(pair("a", "1"), pair("b", "2")),
                inOrder(parse("# a comment\n\n   \n\t# an indented comment\na = 1\n" +
                        " \t\r\n#= not a pair\nb = 2\n")));
    }

    // An empty value is a value, written into the footer as an empty string; only an empty key is
    // refused.
    @Test
    public void anEmptyValueIsKept() {
        assertEquals(pairs(pair("a", ""), pair("b", "")), inOrder(parse("a =\nb = \t\r\n")));
    }

    // A file FragPipe wrote on Windows, or one opened and saved there, ends its lines in \r\n and may
    // begin with a byte order mark. Neither may reach a key or a value.
    @Test
    public void crlfLinesReadLikeLfLines() {
        assertEquals(inOrder(parse(TRYPSIN)), inOrder(parse(TRYPSIN.replace("\n", "\r\n"))));
        assertEquals(pairs(pair("a", "1")), inOrder(parse("a = 1\r")), "a last line's \\r");
    }

    @Test
    public void oneLeadingByteOrderMarkIsSkipped() {
        assertEquals(inOrder(parse(TRYPSIN)), inOrder(parse(BOM + TRYPSIN)));
        // a first line that is a comment is still one behind the mark
        assertEquals(pairs(pair("a", "1")), inOrder(parse(BOM + "# comment\na = 1")));
    }

    // One, and only at the start: a mark is not whitespace, so a second one is part of the key.
    @Test
    public void onlyOneByteOrderMarkIsSkipped() {
        assertEquals(pairs(pair(BOM + "a", "1")), inOrder(parse("" + BOM + BOM + "a = 1")));
        assertEquals(pairs(pair("a", "1"), pair(BOM + "b", "2")), inOrder(parse("a = 1\n" + BOM + "b = 2")));
    }

    // Spaces, tabs and carriage returns are trimmed, and nothing else: every builder trims exactly
    // this set, so a language's own idea of whitespace - Java's String.trim takes every control
    // character - cannot make two builders read one line two ways.
    @Test
    public void onlySpacesTabsAndCarriageReturnsAreTrimmed() {
        String vt = String.valueOf((char) 0x0B);
        String nbsp = String.valueOf((char) 0xA0);
        assertEquals(pairs(pair(vt + "a" + vt, "\f1\f"), pair(nbsp + "b", "2" + nbsp)),
                inOrder(parse(vt + "a" + vt + " = \f1\f\n" + nbsp + "b = 2" + nbsp)));
        assertEquals(pairs(pair("\f", "x")), inOrder(parse("\f = x")), "a form feed is a key, not nothing");
    }

    // --- a file with nothing in it ---------------------------------------------------------------

    // A file of no pairs is as good as no file: the library records no footer at all.
    @Test
    public void aFileOfNoPairsHasNoPairs() {
        assertTrue(parse("").isEmpty());
        assertTrue(parse("\n\n").isEmpty());
        assertTrue(parse("# only a comment\r\n   \r\n").isEmpty());
        assertTrue(parse(String.valueOf(BOM)).isEmpty());
    }

    // --- what is refused, and how it is reported ---------------------------------------------------

    @Test
    public void aLineWithoutEqualsIsRefusedByLine() {
        assertEquals(PATH + ":2: expected 'key = value'", rejection("a = 1\njust text\nb = 2"));
        assertEquals(PATH + ":4: expected 'key = value'",
                rejection("# comment\n\na = 1\nsearch_enzyme_cut_1 KR"), "comment and blank lines count");
    }

    @Test
    public void anEmptyKeyIsRefused() {
        assertEquals(PATH + ":1: empty key", rejection("= value"));
        assertEquals(PATH + ":2: empty key", rejection("a = 1\n \t= value"));
    }

    // A key given twice would leave a builder to choose which one to record. Keys are compared as
    // trimmed, so spacing does not make a second one different.
    @Test
    public void aDuplicateKeyIsRefusedNamingBothLines() {
        assertEquals(PATH + ":3: duplicate key 'a' (first on line 1)", rejection("a = 1\n# x\na = 2"));
        assertEquals(PATH + ":2: duplicate key 'a' (first on line 1)", rejection(" a=1\na =1"),
                "the same value is still a second key");
    }

    // Keys differing only in case are different keys: nothing is folded. (DuckDB cannot record two
    // that differ only in ASCII case, which checkWritable refuses separately; see below.)
    @Test
    public void keysAreCaseSensitive() {
        assertEquals(pairs(pair("clip_nTerm_M", "1"), pair("clip_nterm_m", "0")),
                inOrder(parse("clip_nTerm_M = 1\nclip_nterm_m = 0")));
    }

    @Test
    public void aFileThatIsNotUtf8IsRefused() {
        byte[][] notUtf8 = {
                {'a', ' ', '=', ' ', (byte) 0xFF},
                {'a', '=', (byte) 0xC3}, //a sequence cut short
                {'a', '=', (byte) 0xC0, (byte) 0xAF}, //an overlong '/'
                {'a', '=', (byte) 0xED, (byte) 0xA0, (byte) 0x80}, //a surrogate, encoded
                {'a', '=', 'x', '\n', 'b', '=', (byte) 0xE9}}; //Latin-1, not UTF-8
        for (byte[] bytes : notUtf8) {
            IllegalArgumentException e = assertThrows(IllegalArgumentException.class,
                    () -> MetaFile.parse(PATH, bytes), Arrays.toString(bytes));
            assertEquals("cannot read meta file " + PATH + ": it is not valid UTF-8", e.getMessage());
        }
    }

    // Any Unicode is read as written.
    @Test
    public void utf8IsReadAsWritten() {
        String value = "C:\\Users\\o'brien\\prot" + (char) 0xE9 + "ome " +
                new String(Character.toChars(0x1F9EA)) + ".fasta";
        assertEquals(pairs(pair("database_name", value)), inOrder(parse("database_name = " + value)));
    }

    @Test
    public void aMissingFileIsRefused(@TempDir Path dir) {
        String missing = dir.resolve("absent.meta").toString();
        IllegalArgumentException e = assertThrows(IllegalArgumentException.class,
                () -> MetaFile.parse(missing));
        assertEquals("cannot read meta file " + missing + ": no such file", e.getMessage());
    }

    @Test
    public void aDirectoryIsRefused(@TempDir Path dir) {
        String directory = dir.toString();
        IllegalArgumentException e = assertThrows(IllegalArgumentException.class,
                () -> MetaFile.parse(directory));
        assertEquals("cannot read meta file " + directory + ": it is a directory", e.getMessage());
    }

    // The bytes on disk are what is read: a file is parsed exactly as its contents are.
    @Test
    public void aFileIsReadAsItsBytes(@TempDir Path dir) throws Exception {
        Path file = dir.resolve("fragcast.meta");
        ByteArrayOutputStream bytes = new ByteArrayOutputStream();
        bytes.write(new byte[] {(byte) 0xEF, (byte) 0xBB, (byte) 0xBF}); //the mark as it sits in a file
        bytes.write(TRYPSIN.replace("\n", "\r\n").getBytes(StandardCharsets.UTF_8));
        Files.write(file, bytes.toByteArray());

        assertEquals(inOrder(parse(TRYPSIN)), inOrder(MetaFile.parse(file.toString())));
    }

    // --- the --meta option ---------------------------------------------------------------------

    // No --meta: nothing is read, nothing is forwarded, and the footer records nothing.
    @Test
    public void withoutMetaThereIsNoFooterAndNothingToForward() {
        MetaFile meta = new MetaFile();
        meta.read();
        assertFalse(meta.isSet());
        assertTrue(meta.footer().isEmpty());
        assertTrue(meta.builderArgs().isEmpty());
    }

    // The path is forwarded exactly as given, and the footer is the file's pairs in its order.
    @Test
    public void metaIsForwardedAsGivenAndReadIntoTheFooter(@TempDir Path dir) throws Exception {
        Path file = dir.resolve("fragcast.meta");
        Files.write(file, TRYPSIN.getBytes(StandardCharsets.UTF_8));
        MetaFile meta = new MetaFile();
        meta.accept(file.toString());
        meta.read();

        assertTrue(meta.isSet());
        assertEquals(Arrays.asList("--meta", file.toString()), meta.builderArgs());
        assertEquals(inOrder(parse(TRYPSIN)), inOrder(meta.footer()));
    }

    // A file of no pairs is forwarded all the same - the builder reads it to the same nothing - and
    // the footer this side writes records nothing either.
    @Test
    public void aMetaFileOfNoPairsRecordsNothing(@TempDir Path dir) throws Exception {
        Path file = dir.resolve("empty.meta");
        Files.write(file, "# nothing to record\n".getBytes(StandardCharsets.UTF_8));
        MetaFile meta = new MetaFile();
        meta.accept(file.toString());
        meta.read();

        assertTrue(meta.footer().isEmpty());
        assertEquals(Arrays.asList("--meta", file.toString()), meta.builderArgs());
    }

    @Test
    public void metaGivenTwiceIsAnError() {
        MetaFile meta = new MetaFile();
        meta.accept("first.meta");
        IllegalArgumentException e = assertThrows(IllegalArgumentException.class,
                () -> meta.accept("second.meta"));
        assertTrue(e.getMessage().startsWith("--meta was given more than once"), e.getMessage());
        assertEquals(Arrays.asList("--meta", "first.meta"), meta.builderArgs(),
                "the first one is not replaced");
    }

    // An empty path is a file that cannot be read, not the absence of one.
    @Test
    public void anEmptyMetaPathIsAnError() {
        MetaFile meta = new MetaFile();
        meta.accept("");
        IllegalArgumentException e = assertThrows(IllegalArgumentException.class, meta::read);
        assertEquals("cannot read meta file : no path was given", e.getMessage());
    }

    @Test
    public void aMetaFileThatCannotBeReadIsAnErrorWhenRead(@TempDir Path dir) {
        MetaFile meta = new MetaFile();
        meta.accept(dir.resolve("absent.meta").toString());
        assertThrows(IllegalArgumentException.class, meta::read);
    }

    // --- the footer as DuckDB writes it ------------------------------------------------------------

    // The proof that DuckDB, at the version this build pins, writes KV_METADATA into a Parquet footer
    // in the order given, and that a reader finds exactly those pairs there - nothing more. The keys
    // and values are the awkward ones a file can hold: a quote in a key and in a value, a Windows
    // backslash, a space, '#' and '=', a character outside ASCII, and an empty value.
    @Test
    public void duckDbWritesTheFooterInOrderAndAReaderFindsExactlyIt(@TempDir Path dir) throws Exception {
        String fasta = "C:\\Users\\o'brien\\FASTA files #2\\prot" + (char) 0xE9 + "ome=1.fasta";
        Map<String, String> footer = parse(TRYPSIN.replace("C:\\fasta\\human.fasta", fasta) +
                "it's = a key with a quote\nempty =\n'); DROP TABLE x; -- = {'a': 'b'}\n");

        File parquet = dir.resolve("library.parquet").toFile();
        try (Connection conn = DriverManager.getConnection("jdbc:duckdb:");
             Statement stmt = conn.createStatement()) {
            stmt.execute("COPY (SELECT CAST(1 AS FLOAT) AS PrecursorMz) TO '" +
                    parquet.getAbsolutePath().replace("\\", "/") + "' (FORMAT PARQUET" +
                    MetaFile.kvMetadataOption(footer) + ")");
        }

        Map<String, String> read = ParquetFooters.read(parquet);
        assertEquals(inOrder(footer), inOrder(read));
        assertEquals(fasta, read.get("database_name"), "the path exactly as given");
        assertEquals("", read.get("empty"));
        assertEquals("{'a': 'b'}", read.get("'); DROP TABLE x; --"));
    }

    // With nothing to record the COPY must be the one it always was, not one with an empty clause.
    @Test
    public void anEmptyFooterAddsNothingToTheCopy() {
        assertEquals("", MetaFile.kvMetadataOption(Collections.emptyMap()));
    }

    // DuckDB's parser stops at a NUL, so one inside a literal would cut the statement short. A value's
    // NULs are spliced in instead, and the value is recorded whole.
    @Test
    public void aValueHoldingNulIsRecordedWhole(@TempDir Path dir) throws Exception {
        String nul = String.valueOf((char) 0);
        Map<String, String> footer = parse("a = x" + nul + "y\nb = " + nul + "\nc = " + nul + nul +
                "'it''s'" + nul + "\nd = plain");

        Map<String, String> read = ParquetFooters.read(copyWithFooter(dir, footer));
        assertEquals(inOrder(footer), inOrder(read));
        assertEquals("x" + nul + "y", read.get("a"));
    }

    // --- what DuckDB cannot record -----------------------------------------------------------------

    /** {@code text}, written to a file in {@code dir} and read as {@code --meta}. */
    private static MetaFile metaFile(Path dir, String text) throws Exception {
        Path file = dir.resolve(PATH);
        Files.write(file, text.getBytes(StandardCharsets.UTF_8));
        MetaFile meta = new MetaFile();
        meta.accept(file.toString());
        meta.read();
        return meta;
    }

    /** A one-row Parquet file in {@code dir}, written by DuckDB with {@code footer}. */
    private static File copyWithFooter(Path dir, Map<String, String> footer) throws Exception {
        File parquet = dir.resolve("library.parquet").toFile();
        try (Connection conn = DriverManager.getConnection("jdbc:duckdb:");
             Statement stmt = conn.createStatement()) {
            stmt.execute("COPY (SELECT CAST(1 AS FLOAT) AS PrecursorMz) TO '" +
                    parquet.getAbsolutePath().replace("\\", "/") + "' (FORMAT PARQUET" +
                    MetaFile.kvMetadataOption(footer) + ")");
        }
        return parquet;
    }

    // DuckDB takes keys that differ only in ASCII case for one key, so a footer holding both cannot be
    // written. The format allows them - read takes them as two keys, as FragCast and FragSpecLib do -
    // so checkWritable refuses them, naming both lines, for a predictor to stop on before it predicts
    // anything rather than after. The COPY is run too: when DuckDB records such keys, the check can go.
    @Test
    public void keysDifferingOnlyInAsciiCaseAreReadButRefusedAsUnwritable(@TempDir Path dir) throws Exception {
        MetaFile meta = metaFile(dir, "database_name = /x/a.fasta\n# x\nDatabase_Name = /x/b.fasta\n");
        assertEquals(pairs(pair("database_name", "/x/a.fasta"), pair("Database_Name", "/x/b.fasta")),
                inOrder(meta.footer()));

        IllegalArgumentException e = assertThrows(IllegalArgumentException.class, meta::checkWritable);
        assertEquals(dir.resolve(PATH) + ":3: key 'Database_Name' differs from 'database_name' (line 1) " +
                "only in case, which MSBooster's parquet writer (DuckDB) cannot record as two keys", e.getMessage());

        SQLException duckDb = assertThrows(SQLException.class, () -> copyWithFooter(dir, meta.footer()));
        assertTrue(duckDb.getMessage().contains("Duplicate struct entry name"), duckDb.getMessage());
    }

    // An entry name is written into the SQL as it stands, and DuckDB's parser stops at a NUL.
    @Test
    public void aKeyHoldingNulIsRefusedAsUnwritable(@TempDir Path dir) throws Exception {
        MetaFile meta = metaFile(dir, "a = 1\nk" + (char) 0 + "ey = v\n");
        IllegalArgumentException e = assertThrows(IllegalArgumentException.class, meta::checkWritable);
        assertEquals(dir.resolve(PATH) + ":2: a key holds U+0000, which MSBooster's parquet writer (DuckDB) " +
                "cannot record", e.getMessage());
    }

    // Only ASCII letters are folded, as DuckDB folds them: an accented capital and its small letter,
    // the Kelvin sign and 'k', and the sharp s and "ss" are different keys to it, recorded as given.
    @Test
    public void onlyAsciiCaseMakesKeysUnwritable(@TempDir Path dir) throws Exception {
        MetaFile meta = metaFile(dir, (char) 0xE9 + " = 1\n" + (char) 0xC9 + " = 2\n" +
                "k = 3\n" + (char) 0x212A + " = 4\n" + "ss = 5\n" + (char) 0xDF + " = 6\n" +
                "Key_1 = 7\n");
        meta.checkWritable();
        assertEquals(inOrder(meta.footer()), inOrder(ParquetFooters.read(copyWithFooter(dir, meta.footer()))));
    }

    @Test
    public void withoutMetaThereIsNothingToRefuse() {
        MetaFile meta = new MetaFile();
        meta.read();
        meta.checkWritable();
    }
}
