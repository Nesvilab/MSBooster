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

import java.io.IOException;
import java.nio.ByteBuffer;
import java.nio.charset.CharacterCodingException;
import java.nio.charset.CodingErrorAction;
import java.nio.charset.StandardCharsets;
import java.nio.file.AccessDeniedException;
import java.nio.file.Files;
import java.nio.file.InvalidPathException;
import java.nio.file.NoSuchFileException;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;

/**
 * The {@code --meta} option of both predictors: a text file of key-value pairs that a predicted
 * library records, verbatim and in file order, in its Parquet footer.
 *
 * <p>FragDIA reads the enzyme a library was digested with out of that footer, and FragPipe compares
 * the footers of a library's shards; what goes into it is FragPipe's business. FragPipe writes the
 * file, every library builder - FragCast's {@code build-library}, FragSpecLib, and the two predictors
 * here - takes it as {@code --meta <file>}, and none of them knows a key name. {@link FragCastPredictor}
 * forwards the path to FragCast, which writes the footer itself, and both predictors write the same
 * pairs whenever MSBooster writes a library's Parquet file with DuckDB.
 *
 * <p>The file format, which every builder reads the same way:
 * <ul>
 *   <li>UTF-8; a file that is not is an error. One leading byte order mark is skipped.</li>
 *   <li>Lines end at {@code \n}, and are numbered from 1. Each line is trimmed of spaces, tabs and
 *       carriage returns at both ends - exactly those three, so a {@code \r\n} file reads like a
 *       {@code \n} one and no other character is ever trimmed.</li>
 *   <li>A line that is then empty, or starts with {@code #}, is skipped. There are no inline
 *       comments: a {@code #} anywhere else belongs to the key or value, as it can to a FASTA path.</li>
 *   <li>Any other line is {@code key = value}, split at its first {@code =}, so a value may contain
 *       {@code =}. Key and value are trimmed as the line is. An empty value is a value; a line with
 *       no {@code =}, an empty key and a key already given are errors, each naming the file and line.</li>
 *   <li>A file with no pairs records nothing, exactly as no {@code --meta} does.</li>
 * </ul>
 *
 * <p>The footer this side writes goes through DuckDB ({@link #kvMetadataOption}), which cannot record
 * two of the things the format allows: two keys that differ only in the case of ASCII letters, which
 * it takes for one key, and a key holding U+0000. A predictor that may write the footer itself
 * refuses both with {@link #checkWritable} before anything is predicted. FragCast and FragSpecLib
 * record them; FragPipe writes neither.
 *
 * <p>{@code --fasta} is not part of this. Each predictor still takes it, for the gene names a
 * {@code librarytsv} library reads from the FASTA, but what the footer records about the FASTA
 * ({@code database_name}, when FragPipe writes one) comes from the file like everything else.
 */
public final class MetaFile {
    /** The option that names the file. */
    public static final String FLAG = "--meta";

    /** The option as a usage line lists it. */
    public static final String USAGE = FLAG + " <file> (key = value lines recorded verbatim, in order, " +
            "in a parquet library's footer; given at most once)";

    /** Where the format is described, for a usage message to point at. */
    public static final String FORMAT_USAGE = FLAG + " <file>: one 'key = value' per line, UTF-8, split " +
            "at the first '='; blank lines and lines starting with '#' are skipped. Every pair is " +
            "written into a parquet library's footer as given, in file order; a library.tsv or MGF has " +
            "no footer. Without " + FLAG + ", or with a file of no pairs, the footer records nothing. " +
            "Keys that differ only in ASCII case, and a key holding U+0000, are refused: MSBooster's " +
            "parquet writer (DuckDB) cannot record them.";

    private static final char BOM = 0xFEFF;

    /** The path {@link #accept} took; null until it took one. */
    private String path;

    /** What {@link #read} found in the file, in file order. */
    private Map<String, String> pairs = Collections.emptyMap();

    /** The line {@link #read} found each key on, for {@link #checkWritable} to name. */
    private Map<String, Integer> lines = Collections.emptyMap();

    /**
     * Take {@code --meta}'s value from a predictor's command line. Nothing is read yet, so the
     * predictor's loop sees no difference between this and any other option.
     *
     * @throws IllegalArgumentException when {@code --meta} was already given
     */
    public void accept(String path) {
        if (this.path != null) {
            throw new IllegalArgumentException(FLAG + " was given more than once (" + this.path + ", then " +
                    path + "); a library records one file's pairs.");
        }
        this.path = path;
    }

    /** Whether {@code --meta} was given. */
    public boolean isSet() {
        return path != null;
    }

    /**
     * Read and check the file {@link #accept} took, so a file that cannot be used stops the run before
     * any prediction is paid for. Does nothing when {@code --meta} was not given.
     *
     * @throws IllegalArgumentException naming the file, and the line when a line is at fault
     */
    public void read() {
        if (path != null) {
            final Map<String, Integer> lines = new HashMap<>();
            pairs = Collections.unmodifiableMap(parse(path, lines));
            this.lines = lines;
        }
    }

    /**
     * Refuse the pairs this side's footer writer cannot record as given, so a file that has them stops
     * the run before any prediction is paid for rather than after it. MSBooster writes a footer with
     * DuckDB ({@link #kvMetadataOption}), which takes the pairs as a STRUCT: it compares entry names
     * without regard to the case of ASCII letters, so two keys that differ only there would be one
     * key to it, and an entry name is written into the SQL as it stands, which cannot hold U+0000.
     * A NUL in a value is recorded, so it is not refused. The format allows both, so {@link #read}
     * takes both; a predictor that may write the footer itself calls this right after it. Does nothing
     * when there are no pairs. Only meaningful after {@link #read}.
     *
     * @throws IllegalArgumentException naming the file, and the line of the key at fault
     */
    public void checkWritable() {
        //each key's ASCII lower case, mapped to the first key that had it
        final Map<String, String> folded = new HashMap<>();
        for (String key : pairs.keySet()) {
            final int line = lines.get(key);
            if (key.indexOf('\0') >= 0) {
                throw new IllegalArgumentException(path + ":" + line + ": a key holds U+0000, which " +
                        "MSBooster's parquet writer (DuckDB) cannot record");
            }
            final String first = folded.putIfAbsent(asciiLowerCase(key), key);
            if (first != null) {
                throw new IllegalArgumentException(path + ":" + line + ": key '" + key + "' differs from '" +
                        first + "' (line " + lines.get(first) + ") only in case, which MSBooster's parquet " +
                        "writer (DuckDB) cannot record as two keys");
            }
        }
    }

    /**
     * {@code s} with A-Z made a-z and nothing else changed: the case DuckDB ignores in an entry name,
     * and no more. {@link String#toLowerCase} would fold far more - an accented capital into its small
     * letter, the Kelvin sign into {@code k} - which DuckDB records as different keys.
     */
    private static String asciiLowerCase(String s) {
        final char[] chars = s.toCharArray();
        for (int i = 0; i < chars.length; i++) {
            if (chars[i] >= 'A' && chars[i] <= 'Z') {
                chars[i] += 'a' - 'A';
            }
        }
        return new String(chars);
    }

    /**
     * The library's Parquet footer: the file's pairs in file order. Empty when {@code --meta} was not
     * given or the file holds no pairs, and then no footer is written. Only meaningful after
     * {@link #read}.
     */
    public Map<String, String> footer() {
        return pairs;
    }

    /**
     * What a library builder's command line is given: {@code --meta} and the path exactly as this
     * one took it, or nothing when it took none, so a builder is handed nothing it was not given.
     */
    public List<String> builderArgs() {
        final List<String> args = new ArrayList<>(2);
        if (path != null) {
            args.add(FLAG);
            args.add(path);
        }
        return args;
    }

    /** What the footer will record, in one line for the log. Only meaningful after {@link #read}. */
    public String describe() {
        if (path == null) {
            return "No " + FLAG + " given: the library's footer records nothing";
        }
        if (pairs.isEmpty()) {
            return path + " holds no key-value pairs: the library's footer records nothing";
        }
        final StringBuilder sb = new StringBuilder("Library footer from ").append(path).append(':');
        String separator = " ";
        for (Map.Entry<String, String> pair : pairs.entrySet()) {
            sb.append(separator).append(pair.getKey()).append(" = ").append(pair.getValue());
            separator = "; ";
        }
        return sb.toString();
    }

    /**
     * Read {@code path} as a meta file.
     *
     * @return its pairs, in file order; empty for a file that holds none
     * @throws IllegalArgumentException when the file cannot be read, is not UTF-8, or has a line
     *                                  that is not {@code key = value}
     */
    static LinkedHashMap<String, String> parse(String path) {
        return parse(path, new HashMap<>());
    }

    /** {@link #parse(String)}, giving {@code lines} each key's line. */
    private static LinkedHashMap<String, String> parse(String path, Map<String, Integer> lines) {
        if (path.isEmpty()) {
            throw cannotRead(path, "no path was given");
        }
        final byte[] bytes;
        try {
            final Path file = Paths.get(path);
            if (Files.isDirectory(file)) {
                throw cannotRead(path, "it is a directory");
            }
            bytes = Files.readAllBytes(file);
        } catch (NoSuchFileException e) {
            throw cannotRead(path, "no such file");
        } catch (AccessDeniedException e) {
            throw cannotRead(path, "access denied");
        } catch (IOException | InvalidPathException e) {
            throw cannotRead(path, e.getMessage() == null ? e.toString() : e.getMessage());
        }
        return parse(path, bytes, lines);
    }

    /**
     * Parse a meta file's bytes; {@code path} only names the file in an error.
     *
     * @see #parse(String)
     */
    static LinkedHashMap<String, String> parse(String path, byte[] bytes) {
        return parse(path, bytes, new HashMap<>());
    }

    /** {@link #parse(String, byte[])}, giving {@code firstLine} each key's line. */
    private static LinkedHashMap<String, String> parse(String path, byte[] bytes, Map<String, Integer> firstLine) {
        String text;
        try {
            text = StandardCharsets.UTF_8.newDecoder()
                    .onMalformedInput(CodingErrorAction.REPORT)
                    .onUnmappableCharacter(CodingErrorAction.REPORT)
                    .decode(ByteBuffer.wrap(bytes))
                    .toString();
        } catch (CharacterCodingException e) {
            throw cannotRead(path, "it is not valid UTF-8");
        }
        if (!text.isEmpty() && text.charAt(0) == BOM) {
            text = text.substring(1);
        }

        final LinkedHashMap<String, String> pairs = new LinkedHashMap<>();
        final String[] lines = text.split("\n", -1);
        for (int i = 0; i < lines.length; i++) {
            final int lineNumber = i + 1;
            final String line = trim(lines[i]);
            if (line.isEmpty() || line.charAt(0) == '#') {
                continue;
            }
            final int equals = line.indexOf('=');
            if (equals < 0) {
                throw new IllegalArgumentException(path + ":" + lineNumber + ": expected 'key = value'");
            }
            final String key = trim(line.substring(0, equals));
            final String value = trim(line.substring(equals + 1));
            if (key.isEmpty()) {
                throw new IllegalArgumentException(path + ":" + lineNumber + ": empty key");
            }
            final Integer first = firstLine.putIfAbsent(key, lineNumber);
            if (first != null) {
                throw new IllegalArgumentException(path + ":" + lineNumber + ": duplicate key '" + key +
                        "' (first on line " + first + ")");
            }
            pairs.put(key, value);
        }
        return pairs;
    }

    /**
     * {@code s} without the spaces, tabs and carriage returns at either end. Not {@link String#trim},
     * which takes every control character with it: the set is part of the format, and every builder
     * trims exactly this one.
     */
    static String trim(String s) {
        int start = 0;
        int end = s.length();
        while (start < end && isTrimmed(s.charAt(start))) {
            start++;
        }
        while (end > start && isTrimmed(s.charAt(end - 1))) {
            end--;
        }
        return s.substring(start, end);
    }

    private static boolean isTrimmed(char c) {
        return c == ' ' || c == '\t' || c == '\r';
    }

    private static IllegalArgumentException cannotRead(String path, String reason) {
        return new IllegalArgumentException("cannot read meta file " + path + ": " + reason);
    }

    /**
     * The option that has DuckDB's {@code COPY ... (FORMAT PARQUET ...)} write {@code footer} into the
     * file's key-value metadata, in its order, with its leading comma:
     * {@code , KV_METADATA {'key': 'value', ...}}. Empty when the footer is, so a COPY with nothing to
     * record is exactly the COPY it always was. Keys and values are quoted alike, as DuckDB string
     * literals: the one escape is a doubled single quote, and a backslash is an ordinary character,
     * so a Windows path is recorded as given and no key can end its literal early. DuckDB's parser
     * stops at a NUL, so a value's NULs are spliced in as {@code chr(0)} between literals; a key
     * cannot be, and neither can two keys that differ only in ASCII case, which is why
     * {@link #checkWritable} refuses both first.
     */
    static String kvMetadataOption(Map<String, String> footer) {
        if (footer.isEmpty()) {
            return "";
        }
        final StringBuilder sb = new StringBuilder(", KV_METADATA {");
        boolean first = true;
        for (Map.Entry<String, String> entry : footer.entrySet()) {
            if (!first) {
                sb.append(", ");
            }
            first = false;
            sb.append(literal(entry.getKey())).append(": ").append(value(entry.getValue()));
        }
        return sb.append('}').toString();
    }

    /** {@code s} as a DuckDB string expression: its literal, with each NUL spliced in as {@code chr(0)}. */
    private static String value(String s) {
        final StringBuilder sb = new StringBuilder();
        int start = 0;
        int nul;
        while ((nul = s.indexOf('\0', start)) >= 0) {
            sb.append(literal(s.substring(start, nul))).append(" || chr(0) || ");
            start = nul + 1;
        }
        return sb.append(literal(s.substring(start))).toString();
    }

    private static String literal(String s) {
        return "'" + s.replace("'", "''") + "'";
    }
}
