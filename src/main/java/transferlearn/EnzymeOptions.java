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

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;
import java.util.regex.Pattern;

/**
 * The enzyme a predicted library is digested with, as the library builders take it on the command
 * line, and the Parquet footer that records it.
 *
 * <p>FragDIA no longer takes digestion flags of its own: it reads the enzyme out of the library's
 * Parquet footer. So every library builder - FragCast's {@code build-library}, FragSpecLib, and the
 * two predictors here - takes the same options under the same names and the same rules, and writes
 * the same footer keys. Both predictors hold them here: {@link FragCastPredictor} forwards them to
 * FragCast, which writes the footer itself, and both write the footer whenever MSBooster writes a
 * library's Parquet file with DuckDB.
 *
 * <p>The rules, which no builder relaxes:
 * <ul>
 *   <li>{@code --cut} is given once per enzyme, at most twice; uppercase letters A-Z, or {@code @}
 *       (cuts nothing: nonspecific).</li>
 *   <li>{@code --no-cut} and {@code --sense} are given exactly once per {@code --cut}, paired with it
 *       by order. {@code --no-cut} takes the syntax of {@code --cut} ({@code @}: no restriction);
 *       {@code --sense} is {@code C} (cut after) or {@code N} (cut before).</li>
 *   <li>{@code --enzyme-name} is given not at all, or exactly once per {@code --cut}.</li>
 *   <li>{@code --num-enzyme-termini} (0-3) and {@code --clip-n-term-m} (0 or 1) are each given
 *       exactly once when {@code --cut} is.</li>
 *   <li>Without {@code --cut}, every other one is an error. There are no defaults, and an empty
 *       value is an error.</li>
 * </ul>
 * A value that passes is used exactly as given - no case folding, no trimming - so what the footer
 * says is what the caller said.
 *
 * <p>{@code --fasta} is not one of these options, and none of these rules apply to it. It stays an
 * ordinary option of each predictor's own loop and keeps that loop's leniency, as it always has:
 * given more than once, the last one counts; given empty, it counts as not given, so the footer
 * records no {@code database_name} and FragCast is handed no {@code --fasta}. Where a FASTA is
 * required - for a {@code librarytsv} library, whose gene names come from it - an empty one is
 * refused as a missing one. What the footer records is provenance only: the path exactly as given.
 */
public final class EnzymeOptions {
    static final String ENZYME_NAME = "--enzyme-name";
    static final String CUT = "--cut";
    static final String NO_CUT = "--no-cut";
    static final String SENSE = "--sense";
    static final String NUM_ENZYME_TERMINI = "--num-enzyme-termini";
    static final String CLIP_N_TERM_M = "--clip-n-term-m";

    /** Every flag this class takes, in the order a usage line lists them. */
    public static final List<String> FLAGS = Collections.unmodifiableList(Arrays.asList(
            ENZYME_NAME, CUT, NO_CUT, SENSE, NUM_ENZYME_TERMINI, CLIP_N_TERM_M));

    /** The options as a usage line lists them. */
    public static final String USAGE = "--enzyme-name <name> --cut <residues|@> --no-cut <residues|@> " +
            "--sense <C|N> --num-enzyme-termini <0-3> --clip-n-term-m <0|1> (the enzyme recorded in the " +
            "library's Parquet footer: --cut, --no-cut and --sense once per enzyme, at most 2; " +
            "--enzyme-name once per enzyme or not at all; --num-enzyme-termini and --clip-n-term-m once; " +
            "none of them for a library that records no enzyme)";

    /**
     * What a predictor's usage says of its {@code --fasta}: where the footer records it, and that
     * it keeps the leniency described in the class description.
     */
    public static final String FASTA_USAGE = "--fasta <fasta>: required for librarytsv output, whose " +
            "gene names come from it. A parquet library records it in its footer as database_name, " +
            "exactly as given. Unlike the enzyme options it is not checked: an empty --fasta counts as " +
            "none, so nothing is recorded, and a --fasta given more than once counts its last.";

    /** The footer key a FASTA is recorded under. Provenance only: FragDIA never reads it for behaviour. */
    static final String DATABASE_NAME = "database_name";

    private static final int MAX_ENZYMES = 2;
    private static final Pattern RESIDUES = Pattern.compile("[A-Z]+|@");
    private static final Pattern SENSE_VALUE = Pattern.compile("[CN]");
    private static final Pattern TERMINI_VALUE = Pattern.compile("[0-3]");
    private static final Pattern CLIP_VALUE = Pattern.compile("[01]");

    private final List<String> names = new ArrayList<>();
    private final List<String> cuts = new ArrayList<>();
    private final List<String> noCuts = new ArrayList<>();
    private final List<String> senses = new ArrayList<>();
    private final List<String> numEnzymeTermini = new ArrayList<>();
    private final List<String> clipNTermM = new ArrayList<>();

    /**
     * Takes the enzyme options out of a predictor's command line into {@code enzyme}, and returns
     * the rest of it, in its order, for the predictor's own loop to parse exactly as it did before
     * there were enzyme options. Both predictors call this before their loops, so neither loop ever
     * sees an enzyme option or its value.
     *
     * <p>{@code args} is read as both predictors require it to be: flag-value pairs. Only a flag is
     * ever matched, so a value is never taken for a flag whatever it looks like: in
     * {@code --enzyme-name --cut} the enzyme is named "--cut", and in {@code --basename --cut} the
     * "--cut" stays with {@code --basename} in what is returned. That matters most to
     * {@link Predictor}, whose loop visits values as well as flags. Repeated options are kept in the
     * order given, which is how {@link #validate} pairs a second {@code --no-cut} with the second
     * {@code --cut}. A last flag with no value is returned with the rest, for the predictor to refuse.
     *
     * @return the arguments that are not enzyme options, in their order
     */
    public static String[] extract(String[] args, EnzymeOptions enzyme) {
        final List<String> rest = new ArrayList<>(args.length);
        int i = 0;
        for (; i + 1 < args.length; i += 2) {
            if (!enzyme.accept(args[i], args[i + 1])) {
                rest.add(args[i]);
                rest.add(args[i + 1]);
            }
        }
        if (i < args.length) {
            rest.add(args[i]);
        }
        return rest.toArray(new String[0]);
    }

    /**
     * Take {@code value} if {@code flag} is one of {@link #FLAGS}. Every one is recorded as given,
     * repeats included, so {@link #validate} can say what was wrong with the command line rather than
     * what was left of it.
     *
     * @return whether the flag was one of these
     */
    public boolean accept(String flag, String value) {
        switch (flag) {
            case ENZYME_NAME:
                names.add(value);
                return true;
            case CUT:
                cuts.add(value);
                return true;
            case NO_CUT:
                noCuts.add(value);
                return true;
            case SENSE:
                senses.add(value);
                return true;
            case NUM_ENZYME_TERMINI:
                numEnzymeTermini.add(value);
                return true;
            case CLIP_N_TERM_M:
                clipNTermM.add(value);
                return true;
            default:
                return false;
        }
    }

    /**
     * Check what was given against the rules in the class description.
     *
     * @throws IllegalArgumentException naming the first flag that breaks one, and how
     */
    public void validate() {
        if (cuts.isEmpty()) {
            for (String flag : FLAGS) {
                if (!flag.equals(CUT) && !valuesOf(flag).isEmpty()) {
                    throw new IllegalArgumentException(flag + " was given without " + CUT + ". The " +
                            "enzyme options describe the enzyme a " + CUT + " names, so either give " +
                            "the whole set or none of it.");
                }
            }
            return;
        }
        final int enzymes = cuts.size();
        if (enzymes > MAX_ENZYMES) {
            throw new IllegalArgumentException(CUT + " was given " + times(enzymes) + "; at most " +
                    MAX_ENZYMES + " enzymes are allowed.");
        }
        if (!names.isEmpty() && names.size() != enzymes) {
            throw new IllegalArgumentException(ENZYME_NAME + " was given " + times(names.size()) +
                    "; give it once per " + CUT + ", which was given " + times(enzymes) +
                    ", or not at all.");
        }
        requireCount(NO_CUT, enzymes, "once per " + CUT + ", which was given " + times(enzymes));
        requireCount(SENSE, enzymes, "once per " + CUT + ", which was given " + times(enzymes));
        requireCount(NUM_ENZYME_TERMINI, 1, "once when " + CUT + " is given");
        requireCount(CLIP_N_TERM_M, 1, "once when " + CUT + " is given");

        for (String name : names) {
            requireNonEmpty(ENZYME_NAME, name);
        }
        for (String cut : cuts) {
            requireMatch(CUT, cut, RESIDUES, "uppercase letters A-Z, or @");
        }
        for (String noCut : noCuts) {
            requireMatch(NO_CUT, noCut, RESIDUES, "uppercase letters A-Z, or @");
        }
        for (String sense : senses) {
            requireMatch(SENSE, sense, SENSE_VALUE, "C or N");
        }
        requireMatch(NUM_ENZYME_TERMINI, numEnzymeTermini.get(0), TERMINI_VALUE, "0, 1, 2 or 3");
        requireMatch(CLIP_N_TERM_M, clipNTermM.get(0), CLIP_VALUE, "0 or 1");
    }

    /** Whether an enzyme was given at all - that is, whether {@code --cut} was. */
    public boolean isSet() {
        return !cuts.isEmpty();
    }

    /**
     * The options as a library builder's command line takes them, each value exactly as given:
     * enzyme by enzyme ({@code --enzyme-name} when names were given, then {@code --cut},
     * {@code --no-cut}, {@code --sense}), then {@code --num-enzyme-termini} and
     * {@code --clip-n-term-m}. Empty when no enzyme was given, so a builder is handed nothing it was
     * not given. Only meaningful after {@link #validate}.
     */
    public List<String> builderArgs() {
        final List<String> args = new ArrayList<>();
        if (!isSet()) {
            return args;
        }
        for (int i = 0; i < cuts.size(); i++) {
            if (!names.isEmpty()) {
                args.add(ENZYME_NAME);
                args.add(names.get(i));
            }
            args.add(CUT);
            args.add(cuts.get(i));
            args.add(NO_CUT);
            args.add(noCuts.get(i));
            args.add(SENSE);
            args.add(senses.get(i));
        }
        args.add(NUM_ENZYME_TERMINI);
        args.add(numEnzymeTermini.get(0));
        args.add(CLIP_N_TERM_M);
        args.add(clipNTermM.get(0));
        return args;
    }

    /**
     * The key-value metadata a library's Parquet footer carries, in the keys every library builder
     * writes. With an enzyme: {@code search_enzyme_{name,cut,nocut,sense}_N} for each enzyme N (the
     * name only when names were given), then {@code num_enzyme_termini} and {@code clip_nTerm_M}.
     * With a FASTA: {@code database_name}, its path exactly as given. Nothing else, and nothing at all
     * when neither was given. Only meaningful after {@link #validate}.
     *
     * @param fasta the predictor's {@code --fasta}, the last one when it was given more than once;
     *              empty when it was not given or was given empty, and then no {@code database_name}
     *              is recorded
     */
    public LinkedHashMap<String, String> footer(String fasta) {
        final LinkedHashMap<String, String> footer = new LinkedHashMap<>();
        if (isSet()) {
            for (int i = 0; i < cuts.size(); i++) {
                final int n = i + 1;
                if (!names.isEmpty()) {
                    footer.put("search_enzyme_name_" + n, names.get(i));
                }
                footer.put("search_enzyme_cut_" + n, cuts.get(i));
                footer.put("search_enzyme_nocut_" + n, noCuts.get(i));
                footer.put("search_enzyme_sense_" + n, senses.get(i));
            }
            footer.put("num_enzyme_termini", numEnzymeTermini.get(0));
            footer.put("clip_nTerm_M", clipNTermM.get(0));
        }
        if (fasta != null && !fasta.isEmpty()) {
            footer.put(DATABASE_NAME, fasta);
        }
        return footer;
    }

    /** The enzyme in one line, for the log. Only meaningful after {@link #validate}. */
    public String describe() {
        if (!isSet()) {
            return "No enzyme given: the library's footer records none";
        }
        final StringBuilder sb = new StringBuilder("Library enzyme: ");
        for (int i = 0; i < cuts.size(); i++) {
            if (i > 0) {
                sb.append(" + ");
            }
            if (!names.isEmpty()) {
                sb.append(names.get(i)).append(' ');
            }
            sb.append("(cut ").append(cuts.get(i)).append(", no-cut ").append(noCuts.get(i))
                    .append(", sense ").append(senses.get(i)).append(')');
        }
        return sb.append("; num-enzyme-termini ").append(numEnzymeTermini.get(0))
                .append("; clip-n-term-m ").append(clipNTermM.get(0)).toString();
    }

    /**
     * The option that has DuckDB's {@code COPY ... (FORMAT PARQUET ...)} write {@code footer} into the
     * file's key-value metadata, with its leading comma: {@code , KV_METADATA {'key': 'value', ...}}.
     * Empty when the footer is, so a COPY with nothing to record is exactly the COPY it always was.
     * A DuckDB string literal has one escape - a single quote is doubled - and a backslash is an
     * ordinary character in it, so a Windows path is recorded as given.
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
            sb.append(literal(entry.getKey())).append(": ").append(literal(entry.getValue()));
        }
        return sb.append('}').toString();
    }

    private static String literal(String s) {
        return "'" + s.replace("'", "''") + "'";
    }

    private List<String> valuesOf(String flag) {
        switch (flag) {
            case ENZYME_NAME:
                return names;
            case CUT:
                return cuts;
            case NO_CUT:
                return noCuts;
            case SENSE:
                return senses;
            case NUM_ENZYME_TERMINI:
                return numEnzymeTermini;
            case CLIP_N_TERM_M:
                return clipNTermM;
            default:
                throw new IllegalArgumentException(flag);
        }
    }

    private void requireCount(String flag, int expected, String rule) {
        final int given = valuesOf(flag).size();
        if (given != expected) {
            throw new IllegalArgumentException(flag + " was given " + times(given) +
                    "; it is required exactly " + rule + ".");
        }
    }

    private static String times(int n) {
        return n + (n == 1 ? " time" : " times");
    }

    private static void requireNonEmpty(String flag, String value) {
        if (value.isEmpty()) {
            throw new IllegalArgumentException(flag + " was given an empty value.");
        }
    }

    private static void requireMatch(String flag, String value, Pattern pattern, String allowed) {
        requireNonEmpty(flag, value);
        if (!pattern.matcher(value).matches()) {
            throw new IllegalArgumentException(flag + " " + value + " is not valid: it takes " +
                    allowed + ", exactly as written.");
        }
    }
}
