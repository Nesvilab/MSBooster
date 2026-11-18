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

import java.io.BufferedReader;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.io.UncheckedIOException;
import java.nio.charset.StandardCharsets;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

/**
 * The UniMod delta compositions as Koina reads them. Koina's DeepLC and IM2Deep preprocessing, and the zero-shot
 * encoding of its Prosit-PTM models, turn a [UNIMOD:x] into atom counts from the delta composition of x in the
 * unimod.obo Koina ships. That file is older than the one MSBooster bundles (ptm_resources/unimod.obo), so the codes
 * and compositions here are Koina's: {@value #RESOURCE} holds every term of it that has a delta composition (see the
 * header of that file for its source).
 * <p>
 * A composition is a space-separated list of tokens such as {@code H(2) C(2) O}, {@code 13C(6) N(-2) 15N(2)} or
 * {@code Hex(5) HexNAc(2)}: an optional isotope mass number, an element or monosaccharide symbol, and an optional
 * count in parentheses. The table loads on first use and reads no MSBooster setting.
 */
public final class KoinaUnimod {
    public static final String RESOURCE = "/ptm_resources/koina_unimod_compositions.tsv";
    /** The elements Koina's encodings count atoms of. */
    public static final Set<String> ENCODED_ELEMENTS = Set.of("C", "H", "N", "O", "S", "P");

    //a token as Koina's Prosit zero-shot encoding reads one (zero_shot.py, _TOKEN), matched in full
    private static final Pattern TOKEN = Pattern.compile("(\\d*)([A-Z][A-Za-z]*?)(?:\\((-?\\d+)\\))?");

    private KoinaUnimod() {}

    /**
     * One token of a composition.
     *
     * @param isotope the isotope mass number, such as "13" in 13C, or "" for the element itself
     * @param symbol  an element (C, Na) or monosaccharide (Hex, HexNAc) symbol
     * @param count   how many, 1 if the token gives no count; negative for atoms the mod takes away
     */
    public record Token(String isotope, String symbol, int count) {
        public boolean isIsotope() {
            return !isotope.isEmpty();
        }
    }

    //the table, read on first use
    private static final class Table {
        static final Map<Integer, String> COMPOSITIONS = load();
    }

    /** The delta composition of a UniMod code in Koina's unimod.obo, or null if Koina has none for it. */
    public static String composition(int unimodId) {
        return Table.COMPOSITIONS.get(unimodId);
    }

    /**
     * The delta composition of a UniMod code, given as written after "UNIMOD:", or null if it is not a number or
     * Koina has no composition for it. A mod of the additional mods file without a UniMod code has none.
     */
    public static String composition(String unimodCode) {
        if (unimodCode == null || unimodCode.isEmpty() || unimodCode.length() > 9) {
            return null;
        }
        for (int i = 0; i < unimodCode.length(); i++) {
            if (!Character.isDigit(unimodCode.charAt(i))) {
                return null;
            }
        }
        return composition(Integer.parseInt(unimodCode));
    }

    /**
     * The tokens of a composition, read as the zero-shot encoding of Koina's Prosit-PTM models reads them. Like its
     * pattern, this one cannot read a symbol that starts with a lower-case letter, such as dHex.
     *
     * @throws IllegalArgumentException for a token Koina cannot read
     */
    public static List<Token> tokens(String composition) {
        List<Token> tokens = new ArrayList<>();
        for (String token : composition.trim().split("\\s+")) {
            Matcher m = TOKEN.matcher(token);
            if (!m.matches()) {
                throw new IllegalArgumentException("cannot read composition token '" + token + "' of '"
                        + composition + "'");
            }
            tokens.add(new Token(m.group(1), m.group(2), m.group(3) == null ? 1 : Integer.parseInt(m.group(3))));
        }
        return tokens;
    }

    /**
     * Whether Koina has a composition for this UniMod code made only of the elements C, H, N, O, S and P, and, unless
     * {@code isotopesAllowed}, of none of their isotopes (13C, 15N, 2H, 18O).
     *
     * @param unimodCode the code as written after "UNIMOD:"
     */
    public static boolean hasOnlyEncodedElements(String unimodCode, boolean isotopesAllowed) {
        String composition = composition(unimodCode);
        if (composition == null) {
            return false;
        }
        List<Token> tokens;
        try {
            tokens = tokens(composition);
        } catch (IllegalArgumentException e) {
            return false;
        }
        for (Token token : tokens) {
            if (!ENCODED_ELEMENTS.contains(token.symbol()) || (token.isIsotope() && !isotopesAllowed)) {
                return false;
            }
        }
        return true;
    }

    private static Map<Integer, String> load() {
        InputStream stream = KoinaUnimod.class.getResourceAsStream(RESOURCE);
        if (stream == null) {
            throw new IllegalStateException(RESOURCE + " is missing from the MSBooster jar");
        }
        Map<Integer, String> compositions = new HashMap<>();
        try (BufferedReader reader = new BufferedReader(new InputStreamReader(stream, StandardCharsets.UTF_8))) {
            String line;
            boolean header = true;
            while ((line = reader.readLine()) != null) {
                if (line.startsWith("#") || line.isBlank()) {
                    continue;
                }
                if (header) { //unimod_id, name, delta_composition
                    header = false;
                    continue;
                }
                String[] fields = line.split("\t", -1);
                compositions.put(Integer.parseInt(fields[0]), fields[2]);
            }
        } catch (IOException e) {
            throw new UncheckedIOException("Could not read " + RESOURCE, e);
        }
        return Collections.unmodifiableMap(compositions);
    }
}
