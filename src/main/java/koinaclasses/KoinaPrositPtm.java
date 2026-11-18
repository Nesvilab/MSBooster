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
import java.util.Collections;
import java.util.HashMap;
import java.util.LinkedHashMap;
import java.util.LinkedHashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.concurrent.ConcurrentHashMap;

/**
 * Which mods Koina's Prosit-PTM models (Prosit_2025_intensity_22PTM, Prosit_2025_intensity_40PTM,
 * Prosit_2025_irt_22PTM) can be sent. Their preprocessing encodes each modified position by the atoms (H, C, N, O,
 * P, S) it gains and loses: a pair of site and UniMod code in its hand-written table is encoded as the table says;
 * any other pair goes through its zero-shot encoding, which takes the site's attachment group as the loss and the
 * attachment group plus the code's delta composition as the gain. A pair neither can encode makes Koina answer the
 * whole request with HTTP 400. {@value #RESOURCE} holds the hand-written pairs, the attachment groups and the
 * monosaccharides the zero-shot encoding reads, copied from Koina's code (see the header of that file).
 * <p>
 * A pair is encodable when, as in Koina's preprocessing:
 * <ul>
 *   <li>its code is in the unimod.obo Koina ships ({@link KoinaUnimod}), since the code is looked up there before
 *       anything else; and</li>
 *   <li>the pair is in the hand-written table; or</li>
 *   <li>the site has an attachment group (C, D, H, I, K, L, M, Q, R, S, T, W, Y or the peptide N-terminus), the
 *       code's delta composition is read by the zero-shot encoding (elements H, C, N, O, P, S, isotopes counted as
 *       their element, and the monosaccharides it lists, written by its token pattern: dHex is not), and the
 *       attachment group plus that composition is not negative for any element.</li>
 * </ul>
 * Koina reads a C-terminal mod as one on the N-terminus, so a C-terminal mod is never sent (see
 * {@link KoinaModelSpec.EncodableMods#PROSIT_PTM}). It encodes only the first mod of a residue; the mass notation
 * MSBooster writes peptides in has at most one mod per residue, the C-terminal one apart.
 * <p>
 * The tables load on first use and read no MSBooster setting.
 */
public final class KoinaPrositPtm {
    public static final String RESOURCE = "/ptm_resources/koina_prosit_ptm.tsv";
    /** The site of a mod on the peptide N-terminus in {@value #RESOURCE} and in {@link #refusal}. */
    public static final char NTERM = 'n';
    /** The elements Koina counts, in the order of its encoding. */
    static final List<String> ELEMENTS = List.of("H", "C", "N", "O", "P", "S");
    //Gly, whose mass is that of carbamidomethyl
    private static final String GLY = "1263";

    //residue (or "[") + code, e.g. "S21" or "[737", to whether it is encodable, filled as pairs are asked about
    private static final Map<String, Boolean> KNOWN = new ConcurrentHashMap<>();

    private KoinaPrositPtm() {}

    //the tables, read on first use
    private static final class Tables {
        //site + code of the hand-written pairs, e.g. "S21", "n1"
        static final Set<String> HAND = new LinkedHashSet<>();
        //site to the atoms of its attachment group
        static final Map<Character, int[]> ATTACHMENT = new HashMap<>();
        //monosaccharide symbol to its atoms
        static final Map<String, int[]> SUGARS = new LinkedHashMap<>();

        static {
            load(HAND, ATTACHMENT, SUGARS);
        }

        //the hand-written pairs as residue (or "[") and code, and those a mass is matched against first
        static final Set<String> HAND_MODS = handMods(HAND);
        static final Set<String> PREFERRED = preferred(HAND_MODS);
    }

    /**
     * Whether a Prosit-PTM model can be sent this mod on a residue or the N-terminus.
     *
     * @param aaCode the residue (or "[" for the N-terminus) and UniMod code, as in the mod sets of a
     *               {@link KoinaModelSpec.ModRule}, e.g. "S21" or "[737"
     */
    public static boolean encodes(String aaCode) {
        return KNOWN.computeIfAbsent(aaCode, code -> {
            char site = code.charAt(0) == '[' ? NTERM : code.charAt(0);
            return refusal(site, code.substring(1)) == null;
        });
    }

    /**
     * Why Koina cannot encode this mod, as its preprocessing would say, or null if it can.
     *
     * @param site       the residue, or {@link #NTERM} for the peptide N-terminus
     * @param unimodCode the code as written after "UNIMOD:"
     */
    public static String refusal(char site, String unimodCode) {
        //Koina looks every code up in its unimod.obo first. Its lookup lets one code of the hand-written table that is
        //not in it through (5634, which is not a UniMod code); MSBooster has no mass for it, so it never writes it
        String composition = KoinaUnimod.composition(unimodCode);
        if (composition == null) {
            return "UNIMOD:" + unimodCode + " is not in Koina's unimod.obo";
        }
        if (inHandTable(site, unimodCode)) {
            return null;
        }
        int[] attachment = Tables.ATTACHMENT.get(site);
        if (attachment == null) {
            return "no attachment group for " + siteName(site);
        }
        int[] delta;
        try {
            delta = atoms(composition);
        } catch (IllegalArgumentException e) {
            return e.getMessage();
        }
        for (int i = 0; i < delta.length; i++) {
            if (attachment[i] + delta[i] < 0) {
                return siteName(site) + " plus UNIMOD:" + unimodCode + " (" + composition + ") has a negative "
                        + ELEMENTS.get(i) + " count";
            }
        }
        return null;
    }

    /**
     * Whether the pair is in Koina's hand-written table. A pair of it whose code is not in Koina's unimod.obo is
     * not encodable all the same.
     */
    public static boolean inHandTable(char site, String unimodCode) {
        return Tables.HAND.contains(site + unimodCode);
    }

    /**
     * The pairs of Koina's hand-written table, as residue (or "[" for the N-terminus) and code, as in the mod sets
     * of a {@link KoinaModelSpec.ModRule}.
     */
    public static Set<String> handTableMods() {
        return Tables.HAND_MODS;
    }

    /**
     * The mods a mass is matched against first for the Prosit-PTM models ({@link KoinaModelSpec.ModRule#PROSIT_PTM}):
     * the pairs of Koina's hand-written table, but Gly (UNIMOD:1263), which has the mass and, in that table, the
     * atoms of carbamidomethyl (UNIMOD:4) on C and K. Carbamidomethyl is written instead, as the other models are
     * sent it.
     */
    public static Set<String> preferredMods() {
        return Tables.PREFERRED;
    }

    private static Set<String> handMods(Set<String> hand) {
        Set<String> mods = new LinkedHashSet<>();
        for (String pair : hand) {
            mods.add(pair.charAt(0) == NTERM ? "[" + pair.substring(1) : pair);
        }
        return Collections.unmodifiableSet(mods);
    }

    private static Set<String> preferred(Set<String> handMods) {
        Set<String> mods = new LinkedHashSet<>();
        for (String mod : handMods) {
            if (!mod.substring(1).equals(GLY)) {
                mods.add(mod);
            }
        }
        return Collections.unmodifiableSet(mods);
    }

    /**
     * The H, C, N, O, P and S counts of a delta composition, read as the zero-shot encoding reads one: isotopes
     * count as their element, and a monosaccharide as its atoms.
     *
     * @throws IllegalArgumentException for a token or symbol the zero-shot encoding cannot read
     */
    public static int[] atoms(String composition) {
        List<KoinaUnimod.Token> tokens = KoinaUnimod.tokens(composition);
        int[] counts = new int[ELEMENTS.size()];
        for (KoinaUnimod.Token token : tokens) {
            int[] sugar = Tables.SUGARS.get(token.symbol());
            int element = ELEMENTS.indexOf(token.symbol());
            if (sugar != null) {
                for (int i = 0; i < counts.length; i++) {
                    counts[i] += sugar[i] * token.count();
                }
            } else if (element >= 0) {
                counts[element] += token.count();
            } else {
                throw new IllegalArgumentException("element '" + token.symbol() + "' cannot be encoded (only H C N "
                        + "O P S)");
            }
        }
        return counts;
    }

    private static String siteName(char site) {
        return site == NTERM ? "the N-terminus" : "residue " + site;
    }

    //atoms of a composition the resource writes, all of them elements
    private static int[] elementAtoms(String composition) {
        int[] counts = new int[ELEMENTS.size()];
        if (composition.isBlank()) {
            return counts;
        }
        for (KoinaUnimod.Token token : KoinaUnimod.tokens(composition)) {
            int element = ELEMENTS.indexOf(token.symbol());
            if (element < 0) {
                throw new IllegalStateException(RESOURCE + ": '" + token.symbol() + "' is not an element Koina counts");
            }
            counts[element] += token.count();
        }
        return counts;
    }

    private static void load(Set<String> hand, Map<Character, int[]> attachment, Map<String, int[]> sugars) {
        InputStream stream = KoinaPrositPtm.class.getResourceAsStream(RESOURCE);
        if (stream == null) {
            throw new IllegalStateException(RESOURCE + " is missing from the MSBooster jar");
        }
        try (BufferedReader reader = new BufferedReader(new InputStreamReader(stream, StandardCharsets.UTF_8))) {
            String line;
            boolean header = true;
            while ((line = reader.readLine()) != null) {
                if (line.startsWith("#") || line.isBlank()) {
                    continue;
                }
                if (header) { //kind, name, unimod_id, gain, loss
                    header = false;
                    continue;
                }
                String[] fields = line.split("\t", -1);
                switch (fields[0]) {
                    case "hand":
                        hand.add(fields[1] + fields[2]);
                        break;
                    case "attachment":
                        attachment.put(fields[1].charAt(0), elementAtoms(fields[4]));
                        break;
                    case "sugar":
                        sugars.put(fields[1], elementAtoms(fields[3]));
                        break;
                    default:
                        throw new IllegalStateException(RESOURCE + ": unknown row kind '" + fields[0] + "'");
                }
            }
        } catch (IOException e) {
            throw new UncheckedIOException("Could not read " + RESOURCE, e);
        }
    }
}
