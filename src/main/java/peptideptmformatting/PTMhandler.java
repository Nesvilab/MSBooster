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

package peptideptmformatting;

import allconstants.Constants;
import koinaclasses.KoinaModelSpec.EncodableMods;
import koinaclasses.KoinaModelSpec.ModMassMap;
import koinaclasses.KoinaModelSpec.PeptideFormat;
import org.apache.commons.lang3.math.NumberUtils;
import umich.ms.fileio.filetypes.unimod.UnimodOboReader;
import utils.NumericUtils;

import java.io.*;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.*;
import java.util.stream.Collectors;

import static utils.Print.printError;

public class PTMhandler {
    //handling of PTMs, all in one location.
    public static final double oxidationMass = 15.9949;
    public static final double carbamidomethylationMass = 57.0215;
    public static final double acetylationMass = 42.0106;
    public static final double phosphorylationMass = 79.9663;
    public static final double glyglyMass = 114.042927;
    public static final double pyrogluQMass = -17.0265;
    public static final double pyrogluEMass = -18.01057;
    public static final double pyroCarbamidomethylMass = 39.9950;
    public static final double citrullinationDeamidationMass = 0.984016;
    public static final int carbamidomethylationUnimod = 4;
    public static final int oxidationUnimod = 35;
    public static final int acetylationUnimod = 1;
    public static final int phosphorylationUnimod = 21;
    public static final int glyglyUnimod = 121;
    public static final int pyrogluQUnimod = 28;
    public static final int pyrogluEUnimod = 27;
    public static final int pyroCarbamidomethylUnimod = 26;
    public static final int citrullinationDeamidationUnimod = 7;
    public static final int glyUnimod = 1263; //Gly, which has the mass of carbamidomethyl
    public static final double lactylMass = 72.021129;
    public static final int lactylUnimod = 2114;
    public static final double methylMass = 14.01565;
    public static final double dimethylMass = 28.0313;
    public static final double trimethylMass = 42.04695;
    public static final int methylUnimod = 34;
    public static final int dimethylUnimod = 36;
    public static final int trimethylUnimod = 37;
    public static final double tmtZeroMass = 224.152478; //TMT (TMT zero)
    public static final int tmtZeroUnimod = 739;
    //what Koina's Chronologer reads C[UNIMOD:385] as: a pyro-carbamidomethyl cysteine, carbamidomethyl and the
    //ammonia loss of the cyclization (UNIMOD:385) together
    public static final int chronologerPyroCarbamidomethylUnimod = 385;

    //models trained on tmt6/10/11 extrapolate well to tmt pro and itraq
    public static final double tmt6_10_11Mass = 229.1629;
    public static final double tmtproMass = 304.20715;
    public static final double itraqMass = 144.102063;
    public static final double itraq8Mass = 304.20535;
    public static final int tmt6_10_11Unimod = 737;
    public static final int tmtproUnimod = 2016;
    public static final int itraqUnimod = 214;
    public static final int itraq8Unimod = 730;
    public static final double[] tmtMasses = {tmt6_10_11Mass, tmtproMass, itraqMass, itraq8Mass};
    private static double tmtMass = tmt6_10_11Mass;
    public static int tmtUnimod = tmt6_10_11Unimod;

    public static void setTmtMass(double mass) {
        tmtMass = mass;
        updateUnimodToModMassLimited();
        updateAAUnimodToModMassLimited();
        updateChronologerTables();
    }
    public static double getTmtMass() { return tmtMass; }

    //whether the search in fragger.params has an isobaric label (TMT or iTRAQ), as a fixed or variable mod
    private static boolean isobaricLabel = false;

    /** The search has an isobaric label of this mass (ParameterUtils.readFraggerParams); it becomes the TMT mass. */
    public static void setIsobaricLabel(double mass) {
        setTmtMass(mass);
        isobaricLabel = true;
    }

    /** Whether the search has an isobaric label (TMT or iTRAQ), whose mass is then {@link #getTmtMass()}. */
    public static boolean hasIsobaricLabel() {
        return isobaricLabel;
    }

    /** Forgets the search's isobaric label, and sets the TMT mass back to TMT 6/10/11. For tests. */
    public static void resetIsobaricLabel() {
        isobaricLabel = false;
        setTmtMass(tmt6_10_11Mass);
    }

    public static HashMap<Integer, Double> unimodToModMassLimited = new HashMap<>();
    static {
        updateUnimodToModMassLimited();
    }

    private static void updateUnimodToModMassLimited() {
        unimodToModMassLimited.clear();
        unimodToModMassLimited.put(carbamidomethylationUnimod, carbamidomethylationMass);
        unimodToModMassLimited.put(oxidationUnimod, oxidationMass);
        unimodToModMassLimited.put(acetylationUnimod, acetylationMass);
        unimodToModMassLimited.put(phosphorylationUnimod, phosphorylationMass);
        unimodToModMassLimited.put(glyglyUnimod, glyglyMass);
        unimodToModMassLimited.put(tmtUnimod, tmtMass);
        unimodToModMassLimited.put(pyrogluQUnimod, pyrogluQMass);
        unimodToModMassLimited.put(pyrogluEUnimod, pyrogluEMass);
        unimodToModMassLimited.put(pyroCarbamidomethylUnimod, pyroCarbamidomethylMass);
        unimodToModMassLimited.put(citrullinationDeamidationUnimod, citrullinationDeamidationMass);
    }

    private static void updateAAUnimodToModMassLimited() {
        AAunimodToModMassLimited.clear();
        AAunimodToModMassLimited.put("C" + carbamidomethylationUnimod, carbamidomethylationMass);
        AAunimodToModMassLimited.put("M" + oxidationUnimod, oxidationMass);
        AAunimodToModMassLimited.put("[" + acetylationUnimod, acetylationMass);
        AAunimodToModMassLimited.put("S" + phosphorylationUnimod, phosphorylationMass);
        AAunimodToModMassLimited.put("T" + phosphorylationUnimod, phosphorylationMass);
        AAunimodToModMassLimited.put("Y" + phosphorylationUnimod, phosphorylationMass);
        AAunimodToModMassLimited.put("[" + glyglyUnimod, glyglyMass);
        AAunimodToModMassLimited.put("K" + glyglyUnimod, glyglyMass);
        AAunimodToModMassLimited.put("T" + glyglyUnimod, glyglyMass);
        AAunimodToModMassLimited.put("C" + glyglyUnimod, glyglyMass);
        AAunimodToModMassLimited.put("S" + glyglyUnimod, glyglyMass);
        AAunimodToModMassLimited.put("K" + tmtUnimod, tmtMass);
        AAunimodToModMassLimited.put("S" + tmtUnimod, tmtMass);
        AAunimodToModMassLimited.put("T" + tmtUnimod, tmtMass);
        AAunimodToModMassLimited.put("H" + tmtUnimod, tmtMass);
        AAunimodToModMassLimited.put("[" + tmtUnimod, tmtMass);
        AAunimodToModMassLimited.put("Q" + pyrogluQUnimod, pyrogluQMass);
        AAunimodToModMassLimited.put("E" + pyrogluEUnimod, pyrogluEMass);
        AAunimodToModMassLimited.put("C" + pyroCarbamidomethylUnimod, pyroCarbamidomethylMass);
        AAunimodToModMassLimited.put("Q" + citrullinationDeamidationUnimod, citrullinationDeamidationMass);
        AAunimodToModMassLimited.put("R" + citrullinationDeamidationUnimod, citrullinationDeamidationMass);
        AAunimodToModMassLimited.put("N" + citrullinationDeamidationUnimod, citrullinationDeamidationMass);
    }
    public static final HashMap<String, Double> AAunimodToModMassLimited = new HashMap<>();
    static {
        updateAAUnimodToModMassLimited();
    }

    @SuppressWarnings("unchecked")
    public static <T> Map<T, Double> makeUnimodToModMassAll(boolean includeAA) throws IOException {
        ArrayList<String> modPaths = new ArrayList<>();
        modPaths.add("/ptm_resources/alphapept_koina.csv");
        modPaths.add("/ptm_resources/modification_alphapeptdeep.tsv");
        if (!Constants.additionalModsFile.isEmpty()) {
            modPaths.add(Constants.additionalModsFile);
        }

        Map<T, Double> map = new HashMap<>();
        for (String modPath : modPaths) {
            //get unimod code and monoisotopic mass
            InputStream stream = PTMhandler.class.getResourceAsStream(modPath);
            if (Objects.isNull(stream)) {
                stream = Files.newInputStream(Paths.get(modPath));
            }
            final InputStreamReader reader = new InputStreamReader(stream);
            final BufferedReader ptmFile = new BufferedReader(reader);

            String delimiter = "";
            if (modPath.endsWith(".csv")) {
                delimiter = ",";
            } else if (modPath.endsWith(".tsv")) {
                delimiter = "\t";
            } else {
                printError("Invalid modification file: " + modPath);
                printError("Exiting");
                System.exit(1);
            }

            String line = ptmFile.readLine(); //header
            List<String> fields = Arrays.asList(line.split(delimiter));
            int unimodIdIdx = fields.indexOf("unimod_id");
            int massIdx = Math.max(fields.indexOf("unimod_mass"), fields.indexOf("mod_mass"));
            int nameIdx = fields.indexOf("mod_name");

            if (massIdx == -1) {
                printError("Modification file " + modPath + " is missing 'unimod_mass' or 'mod_mass' column");
                printError("Exiting");
                System.exit(1);
            }
            if (nameIdx == -1) {
                printError("Modification file " + modPath + " is missing 'mod_name' column");
                printError("Exiting");
                System.exit(1);
            }

            while ((line = ptmFile.readLine()) != null) {
                String[] lineSplit = line.split(delimiter, -1);

                String unimod;
                if (unimodIdIdx != -1) {
                    unimod = lineSplit[unimodIdIdx];
                } else { //custom mods
                    unimod = lineSplit[nameIdx].split("@")[0];
                }
                Double mass = Double.parseDouble(lineSplit[massIdx]);
                if (includeAA) {
                    if (modPath.endsWith(".csv")) {
                        map.put((T) (lineSplit[nameIdx].charAt(0) + unimod), mass);
                    } else {
                        //TODO: "Anywhere" should be expanded
                        String localization = lineSplit[nameIdx].split("@")[1];
                        if (localization.contains("term") || localization.equals("n") || localization.equals("c")) {
                            localization = "["; //to match logic for alphapept_koina.csv
                        }
                        map.put((T) (localization + unimod), mass);
                    }
                } else {
                    if (unimodIdIdx != -1) {
                        map.put((T) Integer.valueOf(unimod), mass);
                    }
                    //custom mods without unimod_id are skipped for Integer-keyed map
                }
            }
        }
        return map;
    }
    public static Map<Integer, Double> unimodToModMassAll;

    static {
        try {
            unimodToModMassAll = makeUnimodToModMassAll(false);
        } catch (IOException e) {
            throw new RuntimeException(e);
        }
    }

    public static Map<String, Double> AAunimodToModMassAll;
    public static Set<String> AAunimodToModMassAllKeys;

    static {
        try {
            AAunimodToModMassAll = makeUnimodToModMassAll(true);
            AAunimodToModMassAllKeys = AAunimodToModMassAll.keySet();
        } catch (IOException e) {
            throw new RuntimeException(e);
        }
    }

    //AAunimodToModMassAll as it was when the Prosit-PTM table was last made from it, and that table
    private record DerivedMassTable(Map<String, Double> source, Map<String, Double> table) {}
    private static volatile DerivedMassTable aaUnimodToModMassPrositPtm;

    /**
     * The mass table of Koina's Prosit-PTM models (KoinaModelSpec.ModMassMap.PROSIT_PTM): AAunimodToModMassAll with
     * the two pairs of their hand-written table it lacks, phospho on P (P21) and Gly on C (C1263). It is made again
     * when AAunimodToModMassAll is, as when an additional mods file is read.
     */
    public static Map<String, Double> getAAunimodToModMassPrositPtm() {
        Map<String, Double> source = AAunimodToModMassAll;
        DerivedMassTable derived = aaUnimodToModMassPrositPtm;
        if (derived == null || derived.source() != source) {
            Map<String, Double> table = new HashMap<>(source);
            putIfKnown(table, "P", phosphorylationUnimod);
            putIfKnown(table, "C", glyUnimod);
            derived = new DerivedMassTable(source, Collections.unmodifiableMap(table));
            aaUnimodToModMassPrositPtm = derived;
        }
        return derived.table();
    }

    //the residue with the mass unimodToModMassAll has for the code, unless the table has the pair already
    private static void putIfKnown(Map<String, Double> table, String residue, int unimod) {
        Double mass = unimodToModMassAll.get(unimod);
        if (mass != null) {
            table.putIfAbsent(residue + unimod, mass);
        }
    }

    public static Map<Integer, Double> unimodOboToModMass;

    public static void setUnimodObo() {
        try { //only do this when a library tsv or parquet predicted library is specified
            unimodOboToModMass = new HashMap<>();
            File f = new File(Constants.unimodObo);
            if (! f.exists()) {
                printError("When uploading prediction in spectral library in tsv or parquet format, " +
                        "a unimod.obo file must be provided. " +
                        "In the command line, provide a path at --unimodObo. Exiting.");
                System.exit(1);
            }
            UnimodOboReader uobo = new UnimodOboReader(Paths.get(Constants.unimodObo));
            for (Map.Entry<String, Float> entry : uobo.unimodMassMap.entrySet()) {
                unimodOboToModMass.put(Integer.parseInt(entry.getKey().split(":")[1]), Double.valueOf(entry.getValue()));
            }
        } catch (Exception e) {
            throw new RuntimeException(e);
        }
    }

    //takes start and end coordinates of mod in peptide formatter base field
    //takes model format
    //returns newly formatted peptide string in format of that model
    //if PTM would cause an error or result in the same fragment intensities (PTM not supported), just removes it from string
    //gives some leeway in PTM mass, precision, only needs to be 0.001 away
    //foundUnimods are those PTMs to check first (those that have already been identified in the sample)

    //how a peptide format writes a mod it keeps
    //unimodLabel: the word before the code, as in [UNIMOD:35]
    //ntermPrefix, ntermSuffix: written before and after an N-term mod
    //ctermPrefix: written before a C-term mod
    //customModAsMass: a mod of the additional mods file without a UniMod code is written as its mass
    //parentheses: mods are written in () instead of []
    private record ModNotation(String unimodLabel, String ntermPrefix, String ntermSuffix, String ctermPrefix,
                               boolean customModAsMass, boolean parentheses) {}

    private static final ModNotation DIANN_NOTATION = new ModNotation("UniMod", "", "", "", false, false);
    private static final ModNotation LIBRARYTSV_NOTATION = new ModNotation("UniMod", "_", "", "_", true, true);

    //start is index of [
    //the DIA-NN ("diann") and library.tsv ("librarytsv") formats. A Koina model's peptides are written by
    //formatPeptideBaseToKoina
    public static String[] formatPeptideBaseToSpecific(String peptide, int start, int end, String model,
                                                       Set<String> foundUnimods, boolean cterm) {
        switch (model) {
            case "diann":
                return formatMod(peptide, start, end, DIANN_NOTATION, AAunimodToModMassLimited,
                        foundUnimods, diannAAMods, EncodableMods.ANY, cterm);
            case "librarytsv":
                return formatMod(peptide, start, end, LIBRARYTSV_NOTATION, AAunimodToModMassAll,
                        foundUnimods, AAunimodToModMassAllKeys, EncodableMods.ANY, cterm);
            default:
                throw new IllegalArgumentException(model + " is not a peptide format of formatPeptideBaseToSpecific");
        }
    }

    /**
     * Writes the mod at [start, end] of a base peptide the way a Koina model's peptide format does: as the UniMod
     * code its mass has, looked up first among preferredMods and then among allowedMods, or left out of the
     * peptide if neither has it. A code the model cannot encode is passed over, as one neither set has.
     *
     * @param start     the index of the mod's "["
     * @param end       the index of the mod's "]"
     * @param encodable the codes the model can be sent (KoinaModelSpec.encodableMods)
     * @param cterm     whether a mod on the last residue may be a C-term mod
     * @return the peptide with the mod written, and the residue (or "[") plus code chosen, or "" if it was left out
     */
    public static String[] formatPeptideBaseToKoina(String peptide, int start, int end, PeptideFormat format,
                                                    ModMassMap massMap, Set<String> preferredMods,
                                                    Set<String> allowedMods, EncodableMods encodable,
                                                    boolean cterm) {
        ModNotation notation = new ModNotation("UNIMOD", "", format.ntermSeparator(), format.ctermSeparator(),
                false, false);
        return formatMod(peptide, start, end, notation, massMap.aaUnimodToMass(), preferredMods, allowedMods,
                encodable, cterm);
    }

    private static String[] formatMod(String peptide, int start, int end, ModNotation notation,
                                      Map<String, Double> modMap, Set<String> foundUnimods,
                                      Set<String> modelAllowedUnimods, EncodableMods encodable, boolean cterm) {
        double reportedMass = Double.parseDouble(peptide.substring(start + 1, end));
        String unimod = PTMhandler.findUnimodForMass(foundUnimods, modMap, reportedMass, peptide,
                start - 1, end, true, encodable, cterm);
        if (unimod.isEmpty()) {
            unimod = PTMhandler.findUnimodForMass(modelAllowedUnimods, modMap, reportedMass, peptide,
                    start - 1, end, false, encodable, cterm);
        }

        //need to check if the AA is allowed to hold this PTM
        if (unimod.isEmpty()) {
            //model won't predict this anyway
            peptide = peptide.substring(0, start) + peptide.substring(end + 1);
        } else if (notation.customModAsMass() && !NumberUtils.isDigits(unimod.substring(1))) {
            //custom mod represented by mass
            peptide = peptide.substring(0, start + 1) + reportedMass + peptide.substring(end);
        } else {
            peptide = peptide.substring(0, start + 1) + notation.unimodLabel() + ":" +
                    unimod.substring(1) + peptide.substring(end);
        }

        //nterm
        if (!unimod.isEmpty() && start == 0) {
            int splitpoint = -1;
            if (peptide.startsWith("[")) {
                splitpoint = peptide.indexOf("]");
            } else if (peptide.startsWith("(")) {
                splitpoint = peptide.indexOf(")");
            }
            peptide = notation.ntermPrefix() + peptide.substring(0, splitpoint + 1) + notation.ntermSuffix() +
                    peptide.substring(splitpoint + 1);
        }

        //cterm
        if (cterm && (unimod.startsWith("[") || unimod.startsWith("("))) {
            peptide = peptide.substring(0, start) + notation.ctermPrefix() + peptide.substring(start);
        }

        if (notation.parentheses()) {
            peptide = peptide.replace("[", "(");
            peptide = peptide.replace("]", ")");
        }

        return new String[]{peptide, unimod}; //unimod is accepted unimod, or ""
    }

    //returns unimod for the reportedMass, or empty string if not found. A code the model cannot encode is not found
    private static String findUnimodForMass(Set<String> allowedMods, Map<String, Double> modMap,
                                            Double reportedMass,
                                            String peptide, int start, int end,
                                            boolean removeMods, EncodableMods encodable, boolean cterm) {
        //start is index of amino acid before, or 0 if nterm mod

        //first, remove mods that are not supported
        Set<String> filteredMods;
        if (removeMods) {
            filteredMods = allowedMods.stream()
                    .filter(modMap::containsKey)
                    .collect(Collectors.toSet());
        } else {
            filteredMods = new HashSet<>(allowedMods);  // Always create a copy
        }

        //TODO: could be made faster, like through a binary search
        for (String unimod : filteredMods) {
            Double PTMmass = modMap.get(unimod);
            if (NumericUtils.massesCloseEnough(PTMmass, reportedMass)) {
                char AA = unimod.charAt(0);
                boolean site;
                boolean ctermSite = false;
                if (start == -1) { //nterm
                    site = AA == '[';
                } else if (cterm && end == peptide.length() - 1) { //cterm
                    site = AA == '[';
                    ctermSite = true;
                } else {
                    site = AA == peptide.charAt(start);
                }
                if (site && encodable.encodes(unimod, ctermSite)) {
                    return unimod;
                }
            }
        }
        return "";
    }

    //works for diann
    public static String formatPeptideSpecificToBase(String peptide, Map<Integer, Double> modmap, String brackets) {
        char startBracket = brackets.charAt(0);
        char endBracket = brackets.charAt(1);

        String newpeptide = peptide.toUpperCase();
        newpeptide = newpeptide.replace("UNIMOD:","");

        //for koina
        newpeptide = newpeptide.replace(endBracket + "-", String.valueOf(endBracket));
        newpeptide = newpeptide.replace("-" + startBracket, String.valueOf(startBracket));

        //convert unimod to mass
        ArrayList<Integer> newStarts = new ArrayList<>();
        ArrayList<Integer> newEnds = new ArrayList<>();
        for (int i = 0; i < newpeptide.length(); i++) {
            if (newpeptide.charAt(i) == startBracket) {
                newStarts.add(i);
            } else if (newpeptide.charAt(i) == endBracket) {
                newEnds.add(i);
            }
        }
        for (int i = newStarts.size() - 1; i > -1; i--) {
            String bracketContent = newpeptide.substring(newStarts.get(i) + 1, newEnds.get(i));
            int unimod;
            try {
                unimod = Integer.parseInt(bracketContent);
            } catch (NumberFormatException e) {
                //bracket content is already a mass (e.g., unresolved unimod from library), keep as-is
                //just fix number of sig figs
                String formmatedContent = String.format("%.4f", Double.parseDouble(bracketContent));
                if (!bracketContent.equals(formmatedContent)) {
                    newpeptide = newpeptide.substring(0, newStarts.get(i) + 1) +
                            formmatedContent + newpeptide.substring(newEnds.get(i));
                }
                continue;
            }
            try {
                String modMass = String.format("%.4f", modmap.get(unimod));
                newpeptide = newpeptide.substring(0, newStarts.get(i) + 1) + modMass + newpeptide.substring(newEnds.get(i));
            } catch (Exception ignored) {
                printError("Did not recognize unimod " + unimod + " in peptide " + peptide + ". Exiting");
            }
        }

        return newpeptide;
    }

    /////////////////////////////////////////////////////PROSIT////////////////////////////////////////////////////
    private static HashMap<String, Double> makePrositToModAAmass() {
        HashMap<String, Double> map = new HashMap<>();
        map.put("Oxidation", oxidationMass);
        map.put("Carbamidomethyl", carbamidomethylationMass);
        map.put("TMT_6", tmt6_10_11Mass);
        return map;
    }
    public static final HashMap<String, Double> prositToModAAmass = makePrositToModAAmass();

    /////////////////////////////////////////////////////PDEEP////////////////////////////////////////////////////
    private static HashMap<String, String> makeModAAmassToPdeep() throws IOException {
        HashMap<String, String> map = new HashMap<>();

        //get name and monoisotopic mass
        final InputStream stream = PTMhandler.class.getResourceAsStream("/ptm_resources/modification_alphapeptdeep.tsv");
        final InputStreamReader reader = new InputStreamReader(stream);
        final BufferedReader ptmFile = new BufferedReader(reader);
        String line = ptmFile.readLine(); //header

        while((line = ptmFile.readLine()) != null) {
            String[] lineSplit = line.split("\t");
            String classification = lineSplit[6];
            if (classification.equals("Other") || classification.equals("AA substitution")) { //may need to exclude more
                continue;
            }
            String ptmName = lineSplit[0].split("@")[0];
            double doubleMass = Math.round(Double.parseDouble(lineSplit[1]) * 10000.0) / 10000.0;
            for (double d : new double[]{-0.0001, 0, 0.0001}) { //robust against rounding and small mass differences
                String mass = String.format("%.4f", doubleMass + d);
                map.put(mass, ptmName);
            }
        }
        ptmFile.close();

        return map;
    }
    public static HashMap<String, String> aamassToPdeep = null;

    static {
        try {
            aamassToPdeep = makeModAAmassToPdeep();
        } catch (IOException e) {
            e.printStackTrace();
        }
    }

    private static HashMap<String, String> makePdeepToModAAmass() throws IOException {
        HashMap<String, String> map = new HashMap<>();

        //get name and monoisotopic mass
        final InputStream stream = PTMhandler.class.getResourceAsStream("/ptm_resources/modification_alphapeptdeep.tsv");
        final InputStreamReader reader = new InputStreamReader(stream);
        final BufferedReader ptmFile = new BufferedReader(reader);
        String line = ptmFile.readLine(); //header

        while((line = ptmFile.readLine()) != null) {
            String[] lineSplit = line.split("\t");
            String ptmName = lineSplit[0].split("@")[0];
            String mass = String.format("%.4f", Math.round(Double.parseDouble(lineSplit[1]) * 10000.0) / 10000.0);
            map.put(ptmName, mass);
        }
        ptmFile.close();

        return map;
    }
    public static HashMap<String, String> PdeepToAAmass = null;

    static {
        try {
            PdeepToAAmass = makePdeepToModAAmass();
        } catch (IOException e) {
            e.printStackTrace();
        }
    }

    /////////////////////////////////////////////////////ALPHAPEPTDEEP////////////////////////////////////////////////////
    public static HashSet<String> alphapeptdeepModNames = new HashSet<String>();
    public static HashMap<String, String> writeOutAlphapeptdeepModNames = new HashMap<>();

    public static TreeMap<Double, HashSet<String>> makeModAAmassToAlphapeptdeep() throws IOException {
        TreeMap<Double, HashSet<String>> map = new TreeMap<>();

        ArrayList<String> modPaths = new ArrayList<>();
        modPaths.add("/ptm_resources/modification_alphapeptdeep.tsv");
        if (!Constants.additionalModsFile.isEmpty()) {
            modPaths.add(Constants.additionalModsFile);
        }

        for (String modPath : modPaths) {
            //get name and monoisotopic mass
            InputStream stream = PTMhandler.class.getResourceAsStream(modPath);
            if (Objects.isNull(stream)) {
                stream = Files.newInputStream(Paths.get(modPath));
            }
            final InputStreamReader reader = new InputStreamReader(stream);
            final BufferedReader ptmFile = new BufferedReader(reader);

            String line = ptmFile.readLine(); //header
            //get indices of columns
            String[] lineSplit = line.split("\t", -1);
            List<String> header = Arrays.asList(lineSplit);
            int classificationIdx = header.indexOf("classification");
            int nameIdx = header.indexOf("mod_name");
            int massIdx = Math.max(header.indexOf("unimod_mass"), header.indexOf("mod_mass"));

            if (massIdx == -1) {
                printError("Modification file " + modPath + " is missing 'unimod_mass' or 'mod_mass' column");
                printError("Exiting");
                System.exit(1);
            }
            if (nameIdx == -1) {
                printError("Modification file " + modPath + " is missing 'mod_name' column");
                printError("Exiting");
                System.exit(1);
            }

            while ((line = ptmFile.readLine()) != null) {
                lineSplit = line.split("\t", -1);
                if (classificationIdx != -1) {
                    String classification = lineSplit[classificationIdx];
                    if (classification.equals("Other") || classification.equals("AA substitution")) { //may need to exclude more
                        continue;
                    }
                }
                String nameLocalization =  lineSplit[nameIdx].replace("@n", "@Any_N-term").
                        replace("@c", "@Any_C-term");
                String ptmName = nameLocalization.split("@")[0];
                Double mass = Double.parseDouble(lineSplit[massIdx]);
                HashSet<String> stringMap = map.get(mass);
                if (stringMap == null) {
                    stringMap = new HashSet<>();
                }
                stringMap.add(ptmName);
                map.put(mass, stringMap);
                alphapeptdeepModNames.add(nameLocalization.split("\\^")[0]);
                writeOutAlphapeptdeepModNames.put(nameLocalization.split("\\^")[0], nameLocalization);

            }
            ptmFile.close();
        }

        return map;
    }
    public static TreeMap<Double, HashSet<String>> aamassToAlphapeptdeep = null;

    static {
        try {
            aamassToAlphapeptdeep = makeModAAmassToAlphapeptdeep();
        } catch (IOException e) {
            e.printStackTrace();
        }
    }

    private static HashMap<String, String> makeAlphapeptdeepToModAAmass() throws IOException {
        HashMap<String, String> map = new HashMap<>();

        //get name and monoisotopic mass
        final InputStream stream = PTMhandler.class.getResourceAsStream("/ptm_resources/modification_alphapeptdeep.tsv");
        final InputStreamReader reader = new InputStreamReader(stream);
        final BufferedReader ptmFile = new BufferedReader(reader);
        String line = ptmFile.readLine(); //header

        while((line = ptmFile.readLine()) != null) {
            String[] lineSplit = line.split("\t");
            String ptmName = lineSplit[0].split("@")[0];
            String mass = String.format("%.4f", Math.round(Double.parseDouble(lineSplit[1]) * 10000.0) / 10000.0);
            map.put(ptmName, mass);
        }
        ptmFile.close();

        return map;
    }
    public static HashMap<String, String> AlphapeptdeepToAAmass = null;

    static {
        try {
            AlphapeptdeepToAAmass = makeAlphapeptdeepToModAAmass();
        } catch (IOException e) {
            e.printStackTrace();
        }
    }

    /////////////////////////////////////////////KOINA///////////////////////////////////////////////////////
    public static final HashSet<String> prositAAMods = new HashSet<>(
            Arrays.asList("C4", "M35"));
    public static final HashSet<String> prositTmtAAMods = new HashSet<>(
            Arrays.asList("[737", "K737"));
    public static final HashSet<String> prositCitAAMods = new HashSet<>(
            Arrays.asList("N7", "Q7", "R7"));
    public static final HashSet<String> unispecAAMods = new HashSet<>(
            Arrays.asList("[1", "C4", "Q28", "E27", "M35", "S21", "T21", "Y21", "C26"));
    public static final HashSet<String> diannAAMods = new HashSet<>(
            Arrays.asList("C4", "M35", "[1", "S21", "T21", "Y21", "K121", "T121", "C121", "S121", "[121",
                    "[737", "K737", "S737"));
    public static final HashSet<String> ms2pipAAMods = new HashSet<>(
            Arrays.asList("M35", "C4"));

    //the Prosit lactylation models (KoinaModelSpec.ModMassMap.LAC): lactyl on K and the N-terminus, oxidized M and
    //carbamidomethylated C. Lactyl has the mass of carboxyethyl (UNIMOD:378), which the ALL table has on K
    public static final Map<String, Double> AAunimodToModMassLac = Map.of(
            "K" + lactylUnimod, lactylMass,
            "[" + lactylUnimod, lactylMass,
            "M" + oxidationUnimod, oxidationMass,
            "C" + carbamidomethylationUnimod, carbamidomethylationMass);
    public static final Map<Integer, Double> unimodToModMassLac = Map.of(
            lactylUnimod, lactylMass,
            oxidationUnimod, oxidationMass,
            carbamidomethylationUnimod, carbamidomethylationMass);

    //Chronologer (KoinaModelSpec.ModMassMap.CHRONOLOGER): the pairs of Koina's Chronologer table
    //(models/Chronologer/Chronologer_Preprocess/1/sequence_conversion.py), but succinyl on K, which is left out.
    //C385 is the pyro-carbamidomethyl cysteine as a whole, the token Chronologer has for it. The N-term pyro forms of
    //that table (27, 28 and 385 before the first residue) are left out too: Chronologer reads pyro-Glu and
    //pyro-carbamidomethyl from the residue (Q28, E27, C385) and sets the N-terminus to match, while an N-term code
    //would leave residue 1 unmodified. TMT (UNIMOD:737) has the run's TMT mass, as in the LIMITED table
    public static final HashMap<String, Double> AAunimodToModMassChronologer = new HashMap<>();
    public static final HashMap<Integer, Double> unimodToModMassChronologer = new HashMap<>();
    static {
        updateChronologerTables();
    }

    private static void updateChronologerTables() {
        AAunimodToModMassChronologer.clear();
        AAunimodToModMassChronologer.put("C" + carbamidomethylationUnimod, carbamidomethylationMass);
        AAunimodToModMassChronologer.put("M" + oxidationUnimod, oxidationMass);
        AAunimodToModMassChronologer.put("S" + phosphorylationUnimod, phosphorylationMass);
        AAunimodToModMassChronologer.put("T" + phosphorylationUnimod, phosphorylationMass);
        AAunimodToModMassChronologer.put("Y" + phosphorylationUnimod, phosphorylationMass);
        AAunimodToModMassChronologer.put("K" + acetylationUnimod, acetylationMass);
        AAunimodToModMassChronologer.put("K" + glyglyUnimod, glyglyMass);
        AAunimodToModMassChronologer.put("K" + methylUnimod, methylMass);
        AAunimodToModMassChronologer.put("K" + dimethylUnimod, dimethylMass);
        AAunimodToModMassChronologer.put("K" + trimethylUnimod, trimethylMass);
        AAunimodToModMassChronologer.put("R" + methylUnimod, methylMass);
        AAunimodToModMassChronologer.put("R" + dimethylUnimod, dimethylMass);
        AAunimodToModMassChronologer.put("K" + tmtZeroUnimod, tmtZeroMass);
        AAunimodToModMassChronologer.put("K" + tmt6_10_11Unimod, tmtMass);
        AAunimodToModMassChronologer.put("C" + chronologerPyroCarbamidomethylUnimod, pyroCarbamidomethylMass);
        AAunimodToModMassChronologer.put("E" + pyrogluEUnimod, pyrogluEMass);
        AAunimodToModMassChronologer.put("Q" + pyrogluQUnimod, pyrogluQMass);
        AAunimodToModMassChronologer.put("[" + acetylationUnimod, acetylationMass);
        AAunimodToModMassChronologer.put("[" + tmt6_10_11Unimod, tmtMass);
        AAunimodToModMassChronologer.put("[" + tmtZeroUnimod, tmtZeroMass);

        unimodToModMassChronologer.clear();
        unimodToModMassChronologer.put(carbamidomethylationUnimod, carbamidomethylationMass);
        unimodToModMassChronologer.put(oxidationUnimod, oxidationMass);
        unimodToModMassChronologer.put(phosphorylationUnimod, phosphorylationMass);
        unimodToModMassChronologer.put(acetylationUnimod, acetylationMass);
        unimodToModMassChronologer.put(glyglyUnimod, glyglyMass);
        unimodToModMassChronologer.put(methylUnimod, methylMass);
        unimodToModMassChronologer.put(dimethylUnimod, dimethylMass);
        unimodToModMassChronologer.put(trimethylUnimod, trimethylMass);
        unimodToModMassChronologer.put(tmtZeroUnimod, tmtZeroMass);
        unimodToModMassChronologer.put(tmt6_10_11Unimod, tmtMass);
        unimodToModMassChronologer.put(chronologerPyroCarbamidomethylUnimod, pyroCarbamidomethylMass);
        unimodToModMassChronologer.put(pyrogluEUnimod, pyrogluEMass);
        unimodToModMassChronologer.put(pyrogluQUnimod, pyrogluQMass);
    }
}