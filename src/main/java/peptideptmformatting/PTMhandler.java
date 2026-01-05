/*
 * This file is part of MSBooster.
 *
 * MSBooster is free software: you can redistribute it and/or modify it under
 * the terms of the GNU Lesser General Public License as published by the Free
 * Software Foundation, either version 3 of the License, or (at your
 * option) any later version.
 *
 * MSBooster is distributed in the hope that it will be useful, but WITHOUT
 * ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
 * FITNESS FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public License
 * for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with MSBooster. If not, see <https://www.gnu.org/licenses/>.
 */

package peptideptmformatting;

import allconstants.Constants;
import readers.predictionreaders.LibraryTsvReader;
import umich.ms.fileio.filetypes.unimod.UnimodOboReader;
import utils.NumericUtils;
import utils.Print;

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
    }
    public static double getTmtMass() { return tmtMass; }

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
    private static <T> Map<T, Double> makeUnimodToModMassAll(boolean includeAA) throws IOException {
        ArrayList<String> modPaths = new ArrayList<>();
        modPaths.add("/ptm_resources/alphapept_koina.csv");
        modPaths.add("/ptm_resources/modification_alphapeptdeep.tsv");
        if (!Constants.additionalMods.isEmpty()) {
            modPaths.add(Constants.additionalMods);
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
            String line = ptmFile.readLine(); //header

            while ((line = ptmFile.readLine()) != null) {
                String[] lineSplit = new String[0];
                if (modPath.endsWith(".csv")) {
                    lineSplit = line.split(",", -1);
                } else if (modPath.endsWith(".tsv")) {
                    lineSplit = line.split("\t", -1);
                } else {
                    printError("Invalid modification file: " + modPath);
                    printError("Exiting");
                    System.exit(1);
                }
                String unimod = lineSplit[7];
                Double mass = Double.parseDouble(lineSplit[1]);
                if (includeAA) {
                    if (modPath.endsWith(".csv")) {
                        map.put((T) (lineSplit[0].charAt(0) + unimod), mass);
                    } else {
                        //TODO: "Anywhere" should be expanded
                        String mod = lineSplit[0].split("@")[1];
                        if (mod.contains("term")) {
                            mod = "["; //to match logic for alphapept_koina.csv
                        }
                        map.put((T) (mod + unimod), mass);
                    }
                } else {
                    map.put((T) Integer.valueOf(unimod), mass);
                }
            }
        }
        return map;
    }
    public static final Map<Integer, Double> unimodToModMassAll;

    static {
        try {
            unimodToModMassAll = makeUnimodToModMassAll(false);
        } catch (IOException e) {
            throw new RuntimeException(e);
        }
    }

    public static final Map<String, Double> AAunimodToModMassAll;
    public static final Set<String> AAunimodToModMassAllKeys;

    static {
        try {
            AAunimodToModMassAll = makeUnimodToModMassAll(true);
            AAunimodToModMassAllKeys = AAunimodToModMassAll.keySet();
        } catch (IOException e) {
            throw new RuntimeException(e);
        }
    }

    public static Map<Integer, Double> unimodOboToModMass;

    public static void setUnimodObo() {
        try { //only do this when library tsv or fragpipe speclib predicted library is specified
            unimodOboToModMass = new HashMap<>();
            File f = new File(Constants.unimodObo);
            if (! f.exists()) {
                printError("When uploading prediction in spectral library in tsv or speclib format, " +
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

    //start is index of [
    public static String[] formatPeptideBaseToSpecific(String peptide, int start, int end, String model,
                                                       Set<String> foundUnimods, boolean cterm) {
        //set allowed unimods
        Set<String> modelAllowedUnimods;
        switch(model) {
            case "diann":
                modelAllowedUnimods = diannAAMods;
                break;
            case "unispec":
                modelAllowedUnimods = unispecAAMods;
                break;
            case "prosit":
            case "prosittmt":
                modelAllowedUnimods = prositAAMods;
                break;
            case "ms2pip":
                modelAllowedUnimods = ms2pipAAMods;
                break;
            case "alphapept":
            case "deeplc":
            case "im2deep":
            case "librarytsv":
                modelAllowedUnimods = AAunimodToModMassAllKeys;
                break;
            case "predfull":
                modelAllowedUnimods = predfullAAMods;
                break;
            default:
                modelAllowedUnimods = new HashSet<>();
                break;
        }

        //set unimod string formatting
        String unimodFormat = "";
        switch(model) {
            case "diann":
            case "librarytsv":
                unimodFormat = "UniMod";
                break;
            default: //koina
                unimodFormat = "UNIMOD";
                break;
        }

        //is nterm mod followed by "-"?
        String ntermSuffix = "";
        switch(model) {
            case "diann":
            case "unispec":
            case "deeplc":
                break;
            default:
                ntermSuffix = "-";
                break;
        }
        //is cterm mod preceded by "-"?
        String ctermSuffix = "";
        switch(model) {
            case "alphapept":
                ctermSuffix = "-";
                break;
            default:
                break;
        }

        //all unimod codes allowed
        Map<String, Double> modMap;
        if (model.equals("alphapept") || model.equals("deeplc") ||
                model.equals("im2deep") || model.equals("librarytsv")) {
            modMap = AAunimodToModMassAll;
        } else {
            modMap = AAunimodToModMassLimited;
        }

        double reportedMass = Double.parseDouble(peptide.substring(start + 1, end));
        String unimod = PTMhandler.findUnimodForMass(foundUnimods, modMap, reportedMass, peptide,
                start - 1, end, true, cterm);
        if (unimod.isEmpty()) {
            unimod = PTMhandler.findUnimodForMass(modelAllowedUnimods, modMap, reportedMass, peptide,
                    start - 1, end, false, cterm);
        }

        //need to check if the AA is allowed to hold this PTM
        if (unimod.isEmpty()) {
            //model won't predict this anyway
            peptide = peptide.substring(0, start) + peptide.substring(end + 1);
        } else {
            peptide = peptide.substring(0, start + 1) + unimodFormat + ":" +
                    unimod.substring(1) + peptide.substring(end);
        }

        //nterm
        if (peptide.startsWith("[") && start == 0) {
            int splitpoint = peptide.indexOf("]");
            peptide = peptide.substring(0, splitpoint + 1) + ntermSuffix + peptide.substring(splitpoint + 1);
        }

        //cterm
        if (cterm && unimod.startsWith("[")) {
            peptide = peptide.substring(0, start) + ctermSuffix + peptide.substring(start);
        }

        if (model.equals("librarytsv")) {
            peptide = peptide.replace("[", "(");
            peptide = peptide.replace("]", ")");
        }

        return new String[]{peptide, unimod}; //unimod is accepted unimod, or ""
    }

    //returns unimod for the reportedMass, or empty string if not found
    private static String findUnimodForMass(Set<String> allowedMods, Map<String, Double> modMap,
                                            Double reportedMass,
                                            String peptide, int start, int end,
                                            boolean removeMods, boolean cterm) {
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
                if (start == -1) { //nterm
                    if (AA == '[') {
                        return unimod;
                    }
                } else if (cterm && end == peptide.length() - 1) { //cterm
                    if (AA == '[') {
                        return unimod;
                    }
                } else {
                    if (AA == peptide.charAt(start)) {
                        return unimod;
                    }
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
            int unimod = Integer.parseInt(newpeptide.substring(newStarts.get(i) + 1, newEnds.get(i)));
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

    private static TreeMap<Double, HashSet<String>> makeModAAmassToAlphapeptdeep() throws IOException {
        TreeMap<Double, HashSet<String>> map = new TreeMap<>();

        ArrayList<String> modPaths = new ArrayList<>();
        modPaths.add("/ptm_resources/modification_alphapeptdeep.tsv");
        if (!Constants.additionalMods.isEmpty()) {
            modPaths.add(Constants.additionalMods);
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

            while ((line = ptmFile.readLine()) != null) {
                String[] lineSplit = line.split("\t", -1);
                String classification = lineSplit[6];
                if (classification.equals("Other") || classification.equals("AA substitution")) { //may need to exclude more
                    continue;
                }
                String ptmName = lineSplit[0].split("@")[0];
                Double mass = Double.parseDouble(lineSplit[1]);
                HashSet<String> stringMap = map.get(mass);
                if (stringMap == null) {
                    stringMap = new HashSet<>();
                }
                stringMap.add(ptmName);
                map.put(mass, stringMap);
                alphapeptdeepModNames.add(lineSplit[0].split("\\^")[0]);
                writeOutAlphapeptdeepModNames.put(lineSplit[0].split("\\^")[0], lineSplit[0]);

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
    public static final HashSet<String> predfullAAMods = new HashSet<>(
            Arrays.asList("M35", "C4"));
}