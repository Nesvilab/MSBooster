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

import java.io.*;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.*;

import static utils.Print.printError;

public class PTMhandler {
    //handling of PTMs, all in one location.
    public static final Double oxidationMass = 15.9949;
    public static final Double carbamidomethylationMass = 57.0215;
    public static final Double acetylationMass = 42.0106;
    public static final Double phosphorylationMass = 79.9663;
    public static final Double glyglyMass = 114.042927;
    public static final Double tmtMass = 229.1629;
    public static final Double pyrogluQMass = -17.0265;
    public static final Double pyrogluEMass = -18.01057;
    public static final Double pyroCarbamidomethylMass = 39.9950;
    public static final String carbamidomethylationUnimod = "4";
    public static final String oxidationUnimod = "35";
    public static final String acetylationUnimod = "1";
    public static final String phosphorylationUnimod = "21";
    public static final String glyglyUnimod = "121";
    public static final String tmtUnimod = "737";
    public static final String pyrogluQUnimod = "28";
    public static final String pyrogluEUnimod = "27";
    public static final String pyroCarbamidomethylUnimod = "26";

    private static HashMap<String, Double> makeUnimodToModMass() {
        HashMap<String, Double> map = new HashMap<>();
        map.put(carbamidomethylationUnimod, carbamidomethylationMass);
        map.put(oxidationUnimod, oxidationMass);
        map.put(acetylationUnimod, acetylationMass);
        map.put(phosphorylationUnimod, phosphorylationMass);
        map.put(glyglyUnimod, glyglyMass);
        map.put(tmtUnimod, tmtMass);
        map.put(pyrogluQUnimod, pyrogluQMass);
        map.put(pyrogluEUnimod, pyrogluEMass);
        map.put(pyroCarbamidomethylUnimod, pyroCarbamidomethylMass);
        return map;
    }
    public static final HashMap<String, Double> unimodToModMass = makeUnimodToModMass();

    private static HashMap<String, Double> makeUnimodToModMassAlphaPeptDeep(boolean includeAA) throws IOException {
        ArrayList<String> modPaths = new ArrayList<>();
        modPaths.add("/ptm_resources/alphapept_koina.csv");
        if (!Constants.additionalMods.isEmpty()) {
            modPaths.add(Constants.additionalMods);
        }

        HashMap<String, Double> map = new HashMap<>();
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
                String[] lineSplit = line.split(",", -1);
                String unimod = lineSplit[7];
                Double mass = Double.parseDouble(lineSplit[1]);
                if (includeAA) {
                    String AA = lineSplit[0].substring(0, 1);
                    unimod = AA + unimod;
                }
                map.put(unimod, mass);
            }
        }
        return map;
    }
    public static final HashMap<String, Double> unimodToModMassAlphaPeptDeep;

    static {
        try {
            unimodToModMassAlphaPeptDeep = makeUnimodToModMassAlphaPeptDeep(false);
        } catch (IOException e) {
            throw new RuntimeException(e);
        }
    }

    public static final HashMap<String, Double> AAunimodToModMassAlphaPeptDeep;

    static {
        try {
            AAunimodToModMassAlphaPeptDeep = makeUnimodToModMassAlphaPeptDeep(true);
        } catch (IOException e) {
            throw new RuntimeException(e);
        }
    }

    //takes start and end coordinates of mod in peptide formatter base field
    //takes model format
    //returns newly formatted peptide string in format of that model
    //if PTM would cause an error or result in the same fragment intensities (PTM not supported), just removes it from string
    //gives some leeway in PTM mass, precision, only needs to be 0.001 away

    //start is index of [
    public static String[] formatPeptideBaseToSpecific(String peptide, int start, int end, String model,
            HashSet<String> foundUnimods) {
        //set allowed unimods
        HashSet<String> modelAllowedUnimods = new HashSet<>();
        switch(model) {
            case "diann":
                modelAllowedUnimods = diannMods;
                break;
            case "unispec":
                modelAllowedUnimods = unispecMods;
                break;
            case "prosit":
            case "prosittmt":
                modelAllowedUnimods = prositMods;
                break;
            case "ms2pip":
                modelAllowedUnimods = ms2pipMods;
                break;
            case "deeplc":
                modelAllowedUnimods = deeplcMods;
                break;
            case "alphapept":
                modelAllowedUnimods.addAll(AAunimodToModMassAlphaPeptDeep.keySet());
                break;
        }

        //set unimod string formatting
        String unimodFormat = "";
        switch(model) {
            case "diann":
                unimodFormat = "UniMod";
                break;
            default: //koina
                unimodFormat = "UNIMOD";
                break;
        }

        //are is nterm mod followed by "-"?
        String ntermSuffix = "";
        switch(model) {
            case "diann":
            case "unispec":
                break;
            default:
                ntermSuffix = "-";
                break;
        }

        //all unimod codes allowed
        HashMap<String, Double> modMap;
        if (model.equals("alphapept")) {
            modMap = AAunimodToModMassAlphaPeptDeep;
        } else {
            modMap = unimodToModMass;
        }

        double reportedMass = Double.parseDouble(peptide.substring(start + 1, end));
        String unimod = PTMhandler.findUnimodForMass(foundUnimods, modMap, reportedMass, model, peptide,
                Math.max(0, start - 1), true);
        if (unimod.isEmpty()) {
            unimod = PTMhandler.findUnimodForMass(modelAllowedUnimods, modMap, reportedMass, model, peptide,
                    Math.max(0, start - 1), false);
        }
        if (unimod.isEmpty()) {
            //model won't predict this anyway
            peptide = peptide.substring(0, start) + peptide.substring(end + 1);
        } else {
            if (model.equals("alphapept")) {
                peptide = peptide.substring(0, start + 1) + unimodFormat + ":" +
                        unimod.substring(1) + peptide.substring(end);
            } else {
                peptide = peptide.substring(0, start + 1) + unimodFormat + ":" +
                        unimod + peptide.substring(end);
            }
        }

        if (peptide.startsWith("[") && start == 0) {
            int splitpoint = peptide.indexOf("]");
            peptide = peptide.substring(0, splitpoint + 1) + ntermSuffix + peptide.substring(splitpoint + 1);
        }

        return new String[]{peptide, unimod}; //unimod is accepted unimod, or ""
    }

    //returns unimod for the reportedMass, or empty string if not found
    private static String findUnimodForMass(HashSet<String> allowedMods, HashMap<String, Double> modMap,
                                            Double reportedMass, String model,
                                            String peptide, int start, boolean removeMods) {
        //start is index of amino acid before, or 0 if nterm mod

        //first, remove mods that are not supported
        if (removeMods) {
            ArrayList<String> missingMods = new ArrayList<>();
            for (String unimod : allowedMods) {
                if (!modMap.containsKey(unimod)) {
                    missingMods.add(unimod);
                }
            }
            for (String unimod : missingMods) {
                allowedMods.remove(unimod);
            }
        }

        for (String unimod : allowedMods) {
            Double PTMmass = modMap.get(unimod);
            if (Math.abs(PTMmass - reportedMass) < 0.001) {
                if (model.equals("alphapept")) { //requires right localization
                    String AA = unimod.substring(0, 1);
                    if ((AA.equals("[")) ||
                            (AA.equals(peptide.substring(start, start + 1)))) {
                        return unimod;
                    }
                } else {
                    return unimod;
                }
            }
        }
        return "";
    }

    //works for diann
    public static String formatPeptideSpecificToBase(String peptide, HashMap<String, Double> modmap) {
        String newpeptide = peptide.toUpperCase();
        newpeptide = newpeptide.replace("UNIMOD:","");

        //for koina
        newpeptide = newpeptide.replace("]-", "]");
        newpeptide = newpeptide.replace("-[", "[");

        //convert unimod to mass
        ArrayList<Integer> newStarts = new ArrayList<>();
        ArrayList<Integer> newEnds = new ArrayList<>();
        for (int i = 0; i < newpeptide.length(); i++) {
            if (newpeptide.charAt(i) == '[') {
                newStarts.add(i);
            } else if (newpeptide.charAt(i) == ']') {
                newEnds.add(i);
            }
        }
        for (int i = newStarts.size() - 1; i > -1; i--) {
            String unimod = newpeptide.substring(newStarts.get(i) + 1, newEnds.get(i));
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
        map.put("TMT_6", tmtMass);
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
    public static HashMap<String, ArrayList<String>> sameMass = new HashMap<>();

    private static HashMap<String, String> makeModAAmassToAlphapeptdeep() throws IOException {
        HashMap<String, String> map = new HashMap<>();

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
                String mass = String.format("%.4f", Math.round(Double.parseDouble(lineSplit[1]) * 10000.0) / 10000.0);
                map.put(mass, ptmName);
                alphapeptdeepModNames.add(lineSplit[0].split("\\^")[0]);
                writeOutAlphapeptdeepModNames.put(lineSplit[0].split("\\^")[0], lineSplit[0]);
                ArrayList<String> stringList = sameMass.get(mass);
                if (stringList == null) {
                    stringList = new ArrayList<>();
                }
                stringList.add(ptmName);
                sameMass.put(mass, stringList);
            }
            ptmFile.close();
        }

        return map;
    }
    public static HashMap<String, String> aamassToAlphapeptdeep = null;

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
    public static final HashSet<String> prositMods = new HashSet<>(Arrays.asList("737", "4", "35"));
    public static final HashSet<String> unispecMods = new HashSet<>(Arrays.asList("1", "4", "28", "27", "35", "21", "26"));
    public static final HashSet<String> diannMods = new HashSet<>(Arrays.asList("4", "35", "1", "21", "121", "737"));
    public static final HashSet<String> ms2pipMods = new HashSet<>(Arrays.asList("35", "4"));
    public static final HashSet<String> deeplcMods = new HashSet<>(Arrays.asList("35", "4", "21"));
}
