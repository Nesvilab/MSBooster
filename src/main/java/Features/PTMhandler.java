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

package Features;

import java.io.*;
import java.util.*;

import static utils.Print.printError;

public class PTMhandler {
    //handling of PTMs, all in one location.
    //TODO: might need to import unimod mass file, for alphapeptdeep transfer learning
    public static final Double oxidationMass = 15.9949;
    public static final Double carbamidomethylationMass = 57.0215;
    public static final Double acetylationMass = 42.0106;
    public static final Double phosphorylationMass = 79.9663;
    public static final Double glyglyMass = 114.042927;
    public static final Double tmtMass = 229.1629;
    public static final Double pyrogluQ = -17.0265;
    public static final Double pyrogluE = -18.01057;
    public static final Double pyroCarbamidomethyl = 39.9950;

    private static HashMap<Double, String> makeModMassToUnimod() {
        HashMap<Double, String> map = new HashMap<>();
        map.put(carbamidomethylationMass, "4");
        map.put(oxidationMass, "35");
        map.put(acetylationMass, "1");
        map.put(phosphorylationMass, "21");
        map.put(glyglyMass, "121");
        map.put(tmtMass, "737");
        map.put(pyrogluQ, "28");
        map.put(pyrogluE, "27");
        map.put(pyroCarbamidomethyl, "26");
        return map;
    }
    public static final HashMap<Double, String> modMassToUnimod = makeModMassToUnimod();
    private static HashMap<String, Double> makeUnimodToModMass() {
        HashMap<String, Double> map = new HashMap<>();
        for (Map.Entry<Double, String> entry : modMassToUnimod.entrySet()) {
            map.put(entry.getValue(), entry.getKey( ));
        }
        return map;
    }
    public static final HashMap<String, Double> unimodToModMass = makeUnimodToModMass();

    //takes start and end coordinates of mod in peptide formatter base field
    //takes model format
    //returns newly formatted peptide string in format of that model
    //if PTM would cause an error or result in the same fragment intensities (PTM not supported), just removes it from string
    //gives some leeway in PTM mass, precision, only needs to be 0.001 away
    public static String formatPeptideBaseToSpecific(String peptide, int start, int end, String model) {
        double reportedMass = Double.parseDouble(peptide.substring(start + 1, end));
        boolean foundReplacement = false;
        HashSet<String> modelAllowedUnimods = new HashSet<>();

        //set allowed unimods
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
                modelAllowedUnimods = alphapeptMods;
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

        for (String unimod : modelAllowedUnimods) {
            Double PTMmass = unimodToModMass.get(unimod);

            if (Math.abs(PTMmass - reportedMass) < 0.001) {
                peptide = peptide.substring(0, start + 1) + unimodFormat + ":" + unimod + peptide.substring(end);
                foundReplacement = true;
                break;
            }
        }
        if (! foundReplacement) {
            //model won't predict this anyway
            peptide = peptide.substring(0, start) + peptide.substring(end + 1);
        }

        if (peptide.startsWith("[") && start == 0) {
            int splitpoint = peptide.indexOf("]");
            peptide = peptide.substring(0, splitpoint + 1) + ntermSuffix + peptide.substring(splitpoint + 1);
        }

        return peptide;
    }

    //works for diann
    public static String formatPeptideSpecificToBase(String peptide) {
        String newpeptide = peptide.toUpperCase();
        newpeptide = newpeptide.replace("UNIMOD:","");

        //for koina
        //TODO: figure out when koina format is used
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
                String modMass = String.format("%.4f", unimodToModMass.get(unimod));
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
        if (! Constants.additionalMods.equals("")) {
            modPaths.add(Constants.additionalMods);
        }

        for (String modPath : modPaths) {
            //get name and monoisotopic mass
            InputStream stream = PTMhandler.class.getResourceAsStream(modPath);
            if (Objects.isNull(stream)) {
                stream = new FileInputStream(modPath);
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
    //TODO: alphapeptdeep can take all unimod mods, need a specific method for this
    public static final HashSet<String> alphapeptMods = new HashSet<>(Arrays.asList("1", "4", "28", "27", "35", "21", "26"));
}
