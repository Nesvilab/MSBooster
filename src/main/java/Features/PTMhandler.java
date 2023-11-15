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

public class PTMhandler {
    //handling of PTMs, all in one location. Might need to import unimod mass file
    public static final Double oxidationMass = 15.9949;
    public static final Double carbamidomethylationMass = 57.0215;
    public static final Double acetylationMass = 42.0106;
    public static final Double phosphorylationMass = 79.9663;
    public static final Double glyglyMass = 114.042927;
    public static final Double tmtMass = 229.1629;

    /////////////////////////////////////////////////////DIA-NN////////////////////////////////////////////////////
    private static HashMap<Double, Integer> makeModAAToUnimod() {
        HashMap<Double, Integer> map = new HashMap<>();
        map.put(carbamidomethylationMass, 4);
        map.put(oxidationMass, 35);
        map.put(acetylationMass, 1);
        map.put(phosphorylationMass, 21);
        map.put(glyglyMass, 121);
        map.put(tmtMass, 737);
        return map;
    }
    public static final HashMap<Double, Integer> modAAmassToUnimod = makeModAAToUnimod();
    private static HashMap<String, Double> makeUnimodtoModAA() {
        HashMap<String, Double> map = new HashMap<>();
        for (Map.Entry<Double, Integer> entry : modAAmassToUnimod.entrySet()) {
            map.put(String.valueOf(entry.getValue()), entry.getKey());
        }
        return map;
    }
    public static final HashMap<String, Double> unimodtoModAAmass = makeUnimodtoModAA();

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
        //TODO testing
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
}
