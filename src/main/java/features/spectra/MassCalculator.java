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

package features.spectra;

import allconstants.Constants;

import static utils.Print.printError;

import java.util.*;

public class MassCalculator {
    public final float proton = 1.00727647f;
    private final float H = 1.0078250321f;
    private final float O = 15.9949146221f;
    private final float C = 12.000000f;
    private final float N = 14.0030740052f;
    private final float H2O = H * 2 + O;
    private final float NH3 = H * 3 + N;
    private final String[] neutralLosses = new String[] {"H2O", "NH3"};
    private final String[] internalFragmentTypes = new String[] {"ay", "by"};

    public String fullPeptide;
    public String peptide;
    public int charge;
    public float mass;
    public ArrayList<Double> modMasses = new ArrayList<Double>();
    public MassCalculator[] internalPeptides;
    public SortedMap<Float, String[]> fragmentIons = new TreeMap<>();
    public static final Set<String> allowedFragmentIonTypes = new HashSet<>(Arrays.asList(
            "z", "c", "y", "a", "x", "b", "zdot", "cdot",
            "precursor", "immonium", "internal", "internal-NL",
            "z-NL", "c-NL", "y-NL", "a-NL", "x-NL", "b-NL",
            "precursor-NL", "unknown"));

    //TODO: get immonium ion masses by taking amino acid and subtracting 26.99 Da. This holds for modified AA too (Falick et al 1993)
    //TODO: should we consider related ions, not just immonium?
    public static HashMap<Character, Float> AAmap = new HashMap<Character, Float>()
    {{
        put('A', 71.03711f);
        put('R', 156.10111f);
        put('N', 114.04293f);
        put('D', 115.02694f);
        put('C', 103.00919f);
        put('E', 129.04259f);
        put('Q', 128.05858f);
        put('G', 57.02146f);
        put('H', 137.05891f);
        put('I', 113.08406f);
        put('L', 113.08406f);
        put('K', 128.09496f);
        put('M', 131.04049f);
        put('F', 147.06841f);
        put('P', 97.05276f);
        put('S', 87.03203f);
        put('T', 101.04768f);
        put('W', 186.07931f);
        put('Y', 163.06333f);
        put('V', 99.06841f);
        put('O', 237.14773f);
        put('U', 150.95363f);
        put('Z', 0f);
        put('B', 0f);
        put('X', 0f);
        put('J', 0f);
    }};

    private final HashMap<String, double[]> series = new HashMap<>();
    Set<String> leftIons = new HashSet<>(Arrays.asList("a", "b", "c", "cdot"));
    Set<String> rightIons = new HashSet<>(Arrays.asList("x", "y", "z", "zdot"));

    private static HashMap<Integer, String> makeFlagTOion() {
        HashMap<Integer, String> map = new HashMap<>();
        map.put(0, "b");
        map.put(1, "y");
        map.put(2, "b-NL");
        map.put(3, "y-NL");
        map.put(4, "a");
        map.put(5, "x");
        map.put(6, "a-NL");
        map.put(7, "x-NL");
        map.put(8, "c");
        map.put(9, "z");
        map.put(10, "c-NL");
        map.put(11, "z-NL");
        map.put(12, "imm"); //immonium
        map.put(13, "int"); //internal
        map.put(14, "p"); //precursor
        return map;
    }
    public static HashMap<Integer, String> flagTOion = makeFlagTOion();

    private static HashMap<String, Integer> makeIonToFlag() {
        HashMap<String, Integer> map = new HashMap<>();
        for (Map.Entry<Integer, String> entry : flagTOion.entrySet()) {
            map.put(entry.getValue(), entry.getKey());
        }
        return map;
    }
    public static HashMap<String, Integer> ionTOflag = makeIonToFlag();

    public MassCalculator(String pep, Object charge) {
        fullPeptide = pep + "|" + charge;

        while (true) {
            int ind = pep.indexOf("[");
            if (ind == -1) {
                break;
            }
            int ind2 = pep.indexOf("]");

            String modNum = pep.substring(ind + 1, ind2);
            //add modNum to modMasses
            while (modMasses.size() < ind) {
                modMasses.add(0d);
            }
            modMasses.add(Double.parseDouble(modNum));
//            }

            //replace mod
            pep = pep.substring(0, ind) + pep.substring(ind2 + 1);
        }

        peptide = pep;

        //fill in remaining zeros, accounting for c term mod
        while (modMasses.size() < peptide.length() + 2) {
            modMasses.add(0d);
        }

        //set charge
        if (charge instanceof String) {
            this.charge = Integer.parseInt((String) charge);
        } else {
            this.charge = (int) charge;
        }

        this.mass = calcMass(pep.length(), "y", 1) - H;
        this.internalPeptides = new MassCalculator[peptide.length() - 3];

        //TODO something smarter, like based off fragmentation mode? Or optional parameter to pass flags to initialize.
        // Or fragger.params. Or based on model (DIANN always uses yb)
        String[] flags = {"y", "b"};
        for (String flag : flags) {
            initializeSeries(flag);
        }
    }

    public void initializeSeries(String s) {
        double[] mySeries = new double[peptide.length() - 1];
        mySeries[0] = calcMass(1, s, 1) - proton;

        if (leftIons.contains(s)) {
            for (int i = 1; i < peptide.length() - 1; i++) {
                mySeries[i] = (mySeries[i - 1] + AAmap.get(this.peptide.charAt(i)) + modMasses.get(i + 1));
            }
        } else {
            for (int i = peptide.length() - 3; i > -1; i--) { //TODO this is so f'ing convoluted
                mySeries[peptide.length() - i - 2] = (mySeries[peptide.length() - i - 3] + AAmap.get(this.peptide.charAt(i + 1)) + modMasses.get(i + 2));
            }
        }
        series.put(s, mySeries);
    }

    public float calcMass(int num, String flag, int charge) {
        double mass = 0f;
        if (series.containsKey(flag)) {
            mass = series.get(flag)[num - 1];
        } else {
            switch (flag) {
                case "b":
                    mass += modMasses.get(0);
                    for (int i = 0; i < num; i++) {
                        mass += AAmap.get(this.peptide.charAt(i)); //sum amino acid masses
                        mass += modMasses.get(i + 1); //sum mods
                    }
                    if (num == this.peptide.length()) { //add c term mod
                        mass += modMasses.get(num + 1);
                    }
                    break;
                case "y":
                    mass = H2O;
                    mass += modMasses.get(this.peptide.length() + 1);
                    for (int i = this.peptide.length() - num; i < this.peptide.length(); i++) {
                        mass += AAmap.get(this.peptide.charAt(i)); //sum amino acid masses
                        mass += modMasses.get(i + 1); //sum mods
                    }
                    if (num == this.peptide.length()) { //add n term mod
                        mass += modMasses.get(0);
                    }
                    break;
                case "a":
                    mass -= C;
                    mass -= O;
                    mass += modMasses.get(0);
                    for (int i = 0; i < num; i++) {
                        mass += AAmap.get(this.peptide.charAt(i)); //sum amino acid masses
                        mass += modMasses.get(i + 1); //sum mods
                    }
                    if (num == this.peptide.length()) { //add c term mod
                        mass += modMasses.get(num + 1);
                    }
                    break;
                case "x":
                    mass = 2 * O + C; //+H2O - H2 + CO
                    mass += modMasses.get(this.peptide.length() + 1);
                    for (int i = this.peptide.length() - num; i < this.peptide.length(); i++) {
                        mass += AAmap.get(this.peptide.charAt(i)); //sum amino acid masses
                        mass += modMasses.get(i + 1); //sum mods
                    }
                    if (num == this.peptide.length()) { //add n term mod
                        mass += modMasses.get(0);
                    }
                    break;
                case "c":
                    mass = NH3;
                    mass += modMasses.get(0);
                    for (int i = 0; i < num; i++) {
                        mass += AAmap.get(this.peptide.charAt(i)); //sum amino acid masses
                        mass += modMasses.get(i + 1); //sum mods
                    }
                    if (num == this.peptide.length()) { //add c term mod
                        mass += modMasses.get(num + 1);
                    }
                    break;
                case "z":
                    mass = O - N;
                    mass += modMasses.get(this.peptide.length() + 1);
                    for (int i = this.peptide.length() - num; i < this.peptide.length(); i++) {
                        mass += AAmap.get(this.peptide.charAt(i)); //sum amino acid masses
                        mass += modMasses.get(i + 1); //sum mods
                    }
                    if (num == this.peptide.length()) { //add n term mod
                        mass += modMasses.get(0);
                    }
                    break;
                case "zdot":
                    mass = calcMass(num, "z", 1) - proton - H;
                    break;
                case "cdot":
                    mass = calcMass(num, "c", 1) - proton - H;
                    break;
            }
        }

        //calculate m/z using charge
        return (float) ((mass + charge * proton) / charge);
    }

    //num1 is the y ion the internal fragment is derived from
    //num2 is the a/b ion that is formed from that
    //assume charge 1 only
    public float calcMass(int num1, int num2, String flag) {
        MassCalculator newY = internalPeptides[num1 - 3];
        if (newY == null) {
            newY = makeInternalPeptide(num1);
        }

        switch (flag) {
            case "by":
                return newY.calcMass(num2, "b", 1);
            case "ay":
                return newY.calcMass(num2, "a", 1);
            default:
                printError(flag + " is not supported. Choose ay or by");
                System.exit(-1);
        }
        return 0;
    }

    public float calcMass(int num, String flag, int charge, String neutralLoss) {
        float mass = calcMass(num, flag, 1) - proton;

        //subtract neutral loss
        switch (neutralLoss) {
            case "NH3":
                mass -= NH3;
                break;
            case "H2O":
                mass -= H2O;
                break;
        }

        //calculate m/z using charge
        return (mass + (float) charge * proton) / (float) charge;
    }

    public float calcMass(int num1, int num2, String flag, String neutralLoss) {
        float mass = calcMass(num1, num2, flag);

        //subtract neutral loss
        switch (neutralLoss) {
            case "NH3":
                mass -= NH3;
                break;
            case "H2O":
                mass -= H2O;
                break;
        }

        //calculate m/z using charge
        return mass;
    }

    public MassCalculator() {}
    //y is which y fragment
    private MassCalculator makeInternalPeptide(int y) {
        if (y >= this.peptide.length()) { //only make smaller peptides
            printError("internal fragment must be shorter than initial peptide");
            System.exit(-1);
        }
        String peptide = this.peptide.substring(this.peptide.length() - y);
        //float mass = this.mass;
        ArrayList<Double> modMasses = new ArrayList<Double>();
        modMasses.add(0d);

//        for (int i = 0; i < this.peptide.length() - y + 1; i++) {
//            mass -= modMasses.get(i);
//        }
        for (int i = this.peptide.length() - y + 1; i < this.peptide.length() + 2; i++) {
            modMasses.add(this.modMasses.get(i));
        }

        //don't care about charge for now, since only support charge 1
        MassCalculator mc = new MassCalculator();
        mc.peptide = peptide;
        mc.modMasses = modMasses;
        this.internalPeptides[y - 3] = mc;
        return mc;
    }

    //in case masses overlap
    private void addToFragmentIons(Float mass, String[] info) {
        fragmentIons.merge(mass, info, (a, b) -> new String[]{a[0] + ";" + b[0],
                a[1] + ";" + b[1]});
    }

    //currently up to charge 2
    //format of value is [name, fragment ion type]
    //only support 1 neutral loss
    public void possibleFragmentIons() { //TODO: edit based on allowedFragmentIonTypes? Allow everything above lowest rank allowed
        //calculate precursor isotopic peaks
        for (int iCharge = 1; iCharge < charge + 1; iCharge++) {
            //regular
            if (Constants.lowestFragmentIonType.contains("precursor")) {
                String ionName = "MH" + "+" + iCharge;
                addToFragmentIons(calcMass(peptide.length(), "y", iCharge), new String[]{ionName, "precursor"});
            }

            //neutral loss
            if (Constants.lowestFragmentIonType.contains("precursor-NL")) {
                for (String nl : neutralLosses) {
                    String ionName = "MH-" + nl + "+" + iCharge;
                    addToFragmentIons(calcMass(peptide.length(), "y", iCharge, nl), new String[]{ionName, "precursor-NL"});
                }
            }
        }

        //calculate all abcxyz ions
        String[] ions = new String[]{"a", "b", "c", "x", "y", "z", "zdot", "cdot"};
        int maxCharge = charge; //max fragment ion charge?
        for (int i = 0; i < ions.length; i++) {
            String ionType = ions[i];
            for (int num = 1; num < this.peptide.length(); num++) {
                for (int iCharge = 1; iCharge < maxCharge + 1; iCharge++) {
                    //regular fragment
                    if (Constants.lowestFragmentIonType.contains(ionType)) {
                        String ionName = ionType + num + "+" + iCharge;
                        addToFragmentIons(calcMass(num, ionType, iCharge), new String[]{ionName, ionType});
                    }

                    //neutral loss fragments
                    if (Constants.lowestFragmentIonType.contains(ionType + "-NL")) {
                        if (!ionType.equals("c")) {
                            for (String nl : neutralLosses) {
                                String ionName = ionType + num + "-" + nl + "+" + iCharge;
                                addToFragmentIons(calcMass(num, ionType, iCharge, nl), new String[]{ionName, ionType + "-NL"});
                            }
                        } else { //c-NH3 is same as b
                            String ionName = ionType + num + "-H2O" + "+" + iCharge;
                            addToFragmentIons(calcMass(num, ionType, iCharge, "H2O"), new String[]{ionName, "c-NL"});
                        }
                    }
                }
            }
        }

        //calculate internal fragment ions ay and by
        //what is the proper notation for these?
        for (String ionType : internalFragmentTypes) {
            for (int num1 = 3; num1 < peptide.length(); num1++) {
                for (int num2 = 2; num2 < peptide.length() - 1; num2++) {
                    if (num2 >= num1) {
                        break;
                    }

                    if (Constants.lowestFragmentIonType.contains("internal")) {
                        String ionName = "y" + num1 + ionType.charAt(0) + num2;
                        addToFragmentIons(calcMass(num1, num2, ionType), new String[]{ionName, "internal"});
                    }
                    if (Constants.lowestFragmentIonType.contains("internal-NL")) {
                        for (String nl : neutralLosses) {
                            String ionName = "y" + num1 + ionType.charAt(0) + num2 + "-" + nl;
                            addToFragmentIons(calcMass(num1, num2, ionType, nl), new String[]{ionName, "internal-NL"});
                        }
                    }
                }
            }
        }

        //don't need internal peptides anymore
        //internalPeptides = null;

        //calculate immonium ions
        if (Constants.lowestFragmentIonType.contains("immonium")) {
            for (int i = 0; i < peptide.length(); i++) {
                String ionName;

                //get attached mod
                if (modMasses.get(i + 1) != 0) {
                    ionName = "immonium:" + peptide.charAt(i) + "[" + modMasses.get(i + 1) + "]";
                } else {
                    ionName = "immonium:" + peptide.charAt(i);
                }

                addToFragmentIons((float) (modMasses.get(i + 1) + AAmap.get(peptide.charAt(i)) - 26.99),
                        new String[]{ionName, "immonium"});
            }
        }
    }

    public String[][] annotateMZs(float[] fs) {
        if (fragmentIons.size() == 0) {
            possibleFragmentIons();
        }

        String[] annotations = new String[fs.length];
        String[] fragmentIonTypes = new String[fs.length];

        //convert to arrays to be faster
        float[] mzArray = new float[fragmentIons.size()];
        int k = 0;
        for (float f : fragmentIons.keySet()) {
            mzArray[k] = f;
            k += 1;
        }
        String[][] ionsArray = new String[fragmentIons.size()][];
        k = 0;
        for (String[] s : fragmentIons.values()) {
            ionsArray[k] = s;
            k += 1;
        }

        int q = 0;
        for (int i = 0; i < fs.length; i++) {
            //for predicted peak, +- Da error tolerance
            float minMZ = fs[i] - Constants.DaTolerance;
            float maxMZ = fs[i] + Constants.DaTolerance;

            ArrayList<String[]> consideredFragmentIons = new ArrayList<>();

            for (int p = q; p < mzArray.length; p++) {
                float mz = mzArray[p];

                if (mz >= minMZ) {
                    if (mz <= maxMZ) {
                        String[] ions = ionsArray[p];
                        String[] IDs = ions[0].split(";");
                        String[] ionTypes = ions[1].split(";");

                        for (int j = 0; j < IDs.length; j++) {
                            consideredFragmentIons.add(new String[]{IDs[j], ionTypes[j]});
                        }
                    } else {
                        break;
                    }
                } else { //too low to be matched now
                    q += 1;
                }
            }

            //not able to find match
            if (consideredFragmentIons.size() == 0) {
                annotations[i] = "unknown";
                fragmentIonTypes[i] = "unknown";
                continue;
            }

            //annotate with highest priority fragment ion type
            StringBuilder sb = new StringBuilder();
            for (String ionType : Constants.fragmentIonHierarchy) {
                for (String[] candidates : consideredFragmentIons) {
                    if (candidates[1].equals(ionType)) {
                        if (!sb.toString().equals("")) { //multiple matches
                            sb.append(";");
                        }
                        sb.append(candidates[0]);
                    }
                }

                if (!sb.toString().equals("")) {
                    annotations[i] = sb.toString();
                    fragmentIonTypes[i] = ionType;
                    break;
                }
            }
        }

        return new String[][] {annotations, fragmentIonTypes};
    }
}
