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
import allconstants.FragmentIonConstants;
import predictions.FragmentAnnotationParser;

import java.io.BufferedReader;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.net.URISyntaxException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Map;
import java.util.Set;
import java.util.SortedMap;
import java.util.TreeMap;

import static utils.NumericUtils.sum;
import static utils.Print.printError;

public class MassCalculator {
    public final float proton = 1.00727647f;
    public final float isotopic_massdiff = 1.0033548378f;
    private final float H = 1.0078250321f;
    private final float O = 15.9949146221f;
    private final float C = 12.000000f;
    private final float N = 14.0030740052f;
    private final float H2O = H * 2 + O;
    private final float NH3 = H * 3 + N;
    private final String[] selectNeutralLosses = new String[] {"H2O", "NH3"};
    private final String[] internalFragmentTypes = new String[] {"ay", "by"};

    public static final HashMap<String, Float> makeNlMasses() throws IOException {
        HashMap<String, Float> neutralLossMasses = new HashMap<>();
        final InputStream stream = MassCalculator.class.getClassLoader().getResourceAsStream(
                "fragment_annotation/neutral_losses.txt");
        final InputStreamReader reader = new InputStreamReader(stream);
        final BufferedReader nlFile = new BufferedReader(reader);

        String line;
        while((line = nlFile.readLine()) != null) {
            String[] lineSplit = line.split("\t");
            neutralLossMasses.put(lineSplit[0], Float.valueOf(lineSplit[1]));
        }
        nlFile.close();
        return neutralLossMasses;
    }
    public static final HashMap<String, Float> allNeutralLossMasses;

    static {
        try {
            allNeutralLossMasses = makeNlMasses();
        } catch (IOException e) {
            throw new RuntimeException(e);
        }
    }

    public static final HashMap<String, Float> makeImmoniumMasses() throws IOException {
        HashMap<String, Float> immMasses = new HashMap<>();
        final InputStream stream = MassCalculator.class.getClassLoader().getResourceAsStream(
                "fragment_annotation/immonium_masses.txt");
        final InputStreamReader reader = new InputStreamReader(stream);
        final BufferedReader immFile = new BufferedReader(reader);

        String line;
        while((line = immFile.readLine()) != null) {
            String[] lineSplit = line.split(" ");
            immMasses.put(lineSplit[0], Float.valueOf(lineSplit[1]));
        }
        immFile.close();
        return immMasses;
    }
    public static final HashMap<String, Float> immoniumMasses;

    static {
        try {
            immoniumMasses = makeImmoniumMasses();
        } catch (IOException e) {
            throw new RuntimeException(e);
        }
    }

    public String fullPeptide;
    public String peptide;
    public int charge;
    public float mass;
    public ArrayList<Double> modMasses = new ArrayList<Double>();
    public MassCalculator[] internalPeptides;
    private final int internalAaFromNterm = 3;
    public SortedMap<Float, String[]> fragmentIons = new TreeMap<>();
    public SortedMap<Float, String[]> unispecFragmentIons = new TreeMap<>();
    public HashMap<String, Float> annotationMasses = new HashMap<>();
    public HashMap<String, Float> unispecAnnotationMasses = new HashMap<>();
    public static final Set<String> allowedFragmentIonTypes = new HashSet<>(Arrays.asList(
            "z", "c", "y", "a", "x", "b", "zdot", "cdot",
            "p", "imm", "int", "int-NL",
            "z-NL", "c-NL", "y-NL", "a-NL", "x-NL", "b-NL",
            "p-NL", "unknown"));

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

    private static HashMap<Integer, String> makeFlagTOion() { //TODO: revamp c/zdot and add to this hashmap
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
        map.put(14, "int-NL");
        map.put(15, "p"); //precursor
        map.put(16, "p-NL");
        map.put(17, "unknown");
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

    private MassCalculator() {}

    public MassCalculator(String pep, Object charge) {
        if (pep.endsWith("cterm")) {
            pep = pep.substring(0, pep.length() - 5);
        }
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

        this.mass = calcMass(pep.length(), "y", 1, 0) - H;
        this.internalPeptides = new MassCalculator[peptide.length() - internalAaFromNterm]; //for PEPTIDER, only consider EP,EPT,...EPTIDE, which is 5.

        //TODO something smarter, like based off fragmentation mode? Or optional parameter to pass flags to initialize.
        // Or fragger.params. Or based on model (DIANN always uses yb)
        String[] flags = {"y", "b"};
        for (String flag : flags) {
            initializeSeries(flag);
        }
    }

    public void initializeSeries(String s) {
        double[] mySeries = new double[peptide.length() - 1];
        mySeries[0] = calcMass(1, s, 1, 0) - proton;

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

    public float calcMass(int num, String flag, int charge, int isotope) {
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
                    mass = calcMass(num, "z", 1, isotope) - proton - H;
                    break;
                case "cdot":
                    mass = calcMass(num, "c", 1, isotope) - proton - H;
                    break;
            }
        }

        //calculate m/z using charge
        return (float) ((mass + charge * proton + isotope * isotopic_massdiff) / charge);
    }

    //num1 is the y ion the internal fragment is derived from
    //num2 is the a/b ion that is formed from that
    //assume charge 1 only
    public float calcMass(int num1, int num2, String flag, int isotope) {
        MassCalculator newY = internalPeptides[num1 - internalAaFromNterm];
        if (newY == null) {
            newY = makeInternalPeptide(num1);
        }

        switch (flag) {
            case "by":
                return newY.calcMass(num2, "b", 1, isotope);
            case "ay":
                return newY.calcMass(num2, "a", 1, isotope);
            default:
                printError(flag + " is not supported. Choose ay or by");
                System.exit(1);
        }
        return 0;
    }

    public float calcMass(int num, String flag, int charge, float neutralLoss, int isotope) {
        float mass = calcMass(num, flag.split("-")[0], 1, 0) - proton;

        //subtract neutral loss
        mass -= neutralLoss;

        //calculate m/z using charge
        return (mass + (float) charge * proton + isotope * isotopic_massdiff) / (float) charge;
    }

    public float calcMass(int num1, int num2, String flag, float neutralLoss, int isotope) {
        float mass = calcMass(num1, num2, flag.split("-")[0], isotope);

        //subtract neutral loss
        mass -= neutralLoss;

        //calculate m/z using charge
        return mass;
    }

    public float calcMassImmonium(String annotation, float neutralLoss, int isotope) {
        annotation = annotation.split("[+\\-]")[0];
        return immoniumMasses.get(annotation) - neutralLoss + isotope * isotopic_massdiff;
    }

    public float calcMassPrecursor(int isotope, float neutralLoss, int charge) {
        return (mass + isotope * isotopic_massdiff + proton * charge - neutralLoss) / charge;
    }

    private MassCalculator makeInternalPeptide(int y) {
        if (y >= this.peptide.length()) { //only make smaller peptides
            printError("internal fragment must be shorter than initial peptide");
            System.exit(1);
        }
        String peptide = this.peptide.substring(this.peptide.length() - y);

        ArrayList<Double> modMasses = new ArrayList<Double>();
        modMasses.add(0d);

        for (int i = this.peptide.length() - y + 1; i < this.peptide.length() + 2; i++) {
            modMasses.add(this.modMasses.get(i));
        }

        //don't care about charge for now, since only support charge 1
        MassCalculator mc = new MassCalculator();
        mc.peptide = peptide;
        mc.modMasses = modMasses;
        this.internalPeptides[y - internalAaFromNterm] = mc;
        return mc;
    }

    //in case masses overlap
    private void addToFragmentIons(Float mass, String[] info,
                                   SortedMap<Float, String[]> mzToAnnotationMap,
                                   HashMap<String, Float> annotationToMzMap) {
        mzToAnnotationMap.merge(mass, info, (a, b) -> new String[]{a[0] + ";" + b[0],
                a[1] + ";" + b[1]});
        annotationToMzMap.put(info[0], mass);
    }

    //currently up to charge 2
    //format of value is [name, fragment ion type]
    //only support 1 neutral loss
    private void possibleFragmentIons() { //TODO: edit based on allowedFragmentIonTypes? Allow everything above lowest rank allowed
        //calculate precursor isotopic peaks. Will calculate with lower values of charge
        for (int iCharge = 1; iCharge < charge + 1; iCharge++) {
            //regular
            if (FragmentIonConstants.fragmentIonHierarchySet.contains("p")) {
                String ionName = "p^" + iCharge;
                addToFragmentIons(calcMassPrecursor(0, 0, iCharge),
                        new String[]{ionName, "p"}, fragmentIons, annotationMasses);
            }

            //neutral loss
            if (FragmentIonConstants.fragmentIonHierarchySet.contains("p-NL")) {
                for (String nl : selectNeutralLosses) {
                    String ionName = "p-" + nl + "^" + iCharge;
                    addToFragmentIons(calcMassPrecursor(0, allNeutralLossMasses.get(nl), iCharge),
                            new String[]{ionName, "p-NL"}, fragmentIons, annotationMasses);
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
                    if (FragmentIonConstants.fragmentIonHierarchySet.contains(ionType)) {
                        String ionName = ionType + num + "^" + iCharge;
                        addToFragmentIons(calcMass(num, ionType, iCharge, 0), new String[]{ionName, ionType},
                                fragmentIons, annotationMasses);
                    }

                    //neutral loss fragments
                    if (FragmentIonConstants.fragmentIonHierarchySet.contains(ionType + "-NL")) {
                        if (!ionType.equals("c")) {
                            for (String nl : selectNeutralLosses) {
                                String ionName = ionType + num + "-" + nl + "^" + iCharge;
                                addToFragmentIons(calcMass(num, ionType, iCharge, allNeutralLossMasses.get(nl), 0),
                                        new String[]{ionName, ionType + "-NL"},
                                        fragmentIons, annotationMasses);
                            }
                        } else { //c-NH3 is same as b
                            String ionName = ionType + num + "-H2O^" + iCharge;
                            addToFragmentIons(calcMass(num, ionType, iCharge, allNeutralLossMasses.get("H2O"), 0),
                                    new String[]{ionName, "c-NL"},
                                    fragmentIons, annotationMasses);
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

                    //getting subsequence
                    //start position should be peptide.length - num1
                    int start = peptide.length() - num1;
                    String subsequence = peptide.substring(start, start + num2);

                    if (FragmentIonConstants.fragmentIonHierarchySet.contains("int")) {
                        String ionName = "Int:y" + num1 + ionType.charAt(0) + num2 + "/" + subsequence;
                        addToFragmentIons(calcMass(num1, num2, ionType, 0), new String[]{ionName, "int"},
                                fragmentIons, annotationMasses);
                    }
                    if (FragmentIonConstants.fragmentIonHierarchySet.contains("int-NL")) {
                        for (String nl : selectNeutralLosses) {
                            String ionName = "Int:y" + num1 + ionType.charAt(0) + num2 + "-" + nl + "/" + subsequence;
                            addToFragmentIons(calcMass(num1, num2, ionType,
                                    allNeutralLossMasses.get(nl), 0), new String[]{ionName, "int-NL"},
                                    fragmentIons, annotationMasses);
                        }
                    }
                }
            }
        }

        //calculate immonium ions
        if (FragmentIonConstants.fragmentIonHierarchySet.contains("imm")) {
            HashSet<Character> checkedImmoniumIons = new HashSet<>();
            for (int i = 0; i < peptide.length(); i++) {
                if (! checkedImmoniumIons.contains(peptide.charAt(i))) {
                    String ionName;

                    //get attached mod
                    if (modMasses.get(i + 1) != 0) {
                        ionName = "Imm:" + peptide.charAt(i) + "[" + modMasses.get(i + 1) + "]";
                    } else {
                        ionName = "Imm:" + peptide.charAt(i);
                    }

                    addToFragmentIons((float) (modMasses.get(i + 1) + AAmap.get(peptide.charAt(i)) - 26.99),
                            new String[]{ionName, "imm"},
                            fragmentIons, annotationMasses);
                    checkedImmoniumIons.add(peptide.charAt(i));
                }
            }
        }
    }

    private void possibleUnispecMzs() throws IOException, URISyntaxException {
        final InputStream stream = getClass().getClassLoader().getResourceAsStream(
                "fragment_annotation/unispec_fragments.txt");
        final InputStreamReader reader = new InputStreamReader(stream);
        final BufferedReader fragmentsFile = new BufferedReader(reader);

        String line;
        while((line = fragmentsFile.readLine()) != null) {
            FragmentAnnotationParser fap = new FragmentAnnotationParser(line);
            if (! isFragmentPossible(fap)) {
                continue;
            }

            switch (fap.fragmentIonType) {
                case "y":
                case "b":
                case "a":
                    addToFragmentIons(calcMass(fap.fragnum, fap.fragmentIonType, fap.charge, fap.isotope),
                            new String[]{fap.fullAnnotation, fap.fragmentIonType},
                            unispecFragmentIons, unispecAnnotationMasses);
                    break;
                case "y-NL":
                case "b-NL":
                case "a-NL":
                    addToFragmentIons(calcMass(fap.fragnum, fap.fragmentIonType, fap.charge, fap.neutralLoss, fap.isotope),
                            new String[]{fap.fullAnnotation, fap.fragmentIonType},
                            unispecFragmentIons, unispecAnnotationMasses);
                    break;
                case "p": //p and p-NL are assumed to be the precursor charge
                    int thisCharge = fap.charge;
                    if (thisCharge == 0) { //for p, this means we need to use the precursor charge
                        thisCharge = charge;
                    }
                    addToFragmentIons(calcMassPrecursor(fap.isotope, 0, thisCharge),
                            new String[]{fap.fullAnnotation, fap.fragmentIonType},
                            unispecFragmentIons, unispecAnnotationMasses);
                    break;
                case "p-NL":
                    thisCharge = fap.charge;
                    if (thisCharge == 0) {
                        thisCharge = charge;
                    }
                    addToFragmentIons(calcMassPrecursor(fap.isotope, fap.neutralLoss, thisCharge),
                            new String[]{fap.fullAnnotation, fap.fragmentIonType},
                            unispecFragmentIons, unispecAnnotationMasses);
                    break;
                case "imm":
                    addToFragmentIons(calcMassImmonium(fap.fullAnnotation, fap.neutralLoss, fap.isotope),
                            new String[]{fap.fullAnnotation, fap.fragmentIonType},
                            unispecFragmentIons, unispecAnnotationMasses);
                    break;
                case "int":
                    //num1 is peptide length - start + 1
                    //num2 should be extent
                    //unispec only has b/y internal fragments
                    addToFragmentIons(calcMass(peptide.length() - fap.internalStartPosition + 1, fap.internalExtent,
                            "by", fap.isotope),
                            new String[]{fap.fullAnnotation, fap.fragmentIonType},
                            unispecFragmentIons, unispecAnnotationMasses);
                    break;
            }
        }
        fragmentsFile.close();
    }

    private boolean isFragmentPossible(FragmentAnnotationParser fap) {
        switch (fap.fragmentIonType) {
            case "y":
            case "b":
            case "a":
            case "y-NL":
            case "b-NL":
            case "a-NL":
                if (fap.fragnum >= peptide.length()) {
                    return false;
                }
            break;
            case "p":
            case "p-NL":
                if ((fap.charge != 0) &&
                        (Constants.FragmentationType.equals("HCD") || Constants.FragmentationType.equals("CID"))) {
                    //charge loss from precursor does not happen except in electron based fragmentation
                    return false;
                }
                break;
            case "int":
                if (fap.internalStartPosition < internalAaFromNterm - 1) { //TODO: could adapt
                    return false;
                }
                if (fap.internalStartPosition >= peptide.length()) {
                    return false;
                }
                if (fap.internalStartPosition + fap.internalExtent >= peptide.length() - 1) {
                    return false;
                }
                break;
        }
        return true;
    }

    //takes array of mz values and tries to annotate
    public String[][] annotateMZs(float[] mzs,
                                  String mode,
                                  boolean daltonTolerance) throws IOException, URISyntaxException {
        SortedMap<Float, String[]> mzToAnnotationMap = null;
        if (mode.equals("default")) {
            if (fragmentIons.isEmpty()) {
                possibleFragmentIons();
            }
            mzToAnnotationMap = fragmentIons;
        } else if (mode.equals("unispec")) {
            if (unispecFragmentIons.isEmpty()) {
                possibleUnispecMzs();
            }
            mzToAnnotationMap = unispecFragmentIons;
        } else {
            printError(mode + " not supported for mz annotation. Exiting");
            System.exit(1);
        }

        String[] annotations = new String[mzs.length]; //full annotation with everything
        String[] fragmentIonTypes = new String[mzs.length];

        //convert to arrays to be faster
        float[] mzArray = new float[mzToAnnotationMap.size()];
        int k = 0;
        for (float f : mzToAnnotationMap.keySet()) {
            mzArray[k] = f;
            k += 1;
        }
        String[][] ionsArray = new String[mzToAnnotationMap.size()][];
        k = 0;
        for (String[] s : mzToAnnotationMap.values()) {
            ionsArray[k] = s;
            k += 1;
        }

        int q = 0;
        for (int i = 0; i < mzs.length; i++) {
            //for predicted peak, +- fragment error tolerance
            float minMZ;
            float maxMZ;
            if (daltonTolerance) {
                minMZ = mzs[i] - Constants.DaTolerance;
                maxMZ = mzs[i] + Constants.DaTolerance;
            } else { //ppm tolerance
                minMZ = mzs[i] * (1 - Constants.ppmTolerance);
                maxMZ= mzs[i] * (1 + Constants.ppmTolerance);
            }

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
            if (consideredFragmentIons.isEmpty()) {
                annotations[i] = "unknown";
                fragmentIonTypes[i] = "unknown";
                continue;
            }

            //annotate with highest priority fragment ion type
            StringBuilder sb = new StringBuilder();
            for (String ionType : FragmentIonConstants.fragmentIonHierarchy) {
                for (String[] candidates : consideredFragmentIons) {
                    if (candidates[1].equals(ionType)) {
                        if (!sb.toString().isEmpty()) { //multiple matches
                            sb.append(";");
                        }
                        sb.append(candidates[0]);
                    }
                }

                if (!sb.toString().isEmpty()) {
                    annotations[i] = sb.toString();
                    fragmentIonTypes[i] = ionType;
                    break;
                }
            }
        }

        return new String[][] {annotations, fragmentIonTypes};
    }

    //returns float value of m/z to add onto this fragment
    public float compareModMasses(MassCalculator oldMc,
                                  int fragNum, String ionType, int charge, String fullAnnotation) {
        double newSum;
        double oldSum;
        switch (ionType) {
            case "p":
            case "p-NL":
                newSum = sum(modMasses);
                oldSum = sum(oldMc.modMasses);
                return (float) ((newSum - oldSum) / charge);
            case "imm":
                //if difference is localized on an amino acid, immonium ion derived from amino acid should also change mass
                //will need to have text file specifying immonium abbreviations and localized AA
                //this can be future TODO
                return 0;
            case "unknown":
                return 0;
            case "int":
            case "int-NL":
                //get position within sequence
                //get rid of NL and isotope info
                String annot = "";
                for (int i = 0; i < fullAnnotation.length(); i++) {
                    char c = fullAnnotation.charAt(i);
                    if (Character.isLetter(c)) {
                        annot += fullAnnotation.charAt(i);
                    } else {
                        break;
                    }
                }

                int start = peptide.indexOf(annot);
                if (start == -1) {
                    System.out.println("Internal fragment " + annot +
                            " not in peptide " + peptide + ". Exiting");
                    System.exit(1);
                }

                newSum = sum(modMasses.subList(start + 1, start + 1 + annot.length()));
                oldSum = sum(oldMc.modMasses.subList(start + 1, start + 1 + annot.length()));

                return (float) ((newSum - oldSum) / charge);
            case "a":
            case "b":
            case "c":
            case "a-NL":
            case "b-NL":
            case "c-NL":
                newSum = sum(modMasses.subList(0, fragNum + 1));
                oldSum = sum(oldMc.modMasses.subList(0, fragNum + 1));
                if (this.peptide.length() == fragNum) {
                    newSum += modMasses.get(fragNum + 1);
                    oldSum += oldMc.modMasses.get(fragNum + 1);
                }
                return (float) ((newSum - oldSum) / charge);
            case "x":
            case "y":
            case "z":
            case "x-NL":
            case "y-NL":
            case "z-NL":
                newSum = sum(modMasses.subList(this.peptide.length() - fragNum + 1, this.peptide.length() + 2));
                oldSum = sum(oldMc.modMasses.subList(this.peptide.length() - fragNum + 1, this.peptide.length() + 2));
                if (this.peptide.length() == fragNum) {
                    newSum += modMasses.get(0);
                    oldSum += oldMc.modMasses.get(0);
                }
                return (float) ((newSum - oldSum) / charge);
            default:
                System.out.println(ionType + " not recognized. Exiting");
                System.exit(1);
        }
        return 0;
    }

    public static void main(String[] args) throws IOException, URISyntaxException {
        Constants.FragmentationType = "HCD";
        MassCalculator mc = new MassCalculator("PEPTIDER", 3);
        mc.possibleUnispecMzs();
    }
}
