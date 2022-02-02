package Features;

import java.util.*;

public class MassCalculator {
    private final float proton = 1.00727647f;
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
    private String[] fragmentIonHierarchy;
    public static final Set<String> allowedFragmentIonTypes = new HashSet<>(Arrays.asList(
            "z", "c", "y", "a", "x", "b",
            "precursor", "immonium", "internal", "internal-NL",
            "z-NL", "c-NL", "y-NL", "a-NL", "x-NL", "b-NL",
            "precursor-NL", "unknown"));
    private void makeFragmentIonHierarchy() { //decided by looking at multiple papers
        if (Constants.FragmentationType.equals("ETD")) {
            fragmentIonHierarchy = new String[] {"z", "c", "y", "a", "x", "b",
                    "precursor", "immonium", "internal", "internal-NL",
                    "z-NL", "c-NL", "y-NL", "a-NL", "x-NL", "b-NL",
                    "precursor-NL"};
        } else { //HCD, CID, or not specified
            fragmentIonHierarchy = new String[] {"y", "b", "a", "c", "x", "z",
                    "precursor", "immonium", "internal", "internal-NL",
                    "y-NL", "b-NL", "a-NL", "c-NL", "x-NL", "z-NL",
                    "precursor-NL"};
        }
    }

    //TODO: get immonium ion masses by taking amino acid and subtracting 26.99 Da. This holds for modified AA too (Falick et al 1993)
    //TODO: should we consider related ions, not just immonium?
    public final HashMap<Character, Float> AAmap = new HashMap<Character, Float>()
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
        put('O', 255.15829f);
        put('U', 150.95363f);
        put('Z', 0f);
        put('B', 0f);
        put('X', 0f);
    }};

    //first formatting DIA-NN format to common one
    public MassCalculator(String pep, Object charge) {
        fullPeptide = pep + "|" + charge;

        while (true) {
            int ind = pep.indexOf("[");
            if (ind == -1) {
                break;
            }
            int ind2 = pep.indexOf("]");

            //get mod
//            String modNum;
//            try { //known mod
//                modNum = pep.substring(ind, ind2).split(":")[1];
//                //add modNum to modMasses
//                while (modMasses.size() < ind) {
//                    modMasses.add(0d);
//                }
//                modMasses.add(Constants.unimodtoModAAmass.get(modNum));
//            } catch (Exception e) {
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
    }

    public float calcMass(int num, String flag, int charge) {
        float mass = 0f;

        switch(flag) {
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
        }

        //calculate m/z using charge
        return (mass + (float) charge * proton) / (float) charge;
    }

    //num1 is the y ion the internal fragment is derived from
    //num2 is the a/b ion that is formed from that
    //assume charge 1 only
    public float calcMass(int num1, int num2, String flag) {
        float mass = 0f;

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
                System.out.println(flag + " is not supported. Choose ay or by");
                System.exit(-1);
        }
        return 0;
    }

    public float calcMass(int num, String flag, int charge, String neutralLoss) {
        float mass = 0f;

        switch(flag) {
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
        }

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
        float mass = 0f;

        MassCalculator newY = internalPeptides[num1 - 3];
        if (newY == null) {
            newY = makeInternalPeptide(num1);
        }

        switch (flag) {
            case "by":
                return newY.calcMass(num2, "b", 1, neutralLoss);
            case "ay":
                return newY.calcMass(num2, "a", 1, neutralLoss);
            default:
                System.out.println(flag + " is not supported. Choose ay or by");
                System.exit(-1);
        }
        return 0;
    }

    private MassCalculator() {}
    //y is which y fragment
    private MassCalculator makeInternalPeptide(int y) {
        if (y >= this.peptide.length()) { //only make smaller peptides
            System.out.println("internal fragment must be shorter than initial peptide");
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
        if (! fragmentIons.containsKey(mass)) {
            fragmentIons.put(mass, info);
        } else { //found an overlap
            fragmentIons.put(mass,
                    new String[] {fragmentIons.get(mass)[0] + ";" + info[0],
                            fragmentIons.get(mass)[1] + ";" + info[1]});
        }
    }

    //currently up to charge 2
    //format of value is [name, fragment ion type]
    //only support 1 neutral loss
    public void possibleFragmentIons() {
        //calculate precursor isotopic peaks
        for (int iCharge = 1; iCharge < charge + 1; iCharge++) {
            //regular
            String ionName = "MH" + "+" + iCharge;
            addToFragmentIons(calcMass(peptide.length(), "y", iCharge), new String[] {ionName, "precursor"});

            //neutral loss
            for (String nl : neutralLosses) {
                ionName = "MH-" + nl + "" + "+" + iCharge;
                addToFragmentIons(calcMass(peptide.length(), "y", iCharge, nl), new String[] {ionName, "precursor-NL"});
            }
        }

        //calculate all abcxyz ions
        String ions = "abcxyz";
        int maxCharge = Math.min(2, charge);
        for (int i = 0; i < ions.length(); i++) {
            String ionType = ions.substring(i, i + 1);
            for (int num = 1; num < this.peptide.length(); num++) {
                for (int iCharge = 1; iCharge < maxCharge + 1; iCharge++) {
                    //regular fragment
                    String ionName = ionType + num + "+" + iCharge;
                    addToFragmentIons(calcMass(num, ionType, iCharge), new String[] {ionName, ionType});

                    //neutral loss fragments
                    if (! ionType.equals("c")) {
                        for (String nl : neutralLosses) {
                            ionName = ionType + num + "-" + nl + "+" + iCharge;
                            addToFragmentIons(calcMass(num, ionType, iCharge, nl), new String[]{ionName, ionType + "-NL"});
                        }
                    } else { //c-NH3 is same as b
                        ionName = ionType + num + "-H2O" + "+" + iCharge;
                        addToFragmentIons(calcMass(num, ionType, iCharge, "H2O"), new String[]{ionName, "c-NL"});
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

                    String ionName = "y" + num1 + ionType.charAt(0) + num2;
                    addToFragmentIons(calcMass(num1, num2, ionType), new String[] {ionName, "internal"});

                    for (String nl : neutralLosses) {
                        ionName = "y" + num1 + ionType.charAt(0) + num2 + "-" + nl;
                        addToFragmentIons(calcMass(num1, num2, ionType, nl), new String[] {ionName, "internal-NL"});
                    }
                }
            }
        }

        //don't need internal peptides anymore
        //internalPeptides = null;

        //calculate immonium ions
        for (int i = 0; i < peptide.length(); i++) {
            String ionName;

            //get attached mod
            if (modMasses.get(i + 1) != 0) {
                ionName = "immonium:" + peptide.charAt(i) + "[" + modMasses.get(i + 1) + "]";
            } else {
                ionName = "immonium:" + peptide.charAt(i);
            }

            addToFragmentIons((float) (modMasses.get(i + 1) + AAmap.get(peptide.charAt(i)) - 26.99),
                    new String[] {ionName, "immonium"});
        }
    }

    public String[][] annotateMZs(float[] fs) {
        makeFragmentIonHierarchy();

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
            for (String ionType : fragmentIonHierarchy) {
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

    public static void main(String[] args) {
        MassCalculator mc = new MassCalculator("RTGSVVRQK", 3);
        mc.possibleFragmentIons();
        for (Float mz : mc.fragmentIons.keySet()) {
            System.out.println(mz + " " + mc.fragmentIons.get(mz)[1]);
        }
    }
}
