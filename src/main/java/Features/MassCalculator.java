package Features;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;

public class MassCalculator {
    private final float proton = 1.00727647f;
    private final float H = 1.0078250321f;
    private final float O = 15.9949146221f;
    //private final float N = 14.0030740052f;
    private final float H2O = H * 2 + O;
    //private final float OH = O + H;
    public String fullPeptide;
    private String peptide;
    public int charge;
    public float mass;
    private ArrayList<Float> modMasses = new ArrayList<Float>(); //no inclusion of C-terminal mods?

    //ion series
    public float[] b1;
    public float[] b2;
    public float[] y1;
    public float[] y2;
    private final HashMap<Character, Float> AAmap = new HashMap<Character, Float>()
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
    private final HashMap<String, String> unimodToMods = new HashMap<String, String>()
    {{
        put("1", "Acetyl[AnyN-term]");
        put("4", "Carbamidomethyl[C]");
        put("35", "Oxidation[M]");
    }};
    private final HashMap<String, Float> unimodToMass = new HashMap<String, Float>()
    {{
        put("1", 42.010565f);
        put("4", 57.021464f);
        put("35", 15.994915f);
    }};
    private final HashMap<String, Float> modToMass = new HashMap<String, Float>()
    {{
        put("Acetyl[AnyN-term]", 42.010565f);
        put("Carbamidomethyl[C]", 57.021464f);
        put("Oxidation[M]", 15.994915f);
    }};

    //first formatting DIA-NN format to common one
    public MassCalculator(String pep, Object charge) {
        StringBuilder myMods = new StringBuilder();

        while (true) {
            int ind = pep.indexOf("[");
            if (ind == -1) {
                break;
            }
            int ind2 = pep.indexOf("]");

            //get mod
            String modNum = pep.substring(ind, ind2).split(":")[1];
            String mod = unimodToMods.get(modNum);
            myMods.append(ind).append(",").append(mod).append(";");

            //replace mod
            pep = pep.substring(0, ind) + pep.substring(ind2 + 1);

            //add modNum to modMasses
            while (modMasses.size() < ind) {
                modMasses.add(0f);
            }
            modMasses.add(unimodToMass.get(modNum));
        }

        peptide = pep;
        fullPeptide = pep + "|" + myMods + "|" + charge;

        //fill in remaining zeros
        while (modMasses.size() < peptide.length() + 1) {
            modMasses.add(0f);
        }

        //set charge
        if (charge instanceof String) {
            this.charge = Integer.parseInt((String) charge);
        } else {
            this.charge = (int) charge;
        }

        this.mass = calcMass(pep.length(), 1, 1) - H;
//        this.b1 = new float[peptide.length()];
//        this.y1 = new float[peptide.length()];
//        if (this.charge >= 2) {
//            this.b2 = new float[peptide.length()];
//            this.y2 = new float[peptide.length()];
//        }
    }

    //using common peptide form pep|mod|charge
    //TODO: correct this, but probably won't have to use since no consecutiveFragment feature calculation
    public MassCalculator(String pep) {
        this.fullPeptide = pep;

        String[] pepSplit = pep.split("\\|");
        this.peptide = pepSplit[0];
        this.charge = Integer.parseInt(pepSplit[2]);

        //get mod positions
        float[] modsArray = new float[peptide.length() + 1];
        String modsString = pepSplit[1];
        String[] mods = modsString.split(";");
        for (String mod : mods) {
            if (mod.equals("")) {
                break;
            }
            String[] stringSplit = mod.split(",");
            modsArray[Integer.parseInt(stringSplit[0])] = modToMass.get(stringSplit[1]);
        }

        //add to modMasses
        for (float f : modsArray) {
            modMasses.add(f);
        }

        this.mass = calcMass(peptide.length(), 0, 1);
        this.b1 = new float[peptide.length()];
        this.y1 = new float[peptide.length()];
        if (this.charge >= 2) {
            this.b2 = new float[peptide.length()];
            this.y2 = new float[peptide.length()];
        }
    }

    public float calcMass(int num, int flag, int charge) {
        float mass = 0f;
        if (flag == 0) { //b
            mass += modMasses.get(0);
            for (int i = 0; i < num; i++) {
                mass += AAmap.get(this.peptide.charAt(i)); //sum amino acid masses
                mass += modMasses.get(i + 1); //sum mods
            }
        } else if (flag == 1) { //y
            mass += H2O;
            for (int i = this.peptide.length() - num; i < this.peptide.length(); i++) {
                mass += AAmap.get(this.peptide.charAt(i)); //sum amino acid masses
                mass += modMasses.get(i + 1); //sum mods
            }
            if (num == this.peptide.length()) {
                mass += modMasses.get(0);
            }
        }

        //calculate m/z using charge
        return (mass + (float) charge * proton) / (float) charge;
    }

    public void calcAllMasses() { //currently up to charge 2, does not include precursor peak
        int maxCharge = Math.min(2, charge);

        for (int i = 1; i < peptide.length() + 1; i++) {
            b1[i - 1] = calcMass(i, 0, 1);
        }
        for (int i = 1; i < peptide.length() + 1; i++) {
            y1[i - 1] = calcMass(i, 1, 1);
        }
        if (maxCharge == 2) {
            for (int i = 1; i < peptide.length() + 1; i++) {
                b2[i - 1] = calcMass(i, 0, 2);
            }
            for (int i = 1; i < peptide.length() + 1; i++) {
                y2[i - 1] = calcMass(i, 1, 2);
            }
        }
    }

    //Given an experimental spectrum, see which series has most consecutive fragment ions detected.
    //Returns that int after checking b/y charge 1/2
    //TODO: try with and without normalizing by peptide length
    private int maxConsecutiveIonSeries(float[] expMZs, float[] expIntensities, float[] series) {
        //we don't actually care what the fragment intensities are here
        float[] predIntensities = new float[series.length];
        Arrays.fill(predIntensities, 1f);

        spectrumComparison sc = new spectrumComparison(expMZs, expIntensities, series, predIntensities);
        float[] matchedFrags = sc.matchedIntensities;
        return StatMethods.consecutiveMatches(matchedFrags);
    }

    //TODO: problem with DIANN peptide formatting?
    public float maxConsecutiveIonSeries(float[] expMZs, float[] expIntensities) {
        calcAllMasses();
        int max = 0;
        max = Math.max(max, maxConsecutiveIonSeries(expMZs, expIntensities, y1));
        max = Math.max(max, maxConsecutiveIonSeries(expMZs, expIntensities, b1));
        if (b2 != null) {
            max = Math.max(max, maxConsecutiveIonSeries(expMZs, expIntensities, y2));
            max = Math.max(max, maxConsecutiveIonSeries(expMZs, expIntensities, b2));
        }
        return (float) max / (float) peptide.length();
    }

    public static void main(String[] args) {
        MassCalculator mc = new MassCalculator("[unimod:1]VC[unimod:4]AKIGDFGMAR", 2);
        mc.calcAllMasses();
        System.out.println(mc.mass);
        System.out.println(Arrays.toString(mc.b1));
        System.out.println(Arrays.toString(mc.y1));
        System.out.println(Arrays.toString(mc.b2));
        System.out.println(Arrays.toString(mc.y2));
        System.out.println(mc.calcMass(2,0,2));
    }
}
