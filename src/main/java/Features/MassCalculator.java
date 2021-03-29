package Features;

import java.util.ArrayList;
import java.util.HashMap;

public class MassCalculator {
    private final float proton = 1.00727647f;
    private final float H2O = 18.01056f;
    public String fullPeptide;
    private String peptide;
    public int charge;
    private ArrayList<Float> modMasses = new ArrayList<Float>();
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

    //first formatting DIA-NN format to common one
    public MassCalculator(String pep, String charge) {
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
            pep = pep.substring(0, ind) + pep.substring(ind2 + 1, pep.length());

            //add modNum to modMasses
            while (modMasses.size() < ind - 1) {
                modMasses.add(0f);
            }
            modMasses.add(unimodToMass.get(modNum));
        }

        peptide = pep;
        fullPeptide = pep + "|" + myMods + "|" + charge;

        //fill in remaining zeros
        while (modMasses.size() < peptide.length()) {
            modMasses.add(0f);
        }

        //set charge
        this.charge = Integer.parseInt(charge);
    }

    public float calcMass(int num, int flag, int charge) {
        float mass = 0f;
        if (flag == 0) { //b
            for (int i = 0; i < num; i++) {
                mass += AAmap.get(this.peptide.charAt(i)); //sum amino acid masses
                mass += modMasses.get(i); //sum mods
            }
        } else if (flag == 1) { //y
            mass += H2O;
            for (int i = this.peptide.length() - num; i < this.peptide.length(); i++) {
                mass += AAmap.get(this.peptide.charAt(i)); //sum amino acid masses
                mass += modMasses.get(i); //sum mods
            }
        }

        //calculate m/z using charge
        return (mass + (float) charge * proton) / (float) charge;
    }

    public float[] calcAllMasses() { //currently up to charge 2, does not include precursor peak
        int maxCharge = Math.min(2, charge);

        ArrayList<Float> masses = new ArrayList<Float>();

        for (int i = 1; i < peptide.length(); i++) {
            masses.add(calcMass(i, 0, 1));
        }
        for (int i = 1; i < peptide.length(); i++) {
            masses.add(calcMass(i, 1, 1));
        }
        if (maxCharge == 2) {
            for (int i = 1; i < peptide.length(); i++) {
                masses.add(calcMass(i, 0, 2));
            }
            for (int i = 1; i < peptide.length(); i++) {
                masses.add(calcMass(i, 1, 2));
            }
        }

        float[] results = new float[masses.size()];
        for(int i = 0; i < masses.size(); i++) results[i] = masses.get(i);
        return results;
    }
}
