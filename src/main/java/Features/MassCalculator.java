package Features;

import java.util.ArrayList;
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
    public ArrayList<Double> modMasses = new ArrayList<Double>();

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
            String modNum;
            try { //known mod
                modNum = pep.substring(ind, ind2).split(":")[1];
                //add modNum to modMasses
                while (modMasses.size() < ind) {
                    modMasses.add(0d);
                }
                modMasses.add(Constants.AAmassToUnimod.get(modNum));
            } catch (Exception e) {
                modNum = pep.substring(ind + 1, ind2);
                //add modNum to modMasses
                while (modMasses.size() < ind) {
                    modMasses.add(0d);
                }
                modMasses.add(Double.parseDouble(modNum));
            }

            //replace mod
            pep = pep.substring(0, ind) + pep.substring(ind2 + 1);
        }

        peptide = pep;

        //fill in remaining zeros
        while (modMasses.size() < peptide.length() + 1) {
            modMasses.add(0d);
        }

        //set charge
        if (charge instanceof String) {
            this.charge = Integer.parseInt((String) charge);
        } else {
            this.charge = (int) charge;
        }

        this.mass = calcMass(pep.length(), 1, 1) - H;
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

    public static void main(String[] args) {
        MassCalculator mc = new MassCalculator("Q[-17.0265]LVPALAKV", 1);
        //mc.calcAllMasses();
        System.out.println(mc.mass);
        System.out.println(mc.calcMass(8, 0, 1));
        System.out.println(mc.calcMass(5, 0, 1));
        System.out.println(mc.calcMass(8, 1, 1));
        System.out.println(mc.calcMass(5, 1, 1));
//        System.out.println(Arrays.toString(mc.b1));
//        System.out.println(Arrays.toString(mc.y1));
//        System.out.println(Arrays.toString(mc.b2));
//        System.out.println(Arrays.toString(mc.y2));
    }
}
