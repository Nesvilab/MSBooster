package Features;

import java.util.HashMap;

public class modFormatter {
    int[] pos;
    double[] mass;
    HashMap <Double, String> modifications= new HashMap<Double, String>();

    public modFormatter(int[] positions, double[] masses) {
        pos = positions;
        mass = masses;

        modifications.put(160.03065, "Carbamidomethyl[C]"); //msfragger 3.1
        //modifications.put(160.0307, "Carbamidomethyl[C]"); //msfragger 3.0
        modifications.put(147.0354, "Oxidation[M]");
    }

    public String format() {
        int len = pos.length;
        String formatted = "";
        for (int i = 0; i < len; i++) {
            String toAdd = pos[i] + "," + modifications.get(mass[i]) + ";";
            formatted += toAdd;
        }
        return formatted;
    }
}
