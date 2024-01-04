package External;

public class PeptideSkipper {
    //provide peptide and see if it may be problematic
    public static boolean skipPeptide(String stripped, String charge) {
        //letters
        for (char c : "OUBZJX".toCharArray()) {
            if (stripped.indexOf(c) != -1) {
                return true;
            }
        }
        //length
        if (stripped.length() < 7 || stripped.length() > 20) {
            return true;
        }
        //charge
        if (Integer.parseInt(charge) > 6) {
            return true;
        }
        return false;
    }
}
