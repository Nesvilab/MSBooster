import java.util.ArrayList;
import java.util.Collections;

public class peptidesList {
    //input
    String peptide;
    boolean ntermAcetyl;
    ArrayList<Integer> positions;
    ArrayList<Double> modificationMasses;
    int assumedCharge;

    //created by this class
    String rowInTable; //<peptide (string), modifications (string), charge (int)>

    public peptidesList(String pep, boolean nterm, ArrayList<Integer> pos, ArrayList<Double> modmass, int charge) {
        peptide = pep;
        ntermAcetyl = nterm;
        positions = pos;
        modificationMasses = modmass;
        assumedCharge = charge;

        if (positions.size() > 1) {
            //sorting method
            ArrayList<Integer> clonedPositions = new ArrayList<Integer>(positions);
            Collections.sort(positions);
            ArrayList<Double> modIndex = new ArrayList<Double>();
            for (int i = 0; i < positions.size(); i++) {
                int posIndex = clonedPositions.indexOf(positions.get(i));
                modIndex.add(modmass.get(posIndex));
            }
            modificationMasses = modIndex;
        }
    }

    //convert modification mass to modification string

    //create rowInTable

    //get unique
}
