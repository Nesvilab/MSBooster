package Features;

import java.io.BufferedReader;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Map;

//only interested in unique target peptides
public class FastaReader {
    HashMap<String, ArrayList<String>> protToPep = new HashMap<String, ArrayList<String>>();

    public FastaReader(String fasta) throws FileNotFoundException {
        //load fasta
        BufferedReader reader;
        int prefixLen = Constants.decoyPrefix.length();
        try {
            reader = new BufferedReader(new FileReader(fasta));
            String header = reader.readLine();
            String protein = reader.readLine();

            HashMap<String, ArrayList<String>> pepToProt = new HashMap<String, ArrayList<String>>();
            while (header != null) {
                if (!header.substring(0, prefixLen).equals(Constants.decoyPrefix)) { //only work with target proteins
                    //split by whatever is digestion rules
                    int start = 0;
                    for (int i = 0; i < protein.length(); i++) {
                        boolean doIt = false;
                        if (i == protein.length() - 1) { //end of protein sequence
                            doIt = true;
                        } else if (Constants.cutAfter.contains(protein.substring(i, i + 1))) {
                            if (! Constants.butNotAfter.contains(protein.substring(i + 1, i + 2))) {
                                doIt = true;
                            }
                        }

                        if (doIt) {
                            String protID = header.split("\\|")[1];
                            String pep = protein.substring(start, i + 1);
                            ArrayList<String> value;
                            if (pepToProt.containsKey(pep)) {
                                value = pepToProt.get(pep);
                            } else {
                                value = new ArrayList<String>();    //TODO: is it better to make an empty arrylist object outside loop?
                                                                    //TODO: make it a hashset in the first place?
                            }
                            value.add(protID);
                            pepToProt.put(pep, value);
                            start = i + 1;
                        }
                    }
                }
                header = reader.readLine();
                protein = reader.readLine();
            }
            reader.close();

            //if peptide passes digestion criteria and is unique, add to protToPep
            for (Map.Entry<String, ArrayList<String>> entry : pepToProt.entrySet()) {
                String pep = entry.getKey();
                ArrayList<String> prots = entry.getValue();

                //length
                if (pep.length() < 7 || pep.length() > 50) {
                    continue;
                }

                //mass
                MassCalculator mc = new MassCalculator(pep, 1);
                if (mc.mass < Constants.digestMinMass || mc.mass > Constants.digestMaxMass) {
                    continue;
                }

                //unique
                if (prots.size() != 1) {
                    HashSet<String> protsSet = new HashSet<String>(prots); //in case rev and regular have shared peptide
                    if (protsSet.size() != 1) {
                        continue;
                    }
                }

                //add to protToPep
                ArrayList<String> value;
                String prot = prots.get(0);
                if (protToPep.containsKey(prot)) {
                    value = protToPep.get(prot);
                } else {
                    value = new ArrayList<String>();
                }
                value.add(pep);
                protToPep.put(prot, value);
            }

        } catch (IOException e) {
            e.printStackTrace();
        }
    }

    public static void main(String[] args) throws FileNotFoundException {
        FastaReader fasta = new FastaReader(Constants.fasta);
    }
}
