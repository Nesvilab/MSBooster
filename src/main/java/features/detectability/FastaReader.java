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

package features.detectability;

import allconstants.Constants;
import features.spectra.MassCalculator;

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
    HashMap<String, ProteinEntry> protToPep = new HashMap<>();

    private String[] getHeaderProtein(BufferedReader reader) throws IOException {
        String header = reader.readLine();
        String protein = reader.readLine();
        reader.mark(8192); //default buffer size
        while (true) {
            String s = reader.readLine();
            if (s == null) { //end of file
                break;
            }
            if (s.length() == 0) { //empty line
                reader.mark(8192);
                continue;
            }
            if (s.charAt(0) == '>') {
                reader.reset();
                break;
            } else {
                protein += s;
                reader.mark(8192);
            }
        }
        return new String[]{header, protein};
    }

    public FastaReader(String fasta) throws FileNotFoundException {
        //load fasta
        BufferedReader reader;
        try {
            reader = new BufferedReader(new FileReader(fasta));
            String[] sArray = getHeaderProtein(reader);
            String header = sArray[0];
            String protein = sArray[1];

            HashMap<String, HashSet<String>> pepToProt = new HashMap<String, HashSet<String>>();
            while (header != null) {
                String protID = header.split(" ")[0];

                if (! protID.startsWith(Constants.decoyPrefix)) {
                    protID = protID.substring(1); //remove >

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
                            String pep = protein.substring(start, i + 1);
                            HashSet<String> value;
                            if (pepToProt.containsKey(pep)) {
                                value = pepToProt.get(pep);
                            } else {
                                value = new HashSet<String>();
                            }
                            value.add(protID);
                            pepToProt.put(pep, value);
                            start = i + 1;
                        }
                    }
                }

                sArray = getHeaderProtein(reader);
                header = sArray[0];
                protein = sArray[1];
            }

            reader.close();

            //if peptide passes digestion criteria and is unique, add to protToPep
            for (Map.Entry<String, HashSet<String>> entry : pepToProt.entrySet()) {
                String pep = entry.getKey();
                HashSet<String> prots = entry.getValue();

                //length
                if (pep.length() < Constants.digestMinLength || pep.length() > Constants.digestMaxLength) {
                    continue;
                }

                //mass
                MassCalculator mc = new MassCalculator(pep, 1);
                if (mc.mass < Constants.digestMinMass || mc.mass > Constants.digestMaxMass) {
                    continue;
                }

                //unique
                int protSize = 0;
                for (String prot : prots) {
                    if (! prot.startsWith(Constants.decoyPrefix)) {
                        protSize += 1;
                    }
                }
                if (protSize > 1) {
                    continue;
                }

                //cannot contain ambiguous amino acid
                if (pep.contains("B") || pep.contains("X") || pep.contains("Z")) {
                    continue;
                }

                //add to protToPep
                ArrayList<String> value;
                String prot = prots.iterator().next();
                if (protToPep.containsKey(prot)) {
                    protToPep.get(prot).peptides.add(pep);
                } else {
                    ProteinEntry newProt = new ProteinEntry();
                    newProt.peptides.add(pep);
                    protToPep.put(prot, newProt);
                }
            }
        } catch (IOException e) {
            e.printStackTrace();
            System.exit(1);
        }
    }
}
