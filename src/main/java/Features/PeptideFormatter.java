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

package Features;

import static utils.Print.printError;

import java.util.ArrayList;
import java.util.Map;

//lots of different ways to format peptide string
public class PeptideFormatter {
    String base;
    String diann;
    String predfull;
    String prosit;
    String prositTMT;
    String stripped;
    String baseCharge;
    String dlib;
    String mods = "";
    String alphapeptdeepMods = "";
    String modPositions = "";

    ArrayList<Integer> starts = new ArrayList<>();
    ArrayList<Integer> ends = new ArrayList<>();

    String charge;

    private void findPTMlocations() {
        for (int i = 0; i < base.length(); i++) {
            if (base.charAt(i) == '[') {
                starts.add(i);
            } else if (base.charAt(i) == ']') {
                ends.add(i);
            }
        }
    }

    private void pinTObase(String peptide) {
        //remove AA and period at beginning and end
        peptide = peptide.substring(2, peptide.length() - 2);

        //TODO remove integer part that holds charge
        while (Character.isDigit(peptide.charAt(peptide.length() - 1))) {
            peptide = peptide.substring(0, peptide.length() - 1);
        }

        //remove 0.0000
        peptide = peptide.replace("[0.0000]", "");

        //n and c term mods
        if (peptide.charAt(0) == 'n') {
            peptide = peptide.replace("n", "");
        }
        peptide = peptide.replace("c","");

        base = peptide;
        findPTMlocations();
    }

    private void diannTObase(String peptide) {
        peptide = peptide.toUpperCase();
        peptide = peptide.replace("UNIMOD:","");
        peptide = peptide.replace("[TMT]","[737]");

        //for koina
        peptide = peptide.replace("]-", "]");
        peptide = peptide.replace("-[", "[");
        //convert unimod to mass
        ArrayList<Integer> newStarts = new ArrayList<>();
        ArrayList<Integer> newEnds = new ArrayList<>();
        for (int i = 0; i < peptide.length(); i++) {
            if (peptide.charAt(i) == '[') {
                newStarts.add(i);
            } else if (peptide.charAt(i) == ']') {
                newEnds.add(i);
            }
        }
        for (int i = newStarts.size() - 1; i > -1; i--) {
            try {
                peptide = peptide.substring(0, newStarts.get(i) + 1) +
                        PTMhandler.unimodtoModAAmass.get(peptide.substring(newStarts.get(i) + 1, newEnds.get(i))) +
                        peptide.substring(newEnds.get(i));
            } catch (Exception ignored) {

            }
        }

        base = peptide;
        findPTMlocations();
    }

    private void prositTObase(String peptide) {
        //remove (O) and add [57] to C
        base = peptide.replace("(ox)", "[" + PTMhandler.oxidationMass + "]")
                .replace("C", "C[" + PTMhandler.carbamidomethylationMass + "]");
        findPTMlocations();
    }

    private void predfullTObase(String peptide) {
        //remove (O) and add [57] to C
        base = peptide.replace("(O)", "[" + PTMhandler.oxidationMass + "]")
                .replace("C", "C[" + PTMhandler.carbamidomethylationMass + "]");
        findPTMlocations();
    }

    private void mspTObase(String peptide) {
        String[] pepSplit = peptide.split("\\|");

        base = pepSplit[0];
        if (! pepSplit[1].equals("0")) { //mods
            String[] mods = pepSplit[1].split("/");
            for (int i = mods.length - 1; i > 0; i--) {
                String[] modsSplit = mods[i].split(",");
                int position = Integer.parseInt(modsSplit[0]) + 1;
                base = base.substring(0, position) + "[" +
                        PTMhandler.prositToModAAmass.get(modsSplit[2]) + "]" + base.substring(position);
            }
        }

        findPTMlocations();
    }

    private void pdeep3TObase(String peptide) { //need to make robust for mods not accepted
        //add PTMs to peptide sequence
        String[] pepSplit = peptide.split("\\|");
        String newPeptide = pepSplit[0];
        String[] mods = pepSplit[1].split(";");
        if (mods[0].equals("")) {
            base = newPeptide;
        } else {
            for (int i = mods.length - 1; i > -1; i--) {
                String mod = mods[i];
                String[] commaSplit = mod.split(",");
                int position = Integer.parseInt(commaSplit[0]);
                String PTMtype = commaSplit[1].split("\\[")[0];
                newPeptide = newPeptide.substring(0, position) + "[" + PTMhandler.PdeepToAAmass.get(PTMtype) + "]" +
                        newPeptide.substring(position);
            }
            base = newPeptide;
            findPTMlocations();
        }
    }

    private void baseTOstripped() {
        stripped = base;

        ArrayList<Integer> newStarts = new ArrayList<>(starts);
        newStarts.add(stripped.length());

        ArrayList<Integer> newEnds = new ArrayList<>();
        newEnds.add(-1);
        newEnds.addAll(ends);

        StringBuilder sb = new StringBuilder();
        for (int i = 0; i < newStarts.size(); i++) {
            sb.append(stripped, newEnds.get(i) + 1, newStarts.get(i));
        }
        stripped = sb.toString();
    }

    private void baseTOdiann() {
        diann = base;

        for (int i = starts.size() - 1; i > -1; i--) {
            double reportedMass = Double.parseDouble(diann.substring(starts.get(i) + 1, ends.get(i)));
            boolean foundReplacement = false;
            for (Map.Entry<Double, Integer> entry : PTMhandler.modAAmassToUnimod.entrySet()) {
                Double PTMmass = entry.getKey();

                if (Math.abs(PTMmass - reportedMass) < 0.01) {
                    if (PTMhandler.modAAToLocalization.containsKey(PTMmass)) {
                        if (!PTMhandler.modAAToLocalization.get(PTMmass).contains(
                                diann.substring(starts.get(i) - 1, starts.get(i)))) {
                            break;
                        }
                    }
                    if (PTMmass == 229.1629) {
                        diann = diann.substring(0, starts.get(i) + 1) + "TMT" +
                                diann.substring(ends.get(i));
                    } else {
                        diann = diann.substring(0, starts.get(i) + 1) + "UniMod:" + PTMhandler.modAAmassToUnimod.get(PTMmass) +
                                diann.substring(ends.get(i));
                    }
                    foundReplacement = true;
                    break;
                }
            }
            if (! foundReplacement) {
                //DIANN won't predict this anyway
                diann = diann.substring(0, starts.get(i)) + diann.substring(ends.get(i) + 1);
            }
        }
    }

    private void baseTOpredfull() {
        predfull = base;

        //PredFull only supports oxM and assumes C is carbamidomethylated. Below, we format it this way, but also keep on other PTM masses so m/z values for fragments can be adjusted
        for (int i = starts.size() - 1; i > -1; i--) {
            double reportedMass = Double.parseDouble(predfull.substring(starts.get(i) + 1, ends.get(i)));
            if (starts.get(i) - 1 > -1) { //no changes to nterm mod
                if (Math.abs(PTMhandler.oxidationMass - reportedMass) < 0.01 &&
                        predfull.charAt(starts.get(i) - 1) == 'M') {
                    predfull = predfull.substring(0, starts.get(i)) + "(O)" + predfull.substring(ends.get(i) + 1);
                } else { //mod unsupported, remove for now. Or C carbamidomethylation
                    predfull = predfull.substring(0, starts.get(i)) + predfull.substring(ends.get(i) + 1);
                }
            } else { //deal with nterm mod by deleting
                predfull = predfull.substring(0, starts.get(i)) + predfull.substring(ends.get(i) + 1);
            }
        }
    }

    private void baseTOprosit() {
        prosit = base;

        //PredFull only supports oxM and assumes C is carbamidomethylated. Below, we format it this way, but also keep on other PTM masses so m/z values for fragments can be adjusted
        for (int i = starts.size() - 1; i > -1; i--) {
            double reportedMass = Double.parseDouble(prosit.substring(starts.get(i) + 1, ends.get(i)));
            if (starts.get(i) - 1 > -1) { //no changes to nterm mod
                if (Math.abs(PTMhandler.oxidationMass - reportedMass) < 0.01 &&
                        prosit.charAt(starts.get(i) - 1) == 'M') {
                    prosit = prosit.substring(0, starts.get(i)) + "(ox)" + prosit.substring(ends.get(i) + 1);
                } else if (Math.abs(PTMhandler.tmtMass - reportedMass) < 0.01) {
                    continue;
                } else { //mod unsupported, remove for now. Or C carbamidomethylation
                    prosit = prosit.substring(0, starts.get(i)) + prosit.substring(ends.get(i) + 1);
                }
            } else { //deal with nterm mod by deleting
                if (!prosit.substring(starts.get(i), ends.get(i) + 1).contains(String.valueOf(PTMhandler.tmtMass))) {
                    prosit = prosit.substring(0, starts.get(i)) + prosit.substring(ends.get(i) + 1);
                }
            }
        }
    }

    private void prositTOprositTMT() {
        //Prosit TMT assumes n-term are TMT-labeled
        if (!prosit.startsWith("[")) {
            String tmtLabel = "[" + PTMhandler.tmtMass + "]";
            prositTMT = tmtLabel + prosit;
        } else {
            prositTMT = prosit;
        }

        //S are not labeled
        prositTMT = prositTMT.replace("S[" + PTMhandler.tmtMass + "]", "S");
    }

    private void strippedTOdlib() { dlib = stripped.replace("C", "C[" + PTMhandler.carbamidomethylationMass + "]") + "|" + charge;
    }

    private void baseToMods() {
        baseTOstripped();

        int numMods = starts.size();
        ArrayList<Integer> newEnds = new ArrayList<>();
        newEnds.add(0);
        newEnds.addAll(ends);
        ArrayList<Integer> positions = new ArrayList<>();
        positions.add(0);

        for (int i = 0; i < numMods; i++) {
            String modMass = base.substring(starts.get(i) + 1, ends.get(i));
            double doubleModMass = Double.parseDouble(modMass);
            String modName = PTMhandler.aamassToAlphapeptdeep.get(modMass);
            if (modName == null) {
                doubleModMass -= 0.0001;
                modName = PTMhandler.aamassToAlphapeptdeep.get(String.format("%.4f", doubleModMass));
            }
            if (modName == null) {
                doubleModMass += 0.0002;
                modName = PTMhandler.aamassToAlphapeptdeep.get(String.format("%.4f", doubleModMass));
            }
            if (modName == null) { //can try checking if it's a combo of fixed and var mod
                printError("There is an unknown modification with mass " + doubleModMass +
                        ". Please provide PTM info via additionalMods param in --paramsList.");
                System.exit(1);
            }

            int position = starts.get(i) - newEnds.get(i) + positions.get(i);
            if (i > 0) {
                position -= 1;
            }
            positions.add(position);
            String aa;
            if (position == 0) {
                aa = "Any N-term"; //might need to change
            } else {
                aa = stripped.substring(position - 1, position);
            }

            String modinfo = position + "," + modName + "[" + aa + "]" + ";";
//            if (aa.equals("ProteinN-term")) {
//                aa = "Protein N-term";
//            }

            //check that the ptmName@aa is accepted by alphapeptdeep
            ArrayList<String> ptmSubstitutes = PTMhandler.sameMass.get(String.format("%.4f", doubleModMass));
            while (! PTMhandler.alphapeptdeepModNames.contains(modName + "@" + aa)) {
                if (ptmSubstitutes == null) { //just until the search is right
                    break;
                }
                if (ptmSubstitutes.size() == 1) {
                    printError("This PTM is not supported with mass " + doubleModMass +
                            ": " + modName + "@" + aa + "\n" +
                            "Please provide PTM info via additional_mods param in --paramsList.");
                    System.exit(1);
                }
                ptmSubstitutes.remove(modName);
                modName = ptmSubstitutes.get(0);
                PTMhandler.aamassToAlphapeptdeep.put(String.format("%.4f", doubleModMass), modName);
            }
            String alphapeptdeepModinfo = PTMhandler.writeOutAlphapeptdeepModNames.get(modName + "@" + aa) + ";";

            mods = mods + modinfo;
            alphapeptdeepMods = alphapeptdeepMods + alphapeptdeepModinfo;
            modPositions = modPositions + position + ";";
        }
        if (! mods.equals("")) {
            mods = mods.substring(0, mods.length() - 1);
            alphapeptdeepMods = alphapeptdeepMods.substring(0, alphapeptdeepMods.length() - 1);
            modPositions = modPositions.substring(0, modPositions.length() - 1);
        }
    }

    public PeptideFormatter(String peptide, Object c, String format) {
//        try {
//            charge = (String) c;
//        } catch (Exception e) {
//            charge = c.toString();
//        }
        charge = c.toString();

        if (format.equals("pin")) {
            pinTObase(peptide);
            baseTOdiann();
            baseTOstripped();
            baseTOpredfull();
            baseTOprosit();
            strippedTOdlib();
            if (Constants.spectraModel.contains("pDeep") || Constants.spectraModel.contains("alphapeptdeep") ||
                    Constants.rtModel.contains("pDeep") || Constants.rtModel.contains("alphapeptdeep") ||
                    Constants.imModel.contains("alphapeptdeep")) {
                baseToMods();
            }
            baseCharge = base + "|" + charge;
        }

        if (format.equals("diann")) {
            diann = peptide;
            diannTObase(peptide);
            baseTOstripped();
            baseTOpredfull();
            baseTOprosit();
            baseCharge = base + "|" + charge;
        }

        if (format.equals("base")) {
            base = peptide;

            //find locations of PTMs
            for (int i = 0; i < base.length(); i++) {
                if (base.charAt(i) == '[') {
                    starts.add(i);
                } else if (base.charAt(i) == ']') {
                    ends.add(i);
                }
            }

            baseTOstripped();
            baseTOpredfull();
            baseTOdiann();
            baseTOprosit();
            prositTOprositTMT();
            baseCharge = base + "|" + charge;
        }

        if (format.equals("predfull")) {
            predfull = peptide;
            predfullTObase(peptide);
            baseTOstripped();
            baseTOdiann();
            baseTOprosit();
            baseCharge = base + "|" + charge;
        }

        if (format.equals("msp")) {
            mspTObase(peptide);
            baseTOstripped();
            baseTOpredfull();
            baseTOdiann();
            baseTOprosit();
            baseCharge = base + "|" + charge;
        }

        if (format.equals("prosit")) {
            prosit = peptide;
            prositTObase(peptide);
            baseTOstripped();
            baseTOpredfull();
            baseTOdiann();
            baseCharge = base + "|" + charge;
        }

        if (format.equals("pdeep3")) {
            pdeep3TObase(peptide);
            baseTOstripped();
            baseCharge = base + "|" + charge;
        }
    }

    public String getBaseCharge() {
        return baseCharge;
    }
    public String getBase() {
        return base;
    }
    public String getDiann() {
        return diann;
    }
    public String getProsit() {
        return prosit;
    }
    public String getStripped() {
        return stripped;
    }

    public String getPrositTMT() { return prositTMT; }
}
