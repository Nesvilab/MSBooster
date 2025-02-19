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

package peptideptmformatting;

import java.util.ArrayList;
import java.util.HashSet;

import static peptideptmformatting.PTMhandler.*;
import static utils.Print.printError;

//lots of different ways to format peptide string
public class PeptideFormatter {
    String base;
    private String diann;
    private String predfull;
    private String prosit;
    private String prositTMT;
    private String unispec;
    private String ms2pip;
    private String deeplc;
    private String alphapept;
    private String stripped;
    String baseCharge;
    private String dlib;
    private String mods = "";
    private String alphapeptdeepMods = "";
    private String modPositions = "";

    private ArrayList<Integer> starts = new ArrayList<>();
    private ArrayList<Integer> ends = new ArrayList<>();

    String charge;

    public HashSet<String> foundUnimods = new HashSet<>(); //collection of previously used unimod codes
    //TODO: should this be shared amongst all PFs? Can store in a thread-safe set

    public boolean cterm = false;

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

        while (Character.isDigit(peptide.charAt(peptide.length() - 1))) {
            peptide = peptide.substring(0, peptide.length() - 1);
        }

        //remove 0.0000
        peptide = peptide.replace("[0.0000]", "");

        //n and c term mods
        if (peptide.charAt(0) == 'n') {
            peptide = peptide.substring(1);
        }
        if (peptide.contains("c")) {
            peptide = peptide.replace("c", "");
            cterm = true;
        }

        base = peptide;
    }

    private void diannTObase(String peptide) {
        peptide = peptide.replace("[TMT]", "[" + PTMhandler.tmtUnimod + "]");
        base = PTMhandler.formatPeptideSpecificToBase(peptide, unimodToModMass);
    }
    private void koinaTObase(String peptide) {
        base = PTMhandler.formatPeptideSpecificToBase(peptide, unimodToModMass);
    }

    private void alphapeptTObase(String peptide) {
        base = PTMhandler.formatPeptideSpecificToBase(peptide, unimodToModMassAlphaPeptDeep);
    }

    private void unimodoboTObase(String peptide) {
        base = PTMhandler.formatPeptideSpecificToBase(peptide, unimodOboToModMass);
    }

    private void predfullTObase(String peptide) {
        //remove (O) and add [57] to C
        base = peptide.replace("(O)", "[" + PTMhandler.oxidationMass + "]")
                .replace("C", "C[" + PTMhandler.carbamidomethylationMass + "]");
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

        boolean attemptCterm = cterm;
        for (int i = starts.size() - 1; i > -1; i--) {
            int start = starts.get(i);
            int end = ends.get(i);

            String[] peptideUnimod = PTMhandler.formatPeptideBaseToSpecific(
                    diann, start, end, "diann", foundUnimods, attemptCterm);
            attemptCterm = false;
            diann = peptideUnimod[0];
            if (!peptideUnimod[1].isEmpty()) {
                foundUnimods.add(peptideUnimod[1]);
            }
        }

        //special TMT formatting
        diann = diann.replaceAll("UniMod:" + PTMhandler.tmtUnimod, "TMT");

        //nterm acetyl and TMT on AA1 not allowed
        if (diann.startsWith("[UniMod") && diann.startsWith("TMT", 12)) {
            diann = diann.substring(10);
        }
    }

    private void baseTOunispec() {
        unispec = base;

        boolean attemptCterm = cterm;
        for (int i = starts.size() - 1; i > -1; i--) {
            int start = starts.get(i);
            int end = ends.get(i);

            String[] peptideUnimod = PTMhandler.formatPeptideBaseToSpecific(
                    unispec, start, end, "unispec", foundUnimods, attemptCterm);
            attemptCterm = false;
            unispec = peptideUnimod[0];
            if (!peptideUnimod[1].isEmpty()) {
                foundUnimods.add(peptideUnimod[1]);
            }
        }
    }

    private void baseTOms2pip() {
        ms2pip = base;

        boolean attemptCterm = cterm;
        for (int i = starts.size() - 1; i > -1; i--) {
            int start = starts.get(i);
            int end = ends.get(i);

            String[] peptideUnimod = PTMhandler.formatPeptideBaseToSpecific(
                    ms2pip, start, end, "ms2pip", foundUnimods, attemptCterm);
            attemptCterm = false;
            ms2pip = peptideUnimod[0];
            if (!peptideUnimod[1].isEmpty()) {
                foundUnimods.add(peptideUnimod[1]);
            }
        }
    }

    private void baseTOdeeplc() {
        deeplc = base;

        boolean attemptCterm = cterm;
        for (int i = starts.size() - 1; i > -1; i--) {
            int start = starts.get(i);
            int end = ends.get(i);

            String[] peptideUnimod = PTMhandler.formatPeptideBaseToSpecific(
                    deeplc, start, end, "deeplc", foundUnimods, attemptCterm);
            attemptCterm = false;
            deeplc = peptideUnimod[0];
            if (!peptideUnimod[1].isEmpty()) {
                foundUnimods.add(peptideUnimod[1]);
            }
        }
    }

    private void baseTOalphapept() {
        alphapept = base;

        boolean attemptCterm = cterm;
        for (int i = starts.size() - 1; i > -1; i--) {
            int start = starts.get(i);
            int end = ends.get(i);

            String[] peptideUnimod = PTMhandler.formatPeptideBaseToSpecific(
                    alphapept, start, end, "alphapept", foundUnimods, attemptCterm);
            attemptCterm = false;
            alphapept = peptideUnimod[0];
            if (!peptideUnimod[1].isEmpty()) {
                foundUnimods.add(peptideUnimod[1]);
            }
        }
    }

    private void baseToPredfullKoina() {
        predfull = base;

        boolean attemptCterm = cterm;
        for (int i = starts.size() - 1; i > -1; i--) {
            int start = starts.get(i);
            int end = ends.get(i);

            String[] peptideUnimod = PTMhandler.formatPeptideBaseToSpecific(
                    predfull, start, end, "predfull", foundUnimods, attemptCterm);
            attemptCterm = false;
            predfull = peptideUnimod[0];
            if (!peptideUnimod[1].isEmpty()) {
                foundUnimods.add(peptideUnimod[1]);
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

        boolean attemptCterm = cterm;
        for (int i = starts.size() - 1; i > -1; i--) {
            int start = starts.get(i);
            int end = ends.get(i);

            String[] peptideUnimod = PTMhandler.formatPeptideBaseToSpecific(
                    prosit, start, end, "prosit", foundUnimods, attemptCterm);
            attemptCterm = false;
            prosit = peptideUnimod[0];
            if (!peptideUnimod[1].isEmpty()) {
                foundUnimods.add(peptideUnimod[1]);
            }
        }

        //C is required to have carbamidomethyl
        for (int i = prosit.length() - 1; i > -1; i--) {
            char c = prosit.charAt(i);
            if (c == 'C') {
                boolean substituteC = false;

                if (i == prosit.length() - 1) {
                    substituteC = true;
                } else if (prosit.charAt(i + 1) != '[') {
                    substituteC = true;
                }

                if (substituteC) {
                    prosit = prosit.substring(0, i + 1) + "[UNIMOD:" +
                            PTMhandler.carbamidomethylationUnimod + "]" +
                            prosit.substring(i + 1);
                }
            }
        }
    }

    private void prositTOprositTMT() {
        //Prosit TMT assumes n-term are TMT-labeled
        if (!prosit.startsWith("[")) {
            String tmtLabel = "[UNIMOD:" + PTMhandler.tmtUnimod + "]-";
            prositTMT = tmtLabel + prosit;
        } else {
            prositTMT = prosit;
        }

        //S are not labeled
        prositTMT = prositTMT.replace("S[UNIMOD:" + tmtUnimod + "]", "S");
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

    //only sets base and basecharge
    public PeptideFormatter(String peptide, Object c, String format) {
        if (peptide.endsWith("cterm")) {
            peptide = peptide.substring(0, peptide.length() - 5);
            cterm = true;
        }
        charge = c.toString();

        switch(format) {
            case "pin":
                pinTObase(peptide);
                break;
            case "diann":
                diann = peptide;
                diannTObase(peptide);
                break;
            case "base":
                base = peptide;
                break;
            case "predfull":
                predfull = peptide;
                predfullTObase(peptide); //this version is for standalone
                koinaTObase(peptide); //this version for koina
                break;
            case "unispec":
                unispec = peptide;
                koinaTObase(peptide);
                break;
            case "prosit":
                prosit = peptide;
                koinaTObase(peptide);
                break;
            case "prosittmt":
                prositTMT = peptide;
                koinaTObase(peptide);
                break;
            case "ms2pip":
                ms2pip = peptide;
                koinaTObase(peptide);
                break;
            case "deeplc":
                deeplc = peptide;
                koinaTObase(peptide);
                break;
            case "alphapept":
                alphapept = peptide;
                alphapeptTObase(peptide);
                break;
            case "msp":
                mspTObase(peptide);
                break;
            case "pdeep3":
                pdeep3TObase(peptide);
                break;
            case "unimod.obo":
                unimodoboTObase(peptide);
                break;
        }

        findPTMlocations();
        baseCharge = base + "|" + charge;
    }

    public String getBaseCharge() {
        return baseCharge;
    }

    public String getBase() {
        return base;
    }

    public String getCharge() {
        return charge;
    }

    public String getDiann() {
        if (diann == null) {
            baseTOdiann();
        }
        return diann;
    }

    public String getProsit() {
        if (prosit == null) {
            baseTOprosit();
        }
        return prosit;
    }

    public String getStripped() {
        if (stripped == null) {
            baseTOstripped();
        }
        return stripped;
    }

    public String getPrositTMT() {
        if (prositTMT == null) {
            getProsit();
            prositTOprositTMT();
        }
        return prositTMT;
    }

    public String getPredfull() {
        if (predfull == null) {
            baseTOpredfull();
        }
        return predfull;
    }

    public String getPredfullKoina() {
        if (predfull == null) {
            baseToPredfullKoina();
        }
        return predfull;
    }


    public String getUnispec() {
        if (unispec == null) {
            baseTOunispec();
        }
        return unispec;
    }

    public String getMs2pip() {
        if (ms2pip == null) {
            baseTOms2pip();
        }
        return ms2pip;
    }

    public String getDeeplc() {
        if (deeplc == null) {
            baseTOdeeplc();
        }
        return deeplc;
    }

    public String getAlphapept() {
        if (alphapept == null) {
            baseTOalphapept();
        }
        return alphapept;
    }

    public String getModel(String model) { //model is whole url name
        switch(model.toLowerCase().split("_")[0]) {
            case "diann":
                return getDiann();
            case "prosit":
                if (model.contains("TMT")) {
                    return getPrositTMT();
                }
                return getProsit();
            case "prosittmt":
                return getPrositTMT();
            case "unispec":
                return getUnispec();
            case "ms2pip":
                return getMs2pip();
            case "deeplc":
                return getDeeplc();
            case "alphapept":
                return getAlphapept();
            case "predfull":
                return getPredfullKoina();
            default:
                //TODO implement all
                return "";
        }
    }

    public String getDlib() {
        if (dlib == null) {
            baseTOstripped();
            strippedTOdlib();
        }
        return dlib;
    }

    public String getMods() {
        if (mods == null) {
            baseToMods();
        }
        return mods;
    }

    public String getAlphapeptdeepMods() {
        if (alphapeptdeepMods == null) {
            baseToMods();
        }
        return alphapeptdeepMods;
    }

    public String getModPositions() {
        if (modPositions == null) {
            baseToMods();
        }
        return modPositions;
    }
}
