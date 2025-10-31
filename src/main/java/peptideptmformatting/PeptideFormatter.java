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

import utils.NumericUtils;

import java.util.ArrayList;
import java.util.HashSet;
import java.util.Map;
import java.util.Set;
import java.util.concurrent.ConcurrentHashMap;

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
    private String librarytsv;
    private String stripped;
    String baseCharge;
    private String dlib;
    private String mods = "";
    private String alphapeptdeepMods = "";
    private static ConcurrentHashMap<String, Integer> alphapeptdeepModsSet = new ConcurrentHashMap<>();
    private String modPositions = "";

    private ArrayList<Integer> starts = new ArrayList<>();
    private ArrayList<Integer> ends = new ArrayList<>();

    public String charge;

    public static HashSet<String> foundUnimods = new HashSet<>(); //collection of previously used unimod codes
    //TODO: Can store in a thread-safe set

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

    private void pinTObase(String peptide, boolean adjacentAA) {
        //remove AA and period at beginning and end
        if (adjacentAA) {
            peptide = peptide.substring(2, peptide.length() - 2);
        }

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
        base = PTMhandler.formatPeptideSpecificToBase(peptide, unimodToModMassLimited);
    }
    private void koinaTObase(String peptide, Map<Integer, Double> unimodToMassMap) {
        base = PTMhandler.formatPeptideSpecificToBase(peptide, unimodToMassMap);
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
                    diann, start, end, "diann", diannAAMods, attemptCterm);
            attemptCterm = false;
            diann = peptideUnimod[0];
        }

        //special TMT formatting
        diann = diann.replaceAll("UniMod:" + PTMhandler.tmtUnimod, "TMT");

        //nterm mod and TMT on AA1 not allowed
        if (diann.startsWith("[")) {
            String[] diannsplit = diann.split("]");
            if (diannsplit.length > 2) {
                if (diannsplit[1].substring(1).startsWith("[TMT")) {
                    diann = diann.substring(diann.indexOf("]") + 1);
                }
            }
        }
    }

    private void baseTOunispec() {
        unispec = base;

        boolean attemptCterm = cterm;
        for (int i = starts.size() - 1; i > -1; i--) {
            int start = starts.get(i);
            int end = ends.get(i);

            String[] peptideUnimod = PTMhandler.formatPeptideBaseToSpecific(
                    unispec, start, end, "unispec", unispecAAMods, attemptCterm);
            attemptCterm = false;
            unispec = peptideUnimod[0];
        }
    }

    private void baseTOms2pip() {
        ms2pip = base;

        boolean attemptCterm = cterm;
        for (int i = starts.size() - 1; i > -1; i--) {
            int start = starts.get(i);
            int end = ends.get(i);

            String[] peptideUnimod = PTMhandler.formatPeptideBaseToSpecific(
                    ms2pip, start, end, "ms2pip", ms2pipAAMods, attemptCterm);
            attemptCterm = false;
            ms2pip = peptideUnimod[0];
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

    private void baseTOlibrarytsv() {
        librarytsv = base;

        boolean attemptCterm = cterm;
        for (int i = starts.size() - 1; i > -1; i--) {
            int start = starts.get(i);
            int end = ends.get(i);

            String[] peptideUnimod = PTMhandler.formatPeptideBaseToSpecific(
                    librarytsv, start, end, "librarytsv", foundUnimods, attemptCterm);
            attemptCterm = false;
            librarytsv = peptideUnimod[0];
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
                    predfull, start, end, "predfull", predfullAAMods, attemptCterm);
            attemptCterm = false;
            predfull = peptideUnimod[0].replace("[UNIMOD:4]", "");
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

    private void baseTOprosit(HashSet<String> uniqMods) {
        prosit = base;

        boolean attemptCterm = cterm;
        for (int i = starts.size() - 1; i > -1; i--) {
            int start = starts.get(i);
            int end = ends.get(i);

            String[] peptideUnimod = PTMhandler.formatPeptideBaseToSpecific(
                    prosit, start, end, "prosit", uniqMods, attemptCterm);
            attemptCterm = false;
            prosit = peptideUnimod[0];
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
            HashSet<String> possibleMods = new HashSet<>();
            Map<Double, HashSet<String>> tt = PTMhandler.aamassToAlphapeptdeep.subMap(doubleModMass - 0.0001, false, doubleModMass + 0.0001, false);

            if (tt.isEmpty()) {
                printError("There is an unknown modification with mass " + doubleModMass +
                        ". Please provide PTM info via additionalMods param in --paramsList.");
                System.exit(1);
            }

            for (Map.Entry<Double, HashSet<String>> e : tt.entrySet()) {
                possibleMods.addAll(e.getValue());
            }

            int position = starts.get(i) - newEnds.get(i) + positions.get(i);
            if (i > 0) {
                position -= 1;
            }
            positions.add(position);
            String aa;
            if (position == 0) {
                aa = "Any_N-term"; //might need to change
            } else {
                aa = stripped.substring(position - 1, position);
            }

            //check that the ptmName@aa is accepted by alphapeptdeep
            for (String name : possibleMods) {
                if (PTMhandler.alphapeptdeepModNames.contains(name + "@" + aa)) {
                    String modinfo = position + "," + name + "[" + aa + "]" + ";";
                    String alphapeptdeepModinfo = PTMhandler.writeOutAlphapeptdeepModNames.get(name + "@" + aa) + ";";

                    mods = mods + modinfo;
                    alphapeptdeepMods = alphapeptdeepMods + alphapeptdeepModinfo;
                    alphapeptdeepModsSet.put(alphapeptdeepModinfo, 0);
                    modPositions = modPositions + position + ";";
                    break;
                }
            }
        }
        if (!mods.isEmpty()) {
            mods = mods.substring(0, mods.length() - 1);
            alphapeptdeepMods = alphapeptdeepMods.substring(0, alphapeptdeepMods.length() - 1);
            //modPositions = modPositions.substring(0, modPositions.length() - 1);
        }
    }

    public HashSet<String> getAlphapeptdeepModsSet() {
        HashSet<String> keys = new HashSet<>(alphapeptdeepModsSet.keySet());
        return keys;
    }

    //this version just to get access to fields
    public PeptideFormatter() {}

    //only sets base and basecharge
    public PeptideFormatter(String peptide, Object c, String format) {
        if (peptide.endsWith("cterm")) {
            peptide = peptide.substring(0, peptide.length() - 5);
            cterm = true;
        }
        charge = c.toString();

        switch(format) {
            case "pin":
                pinTObase(peptide, true);
                break;
            case "apdpred":
                pinTObase(peptide, false);
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
                koinaTObase(peptide, unimodToModMassLimited); //this version for koina
                break;
            case "unispec":
                unispec = peptide;
                koinaTObase(peptide, unimodToModMassLimited);
                break;
            case "prosit":
            case "prosit_cit":
                prosit = peptide;
                koinaTObase(peptide, unimodToModMassLimited);
                break;
            case "prosittmt":
                prositTMT = peptide;
                koinaTObase(peptide, unimodToModMassLimited);
                break;
            case "ms2pip":
                ms2pip = peptide;
                koinaTObase(peptide, unimodToModMassLimited);
                break;
            case "deeplc":
            case "im2deep":
                deeplc = peptide;
                koinaTObase(peptide, unimodToModMassAll);
                break;
            case "alphapept":
                alphapept = peptide;
                koinaTObase(peptide, unimodToModMassAll);
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
        baseTOdiann();
        return diann;
    }

    public String getProsit(HashSet<String> uniqMods) {
        baseTOprosit(uniqMods);
        return prosit;
    }

    public String getStripped() {
        baseTOstripped();
        return stripped;
    }

    public String getPrositTMT() {
        getProsit(prositTmtAAMods);
        prositTOprositTMT();
        return prositTMT;
    }

    public String getPredfull() {
        baseTOpredfull();
        return predfull;
    }

    public String getPredfullKoina() {
        baseToPredfullKoina();
        return predfull;
    }

    public String getUnispec() {
        baseTOunispec();
        return unispec;
    }

    public String getMs2pip() {
        baseTOms2pip();
        return ms2pip;
    }

    public String getDeeplc() {
        baseTOdeeplc();
        return deeplc;
    }

    public String getAlphapept() {
        baseTOalphapept();
        return alphapept;
    }

    public String getLibrarytsv() {
        if (librarytsv == null) {
            baseTOlibrarytsv();
        }
        return librarytsv;
    }

    public String getModel(String model) { //model is whole url name
        switch(model.toLowerCase().split("_")[0]) {
            case "diann":
                return getDiann();
            case "prosit":
                if (model.contains("TMT")) {
                    return getPrositTMT();
                } else if (model.contains("_cit")) {
                    getProsit(prositCitAAMods);
                }
                return getProsit(prositAAMods);
            case "prosittmt":
                return getPrositTMT();
            case "unispec":
                return getUnispec();
            case "ms2pip":
                return getMs2pip();
            case "deeplc":
            case "im2deep":
                return getDeeplc();
            case "alphapept":
                return getAlphapept();
            case "predfull":
                return getPredfullKoina();
            default:
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
        if (alphapeptdeepMods.isEmpty()) {
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
