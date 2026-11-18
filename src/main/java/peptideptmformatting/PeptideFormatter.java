/*
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or
 * implied. See the License for the specific language governing
 * permissions and limitations under the License.
 */

package peptideptmformatting;

import koinaclasses.KoinaModelSpec;
import koinaclasses.KoinaModelSpec.EncodableMods;
import koinaclasses.KoinaModelSpec.ModMassMap;
import koinaclasses.KoinaModelSpec.ModRule;
import koinaclasses.KoinaModelSpec.PeptideFormat;
import utils.NumericUtils;

import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.concurrent.ConcurrentHashMap;

import static peptideptmformatting.PTMhandler.*;
import static utils.Print.printError;

//lots of different ways to format peptide string
public class PeptideFormatter {
    String base;
    private String diann;
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

    public static Set<String> foundUnimods = ConcurrentHashMap.newKeySet(); //collection of previously used unimod codes

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
        base = PTMhandler.formatPeptideSpecificToBase(peptide, unimodToModMassLimited, "[]");
    }
    private void koinaTObase(String peptide, Map<Integer, Double> unimodToMassMap) {
        base = PTMhandler.formatPeptideSpecificToBase(peptide, unimodToMassMap, "[]");
    }

    private void unimodoboTObase(String peptide) {
        base = PTMhandler.formatPeptideSpecificToBase(peptide, unimodOboToModMass, "[]");
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

    private void baseTOlibrarytsv() {
        librarytsv = base;

        boolean attemptCterm = cterm;
        for (int i = starts.size() - 1; i > -1; i--) {
            int start = starts.get(i);
            int end = ends.get(i);

            //TODO: extra logic for librarytsv periods for termini
            String[] peptideUnimod = PTMhandler.formatPeptideBaseToSpecific(
                    librarytsv, start, end, "librarytsv", foundUnimods, attemptCterm);
            attemptCterm = false;
            librarytsv = peptideUnimod[0];
            if (!peptideUnimod[1].isEmpty()) {
                foundUnimods.add(peptideUnimod[1]);
            }
        }
    }

    //writes the base peptide as a Koina model with this format, mass table and mods is sent it
    //encodable: the codes the model can be sent; a mod with no such code is left out
    //rememberFoundMods: each code found goes into foundUnimods (ModRule.remembersFoundMods)
    //droppedMods: if not null, gets each mod left out of the peptide, as modLabel writes it
    private String baseToKoina(PeptideFormat format, ModMassMap massMap, Set<String> preferredMods,
                               Set<String> allowedMods, EncodableMods encodable, boolean rememberFoundMods,
                               List<String> droppedMods) {
        String peptide = base;

        boolean attemptCterm = cterm;
        boolean ntermKept = false;
        boolean ctermKept = false;
        for (int i = starts.size() - 1; i > -1; i--) {
            int start = starts.get(i);
            int end = ends.get(i);

            String[] peptideUnimod = PTMhandler.formatPeptideBaseToKoina(
                    peptide, start, end, format, massMap, preferredMods, allowedMods, encodable, attemptCterm);
            peptide = peptideUnimod[0];
            if (peptideUnimod[1].isEmpty()) {
                if (droppedMods != null) {
                    droppedMods.add(modLabel(start, end, attemptCterm));
                }
            } else {
                if (rememberFoundMods) {
                    foundUnimods.add(peptideUnimod[1]);
                }
                if (start == 0) {
                    ntermKept = true;
                } else if (attemptCterm && peptideUnimod[1].charAt(0) == '[') {
                    ctermKept = true;
                }
            }
            attemptCterm = false;
        }

        if (format.termModsOnResidues()) {
            peptide = termModsOnResidues(peptide, ntermKept, ctermKept, droppedMods);
        }

        if (format.camOnBareCys()) {
            //C is required to have carbamidomethyl
            for (int i = peptide.length() - 1; i > -1; i--) {
                char c = peptide.charAt(i);
                if (c == 'C') {
                    boolean substituteC = false;

                    if (i == peptide.length() - 1) {
                        substituteC = true;
                    } else if (peptide.charAt(i + 1) != '[') {
                        substituteC = true;
                    }

                    if (substituteC) {
                        peptide = peptide.substring(0, i + 1) + "[UNIMOD:" +
                                PTMhandler.carbamidomethylationUnimod + "]" +
                                peptide.substring(i + 1);
                    }
                }
            }
        }

        if (format.tmtOnNterm()) {
            //Prosit TMT assumes n-term are TMT-labeled
            if (!peptide.startsWith("[")) {
                String tmtLabel = "[UNIMOD:" + PTMhandler.tmtUnimod + "]-";
                peptide = tmtLabel + peptide;
            }

            //S are not labeled
            peptide = peptide.replace("S[UNIMOD:" + tmtUnimod + "]", "S");
        }
        return peptide;
    }

    /**
     * Moves the terminal mods of a peptide written for a Koina model onto the terminal residues
     * (PeptideFormat.termModsOnResidues): "[UNIMOD:1]PEP" becomes "P[UNIMOD:1]EP". The C-term mod is already
     * written after the last residue. A terminal residue that keeps a mod of its own keeps only that one: the
     * terminal mod is taken out and added to droppedMods (D24).
     *
     * @param peptide   the peptide as written, with its N-term mod first if ntermKept, and its C-term mod last if
     *                  ctermKept
     * @param ntermKept whether the base peptide's N-term mod was written
     * @param ctermKept whether the base peptide's C-term mod was written
     */
    private String termModsOnResidues(String peptide, boolean ntermKept, boolean ctermKept,
                                      List<String> droppedMods) {
        if (ctermKept) {
            int open = peptide.lastIndexOf('[');
            if (open > 0 && peptide.charAt(open - 1) == ']') { //the last residue has a mod of its own
                peptide = peptide.substring(0, open);
                if (droppedMods != null) {
                    int last = starts.size() - 1;
                    droppedMods.add(modLabel(starts.get(last), ends.get(last), true));
                }
            }
        }
        if (ntermKept) {
            int close = peptide.indexOf(']');
            int residue = close + 1;
            if (residue + 1 < peptide.length() && peptide.charAt(residue + 1) == '[') { //residue 1 has a mod
                peptide = peptide.substring(residue);
                if (droppedMods != null) {
                    droppedMods.add(modLabel(starts.get(0), ends.get(0), false));
                }
            } else {
                peptide = peptide.charAt(residue) + peptide.substring(0, residue) + peptide.substring(residue + 1);
            }
        }
        return peptide;
    }

    /**
     * The mod at [start, end] of the base peptide as MSFragger writes it in a pin: its residue, or n for an N-term
     * and c for a C-term mod, then its mass in brackets, e.g. S[79.9663] or n[42.0106].
     *
     * @param cterm whether a mod on the last residue is read as a C-term mod
     */
    private String modLabel(int start, int end, boolean cterm) {
        //the residue before the mod, past any mods written before it
        int i = start - 1;
        while (i >= 0 && base.charAt(i) == ']') {
            i = base.lastIndexOf('[', i) - 1;
        }
        String site;
        if (i < 0) {
            site = "n";
        } else if (cterm && end == base.length() - 1) {
            site = "c";
        } else {
            site = String.valueOf(base.charAt(i));
        }
        return site + base.substring(start, end + 1);
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
            Map<Double, HashSet<String>> tt = PTMhandler.aamassToAlphapeptdeep.subMap(doubleModMass - 0.001, false, doubleModMass + 0.001, false);

            if (tt.isEmpty()) {
                printError("There is an unknown modification with mass " + doubleModMass);
                System.exit(1);
            }

            for (Map.Entry<Double, HashSet<String>> e : tt.entrySet()) {
                possibleMods.addAll(e.getValue());
            }

            int position = starts.get(i) - newEnds.get(i) + positions.get(i);
            if (i > 0) {
                position -= 1;
            }
            String aa;
            if (position == 0) {
                aa = "Any_N-term"; //might need to change
            } else {
                if ((cterm) && (i == numMods - 1)) {
                    aa = "Any_C-term";
                    position = -1;
                } else {
                    aa = stripped.substring(position - 1, position);
                }
            }
            positions.add(position);

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
        peptide = readCterm(peptide);
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
            case "prosit": //the local Prosit predictions (MspReader); a Koina model's peptides take its spec
                koinaTObase(peptide, unimodToModMassLimited);
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

    //a peptide as a Koina model is sent it, read back with the model's mass table. Only sets base and basecharge
    public PeptideFormatter(String peptide, Object c, KoinaModelSpec spec) {
        peptide = readCterm(peptide);
        charge = c.toString();
        if (spec.peptideFormat().termModsOnResidues()) {
            peptide = ntermModOffFirstResidue(peptide, spec.modMassMap().aaUnimodToMass());
        }
        koinaTObase(peptide, spec.modMassMap().unimodToMass());

        findPTMlocations();
        baseCharge = base + "|" + charge;
    }

    /**
     * A peptide written with its N-term mod on the first residue (PeptideFormat.termModsOnResidues), with that mod
     * put back before the residue: "A[UNIMOD:1]DK" becomes "[UNIMOD:1]ADK". A mod on the first residue is read as
     * the N-term mod when the mass table has its code as an N-term mod and not as a mod of that residue; "K[UNIMOD:1]"
     * stays acetyl on K. Koina reads both the same, so either reading names the same prediction.
     */
    static String ntermModOffFirstResidue(String peptide, Map<String, Double> aaUnimodToMass) {
        if (peptide.length() < 3 || peptide.charAt(1) != '[') {
            return peptide;
        }
        int close = peptide.indexOf(']');
        String mod = peptide.substring(2, close); //UNIMOD:1
        String code = mod.substring(mod.indexOf(':') + 1);
        char residue = peptide.charAt(0);
        if (aaUnimodToMass.containsKey(residue + code) || !aaUnimodToMass.containsKey("[" + code)) {
            return peptide;
        }
        return "[" + mod + "]" + residue + peptide.substring(close + 1);
    }

    //the peptide without the "cterm" a C-term mod is marked with, if it has one
    private String readCterm(String peptide) {
        if (peptide.endsWith("cterm")) {
            peptide = peptide.substring(0, peptide.length() - 5);
            cterm = true;
        }
        return peptide;
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

    //the local Prosit format (PinReader.createPrositList and createPrositTMTList, MspReader), with these mods
    //preferred. A Koina model's peptides are written by getModel
    public String getProsit(Set<String> uniqMods) {
        return baseToKoina(PeptideFormat.PROSIT, ModMassMap.LIMITED, uniqMods, prositAAMods, EncodableMods.ANY, false,
                null);
    }

    public String getStripped() {
        baseTOstripped();
        return stripped;
    }

    //the local PrositTMT format (MspReader)
    public String getPrositTMT() {
        return baseToKoina(PeptideFormat.PROSIT_TMT, ModMassMap.LIMITED, prositTmtAAMods, prositAAMods,
                EncodableMods.ANY, false, null);
    }

    public String getLibrarytsv() {
        if (librarytsv == null) {
            baseTOlibrarytsv();
        }
        return librarytsv;
    }

    /**
     * The peptide as a Koina model is sent it: in the model's peptide format, with the mods its rule allows written
     * as UniMod codes and the others left out. For a rule that remembers the codes it finds, this adds them to
     * {@link #foundUnimods}.
     */
    public String getModel(KoinaModelSpec spec) {
        return getModel(spec, null);
    }

    /**
     * {@link #getModel(KoinaModelSpec)}, also adding to {@code droppedMods} each mod left out of the peptide, as
     * its residue (n for an N-term mod, c for a C-term mod) and its mass in brackets, e.g. S[79.9663]. Its mass is
     * missing from the peptide sent, and transferKoinaPreds adds it back to the fragments that hold it.
     *
     * @param droppedMods the list the mods left out are added to; null if they are not wanted
     */
    public String getModel(KoinaModelSpec spec, List<String> droppedMods) {
        ModRule rule = spec.modRule();
        return baseToKoina(spec.peptideFormat(), spec.modMassMap(), rule.preferredMods(), rule.allowedMods(),
                spec.encodableMods(), rule.remembersFoundMods(), droppedMods);
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
