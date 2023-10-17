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

import org.apache.commons.lang3.ArrayUtils;
import org.checkerframework.checker.units.qual.A;
import org.checkerframework.checker.units.qual.N;
import umich.ms.fileio.exceptions.FileParsingException;
import umontreal.ssj.probdist.EmpiricalDist;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.*;
import java.util.concurrent.ConcurrentHashMap;
import java.util.concurrent.ExecutionException;

//TODO: which of these tools allows O abd U amino acids?
public class PinReader {
    String name; //used for resetting
    BufferedReader in;
    String[] header;
    private String[] row;

    int scanNumIdx;
    int labelIdx;
    int rankIdx;
    int specIdx;
    int pepIdx;
    int eScoreIdx;
    int rtIdx;
    private boolean calcEvalue = false;

    MzmlReader mzml;
    int length;
    Double rtCutoff = Double.NaN;

    public PinReader(String pin) throws IOException {
        name = pin;
        in = new BufferedReader(new FileReader(name));
        String line = in.readLine();
        header = line.split("\t");

        //set column indices
        scanNumIdx = ArrayUtils.indexOf(header, "ScanNr");
        labelIdx = ArrayUtils.indexOf(header, "Label");
        rankIdx = ArrayUtils.indexOf(header, "rank");
        specIdx = ArrayUtils.indexOf(header, "SpecId");
        pepIdx = ArrayUtils.indexOf(header, "Peptide");
        rtIdx = ArrayUtils.indexOf(header, "retentiontime");
        if (Arrays.asList(header).contains("log10_evalue")) {
            eScoreIdx = ArrayUtils.indexOf(header, "log10_evalue"); //DDA
        } else {
            eScoreIdx = ArrayUtils.indexOf(header, "hyperscore"); //DIA
            calcEvalue = true;
        }

        getLength();
    }

    //reload from start
    public void reset() throws IOException {
        in = new BufferedReader(new FileReader(name));
        String line = in.readLine();
    }

    //get next row ready
    public boolean next() throws IOException {
        String line = in.readLine();
        if (line != null) {
            row = line.split("\t");
            if (! rtCutoff.isNaN()) {
                return !(getRT() > rtCutoff);
            }
            return true;
        }
        //in.close();
        return false;
    }

    public void close() throws IOException {
        in.close();
    }

    public void attachMzml(MzmlReader mzml) {
        this.mzml = mzml;
    }

    public void getLength() throws IOException {
        length = 0;
        while (next()) {
            length += 1;
        }
        reset();
    }

    public String[] getRow() {return row;}

    public PeptideFormatter getPep() {
        String[] periodSplit = row[specIdx].split("\\.");
        return new PeptideFormatter(row[pepIdx], periodSplit[periodSplit.length - 1].split("_")[0], "pin");
    }

    public int getTD() {return Math.max(0, Integer.parseInt(row[labelIdx]));} //just leave as -1?

    public int getScanNum() {return Integer.parseInt(row[scanNumIdx]);}

    public int getRank() {
        try {
            return Integer.parseInt(row[rankIdx]);
        } catch (Exception e) {
            String[] specIdxSplit = row[specIdx].split("_");
            return Integer.parseInt(specIdxSplit[specIdxSplit.length - 1]);
        }
    }

    //public String getEScore() {return String.valueOf(Math.pow(10, Double.parseDouble(row[eScoreIdx])));}
    public String getEScore() {
        if (calcEvalue) {
            return String.valueOf(Math.exp(15.0 - Double.parseDouble(row[eScoreIdx])));
        } else {
            return row[eScoreIdx];
        }
    }

    public float getRT() {return Float.parseFloat(row[rtIdx]);}

    public void findWashGradientCutoff(ConcurrentHashMap<String, PredictionEntry> preds) throws IOException {
        int bins = Constants.washGradientBins;
        //TODO read in RT directly from pin
        next();
        float minRT = getRT();
        float maxRT = 0;
        while (next()) {
            maxRT = getRT();
        }
        reset();
        float RTrange = maxRT - minRT;

        int[] counts = new int[bins];
        while (next()) {
            int bin = Math.min(bins-1,Math.round(((mzml.scanNumberObjects.get(getScanNum()).RT - minRT) / RTrange * bins)));
            counts[bin] += 1;
        }
        reset();
        int maxAt = 0;
        for (int i = 0; i < counts.length; i++) {
            maxAt = counts[i] > counts[maxAt] ? i : maxAt;
        }
        int zeroAt = -1;
        for (int i = maxAt + 1; i < counts.length; i++) {
            if (counts[i] == 0) {
                zeroAt = i;
                break;
            }
        }

        float cutoffRT = maxRT * 2;
        if (zeroAt != -1) {
            cutoffRT = minRT + zeroAt * RTrange / bins;
        } else {
            System.out.println("Washing gradient not detected. If there is one, you can try increasing the parameter " +
                    "washGradientBins (default 100) or manually setting it with the parameter rtCutoff.");
        }

        rtCutoff = (double) cutoffRT;
        System.out.println("Setting RT cutoff at " + rtCutoff);

        //edit length
        getLength();
    }

    public String[] createPDeep2List() throws IOException {
        ArrayList<String> peps = new ArrayList<String>();
        while (next()) {
            PeptideFormatter pf = getPep();
            peps.add(pf.stripped + "\t" + pf.mods + "\t" + pf.charge);
        }
        return peps.toArray(new String[0]);
    }

    public String[] createPDeep3List() throws IOException {
        ArrayList<String> peps = new ArrayList<String>();
        while (next()) {
            PeptideFormatter pf = getPep();
            peps.add("." + "\t" + "." + "\t" + pf.stripped + "\t" + pf.mods + "\t" + pf.charge);
        }
        return peps.toArray(new String[0]);
    }

    public String[] createDeepMSPeptideList() throws IOException {
        ArrayList<String> peps = new ArrayList<String>();
        while (next()) {
            peps.add(getPep().stripped);
        }
        return peps.toArray(new String[0]);
    }

    public String[] createDiannList() throws IOException {
        ArrayList<String> peps = new ArrayList<String>();
        //TreeMap<Integer, Integer> modMap = new TreeMap<>(); //sorted for future use
        while (next()) {
            PeptideFormatter pf = getPep();
            peps.add(pf.diann + "\t" + pf.charge);
        }
        return peps.toArray(new String[0]);
    }

    public String[] createPredFullList(File mzmlFile) throws IOException, InterruptedException, ExecutionException, FileParsingException {
        ArrayList<String> peps = new ArrayList<String>();
        if (Constants.NCE.equals("")) {
            mzml = new MzmlReader(mzmlFile.getCanonicalPath());
        }
        while (next()) {
            PeptideFormatter pf = getPep();
            if (! pf.stripped.contains("O") && ! pf.stripped.contains("U") &&
                    ! pf.stripped.contains("Z") && ! pf.stripped.contains("B") &&
                    ! pf.stripped.contains("X")) {
                String NCE = getNCE();
                peps.add(pf.predfull + "\t" + pf.charge + "\t" + Constants.FragmentationType + "\t" + NCE);
            }
        }
        return peps.toArray(new String[0]);
    }

    public String[] createPrositList(File mzmlFile) throws IOException, InterruptedException, ExecutionException, FileParsingException {
        ArrayList<String> peps = new ArrayList<String>();
        if (Constants.NCE.equals("")) {
            mzml = new MzmlReader(mzmlFile.getCanonicalPath());
        }
        while (next()) {
            PeptideFormatter pf = getPep();
            String NCE = getNCE();
            peps.add(pf.prosit + "," + NCE + "," + pf.charge);
        }
        return peps.toArray(new String[0]);
    }

    public String[] createPrositTMTList(File mzmlFile) throws IOException, InterruptedException, ExecutionException, FileParsingException {
        ArrayList<String> peps = new ArrayList<String>();
        if (Constants.NCE.equals("")) {
            mzml = new MzmlReader(mzmlFile.getCanonicalPath());
        }
        while (next()) {
            PeptideFormatter pf = getPep();
            String NCE = getNCE();
            peps.add(pf.prosit + "," + NCE + "," + pf.charge + "," + Constants.FragmentationType);
        }
        return peps.toArray(new String[0]);
    }

    public String[] createAlphapeptdeepList(File mzmlFile) throws IOException, InterruptedException, ExecutionException, FileParsingException {
        ArrayList<String> peps = new ArrayList<String>();
        if (Constants.NCE.equals("")) {
            mzml = new MzmlReader(mzmlFile.getCanonicalPath());
        }
        Constants.instrument = getInstrument();
        while (next()) {
            PeptideFormatter pf = getPep();
            String NCE = getNCE();
            peps.add(pf.stripped + "," + pf.alphapeptdeepMods + "," + pf.modPositions + "," + pf.charge + "," +
                    NCE + "," + Constants.instrument + "," + pf.base);
        }
        return peps.toArray(new String[0]);
    }

    public String[] createFull() throws IOException {
        ArrayList<String> peps = new ArrayList<String>();
        //TreeMap<Integer, Integer> modMap = new TreeMap<>(); //sorted for future use
        while (next()) {
            PeptideFormatter pf = getPep();
            peps.add(pf.base + "\t" + pf.charge);
        }
        return peps.toArray(new String[0]);
    }

    public String[] createJSON(File mzmlFile, String modelFormat)
            throws IOException, InterruptedException, ExecutionException, FileParsingException {
        ArrayList<String> peps = new ArrayList<String>();
        if (Constants.NCE.equals("")) {
            mzml = new MzmlReader(mzmlFile.getCanonicalPath());
        }
        Constants.instrument = getInstrument();
        while (next()) {
            PeptideFormatter pf = getPep();
            if ((modelFormat.contains("Prosit") || modelFormat.contains("ms2pip"))
                    && pf.stripped.contains("U")) { // no peptides with U
                continue;
            }
            if (modelFormat.contains("ms2pip") && pf.stripped.length() > 30) { //peptide has length limit
                continue;
            }
            String fragmentation = "";
            Set<String> fragTypes = mzml.scanNumberObjects.get(getScanNum()).NCEs.keySet();
            if (fragTypes.contains("CID")) {
                fragmentation = "CID";
            } else {
                fragmentation = "HCD";
            }
            String NCE = getNCE(fragmentation);
            String pep = pf.diann.replace("UniMod", "UNIMOD");
            if (pep.contains("[TMT]")) {
                pep = pep.replace("[TMT]", "[UNIMOD:737]");
            }

            if (pep.startsWith("[")) { //this is the case for all n term mods //TODO deal with c term mods
                int splitpoint = pep.indexOf("]");
                if (modelFormat.contains("Prosit")) {
                    pep = pep.substring(splitpoint + 1);
                } else {
                    pep = pep.substring(0, splitpoint + 1) + "-" + pep.substring(splitpoint + 1);
                }
            }
            StringBuilder sb = new StringBuilder();
            sb.append(pep).append(",").append(pf.charge).append(",").append(NCE).append(",").
                    append(Constants.instrument).append(",").append(fragmentation);
            peps.add(sb.toString());
        }
        return peps.toArray(new String[0]);
    }

    private String getNCE() {
        if (Constants.NCE.equals("")) {
            return String.valueOf(mzml.scanNumberObjects.get(getScanNum()).NCEs.get(Constants.FragmentationType));
        } else {
            return Constants.NCE;
        }

    }

    private String getNCE(String frag) {
        if (Constants.NCE.equals("")) {
            return String.valueOf(mzml.scanNumberObjects.get(getScanNum()).NCEs.get(frag));
        } else {
            return Constants.NCE;
        }

    }

    //TODO: support for astral model?
    private String getInstrument() {
        HashSet<String> LumosKeys = new HashSet<>(Arrays.asList("LTQ", "Lumos", "Fusion", "Elite", "Velos", "Eclipse", "Tribrid"));
        HashSet<String> QEKeys = new HashSet<>(Arrays.asList("QE", "Exactive", "Exploris"));
        HashSet<String> SciexTOFKeys = new HashSet<>(Arrays.asList("Sciex", "TripleTOF"));
        HashSet<String> timsTOFKeys = new HashSet<>(Arrays.asList("flight"));

        if (Constants.instrument.equals("")) {
            String model = mzml.scans.getRunInfo().getDefaultInstrument().getModel();
            String analyzer = mzml.scans.getRunInfo().getDefaultInstrument().getAnalyzer();
            for (String k : LumosKeys) {
                if (model.contains(k) || analyzer.contains(k)) {
                    return "Lumos";
                }
            }
            for (String k : QEKeys) {
                if (model.contains(k) || analyzer.contains(k)) {
                    return "QE";
                }
            }
            for (String k : SciexTOFKeys) {
                if (model.contains(k) || analyzer.contains(k)) {
                    return "SciexTOF";
                }
            }
            for (String k : timsTOFKeys) {
                if (model.contains(k) || analyzer.contains(k)) {
                    return "timsTOF";
                }
            }
            System.out.println("Could not detect instrument type. Setting to Lumos");
            return "Lumos"; //default if nothing found
        } else {
            return Constants.instrument;
        }
    }
}