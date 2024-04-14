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

import static utils.Print.printInfo;

import org.apache.commons.lang3.ArrayUtils;
import umich.ms.fileio.exceptions.FileParsingException;

import java.io.*;
import java.util.*;
import java.util.concurrent.ExecutionException;

//TODO: which of these tools allows O abd U amino acids?
public class PinReader {
    String name; //used for resetting
    BufferedReader in;
    String[] header;
    String line;
    private String[] row;

    final int scanNumIdx;
    final int labelIdx;
    final int rankIdx;
    final int specIdx;
    final int pepIdx;
    final int eScoreIdx;
    final int rtIdx;

    HashMap<String, Integer> idxMap = new HashMap<>();
    boolean calcEvalue = false;

    MzmlReader mzml;
    private int length = -1;
    Double rtCutoff = Double.NaN;

    private static final Object lock = new Object(); // Create a lock object for synchronization

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
    }

    //reload from start
    public void reset() throws IOException {
        in = new BufferedReader(new FileReader(name));
        String line = in.readLine();
    }

    //get next row ready
    //do not split line for parallel processes
    public boolean next(boolean splitLine) throws IOException {
        String l = in.readLine();
        if (l != null) {
            line = l;
            if (splitLine) {
                row = line.split("\t");
            }
            if (! rtCutoff.isNaN()) {
                return !(getRT() > rtCutoff);
            }
            return true;
        }
        row = line.split("\t");
        return false;
    }

    public void close() throws IOException {
        in.close();
    }

    public void attachMzml(MzmlReader mzml) {
        this.mzml = mzml;
    }

    public String getColumn(String col) {
        int colNum;
        if (idxMap.containsKey(col)) {
            colNum = idxMap.get(col);
        } else {
            colNum = ArrayUtils.indexOf(header, col);
            idxMap.put(col, colNum);
        }
        return getRow()[colNum];
    }

    public int getLength() throws IOException {
        if (length == -1) {
                //https://stackoverflow.com/questions/453018/number-of-lines-in-a-file-in-java
            InputStream is = new BufferedInputStream(new FileInputStream(name));
            try {
                byte[] c = new byte[1024];

                int readChars = is.read(c);
                if (readChars == -1) {
                    // bail out if nothing to read
                    length = 0;
                }

                // make it easy for the optimizer to tune this loop
                int count = -1;
                while (readChars == 1024) {
                    for (int i = 0; i < 1024; ) {
                        if (c[i++] == '\n') {
                            ++count;
                        }
                    }
                    readChars = is.read(c);
                }

                // count remaining characters
                while (readChars != -1) {
                    for (int i = 0; i < readChars; ++i) {
                        if (c[i] == '\n') {
                            ++count;
                        }
                    }
                    readChars = is.read(c);
                }

                if (count <= 0) {
                    length = 1;
                } else {
                    length = count;
                }
            } finally {
                is.close();
            }
            printInfo(name + " has " + length + " PSMs");
        }
        return length;
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
            return String.valueOf(Math.pow(10, Double.parseDouble(row[eScoreIdx])));
        }
    }

    public float getRT() {return Float.parseFloat(row[rtIdx]);}

    //TODO: new algorithm?
    public void findWashGradientCutoff() throws IOException, FileParsingException {
        int bins = Constants.washGradientBins;
        //TODO read in RT directly from pin
        next(true);
        float minRT = getRT();
        float maxRT = 0;
        while (next(true)) {
            if (getRT() > maxRT) {
                maxRT = getRT();
            }
            if (getRT() < minRT) {
                minRT = getRT();
            }
        }
        reset();
        float RTrange = maxRT - minRT;

        int[] counts = new int[bins];
        while (next(true)) {
            int bin = Math.min(bins-1,Math.round(((mzml.getScanNumObject(getScanNum()).RT - minRT) / RTrange * bins)));
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
            printInfo("Washing gradient not detected. If there is one, you can try increasing the parameter " +
                    "washGradientBins (default 100) or manually setting it with the parameter rtCutoff.");
        }

        rtCutoff = (double) cutoffRT;
        printInfo("Setting RT cutoff at " + rtCutoff);

        //edit length
        length = getLength();
    }

    public String[] createPDeep2List() throws IOException {
        ArrayList<String> peps = new ArrayList<String>();
        while (next(true)) {
            PeptideFormatter pf = getPep();
            peps.add(pf.stripped + "\t" + pf.mods + "\t" + pf.charge);
        }
        return peps.toArray(new String[0]);
    }

    public String[] createPDeep3List() throws IOException {
        ArrayList<String> peps = new ArrayList<String>();
        while (next(true)) {
            PeptideFormatter pf = getPep();
            peps.add("." + "\t" + "." + "\t" + pf.stripped + "\t" + pf.mods + "\t" + pf.charge);
        }
        return peps.toArray(new String[0]);
    }

    public String[] createDeepMSPeptideList() throws IOException {
        ArrayList<String> peps = new ArrayList<String>();
        while (next(true)) {
            peps.add(getPep().stripped);
        }
        return peps.toArray(new String[0]);
    }

    public String[] createDiannList() throws IOException {
        ArrayList<String> peps = new ArrayList<String>();
        //TreeMap<Integer, Integer> modMap = new TreeMap<>(); //sorted for future use
        while (next(true)) {
            PeptideFormatter pf = getPep();
            peps.add(pf.diann + "\t" + pf.charge);
        }
        return peps.toArray(new String[0]);
    }

    public String[] createPredFullList(File mzmlFile, PinMzmlMatcher pmm) throws IOException, InterruptedException, ExecutionException, FileParsingException {
        ArrayList<String> peps = new ArrayList<String>();
        int fileI = 0;
        if (Constants.NCE.equals("") && pmm.mzmlReaders[fileI] == null) {
            mzml = new MzmlReader(mzmlFile.getCanonicalPath());
            pmm.mzmlReaders[fileI] = mzml;
        } else if (pmm.mzmlReaders[fileI] != null) {
            mzml = pmm.mzmlReaders[fileI];
        }
        while (next(true)) {
            PeptideFormatter pf = getPep();
            if (! pf.stripped.contains("O") && ! pf.stripped.contains("U") &&
                    ! pf.stripped.contains("Z") && ! pf.stripped.contains("B") &&
                    ! pf.stripped.contains("X")) {
                peps.add(pf.predfull + "\t" + pf.charge + "\t" + Constants.FragmentationType + "\t" + Constants.NCE);
            }
        }
        return peps.toArray(new String[0]);
    }

    public String[] createPrositList(File mzmlFile, PinMzmlMatcher pmm) throws IOException, InterruptedException, ExecutionException, FileParsingException {
        ArrayList<String> peps = new ArrayList<String>();
        int fileI = 0;
        if (Constants.NCE.equals("") && pmm.mzmlReaders[fileI] == null) {
            mzml = new MzmlReader(mzmlFile.getCanonicalPath());
            pmm.mzmlReaders[fileI] = mzml;
        } else if (pmm.mzmlReaders[fileI] != null) {
            mzml = pmm.mzmlReaders[fileI];
        }
        while (next(true)) {
            PeptideFormatter pf = getPep();
            peps.add(pf.prosit + "," + Constants.NCE + "," + pf.charge);
        }
        return peps.toArray(new String[0]);
    }

    public String[] createPrositTMTList(File mzmlFile, PinMzmlMatcher pmm) throws IOException, InterruptedException, ExecutionException, FileParsingException {
        ArrayList<String> peps = new ArrayList<String>();
        int fileI = 0;
        if (Constants.NCE.equals("") && pmm.mzmlReaders[fileI] == null) {
            mzml = new MzmlReader(mzmlFile.getCanonicalPath());
            pmm.mzmlReaders[fileI] = mzml;
        } else if (pmm.mzmlReaders[fileI] != null) {
            mzml = pmm.mzmlReaders[fileI];
        }
        while (next(true)) {
            PeptideFormatter pf = getPep();
            peps.add(pf.prosit + "," + Constants.NCE + "," + pf.charge + "," + Constants.FragmentationType);
        }
        return peps.toArray(new String[0]);
    }

    public String[] createAlphapeptdeepList(File mzmlFile, PinMzmlMatcher pmm) throws IOException, InterruptedException, ExecutionException, FileParsingException {
        ArrayList<String> peps = new ArrayList<String>();
        int fileI = 0;
        if (Constants.NCE.equals("") && pmm.mzmlReaders[fileI] == null) {
            mzml = new MzmlReader(mzmlFile.getCanonicalPath());
            pmm.mzmlReaders[fileI] = mzml;
        } else if (pmm.mzmlReaders[fileI] != null) {
            mzml = pmm.mzmlReaders[fileI];
        }
        while (next(true)) {
            PeptideFormatter pf = getPep();
            peps.add(pf.stripped + "," + pf.alphapeptdeepMods + "," + pf.modPositions + "," + pf.charge + "," +
                    Constants.NCE + "," + Constants.instrument + "," + pf.base);
        }
        return peps.toArray(new String[0]);
    }

    public String[] createFull() throws IOException {
        ArrayList<String> peps = new ArrayList<String>();
        while (next(true)) {
            PeptideFormatter pf = getPep();
            peps.add(pf.base + "\t" + pf.charge);

            if (Constants.features.contains("peptideCounts")) {
                if (Constants.peptideCounter.containsKey(pf.stripped)) {
                    Constants.peptideCounter.put(pf.stripped,
                            Constants.peptideCounter.get(pf.stripped) + 1);
                } else {
                    Constants.peptideCounter.put(pf.stripped, 1);
                }
            }
        }
        return peps.toArray(new String[0]);
    }

    public String[] createJSON(File mzmlFile, PinMzmlMatcher pmm, String modelFormat)
            throws IOException, InterruptedException, ExecutionException, FileParsingException {
        ArrayList<String> peps = new ArrayList<String>();
        int fileI = 0;
        for (File f : pmm.mzmlFiles) {
            if (f.toString().equals(mzmlFile.toString())) {
                break;
            }
            fileI++;
        }
        if (pmm.mzmlReaders[fileI] == null) {
            mzml = new MzmlReader(mzmlFile.getCanonicalPath());
            pmm.mzmlReaders[fileI] = mzml;
        } else if (pmm.mzmlReaders[fileI] != null) {
            mzml = pmm.mzmlReaders[fileI];
        }
        while (next(true)) {
            PeptideFormatter pf = getPep();
            if ((modelFormat.contains("Prosit") || modelFormat.contains("ms2pip") || modelFormat.contains("Deeplc"))
                    && pf.stripped.contains("U")) { // no peptides with U
                continue;
            }
            if (modelFormat.contains("ms2pip") && pf.stripped.length() > 30) { //peptide has length limit
                continue;
            }
            if (modelFormat.contains("Prosit") && modelFormat.contains("TMT") && pf.stripped.length() > 30) {
                continue;
            }

            String pep = pf.diann.replace("UniMod", "UNIMOD");
            if (pep.contains("[TMT]")) {
                pep = pep.replace("[TMT]", "[UNIMOD:737]");
            }

            if (pep.startsWith("[")) { //this is the case for all n term mods //TODO deal with c term mods
                int splitpoint = pep.indexOf("]");
                if (modelFormat.contains("Prosit")) {
                    if (modelFormat.contains("TMT") && pep.startsWith("[UNIMOD:737]")) {
                        pep = pep.substring(0, splitpoint + 1) + "-" + pep.substring(splitpoint + 1);
                    } else {
                        pep = pep.substring(splitpoint + 1);
                    }
                } else {
                    pep = pep.substring(0, splitpoint + 1) + "-" + pep.substring(splitpoint + 1);
                }
            }
            if (modelFormat.contains("Prosit") && modelFormat.contains("TMT")) {
                pep = pep.replace("S[UNIMOD:737]", "S");
                if (!pep.startsWith("[")) {
                    pep = "[UNIMOD:737]-" + pep;
                }
            }
            String sb = pep + "," + pf.charge + "," +  Constants.NCE + "," +
                    Constants.instrument + "," + Constants.FragmentationType;
            peps.add(sb);
        }
        return peps.toArray(new String[0]);
    }

    public LinkedList[] getTopPSMs(int num) throws IOException {
        LinkedList<String> PSMs = new LinkedList<>();
        LinkedList<Integer> scanNums = new LinkedList<>();

        //quicker way of getting top NCE
        PriorityQueue<Float> topEscores = new PriorityQueue<>(num, (a, b) -> Float.compare(b, a));
        while (next(true)) {
            float escore = Float.parseFloat(getEScore());
            if (topEscores.size() < num) {
                topEscores.offer(escore);
            } else if (escore < topEscores.peek()) {
                topEscores.poll();
                topEscores.offer(escore);
            }
        }
        reset();
        float eScoreCutoff = topEscores.peek();

        while (next(true) && PSMs.size() < num) {
            float escore = Float.parseFloat(getEScore());
            if (escore <= eScoreCutoff) {
                PeptideFormatter pf = getPep();
                PSMs.add(pf.getBaseCharge() + "," + pf.getDiann() + "," + pf.getStripped());
                scanNums.add(getScanNum());
            }
        }
        reset();
        return new LinkedList[]{PSMs, scanNums};
    }

    public LinkedList[] getDecoyPSMs(int num) throws IOException {
        LinkedList<String> PSMs = new LinkedList<>();
        LinkedList<Integer> scanNums = new LinkedList<>();

        //quicker way of getting top NCE
        PriorityQueue<Float> topEscores = new PriorityQueue<>(num, (a, b) -> Float.compare(b, a));
        while (next(true)) {
            if (getTD() == 0) {
                float escore = Float.parseFloat(getEScore());
                if (topEscores.size() < num) {
                    topEscores.offer(escore);
                } else if (escore < topEscores.peek()) {
                    topEscores.poll();
                    topEscores.offer(escore);
                }
            }
        }
        reset();
        float eScoreCutoff = topEscores.peek();

        while (next(true) && PSMs.size() < num) {
            float escore = Float.parseFloat(getEScore());
            if ((escore <= eScoreCutoff) && (getTD() == 0)) {
                PeptideFormatter pf = getPep();
                PSMs.add(pf.getBaseCharge() + "," + pf.getDiann() + "," + pf.getStripped());
                scanNums.add(getScanNum());
            }
        }
        reset();
        return new LinkedList[]{PSMs, scanNums};
    }
}