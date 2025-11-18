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

package readers.datareaders;

import allconstants.Constants;
import allconstants.NceConstants;
import org.apache.commons.lang3.ArrayUtils;
import peptideptmformatting.PTMhandler;
import peptideptmformatting.PeptideFormatter;
import peptideptmformatting.PeptideSkipper;
import umich.ms.fileio.exceptions.FileParsingException;
import utils.Print;

import java.io.BufferedInputStream;
import java.io.BufferedReader;
import java.io.FileInputStream;
import java.io.FileReader;
import java.io.IOException;
import java.io.InputStream;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;
import java.util.LinkedList;
import java.util.PriorityQueue;
import java.util.concurrent.ExecutionException;

import static utils.InstrumentUtils.mapInstrumentToModelSpecific;
import static utils.Print.printInfo;

public class PinReader {
    public String name; //used for resetting
    BufferedReader in;
    public String[] header;
    public String line;
    private String[] row;

    public final int scanNumIdx;
    final int labelIdx;
    final int rankIdx;
    public final int specIdx;
    public final int pepIdx;
    final int eScoreIdx;
    final int rtIdx;

    HashMap<String, Integer> idxMap = new HashMap<>();
    boolean calcEvalue = false;

    MzmlReader mzml;
    private int length = -1;
    public Double rtCutoff = Double.NaN;

    private static final Object lock = new Object(); // Create a lock object for synchronization

    public PinReader(String pin) throws IOException {
        name = pin;
        in = new BufferedReader(new FileReader(name));
        line = in.readLine();
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
                    printInfo(name + " has 0 PSMs");
                } else {
                    length = count;
                    printInfo(name + " has " + length + " PSMs");
                }
            } finally {
                is.close();
            }
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
            peps.add(pf.getStripped() + "\t" + pf.getMods() + "\t" + pf.getCharge());
        }
        return peps.toArray(new String[0]);
    }

    public String[] createPDeep3List() throws IOException {
        ArrayList<String> peps = new ArrayList<String>();
        while (next(true)) {
            PeptideFormatter pf = getPep();
            peps.add("." + "\t" + "." + "\t" + pf.getStripped() + "\t" + pf.getMods() + "\t" + pf.getCharge());
        }
        return peps.toArray(new String[0]);
    }

    public String[] createDeepMSPeptideList() throws IOException {
        ArrayList<String> peps = new ArrayList<String>();
        while (next(true)) {
            peps.add(getPep().getStripped());
        }
        return peps.toArray(new String[0]);
    }

    public String[] createDiannList() throws IOException {
        ArrayList<String> peps = new ArrayList<String>();
        //TreeMap<Integer, Integer> modMap = new TreeMap<>(); //sorted for future use
        while (next(true)) {
            PeptideFormatter pf = getPep();
            peps.add(pf.getDiann() + "\t" + pf.getCharge());
        }
        return peps.toArray(new String[0]);
    }

    public String[] createPrositList()
            throws IOException, InterruptedException, ExecutionException, FileParsingException {
        ArrayList<String> peps = new ArrayList<String>();
        while (next(true)) {
            PeptideFormatter pf = getPep();
            peps.add(pf.getProsit(PTMhandler.prositAAMods) + "," + NceConstants.getNCE() + "," + pf.getCharge());
        }
        return peps.toArray(new String[0]);
    }

    public String[] createPrositTMTList()
            throws IOException, InterruptedException, ExecutionException, FileParsingException {
        ArrayList<String> peps = new ArrayList<String>();
        while (next(true)) {
            PeptideFormatter pf = getPep();
            peps.add(pf.getProsit(PTMhandler.prositTmtAAMods) + "," + NceConstants.getNCE() + "," + pf.getCharge() + "," + Constants.FragmentationType);
        }
        return peps.toArray(new String[0]);
    }

    public String[] createAlphapeptdeepList()
            throws IOException {
        if (Constants.keepDecoys == 0) {
            Print.printInfo("Filtering out decoys");
        }
        ArrayList<String> peps = new ArrayList<String>();
        while (next(true)) {
            PeptideFormatter pf = getPep();
            if (Constants.keepDecoys == 0 && getTD() == 0) {
                continue;
            }
            peps.add(pf.getStripped() + "," + pf.getAlphapeptdeepMods() + "," + pf.getModPositions() + "," +
                    pf.getCharge() + "," + NceConstants.getNCE() + "," + Constants.instrument + "," +
                    pf.getLibrarytsv() + "," + getColumn("Proteins") + "," + -1 * (getTD() - 1));
        }
        return peps.toArray(new String[0]);
    }

    public String[] createFull() throws IOException {
        ArrayList<String> peps = new ArrayList<String>();
        while (next(true)) {
            PeptideFormatter pf = getPep();
            if (pf.cterm) {
                peps.add(pf.getBase() + "cterm" + "\t" + pf.getCharge());
            } else {
                peps.add(pf.getBase() + "\t" + pf.getCharge());
            }

            if (Constants.features.contains("peptideCounts")) {
                if (Float.valueOf(getColumn("hyperscore")) > 10) {
                    if (Constants.peptideCounter.containsKey(pf.getStripped())) {
                        HashSet<String> peptideSet = Constants.peptideCounter.get(pf.getStripped());
                        peptideSet.add(name);
                        Constants.peptideCounter.put(pf.getStripped(), peptideSet);
                    } else {
                        HashSet<String> peptideSet = new HashSet<>();
                        peptideSet.add(name);
                        Constants.peptideCounter.put(pf.getStripped(), peptideSet);
                    }
                } else {
                    if (! Constants.peptideCounter.containsKey(pf.getStripped())) {
                        HashSet<String> peptideSet = new HashSet<>();
                        Constants.peptideCounter.put(pf.getStripped(), peptideSet);
                    }
                }
            }
        }
        return peps.toArray(new String[0]);
    }

    public String[] createJSON(String modelFormat) throws IOException {
        ArrayList<String> peps = new ArrayList<String>();

        while (next(true)) {
            PeptideFormatter pf = getPep();
            if (PeptideSkipper.skipPeptide(pf, modelFormat)) {
                continue;
            }

            String instrument = mapInstrumentToModelSpecific(modelFormat, pf.getCharge());
            String pep = pf.getModel(modelFormat);
            String sb = pep + "," + pf.getCharge() + "," + NceConstants.getCalibratedNCE(modelFormat) + "," +
                    instrument + "," + Constants.FragmentationType + "," + pf.getStripped().length();
            peps.add(sb);
        }
        return peps.toArray(new String[0]);
    }

    public LinkedList[] getTopPSMs(int num, boolean charge2forIM) throws IOException {
        LinkedList<PeptideFormatter> PSMs = new LinkedList<>();
        LinkedList<Integer> scanNums = new LinkedList<>();

        //quicker way of getting top NCE
        PriorityQueue<Float> topEscores = new PriorityQueue<>(num, (a, b) -> Float.compare(b, a));
        while (next(true)) {
            if (charge2forIM) {
                if (getColumn("charge_2").equals("0")) {
                    continue;
                }
            }
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
            if (charge2forIM) {
                if (getColumn("charge_2").equals("0")) {
                    continue;
                }
            }
            float escore = Float.parseFloat(getEScore());
            if (escore <= eScoreCutoff) {
                PeptideFormatter pf = getPep(); //pf.getBaseCharge() + "," + pf.getDiann() + "," + pf.getStripped()
                PSMs.add(pf);
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