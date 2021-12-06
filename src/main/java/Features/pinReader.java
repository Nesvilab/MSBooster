package Features;

import org.apache.commons.lang3.ArrayUtils;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.TreeMap;

public class pinReader {
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

    public pinReader(String pin) throws IOException {
        //TODO: use log10eval if available
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
        eScoreIdx = ArrayUtils.indexOf(header, "hyperscore");
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
            return true;
        }
        //in.close();
        return false;
    }

    public void close() throws IOException {
        in.close();
    }

    public String[] getRow() {return row;}

    public String getPep() {return percolatorFormatter.percolatorPepFormat(row, pepIdx, specIdx);}

    public String getStrippedPep() {
        String peptide = percolatorFormatter.percolatorPepFormat(row, pepIdx, specIdx).split("\\|")[0];
        ArrayList<Integer> starts = new ArrayList<>();
        ArrayList<Integer> ends = new ArrayList<>();
        ends.add(0);
        for (int i = 0; i < peptide.length(); i++) {
            String myChar = peptide.substring(i, i + 1);
            if (myChar.equals("[")) {
                starts.add(i);
            } else if (myChar.equals("]")) {
                ends.add(i + 1);
            }
        }
        starts.add(peptide.length());

        StringBuilder sb = new StringBuilder();
        for (int i = 0; i < starts.size(); i++) {
            sb.append(peptide.substring(ends.get(i), starts.get(i)));
        }
        return sb.toString();
    }

    public String getFullPep() {return percolatorFormatter.percolatorPepFormatFull(row, pepIdx, specIdx);}

    public int getTD() {return Math.max(0, Integer.parseInt(row[labelIdx]));}

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
    public String getEScore() {return String.valueOf(Math.exp(15.0 - Double.parseDouble(row[eScoreIdx])));}

//    public HashSet<String> getAllPep() throws IOException {
//        HashSet<String> peps = new HashSet<String>();
//        while (next()) {
//            peps.add(getPep().split("\\|")[0]);
//        }
//        reset();
//        return peps;
//    }

    public String[] createPDeep3List() throws IOException {
        ArrayList<String> peps = new ArrayList<String>();
        while (next()) {
            String[] pepSplit = getPep().split("\\|");
            peps.add("." + "\t" + "." + "\t" + pepSplit[0] + "\t" + pepSplit[1] + "\t" + pepSplit[2]);
        }
        return peps.toArray(new String[0]);
    }

    public String[] createDeepMSPeptideList() throws IOException {
        ArrayList<String> peps = new ArrayList<String>();
        while (next()) {
            peps.add(getStrippedPep());
        }
        return peps.toArray(new String[0]);
    }
    //TODO: how to handle shifting of uncommon PTMS?
    public String[] createDiannList() throws IOException {
        ArrayList<String> peps = new ArrayList<String>();
        TreeMap<Integer, Integer> modMap = new TreeMap<>(); //sorted for future use
        while (next()) {
            String[] pepSplit = getPep().split("\\|");
            peps.add(pepSplit[0] + "\t" + pepSplit[1]);
        }
        return peps.toArray(new String[0]);
    }

    public String[] createDiannListFull() throws IOException {
        ArrayList<String> peps = new ArrayList<String>();
        TreeMap<Integer, Integer> modMap = new TreeMap<>(); //sorted for future use
        while (next()) {
            String[] pepSplit = getFullPep().split("\\|");
            peps.add(pepSplit[0] + "\t" + pepSplit[1]);
        }
        return peps.toArray(new String[0]);
    }

    public static void main(String[] args) throws IOException {
//        pinReader p = new pinReader("C:/Users/kevin/Downloads/proteomics/wideWindow/1-18/perc/combined.pin");
//        while (p.next()) {
//            if (p.getRow().length > 33) {
//                System.out.println(p.getRow().length);
//            }
//        }
    }
}