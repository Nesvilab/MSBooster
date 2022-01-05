package Features;

import org.apache.commons.lang3.ArrayUtils;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;

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
    private boolean calcEvalue = false;

    public pinReader(String pin) throws IOException {
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
        if (Arrays.stream(header).anyMatch("log10_evalue"::equals)) {
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

    public String getPep() {
        if (Constants.spectraRTPredModel.equals("DIA-NN")) {
            return percolatorFormatter.DiannPepFormat(row, pepIdx, specIdx);
        } else { //PredFull
            return percolatorFormatter.PredfullPepFormat(row, pepIdx, specIdx);
        }
    }

    public String getStrippedPep() {
        String peptide = percolatorFormatter.DiannPepFormat(row, pepIdx, specIdx).split("\\|")[0];
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
    public String getEScore() {
        if (calcEvalue) {
            return String.valueOf(Math.exp(15.0 - Double.parseDouble(row[eScoreIdx])));
        } else {
            return row[eScoreIdx];
        }
    }

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

    public String[] createDiannList() throws IOException {
        ArrayList<String> peps = new ArrayList<String>();
        //TreeMap<Integer, Integer> modMap = new TreeMap<>(); //sorted for future use
        while (next()) {
            String[] pepSplit = getPep().split("\\|");
            peps.add(pepSplit[0] + "\t" + pepSplit[1]);
        }
        return peps.toArray(new String[0]);
    }

    public String[] createDiannListFull() throws IOException {
        ArrayList<String> peps = new ArrayList<String>();
        //TreeMap<Integer, Integer> modMap = new TreeMap<>(); //sorted for future use
        while (next()) {
            String[] pepSplit = getFullPep().split("\\|");
            peps.add(pepSplit[0] + "\t" + pepSplit[1]);
        }
        return peps.toArray(new String[0]);
    }

    public String[] createPredFullList() throws IOException {
        ArrayList<String> peps = new ArrayList<String>();
        while (next()) {
            //first is peptide, then missed masses
            String pep = row[pepIdx];
            pep = pep.substring(2, pep.length() - 2);

            //n term acetylation
            if (pep.charAt(0) == 'n') {
                pep = pep.replace("n", "");
            }
            pep = pep.replace("c","");

            //replace oxidized M
            pep = pep.replace("M[15.9949]", "M(O)");


            //find locations of PTMs
            ArrayList<Integer> starts = new ArrayList<>();
            ArrayList<Integer> ends = new ArrayList<>();
            for (int i = 0; i < pep.length(); i++) {
                if (pep.charAt(i) == '[') {
                    starts.add(i);
                } else if (pep.charAt(i) == ']') {
                    ends.add(i);
                }
            }

            for (int i = starts.size() - 1; i > -1; i--) {
                pep = pep.substring(0, starts.get(i)) + pep.substring(ends.get(i) + 1);
            }

            //charge
            String[] charges = row[specIdx].split("_");
            String chargeStr = charges[charges.length - 2];
            int charge = Integer.parseInt(chargeStr.substring(chargeStr.length() - 1));

            peps.add(pep + "\t" + charge + "\t" + Constants.FragmentationType + "\t" + Constants.NCE);
        }
        return peps.toArray(new String[0]);
    }

    public String[] createPredFullListFull() throws IOException {
        ArrayList<String> peps = new ArrayList<String>();
        while (next()) {
//            //first is peptide, then missed masses
//            String pep = row[pepIdx];
//            pep = pep.substring(2, pep.length() - 2);
//
//            //n term acetylation
//            if (pep.charAt(0) == 'n') {
//                pep = pep.replace("n", "");
//            }
//            pep = pep.replace("c","");
//
//            //charge
//            String[] charges = row[specIdx].split("_");
//            String chargeStr = charges[charges.length - 2];
//            int charge = Integer.parseInt(chargeStr.substring(chargeStr.length() - 1));
//
//            peps.add(pep + "\t" + charge);
            peps.add(getPep());
        }
        return peps.toArray(new String[0]);
    }

    public static void main(String[] args) throws IOException {
        pinReader p = new pinReader("C:/Users/kevin/Downloads/proteomics/melanoma/201905024_F_7951_pro_1.pin");
        p.next();
        System.out.println(p.getEScore());
    }
}