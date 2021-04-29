package Features;

import org.apache.commons.lang3.ArrayUtils;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.util.HashSet;

public class pinReader {
    String name; //used for resetting
    BufferedReader in;
    String[] header;
    //List<String[]> rows = new ArrayList<>();
    private String[] row;
    int rowNum = -1;

    int scanNumIdx;
    int labelIdx;
    int rankIdx;
    int specIdx;
    int pepIdx;
    int eScoreIdx;

    public pinReader(String pin) throws IOException {
        name = pin;
        in = new BufferedReader(new FileReader(name));
        String line = in.readLine();
        header = line.split("\t");
//        while ((line = in.readLine()) != null) {
//            rows.add(line.split("\t"));
//        }
//        in.close();

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
//        rowNum += 1;
//        try {
//            row = rows.get(rowNum);
//            return true;
//        } catch (Exception e) {
//            return false;
//        }
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

    public int getTD() {return Math.max(0, Integer.parseInt(row[labelIdx]));}

    public int getScanNum() {return Integer.parseInt(row[scanNumIdx]);}

    public int getRank() {return Integer.parseInt(row[rankIdx]);}

    //public String getEScore() {return String.valueOf(Math.pow(10, Double.parseDouble(row[eScoreIdx])));}
    public String getEScore() {return String.valueOf(Math.exp(15.0 - Double.parseDouble(row[eScoreIdx])));}

    public HashSet<String> getAllPep() throws IOException {
        HashSet<String> peps = new HashSet<String>();
        while (next()) {
            peps.add(getPep().split("\\|")[0]);
        }
        reset();
        return peps;
    }

    public static void main(String[] args) throws IOException {
        pinReader p = new pinReader("C:/Users/kevin/Downloads/proteomics/wideWindow/1-18/perc/combined.pin");
        while (p.next()) {
            if (p.getRow().length > 33) {
                System.out.println(p.getRow().length);
            }
        }
    }
}