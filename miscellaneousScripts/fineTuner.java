package Features;

import java.util.ArrayList;
import java.util.HashSet;

public class fineTuner {

    pepXMLReader rank1;
    pepXMLReader rank2;
    fineTuneHit[] hits = new fineTuneHit[Constants.fineTuneSize];

    //method for getting fine tuning set for pDeep3
    //consider own class for this, so don't need to read twice
    public void addPepXMLReader(pepXMLReader x) {
        if (rank1 == null) {
            rank1 = x;
        } else if (rank2 == null) {
            rank2 = x;
        }
    }

    public class fineTuneHit {
        String pep;
        double e;
        double rt;
        int scanNum;

        fineTuneHit(String pep, String e, double rt, int scanNum) {
            this.pep = pep;
            this.e = Double.parseDouble(e);
            this.rt = rt;
            this.scanNum = scanNum; //in case it needs scannum below
        }
    }

    public void fineTune(String outfile) {

        //find unique scan numbers between rank 1 and 2
        int[] rank1nums = rank1.getScanNumbers();
        int[] rank2nums = rank2.getScanNumbers();

        String[] pep = rank1.getXMLpeptides();
        int[] td = rank1.getTargetOrDecoy();
        String[] e = rank1.getEScore();
        double[] rt = rank1.getRT();

        //get scans that are target peptides and only contain that one peptide
        ArrayList<fineTuneHit> fineTuneHits = new ArrayList<>();
        int i = 0; // for rank2nums
        for (int j = 0; j < rank1nums.length; j++) {
            if (rank1nums[j] != rank2nums[i]) {
                if (td[j] == 1) { //target
                    fineTuneHits.add(new fineTuneHit(pep[j], e[j], rt[j], rank1nums[j]));
                }
            } else { //scan with multiple peptides identified
                if (i < rank2nums.length - 1) { //no issue when run out
                    i++;
                }
            }
        }

        //get top scorers
        fineTuneHits.sort((a, b) -> Double.compare(a.e, b.e));
        HashSet<String> pepNames = new HashSet<>();
        String[] finalHits = new String[Constants.fineTuneSize];
        int finalHitsIdx = 0;

        int fineTuneHitsIdx = -1;
        while (finalHitsIdx < Constants.fineTuneSize) {
            try {
                fineTuneHitsIdx += 1;
                fineTuneHit fth = fineTuneHits.get(fineTuneHitsIdx);
                //System.out.println(fineTuneHitsIdx);
                String xmlpep = fth.pep;

                if (!pepNames.contains(xmlpep)) { //skips XMLpeptide if detected before,
                    // since new instance will have  worse expectation
                    String[] pepSplit = xmlpep.split("\\|");

                    String hitToAdd = rank1.baseName + "\t" + fth.scanNum + "\t" + pepSplit[0] + "\t" + pepSplit[1] +
                            "\t" + pepSplit[2] + "\t" + fth.rt;
                    finalHits[finalHitsIdx] = hitToAdd;
                    hits[finalHitsIdx] = fth; //for creating mgf
                    finalHitsIdx += 1;
                    pepNames.add(xmlpep);
                }
            } catch (Exception exc) { //no more fineTune hits
                System.out.println("using only " + finalHitsIdx + " peptides for finetuning");
            }
        }

        //create fineTune file
//        peptideFileCreator.createPeptideFile(finalHits, constants.outputFolder + outfile,
//                "fineTune");
    }


    public static void main(String[] args) {
        pepXMLReader x = new pepXMLReader("C:/Users/kevin/Downloads/proteomics/wideWindow/wide1/" +
                "23aug2017_hela_serum_timecourse_pool_wide_001_170829031834_rank1.pepXML");
        pepXMLReader y = new pepXMLReader("C:/Users/kevin/Downloads/proteomics/wideWindow/wide1/" +
                "23aug2017_hela_serum_timecourse_pool_wide_001_170829031834_rank2.pepXML");

        fineTuner ft = new fineTuner();
        ft.addPepXMLReader(x);
        ft.addPepXMLReader(y);
        ft.fineTune("fineTune_wide1.txt");
    }
}
