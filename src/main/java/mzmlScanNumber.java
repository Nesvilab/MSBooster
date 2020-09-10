import umich.ms.fileio.exceptions.FileParsingException;

import java.io.IOException;
import java.lang.reflect.InvocationTargetException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.TreeMap;

public class mzmlScanNumber {
    final int scanNum;
    final mzMLReader mzmlScans;
    double[] expMZs;
    double[] expIntensities;
    TreeMap<Integer, peptideObj> peptideObjects = new TreeMap<>();

    public mzmlScanNumber(mzMLReader mzmlScans, int scanNum) throws FileParsingException {
        this.scanNum = scanNum;
        this.mzmlScans = mzmlScans;
        this.expMZs = mzmlScans.getMZ(scanNum);
        this.expIntensities = mzmlScans.getIntensity(scanNum);
    }

    public void setPeptideObject(String name, int rank, int targetORdecoy,
                                 HashMap<String, double[]> allPredMZs, HashMap<String, double[]> allPredIntensities) {
        double[] predMZs = allPredMZs.get(name);
        double[] predIntensities = allPredIntensities.get(name);

        peptideObjects.put(rank, new peptideObj(this, name, rank, targetORdecoy, predMZs, predIntensities));
    }

    public peptideObj getPeptideObject(int rank) {
        return peptideObjects.get(rank);
    }

    public peptideObj getPeptideObject(String name) {
        peptideObj returnedP = null;

        for (peptideObj p : peptideObjects.values()) {
            if (p.name == name) {
                return p;
            }
        }
        System.out.println("no peptideObj with name " + name + ", returning null");
        return returnedP; //only if peptide name not there
    }

    public int[] targetDecoyOrder() { // should work because treemap is sorted
        //good for control when comparing how similarity measures rerank
        //need to account for when peptide isn't added, like those with U
        int[] tdList = new int[peptideObjects.size()];

        int i = 0;
        for (peptideObj pObj : peptideObjects.values()) {
            int rank = pObj.rank;
            int targetORdecoy = pObj.targetORdecoy;
            tdList[i] = targetORdecoy;
            i++;
        }

        return tdList;
    }

    //this version provides new ranks based on similarity rescoring
    public ArrayList<Integer> targetDecoyOrder(String similarityMeasure) throws InvocationTargetException, NoSuchMethodException,
            FileParsingException, IllegalAccessException {
        ArrayList<Integer> tdList = new ArrayList<>();
        ArrayList<Double> scoreList = new ArrayList<>();

        for (peptideObj pObj : peptideObjects.values()) {
            int targetORdecoy = pObj.targetORdecoy;
            double score = pObj.getScore(similarityMeasure);

            if (tdList.isEmpty()) {
                tdList.add(targetORdecoy);
                scoreList.add(score);
            } else {
                boolean added = false;

                for (int i = 0; i < tdList.size(); i++) {
                    double prevScore = scoreList.get(i);

                    if (score > prevScore) {
                        tdList.add(i, targetORdecoy);
                        scoreList.add(i, score);
                        added = true;
                        break;
                    }
                }

                if (! added) { //lowest similarity yet
                    tdList.add(targetORdecoy);
                    scoreList.add(score);
                }
            }
        }
        return tdList;
    }

    public static void main(String[] args) throws IOException {
        mgfFileReader mgf = new mgfFileReader("preds/rank1_preds.txt");
        HashMap<String, double[]> x = mgf.getMzDict();
        System.out.println(Arrays.toString(x.get("IVNPDPVVK||2")));
    }
}
