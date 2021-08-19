package Features;

import org.apache.commons.lang3.ArrayUtils;
import umich.ms.datatypes.scan.IScan;
import umich.ms.datatypes.spectrum.ISpectrum;
import umich.ms.fileio.exceptions.FileParsingException;

import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;

import static Features.floatUtils.doubleToFloat;

public class mzmlScanNumber {
    final int scanNum;
    //mzMLReader mzmlScans;
    //final String mzmlPath;
    private float[] expMZs; //made private so don't accidentally access it, given that it may change
    private float[] expIntensities;
    public float[] OGexpMZs;
    public float[] OGexpIntensities;
    float RT;
    int RTbinSize;
    float normalizedRT;
    Float IM;
    int IMbinSize;
    ArrayList<peptideObj> peptideObjects = new ArrayList<>();
    //double[] mzFreqs;
    public static float[] zeroFloatArray = new float[]{0};

//    public mzmlScanNumber(mzMLReader mzmlScans, int scanNum, float RT, Float IM) throws FileParsingException {
//        this.scanNum = scanNum;
//        //this.mzmlScans = mzmlScans;
//        this.mzmlPath = mzmlScans.pathStr;
//        this.expMZs = mzmlScans.getMZ(scanNum);
//        this.expIntensities = mzmlScans.getIntensity(scanNum);
//        this.RT = RT;
//        this.IM = IM;
//        //this.mzFreqs = mzmlScans.getMzFreq();
//        //this.windowStart = windowStart;
//    }

    public mzmlScanNumber(IScan scan) throws FileParsingException {
        this.scanNum = scan.getNum();
        ISpectrum spectrum = scan.fetchSpectrum();
        this.expMZs = doubleToFloat(spectrum.getMZs());
        this.expIntensities = doubleToFloat(spectrum.getIntensities());
        this.OGexpMZs = new float[expMZs.length];
        this.OGexpIntensities = new float[expIntensities.length];
        for (int i = 0; i < expMZs.length; i++) { //for maxConsecutive score
            OGexpMZs[i] = expMZs[i];
            OGexpIntensities[i] = expIntensities[i];
        }
        this.RT = scan.getRt().floatValue();
        if (Constants.useIM) {
            this.IM = scan.getIm().floatValue();
        }
    }

    public mzmlScanNumber(int scanNum, float[] expMZs, float[] expInts, float RT, float IM) throws FileParsingException {
        this.scanNum = scanNum;
        this.expMZs = expMZs;
        this.expIntensities = expInts;
        this.OGexpMZs = new float[expMZs.length];
        this.OGexpIntensities = new float[expIntensities.length];
        for (int i = 0; i < expMZs.length; i++) { //for maxConsecutive score
            OGexpMZs[i] = expMZs[i];
            OGexpIntensities[i] = expIntensities[i];
        }
        this.RT = RT;
        this.IM = IM;
    }
    public mzmlScanNumber(int scanNum, float[] expMZs, float[] expInts, float RT) throws FileParsingException {
        this.scanNum = scanNum;
        this.expMZs = expMZs;
        this.expIntensities = expInts;
        this.OGexpMZs = new float[expMZs.length];
        this.OGexpIntensities = new float[expIntensities.length];
        for (int i = 0; i < expMZs.length; i++) { //for maxConsecutive score
            OGexpMZs[i] = expMZs[i];
            OGexpIntensities[i] = expIntensities[i];
        }
        this.RT = RT;
    }

//    public void setMZMLScans(mzMLReader mr){ //we don't always want the scans to stick around
//        this.mzmlScans = mr;
//    }

    public float[] getExpMZs() { return expMZs; }
    public float[] getExpIntensities() { return expIntensities; }

    public void setPeptideObject(String name, int rank, int targetORdecoy, String escore,
                                 HashMap<String, float[]> allPredMZs, HashMap<String, float[]> allPredIntensities,
                                 HashMap<String, Float> allPredRTs, HashMap<String, Float> allPredIMs) throws Exception {

        if (rank != peptideObjects.size() + 1) { //need to add entries in order
            throw new AssertionError("must add next rank");
        }
        try {
            float[] predMZs = allPredMZs.get(name);
            float[] predIntensities = allPredIntensities.get(name);
            float predRT = allPredRTs.get(name);
            float predIM = allPredIMs.get(name);

            peptideObj newPepObj;
            if (predMZs.length > 1) {
                newPepObj = new peptideObj(this, name, rank, targetORdecoy, escore, predMZs,
                        predIntensities, predRT, predIM);
            } else { //only 1 frag to match
                newPepObj = new peptideObj(this, name, rank, targetORdecoy, escore, zeroFloatArray,
                        zeroFloatArray, predRT, predIM);
            }
            peptideObjects.add(rank - 1, newPepObj);

            //remove higher ranked peaks
            if (Constants.removeRankPeaks) {
                for (int i = newPepObj.spectralSimObj.matchedIdx.size() - 1; i > -1; i--) {
                    expMZs = ArrayUtils.remove(expMZs, i);
                    expIntensities = ArrayUtils.remove(expIntensities, i);
                }
            }
        } catch (Exception e) {
            //when peptide isn't in predictions, like U peptides.
            //Set to arbitrary 0 vectors so nothing matches, similarity 0
            peptideObjects.add(rank - 1, new peptideObj(this, name, rank, targetORdecoy, escore,
                    zeroFloatArray, zeroFloatArray, 0.0f, null));
        }

//        //recording experimental IMs for peptides
//        if (! mzmlScans.peptideIMs.containsKey(name)) {
//            mzmlScans.peptideIMs.put(name, new ArrayList<Float>());
//        }
//        mzmlScans.peptideIMs.get(name).add(IM);
//        //for weighted
//        for (int i = 0; i < Math.max(Math.ceil(-1 * Math.log10(Double.parseDouble(peptideObjects.get(rank - 1).escore))), 1); i++) {
//            mzmlScans.peptideIMs.get(name).add(IM);
//        }
    }

    public peptideObj getPeptideObject(int rank) {
        return peptideObjects.get(rank - 1);
    } //arraylist now, not hashmap

    public peptideObj getPeptideObject(String name) {
        peptideObj returnedP = null;

        for (peptideObj p : peptideObjects) {
            if (p.name.equals(name)) {
                return p;
            }
        }
        System.out.println("no peptideObj with name " + name + ", returning null");
        return returnedP; //only if peptide name not there
    }

//    public int[] targetDecoyOrder() { // should work because treemap is sorted
//        //good for control when comparing how similarity measures rerank
//        //need to account for when peptide isn't added, like those with U
//        int[] tdList = new int[peptideObjects.size()];
//        for (peptideObj pObj : peptideObjects) {tdList[pObj.rank - 1] = pObj.targetORdecoy;}
//        return tdList;
//    }
//
//    //this version provides new ranks based on similarity rescoring
//    public ArrayList<Integer> targetDecoyOrder(String similarityMeasure) throws InvocationTargetException, NoSuchMethodException,
//            FileParsingException, IllegalAccessException {
//        ArrayList<Integer> tdList = new ArrayList<>();
//        ArrayList<Double> scoreList = new ArrayList<>();
//
//        for (peptideObj pObj : peptideObjects) {
//            int targetORdecoy = pObj.targetORdecoy;
//            double score = pObj.getScore(similarityMeasure); //may fail for U peptides
//
//            if (tdList.isEmpty()) {
//                tdList.add(targetORdecoy);
//                scoreList.add(score);
//            } else {
//                boolean added = false;
//
//                //add score in place in sorted list, biggest to smallest
//                for (int i = 0; i < tdList.size(); i++) {
//                    double prevScore = scoreList.get(i);
//
//                    if (score >= prevScore) {
//                        tdList.add(i, targetORdecoy);
//                        scoreList.add(i, score);
//                        added = true;
//                        break;
//                    }
//                }
//
//                if (! added) { //lowest similarity yet
//                    tdList.add(targetORdecoy);
//                    scoreList.add(score);
//                }
//            }
//        }
//        return tdList;
//    }
//
//    //use in method above
//    public Object[] reRankOrder(String similarityMeasure) throws InvocationTargetException,
//            NoSuchMethodException, FileParsingException, IllegalAccessException {
//        //return 1) reranked list; 2) scores in decreasing order
//
//        ArrayList<Integer> reRankList = new ArrayList<>();
//        ArrayList<Double> scoreList = new ArrayList<>();
//
//        for (peptideObj pObj : peptideObjects) {
//            int rank = pObj.rank;
//            double score = pObj.getScore(similarityMeasure); //may fail for U peptides, need to put at end of list
//
//            if (reRankList.isEmpty()) {
//                reRankList.add(rank);
//                scoreList.add(score);
//            } else {
//                boolean added = false;
//
//                //add score in place in sorted list, biggest to smallest
//                for (int i = 0; i < reRankList.size(); i++) {
//                    double prevScore = scoreList.get(i);
//
//                    if (score >= prevScore) {
//                        reRankList.add(i, rank);
//                        scoreList.add(i, score);
//                        added = true;
//                        break;
//                    }
//                }
//
//                if (! added) { //lowest similarity yet
//                    reRankList.add(rank);
//                    scoreList.add(score);
//                }
//            }
//        }
//
//        return new Object[] {reRankList, scoreList};
//    }
//
//    public String[] getNames(boolean targetOnly){
//        if (targetOnly) {
//            ArrayList<String> namesList = new ArrayList<>();
//            for (peptideObj p : peptideObjects) {
//                if (p.targetORdecoy == 1){
//                    namesList.add(p.name);
//                }
//            }
//            String[] names = new String[namesList.size()];
//            names = namesList.toArray(names);
//            return names;
//        } else {
//            int i = 0;
//            String[] names = new String[peptideObjects.size()];
//            for (peptideObj p : peptideObjects) {
//                names[i] = p.name;
//                i++;
//            }
//            return names;
//        }
//    }

    public static void main(String[] args) throws IOException {
        System.out.println(Arrays.toString("LAHYNKR||2".split("\\|", 2)));
    }
}
