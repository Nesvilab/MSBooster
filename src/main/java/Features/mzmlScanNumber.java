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
    private float[] expMZs; //made private so don't accidentally access it, given that it may change
    private float[] expIntensities;
    float RT;
    int RTbinSize;
    float normalizedRT;
    Float IM;
    int IMbinSize;
    ArrayList<peptideObj> peptideObjects = new ArrayList<>();
    //double[] mzFreqs;
    public static float[] zeroFloatArray = new float[]{0};

    public mzmlScanNumber(IScan scan) throws FileParsingException {
        this.scanNum = scan.getNum();
        ISpectrum spectrum = scan.fetchSpectrum();
        this.expMZs = doubleToFloat(spectrum.getMZs());
        if (expMZs.length == 0) {
            System.out.println("scan with scan number " + scanNum + " is empty");
        }
        this.expIntensities = doubleToFloat(spectrum.getIntensities());
        this.RT = scan.getRt().floatValue();
        if (Constants.useIM) {
            this.IM = scan.getIm().floatValue();
        }
    }

    public mzmlScanNumber(int scanNum, float[] expMZs, float[] expInts, float RT, float IM) {
        this.scanNum = scanNum;
        this.expMZs = expMZs;
        this.expIntensities = expInts;
        this.RT = RT;
        this.IM = IM;
    }
    public mzmlScanNumber(int scanNum, float[] expMZs, float[] expInts, float RT) {
        this.scanNum = scanNum;
        this.expMZs = expMZs;
        this.expIntensities = expInts;
        this.RT = RT;
    }

    public float[] getExpMZs() { return expMZs; }
    public float[] getExpIntensities() { return expIntensities; }

    public void setPeptideObject(PeptideFormatter name, int rank, int targetORdecoy, String escore,
                                 HashMap<String, PredictionEntry> allPreds) throws Exception {

        if (rank != peptideObjects.size() + 1) { //need to add entries in order
            throw new AssertionError("must add next rank");
        }
        try {
            PredictionEntry predictionEntry = allPreds.get(name.baseCharge);
            float[] predMZs = predictionEntry.mzs;
            float[] predIntensities = predictionEntry.intensities;
            float predRT = predictionEntry.RT;
            float predIM = predictionEntry.IM;
            String[] fragmentIonTypes = predictionEntry.fragmentIonTypes;

            peptideObj newPepObj;
            if (predMZs.length > 1) {
                newPepObj = new peptideObj(this, name.baseCharge, rank, targetORdecoy, escore, predMZs,
                        predIntensities, predRT, predIM, fragmentIonTypes);
            } else { //only 1 frag to match
                newPepObj = new peptideObj(this, name.baseCharge, rank, targetORdecoy, escore, zeroFloatArray,
                        zeroFloatArray, predRT, predIM, fragmentIonTypes);
            }
            peptideObjects.add(rank - 1, newPepObj);

            //remove higher ranked peaks
            //TODO: base off fragger.params topN?
            if (Constants.removeRankPeaks) {
                for (int i = newPepObj.spectralSimObj.matchedIdx.size() - 1; i > -1; i--) {
                    expMZs = ArrayUtils.remove(expMZs, i);
                    expIntensities = ArrayUtils.remove(expIntensities, i);
                }
            }
        } catch (Exception e) { //TODO: percolator imputation
            //when peptide isn't in predictions, like unsupported amino acids
            //Set to arbitrary 0 vectors so nothing matches, similarity 0
            //may need to adapt this if using percolator imputation
            if (name.stripped.contains("U") || name.stripped.contains("O") || name.stripped.contains("X") ||
                    name.stripped.contains("B") || name.stripped.contains("Z") ) {
                peptideObjects.add(rank - 1, new peptideObj(this, name.baseCharge, rank, targetORdecoy, escore,
                        zeroFloatArray, zeroFloatArray, 0.0f, null, null));
            } else if (name.stripped.length() > 20) { //TODO: update this when longer ones are supported
                peptideObjects.add(rank - 1, new peptideObj(this, name.baseCharge, rank, targetORdecoy, escore,
                        zeroFloatArray, zeroFloatArray, 0.0f, null, null));
            } else if (Integer.parseInt(name.charge) > 6) { //TODO: update this for different tools
                peptideObjects.add(rank - 1, new peptideObj(this, name.baseCharge, rank, targetORdecoy, escore,
                        zeroFloatArray, zeroFloatArray, 0.0f, null, null));
            } else {
                String[] periodSplit = Constants.spectraRTPredFile.split("\\.");
                if (periodSplit[periodSplit.length - 1].equals("dlib")) { //won't always include every entry
                    peptideObjects.add(rank - 1, new peptideObj(this, name.baseCharge, rank, targetORdecoy, escore,
                            zeroFloatArray, zeroFloatArray, 0.0f, null, null));
                } else {
                    System.out.println("Prediction missing in file for " + name.baseCharge);
                    e.printStackTrace();
                    System.exit(-1);
                }
            }
        }
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

    public static void main(String[] args) throws IOException {
        System.out.println(Arrays.toString("LAHYNKR||2".split("\\|", 2)));
    }
}
