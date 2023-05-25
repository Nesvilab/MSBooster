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
import umich.ms.datatypes.scan.IScan;
import umich.ms.datatypes.spectrum.ISpectrum;
import umich.ms.fileio.exceptions.FileParsingException;

import java.io.IOException;
import java.util.*;

import static Features.floatUtils.doubleToFloat;

public class mzmlScanNumber {
    final int scanNum;
    private float[] expMZs;
    private float[] expIntensities;
    float RT;
    int RTbinSize;
    float normalizedRT;
    Float IM;
    int IMbinSize;
    ArrayList<Float> NCEs = new ArrayList<>(); //TODO use these to write peptide input files for prediction
    ArrayList<peptideObj> peptideObjects = new ArrayList<>();
    //double[] mzFreqs;
    public static float[] zeroFloatArray = new float[]{0};
    private static final HashSet<String> nceModels =
            new HashSet<>(Arrays.asList("PredFull", "Prosit", "PrositTMT", "alphapeptdeep"));

    //this version if creating from mzml scan number
    public mzmlScanNumber(IScan scan) throws FileParsingException {
        this.scanNum = scan.getNum();
        ISpectrum spectrum = scan.fetchSpectrum();
        this.expMZs = doubleToFloat(spectrum.getMZs());
        this.expIntensities = doubleToFloat(spectrum.getIntensities());
        this.RT = scan.getRt().floatValue();
        if (Constants.useIM) {
            this.IM = scan.getIm().floatValue();
        }
        if (nceModels.contains(Constants.spectraRTPredModel)) {
            String[] nceInfo = scan.getFilterString().split("@");
            if (nceInfo.length > 1) {
                for (String s : Arrays.copyOfRange(nceInfo, 1, nceInfo.length)) {
                    String fragmentationInfo = s.split(" ")[0];
                    for (int i = 0; i < fragmentationInfo.length(); i++) {
                        if (Character.isDigit(fragmentationInfo.charAt(i))) {
                            NCEs.add(Float.parseFloat(fragmentationInfo.substring(i)));
                            break;
                        }
                    }
                }
            }
        }
    }

    //this version if creating from mgf file
    public mzmlScanNumber(int scanNum, float[] expMZs, float[] expInts, float RT, float IM) {
        this.scanNum = scanNum;
        this.expMZs = expMZs;
        this.expIntensities = expInts;
        this.RT = RT;
        this.IM = IM;
    }

    public float[] getExpMZs() { return expMZs; }
    public float[] getExpIntensities() { return expIntensities; }

    public void setPeptideObject(PeptideFormatter name, int rank, int targetORdecoy, String escore,
                                 HashMap<String, PredictionEntry> allPreds) throws Exception {

        if (rank != peptideObjects.size() + 1) { //need to add entries in order
            throw new AssertionError("must add next rank");
        }
        PredictionEntry predictionEntry = allPreds.get(name.baseCharge);
        try {
            float[] predMZs = predictionEntry.mzs;
            float[] predIntensities = predictionEntry.intensities;
            float predRT = predictionEntry.RT;
            float predIM = predictionEntry.IM;

            peptideObj newPepObj;
            if (predMZs.length > 1) {
                if (Constants.divideFragments.equals("0")) {
                    newPepObj = new peptideObj(this, name.baseCharge, rank, targetORdecoy, escore, predMZs,
                            predIntensities, predRT, predIM);
                } else {
                    newPepObj = new peptideObj(this, name.baseCharge, rank, targetORdecoy, escore, predMZs,
                            predIntensities, predRT, predIM, predictionEntry.fragmentIonTypes);
                }
            } else { //only 1 frag to match
                newPepObj = new peptideObj(this, name.baseCharge, rank, targetORdecoy, escore, zeroFloatArray,
                        zeroFloatArray, predRT, predIM);
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
            float predRT = 0f;
            float predIM = 0f;
            if (predictionEntry != null) {
                predRT = predictionEntry.RT; //best option?
                predIM = predictionEntry.IM; //best option?
            }

            //when peptide isn't in predictions, like unsupported amino acids
            //Set to arbitrary 0 vectors so nothing matches, similarity 0
            //may need to adapt this if using percolator imputation
            if (name.stripped.contains("U") || name.stripped.contains("O") || name.stripped.contains("X") ||
                    name.stripped.contains("B") || name.stripped.contains("Z") ) {
                peptideObjects.add(rank - 1, new peptideObj(this, name.baseCharge, rank, targetORdecoy, escore,
                        zeroFloatArray, zeroFloatArray, predRT, predIM));
            } else if (name.stripped.length() > 15) { //TODO: update this when longer ones are supported
                peptideObjects.add(rank - 1, new peptideObj(this, name.baseCharge, rank, targetORdecoy, escore,
                        zeroFloatArray, zeroFloatArray, predRT, predIM));
            } else if (name.stripped.length() < 7) { //TODO: update this when shorter ones are supported
                peptideObjects.add(rank - 1, new peptideObj(this, name.baseCharge, rank, targetORdecoy, escore,
                        zeroFloatArray, zeroFloatArray, predRT, predIM));
            } else if (Integer.parseInt(name.charge) > 6) { //TODO: update this for different tools
                peptideObjects.add(rank - 1, new peptideObj(this, name.baseCharge, rank, targetORdecoy, escore,
                        zeroFloatArray, zeroFloatArray, predRT, predIM));
            } else {
                String[] periodSplit = Constants.spectraRTPredFile.split("\\.");
                if (periodSplit[periodSplit.length - 1].equals("dlib")) { //won't always include every entry
                    peptideObjects.add(rank - 1, new peptideObj(this, name.baseCharge, rank, targetORdecoy, escore,
                            zeroFloatArray, zeroFloatArray, predRT, predIM));
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
}
