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

import static Features.FloatUtils.doubleToFloat;
import static utils.Print.printError;
import static utils.Print.printInfo;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.concurrent.ConcurrentHashMap;
import org.apache.commons.lang3.ArrayUtils;
import umich.ms.datatypes.scan.IScan;
import umich.ms.datatypes.scan.props.PrecursorInfo;
import umich.ms.datatypes.spectrum.ISpectrum;
import umich.ms.fileio.exceptions.FileParsingException;

public class MzmlScanNumber {
    final int scanNum;
    public double isolationLower;
    public double isolationUpper;
    public float[] expMZs;
    public float[] expIntensities;
    public float[] savedExpMZs;
    public float[] savedExpIntensities;
    float RT;
    int RTbinSize;
    float normalizedRT;
    Float IM;
    int IMbinSize;
    Double lowerLimit;
    Double upperLimit;
    public HashMap<String, Float> NCEs = new HashMap<>();
    ArrayList<PeptideObj> peptideObjects = new ArrayList<>();
    //double[] mzFreqs;
    public static float[] zeroFloatArray = new float[]{0};

    //this version if creating from mzml scan number
    public MzmlScanNumber(IScan scan) throws FileParsingException {
        this.scanNum = scan.getNum();
        PrecursorInfo pi = scan.getPrecursor();
        if (pi.getMzRangeStart() != null) {
            this.isolationLower = pi.getMzRangeStart() == 0 ? Double.NaN : pi.getMzRangeStart();
            this.isolationUpper = pi.getMzRangeEnd() == 0 ? Double.NaN : pi.getMzRangeEnd();
        }
        ISpectrum spectrum = scan.fetchSpectrum();
        this.expMZs = doubleToFloat(spectrum.getMZs());
        this.expIntensities = doubleToFloat(spectrum.getIntensities());
        sortArrays(this.expMZs, this.expIntensities);
        if (Constants.removeRankPeaks && (Constants.features.contains("hypergeometricProbability") ||
                Constants.features.contains("intersection") ||
                Constants.features.contains("adjacentSimilarity"))) {
            this.savedExpMZs = this.expMZs;
            this.savedExpIntensities = this.expIntensities;
        }
        this.RT = scan.getRt().floatValue();
        if (Constants.useIM) {
            this.IM = scan.getIm().floatValue();
        }
        lowerLimit = scan.getScanMzWindowLower();
        upperLimit = scan.getScanMzWindowUpper();

        if (!Constants.NCE.isEmpty() && !Constants.FragmentationType.isEmpty()) {
            NCEs.put(Constants.FragmentationType, Float.parseFloat(Constants.NCE));
        } else {
            try {
                String[] nceInfo = scan.getFilterString().split("@");
                if (nceInfo.length > 1) {
                    for (String s : Arrays.copyOfRange(nceInfo, 1, nceInfo.length)) {
                        String fragmentationInfo = s.split(" ")[0];
                        StringBuilder fragmentationType = new StringBuilder();
                        for (int i = 0; i < fragmentationInfo.length(); i++) {
                            char myChar = fragmentationInfo.charAt(i);
                            if (Character.isDigit(myChar)) {
                                NCEs.put(fragmentationType.toString().toUpperCase(),
                                        Float.parseFloat(fragmentationInfo.substring(i)));
                                break;
                            } else {
                                fragmentationType.append(myChar);
                            }
                        }
                    }
                }
            } catch (NullPointerException e) { //like in DIA-Umpire, there is no filter string
                printInfo("mzml file does not contain filter string. Setting NCE to 25 and " +
                        "FragmentationType to HCD. If other values are desired, please set them in " +
                        "the parameter file with NCE and FragmentationType.");
                Constants.NCE = "25";
                Constants.FragmentationType = "HCD";
            }
//            }
        }
    }

    //this version if creating from mgf file
    public MzmlScanNumber(int scanNum, float[] expMZs, float[] expInts, float RT, float IM) {
        this.scanNum = scanNum;
        this.expMZs = expMZs;
        this.expIntensities = expInts;
        sortArrays(this.expMZs, this.expIntensities);
        this.RT = RT;
        this.IM = IM;
    }

    public static void sortArrays(float[] eMZs, float[] eIntensities) {
        int n = eMZs.length;
        float[][] pairs = new float[n][2];

        // Create pairs of eMZs and corresponding eIntensities
        for (int i = 0; i < n; i++) {
            pairs[i][0] = eMZs[i];
            pairs[i][1] = eIntensities[i];
        }

        // Sort the pairs based on eMZs
        Arrays.sort(pairs, (a, b) -> Float.compare(a[0], b[0]));

        // Update eMZs and eIntensities based on sorted pairs
        for (int i = 0; i < n; i++) {
            eMZs[i] = pairs[i][0];
            eIntensities[i] = pairs[i][1];
        }
    }

    public float[] getExpMZs() { return expMZs; }
    public float[] getExpIntensities() { return expIntensities; }

    public PeptideObj setPeptideObject(PeptideFormatter name, int rank, int targetORdecoy, String escore,
                                       PredictionEntryHashMap allPreds, boolean set) {
        PredictionEntry predictionEntry = allPreds.get(name.baseCharge);
        PeptideObj newPepObj = null;
        try {
            float[] predMZs = predictionEntry.mzs;
            float[] predIntensities = predictionEntry.intensities;
            float predRT = predictionEntry.RT;
            float predIM = predictionEntry.IM;

            //filter out predicted fragments if not in range
            //TODO can make this faster by only comparing beginning and end
            ArrayList<Float> tmpMZs = new ArrayList<>();
            ArrayList<Float> tmpIntensities = new ArrayList<>();
            for (int i = 0; i < predMZs.length; i++) {
                if ((predMZs[i] >= lowerLimit) & (predMZs[i] <= upperLimit)) {
                    tmpMZs.add(predMZs[i]);
                    tmpIntensities.add(predIntensities[i]);
                }
            }
            predMZs = new float[tmpMZs.size()];
            predIntensities = new float[tmpIntensities.size()];
            for (int i = 0; i < tmpMZs.size(); i++) {
                predMZs[i] = tmpMZs.get(i);
                predIntensities[i] = tmpIntensities.get(i);
            }

            if (predMZs.length > 1) {
                if (Constants.divideFragments.equals("0")) {
                    newPepObj = new PeptideObj(this, name.baseCharge, rank, targetORdecoy, escore, predMZs,
                            predIntensities, predRT, predIM);
                } else {
                    newPepObj = new PeptideObj(this, name.baseCharge, rank, targetORdecoy, escore, predMZs,
                            predIntensities, predRT, predIM, predictionEntry.fragmentIonTypes);
                }
            } else { //only 1 frag to match
                newPepObj = new PeptideObj(this, name.baseCharge, rank, targetORdecoy, escore, zeroFloatArray,
                        zeroFloatArray, predRT, predIM);
            }
            if (set) {
                peptideObjects.add(newPepObj);
            }

            //remove higher ranked peaks
            //TODO: base off fragger.params topN?
            if (Constants.removeRankPeaks) {
                for (int i = newPepObj.spectralSimObj.matchedIdx.size() - 1; i > -1; i--) {
                    expMZs = ArrayUtils.remove(expMZs, i);
                    expIntensities = ArrayUtils.remove(expIntensities, i);
                }
            }
            newPepObj.spectralSimObj.matchedIdx.clear();
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
            if (PeptideSkipper.skipPeptide(name.stripped, name.charge)) {
                newPepObj = new PeptideObj(this, name.baseCharge, rank, targetORdecoy, escore,
                        zeroFloatArray, zeroFloatArray, predRT, predIM);
                if (set) {
                    peptideObjects.add(newPepObj);
                }
            } else {
                String[] periodSplit = Constants.spectraRTPredFile.split("\\.");
                if (periodSplit[periodSplit.length - 1].equals("dlib")) { //won't always include every entry
                    newPepObj = new PeptideObj(this, name.baseCharge, rank, targetORdecoy, escore,
                            zeroFloatArray, zeroFloatArray, predRT, predIM);
                    peptideObjects.add(newPepObj);
                } else {
                    printError("Prediction missing in file for " + name.baseCharge);
                    e.printStackTrace();
                    System.exit(-1);
                }
            }
        }
        return newPepObj;
    }

    public PeptideObj getPeptideObject(int rank) {
        return peptideObjects.get(rank - 1);
    } //arraylist now, not hashmap

    public PeptideObj getPeptideObject(String name) {
        PeptideObj returnedP = null;

        for (PeptideObj p : peptideObjects) {
            if (p.name.equals(name)) {
                return p;
            }
        }
        printInfo("no peptideObj with name " + name + ", returning null");
        return returnedP; //only if peptide name not there
    }
}
