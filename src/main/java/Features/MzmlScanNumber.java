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

import External.PeptideSkipper;
import org.apache.commons.lang3.ArrayUtils;
import umich.ms.datatypes.scan.IScan;
import umich.ms.datatypes.scan.props.PrecursorInfo;
import umich.ms.datatypes.spectrum.ISpectrum;
import umich.ms.fileio.exceptions.FileParsingException;

import java.util.*;
import java.util.concurrent.ConcurrentHashMap;

import static Features.FloatUtils.doubleToFloat;

public class MzmlScanNumber {
    final int scanNum;
    public double isolationLower;
    public double isolationUpper;
    public float[] expMZs;
    public float[] expIntensities;
    public float[] savedExpMZs;
    public float[] savedExpIntensities;
    float RT;
    Double calibratedRT;
    int RTbinSize;
    float normalizedRT;
    Float IM;
    int IMbinSize;
    HashMap<String, Float> NCEs = new HashMap<>();
    ArrayList<PeptideObj> peptideObjects = new ArrayList<>();
    //double[] mzFreqs;
    public static float[] zeroFloatArray = new float[]{0};
    private static final HashSet<String> nceModels =
            new HashSet<>(Arrays.asList("PredFull", "Prosit", "PrositTMT", "alphapeptdeep"));

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
        this.savedExpMZs = this.expMZs;
        this.savedExpIntensities = this.expIntensities;
        this.RT = scan.getRt().floatValue();
        if (Constants.useIM) {
            this.IM = scan.getIm().floatValue();
        }
        if (! Constants.NCE.equals("") && ! Constants.FragmentationType.equals("")) {
            NCEs.put(Constants.FragmentationType, Float.parseFloat(Constants.NCE));
        } else {
            //decide if we read in NCE and fragment info
            boolean read = false;
            for (String substring : nceModels) {
                if (Constants.spectraRTPredModel.contains(substring)) {
                    read = true;
                    break;
                }
            }
            if (read || Constants.useKoina) {
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
                    System.out.println("mzml file does not contain filter string. Setting NCE to 25 and " +
                            "FragmentationType to HCD. If other values are desired, please set them in " +
                            "the parameter file with NCE and FragmentationType.");
                    Constants.NCE = "25";
                    Constants.FragmentationType = "HCD";
                }
            }
        }
    }

    //this version if creating from mgf file
    public MzmlScanNumber(int scanNum, float[] expMZs, float[] expInts, float RT, float IM) {
        this.scanNum = scanNum;
        this.expMZs = expMZs;
        this.expIntensities = expInts;
        this.RT = RT;
        this.IM = IM;
    }

    public float[] getExpMZs() { return expMZs; }
    public float[] getExpIntensities() { return expIntensities; }

    public void setPeptideObject(PeptideFormatter name, int rank, int targetORdecoy, String escore,
                                 ConcurrentHashMap<String, PredictionEntry> allPreds) {

        if (rank != peptideObjects.size() + 1) { //need to add entries in order
            throw new AssertionError("must add next rank");
        }
        PredictionEntry predictionEntry = allPreds.get(name.baseCharge);
        try {
            float[] predMZs = predictionEntry.mzs;
            float[] predIntensities = predictionEntry.intensities;
            float predRT = predictionEntry.RT;
            float predIM = predictionEntry.IM;

            PeptideObj newPepObj;
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
            if (PeptideSkipper.skipPeptide(name.stripped, name.charge)) {
                peptideObjects.add(rank - 1, new PeptideObj(this, name.baseCharge, rank, targetORdecoy, escore,
                        zeroFloatArray, zeroFloatArray, predRT, predIM));
            } else {
                String[] periodSplit = Constants.spectraRTPredFile.split("\\.");
                if (periodSplit[periodSplit.length - 1].equals("dlib")) { //won't always include every entry
                    peptideObjects.add(rank - 1, new PeptideObj(this, name.baseCharge, rank, targetORdecoy, escore,
                            zeroFloatArray, zeroFloatArray, predRT, predIM));
                } else {
                    System.out.println("Prediction missing in file for " + name.baseCharge);
                    e.printStackTrace();
                    System.exit(-1);
                }
            }
        }
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
        System.out.println("no peptideObj with name " + name + ", returning null");
        return returnedP; //only if peptide name not there
    }
}
