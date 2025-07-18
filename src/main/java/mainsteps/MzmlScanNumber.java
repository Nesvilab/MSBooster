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

package mainsteps;

import allconstants.Constants;
import allconstants.NceConstants;
import org.apache.commons.lang3.ArrayUtils;
import peptideptmformatting.PeptideFormatter;
import peptideptmformatting.PeptideSkipper;
import predictions.PredictionEntry;
import predictions.PredictionEntryHashMap;
import umich.ms.datatypes.scan.IScan;
import umich.ms.datatypes.scan.props.PrecursorInfo;
import umich.ms.datatypes.spectrum.ISpectrum;
import umich.ms.fileio.exceptions.FileParsingException;

import java.io.IOException;
import java.net.URISyntaxException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.concurrent.atomic.AtomicInteger;

import static utils.NumericUtils.doubleToFloat;
import static utils.NumericUtils.getRanks;
import static utils.Print.printError;
import static utils.Print.printInfo;

public class MzmlScanNumber {
    final int scanNum;
    public double ms1IsolationLowerLimit;
    public double ms1IsolationUpperLimit;
    public float[] expMZs;
    public float[] expIntensities;
    public float RT;
    public int RTbinSize;
    public float normalizedRT;
    public Float IM;
    public int IMbinSize;
    public Double ms2LowerLimit;
    public Double ms2UpperLimit;
    public ArrayList<PeptideObj> peptideObjects = new ArrayList<>();
    //double[] mzFreqs;
    public static float[] zeroFloatArray = new float[]{0};
    public static String[] zeroStringArray = new String[]{""};
    public AtomicInteger skippedPSMs = new AtomicInteger(0);

    //this version if creating from mzml scan number
    public MzmlScanNumber(IScan scan) throws FileParsingException {
        this.scanNum = scan.getNum();
        PrecursorInfo pi = scan.getPrecursor();
        if (pi.getMzRangeStart() != null) {
            this.ms1IsolationLowerLimit = pi.getMzRangeStart() == 0 ? Double.NaN : pi.getMzRangeStart();
            this.ms1IsolationUpperLimit = pi.getMzRangeEnd() == 0 ? Double.NaN : pi.getMzRangeEnd();
        }
        ISpectrum spectrum = scan.fetchSpectrum();
        this.expMZs = doubleToFloat(spectrum.getMZs());
        this.expIntensities = doubleToFloat(spectrum.getIntensities());
        sortArrays(this.expMZs, this.expIntensities);
        this.RT = scan.getRt().floatValue();
        if (Constants.useIM) {
            try {
                this.IM = scan.getIm().floatValue();
            } catch (Exception e) {
                printInfo("Ion mobility values not found in data files, setting useIM to false.");
                Constants.useIM = false;
            }
        }
        ms2LowerLimit = scan.getScanMzWindowLower();
        ms2UpperLimit = scan.getScanMzWindowUpper();
    }

    //this version if creating from mgf file
    //TODO: update support
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
                                       PredictionEntryHashMap allPreds, boolean set) throws IOException, URISyntaxException {
        PredictionEntry predictionEntry = allPreds.get(name.getBaseCharge());
        PeptideObj newPepObj = null;
        try {
            //merge spectra and aux spectra
            float[] predMZs;
            float[] predIntensities;
            String[] predFragmentIonTypes;

            if (predictionEntry.auxSpectra != null) {
                List<Integer> rankListMain = new ArrayList<>();
                List<Integer> rankListAux = new ArrayList<>();

                getRanks(predictionEntry.mzs, predictionEntry.auxSpectra.mzs, rankListMain, rankListAux);

                predMZs = new float[rankListMain.size() + rankListAux.size()];
                predIntensities = new float[rankListMain.size() + rankListAux.size()];
                predFragmentIonTypes = new String[rankListMain.size() + rankListAux.size()];

                for (int i = 0; i < rankListMain.size(); i++) {
                    predMZs[rankListMain.get(i)] = predictionEntry.mzs[i];
                    predIntensities[rankListMain.get(i)] = predictionEntry.intensities[i];
                    predFragmentIonTypes[rankListMain.get(i)] = predictionEntry.fragmentIonTypes[i];
                }
                for (int i = 0; i < rankListAux.size(); i++) {
                    predMZs[rankListAux.get(i)] = predictionEntry.auxSpectra.mzs[i];
                    predIntensities[rankListAux.get(i)] = predictionEntry.auxSpectra.intensities[i];
                    predFragmentIonTypes[rankListAux.get(i)] = predictionEntry.auxSpectra.fragmentIonTypes[i];
                }
            } else {
                predMZs = predictionEntry.mzs;
                predIntensities = predictionEntry.intensities;
                predFragmentIonTypes = predictionEntry.fragmentIonTypes;
            }

            float predRT = predictionEntry.RT;
            float predIM = predictionEntry.IM;

            //filter out predicted fragments if not in range
            //TODO can make this faster by only comparing beginning and end
            ArrayList<Float> tmpMZs = new ArrayList<>();
            ArrayList<Float> tmpIntensities = new ArrayList<>();
            ArrayList<String> tmpFragmentIonTypes = new ArrayList<>();
            for (int i = 0; i < predMZs.length; i++) {
                if ((predMZs[i] >= ms2LowerLimit) & (predMZs[i] <= ms2UpperLimit)) {
                    tmpMZs.add(predMZs[i]);
                    tmpIntensities.add(predIntensities[i]);
                    tmpFragmentIonTypes.add(predFragmentIonTypes[i]);
                }
            }
            predMZs = new float[tmpMZs.size()];
            predIntensities = new float[tmpIntensities.size()];
            predFragmentIonTypes = new String[tmpFragmentIonTypes.size()];
            for (int i = 0; i < tmpMZs.size(); i++) {
                predMZs[i] = tmpMZs.get(i);
                predIntensities[i] = tmpIntensities.get(i);
                predFragmentIonTypes[i] = tmpFragmentIonTypes.get(i);
            }

            if (predMZs.length > 1) {
                newPepObj = new PeptideObj(this, name.getBaseCharge(), rank, targetORdecoy, escore,
                        predMZs, predIntensities, predFragmentIonTypes, predRT, predIM, predictionEntry.daltonMatching);
            } else { //only 1 frag to match
                newPepObj = new PeptideObj(this, name.getBaseCharge(), rank, targetORdecoy, escore,
                        zeroFloatArray, zeroFloatArray, zeroStringArray, predRT, predIM, predictionEntry.daltonMatching);
            }
            if (set) {
                peptideObjects.add(newPepObj);
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
            if (PeptideSkipper.skipPeptide(name.getStripped(), name.getCharge(),
                    Constants.spectraModel + Constants.rtModel + Constants.imModel)) {
                newPepObj = new PeptideObj(this, name.getBaseCharge(), rank, targetORdecoy, escore,
                        zeroFloatArray, zeroFloatArray, zeroStringArray, predRT, predIM, predictionEntry.daltonMatching);
                if (set) {
                    peptideObjects.add(newPepObj);
                }
                skippedPSMs.incrementAndGet();
            } else {
                //TODO: if we start using dlib again
//                String[] periodSplit = Constants.spectraRTPredFile.split("\\.");
//                if (periodSplit[periodSplit.length - 1].equals("dlib")) { //won't always include every entry
//                    newPepObj = new PeptideObj(this, name.baseCharge, rank, targetORdecoy, escore,
//                            zeroFloatArray, zeroFloatArray, predRT, predIM);
//                    peptideObjects.add(newPepObj);
//                } else {
//                    printError("Prediction missing in file for " + name.baseCharge);
//                    e.printStackTrace();
//                    System.exit(1);
//                }
                printError("Prediction missing in file for " + name.getBaseCharge());
                e.printStackTrace();
                System.exit(1);
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
