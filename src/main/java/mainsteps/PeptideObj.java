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
import allconstants.FragmentIonConstants;
import features.spectra.MassCalculator;
import features.spectra.SpectrumComparison;

import java.io.IOException;
import java.net.URISyntaxException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Map;


public class PeptideObj {
    public String name;
    public int charge;
    int rank;
    public double precursorMz;
    public MzmlScanNumber scanNumObj;
    int scanNum;
    int previousScan;
    int nextScan;
    int targetORdecoy;
    public int length;
    public String escore;
    public float RT;
    public float deltaRT;
    public float deltaRTbin;
    public float RTzscore;
    public double RTprob;
    public double deltaRTLOESS;
    public double deltaRTLOESS_real;
    public double deltaRTLOESSnormalized;
    public double calibratedRT;
    public double predRTrealUnits;
    public double deltaIMLOESS;
    public double deltaIMLOESSnormalized;
    public double IMprobabilityUnifPrior;
    public double RTprobabilityUnifPrior;
    public Float IM;
    public double IMprob;
    String[] fragmentIonTypes;
    public int chromatogramWindowQuery;
    public int peptideCounts;

    static final HashMap<String, Float> baseMap = makeBaseMap();
    static private HashMap<String, Float> makeBaseMap() {
        HashMap<String, Float> map = new HashMap<>();
        for (String s : FragmentIonConstants.fragmentIonHierarchy) {
            map.put(s, 0f);
        }
        return map;
    }

    public HashMap<String, Float> matchedIntensities = baseMap;
    //HashMap<String, Float> intensitiesDifference = baseMap;
    public HashMap<String, Float> predIntensities = baseMap;
    public HashMap<String, Float> peakCounts = baseMap;
    public HashMap<String, Float> individualSpectralSimilarities = baseMap;

    public SpectrumComparison spectralSimObj;
    public double intensity_distribution_similarity;

    public PeptideObj(MzmlScanNumber scanNumObj, String name, int rank, int targetORdecoy, String escore,
                      float[] predMZs, float[] predIntensities, float predRT, Float predIM) throws IOException, URISyntaxException {
        this.name = name;
        this.charge = Integer.parseInt(name.split("\\|")[1]);
        this.rank = rank;
        this.scanNumObj = scanNumObj;
        this.scanNum = scanNumObj.scanNum;
        this.targetORdecoy = targetORdecoy;
        int length = 0;
        for (int i = 0; i < name.length() - 2; i++) {
            if (Character.isAlphabetic(name.charAt(i))) {
                length += 1;
            }
        }
        this.length = length;
        this.escore = escore;
        //this.predMZs = predMZs;
        //this.predInts = predIntensities;
        this.spectralSimObj = new SpectrumComparison(this, scanNumObj.getExpMZs(), scanNumObj.getExpIntensities(),
                predMZs, predIntensities, length); //calculate similarity with subset of fragments
        this.RT = predRT;
        this.IM = predIM;
        if (Constants.useMatchedIntensities || Constants.usePeakCounts || Constants.useIntensitiesDifference ||
                Constants.usePredIntensities || Constants.useIndividualSpectralSimilarities ||
                Constants.useIntensityDistributionSimilarity) {
            makeFragmentAnnotationFeatures(predMZs, predIntensities);
        }
    }

    public PeptideObj(MzmlScanNumber scanNumObj, String name, int rank, int targetORdecoy, String escore,
                      float[] predMZs, float[] predIntensities, float predRT, Float predIM, String[] fragmentIonTypes)
            throws IOException, URISyntaxException {
        this.name = name;
        this.charge = Integer.parseInt(name.split("\\|")[1]);
        this.rank = rank;
        this.scanNumObj = scanNumObj;
        this.scanNum = scanNumObj.scanNum;
        this.targetORdecoy = targetORdecoy;
        int length = 0;
        for (int i = 0; i < name.length() - 2; i++) {
            if (Character.isAlphabetic(name.charAt(i))) {
                length += 1;
            }
        }
        this.length = length;
        this.escore = escore;
        //this.predMZs = predMZs;
        //this.predInts = predIntensities;
        this.fragmentIonTypes = fragmentIonTypes;
        this.spectralSimObj = new SpectrumComparison(this, scanNumObj.getExpMZs(), scanNumObj.getExpIntensities(),
                predMZs, predIntensities, length, fragmentIonTypes); //calculate similarity with subset of fragments
        this.RT = predRT;
        this.IM = predIM;
        if (Constants.useMatchedIntensities || Constants.usePeakCounts || Constants.useIntensitiesDifference ||
                Constants.usePredIntensities || Constants.useIndividualSpectralSimilarities ||
                Constants.useIntensityDistributionSimilarity) {
            makeFragmentAnnotationFeatures(predMZs, predIntensities);
        }
    }

    public PeptideObj(){}

    //how to deal with this if ignored fragment ions types, so matchedIntensities and fragmentIonTypes not same length?
    //save masscalculator and annotateMZs
    private void makeFragmentAnnotationFeatures(float[] predMZs, float[] predInts) throws IOException, URISyntaxException { //TODO what if we want to use unispec?
        //filter for top fragments for all experimental and pred vectors
        ArrayList<Float> expIntensitiesList = new ArrayList<>();
        ArrayList<Float> expMZsList = new ArrayList<>();
        float maxInt = 0f;
        for (float f : scanNumObj.getExpIntensities()) {
            if (f > maxInt) {
                maxInt = f;
            }
        }
        float minInt = maxInt / 100f * Constants.percentBasePeak;
        float expTotalIntensity = 0f;
        for (int i = 0; i < scanNumObj.getExpIntensities().length; i++) {
            float f = scanNumObj.getExpIntensities()[i];
            if (f > minInt) {
                expIntensitiesList.add(f);
                expMZsList.add(scanNumObj.getExpMZs()[i]);
                expTotalIntensity += f;
            }
        }

        float[] expIntensities = new float[expIntensitiesList.size()];
        float[] expMZs = new float[expMZsList.size()];
        for (int i = 0; i < expIntensities.length; i++) {
            expIntensities[i] = expIntensitiesList.get(i) / expTotalIntensity;
            expMZs[i] = expMZsList.get(i);
        }

        //pred
        ArrayList<Float> predIntensitiesList = new ArrayList<>();
        ArrayList<Float> predMZsList = new ArrayList<>();
        ArrayList<String> predTypesList = new ArrayList<>();
        maxInt = 0f;
        for (float f : predInts) {
            if (f > maxInt) {
                maxInt = f;
            }
        }
        minInt = maxInt / 100f * Constants.percentBasePeak;

        float predTotalIntensity = 0f;
        for (int i = 0; i < predInts.length; i++) {
            float f = predInts[i];
            if (f > minInt) {
                predIntensitiesList.add(f);
                predMZsList.add(predMZs[i]);
                predTypesList.add(fragmentIonTypes[i]);
                predTotalIntensity += f;
            }
        }

        float[] predIntensities1 = new float[predIntensitiesList.size()];
        float[] predMZs1 = new float[predMZsList.size()];
        String[] predTypes1 = new String[predTypesList.size()];
        for (int i = 0; i < predIntensities1.length; i++) {
            predIntensities1[i] = predIntensitiesList.get(i) / predTotalIntensity;
            predMZs1[i] = predMZsList.get(i);
            predTypes1[i] = predTypesList.get(i);
        }

        MassCalculator mc = new MassCalculator(name.split("\\|")[0], charge);
        String[] expFragmentIonTypes = mc.annotateMZs(expMZs, "default", true)[1];
        String[] predFragmentIonTypes = predTypes1;

        for (int i = 0; i < expFragmentIonTypes.length; i++) {
            String presentType = expFragmentIonTypes[i];
            if (FragmentIonConstants.fragmentIonHierarchySet.contains(presentType)) {
                if (Constants.useMatchedIntensities || Constants.useIntensityDistributionSimilarity
                        || Constants.useIntensitiesDifference) {
                    matchedIntensities.put(presentType,
                            matchedIntensities.get(presentType) + (expIntensities[i]));
                }

                if (Constants.usePeakCounts) {
                    peakCounts.put(presentType,
                            peakCounts.get(presentType) + (1f / expIntensities.length));
                }
            }
        }

        for (int i = 0; i < predFragmentIonTypes.length; i++) {
            String presentType = predFragmentIonTypes[i];
            if (FragmentIonConstants.fragmentIonHierarchySet.contains(presentType)) {
                if (Constants.useIntensityDistributionSimilarity || Constants.useIntensitiesDifference) {
                    predIntensities.put(presentType,
                            predIntensities.get(presentType) + predIntensities1[i]);
                }
            }
        }

        if (Constants.useIndividualSpectralSimilarities) {
            for (Map.Entry<String, Float> entry : individualSpectralSimilarities.entrySet()) {
                //get just the experimental intensities belonging to appropriate fragment ion type
                ArrayList<Float> subsetPredMZs = new ArrayList<>();
                ArrayList<Float> subsetPredInts = new ArrayList<>();
                for (int i = 0; i < predFragmentIonTypes.length; i++) {
                    if (predFragmentIonTypes[i].equals(entry.getKey())) {
                        subsetPredMZs.add(predMZs1[i]);
                        subsetPredInts.add(predIntensities1[i]);
                    }
                }

                //convert arraylist to array
                float[] subsetPredMZsArray = new float[subsetPredMZs.size()];
                float[] subsetPredIntsArray = new float[subsetPredInts.size()];
                for (int i = 0; i < subsetPredMZs.size(); i++) {
                    subsetPredMZsArray[i] = subsetPredMZs.get(i);
                    subsetPredIntsArray[i] = subsetPredInts.get(i);
                }

                individualSpectralSimilarities.put(entry.getKey(),
                        (float) new SpectrumComparison(this, expMZs, expIntensities,
                                subsetPredMZsArray, subsetPredIntsArray, this.length).unweightedSpectralEntropy());
            }
        }
    }
}
