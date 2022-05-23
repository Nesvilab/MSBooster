package Features;

import java.util.*;

public class peptideObj {
    final String name;
    final int charge;
    final int rank;
    final mzmlScanNumber scanNumObj;
    final int scanNum;
    final int targetORdecoy;
    final String escore;
    final float RT;
    float deltaRT;
    float deltaRTbin;
    float RTzscore;
    double RTprob;
    double deltaRTLOESS;
    double deltaRTLOESSnormalized;
    double deltaIMLOESS;
    double deltaIMLOESSnormalized;
    Float IM;
    double IMprob;
    float[] predMZs;
    float[] predInts;

    HashMap<String, Float> matchedIntensities = makeBaseMap();
    //HashMap<String, Float> intensitiesDifference = makeBaseMap();
    HashMap<String, Float> predIntensities = makeBaseMap();
    HashMap<String, Float> peakCounts = makeBaseMap();
    HashMap<String, Float> individualSpectralSimilarities = makeBaseMap();

    spectrumComparison spectralSimObj;

    public peptideObj(mzmlScanNumber scanNumObj, String name, int rank, int targetORdecoy, String escore,
                      float[] predMZs, float[] predIntensities, float predRT, Float predIM) {
        this.name = name;
        this.charge = Integer.parseInt(name.split("\\|")[1]);
        this.rank = rank;
        this.scanNumObj = scanNumObj;
        this.scanNum = scanNumObj.scanNum;
        this.targetORdecoy = targetORdecoy;
        this.escore = escore;
        this.predMZs = predMZs;
        this.predInts = predIntensities;
        this.spectralSimObj = new spectrumComparison(scanNumObj.getExpMZs(), scanNumObj.getExpIntensities(),
                predMZs, predIntensities, Constants.useTopFragments, Constants.useBasePeak);
        this.RT = predRT;
        this.IM = predIM;
        if (Constants.useMatchedIntensities || Constants.usePeakCounts || Constants.useIntensitiesDifference ||
                Constants.usePredIntensities || Constants.useIndividualSpectralSimilarities ||
                Constants.useIntensityDistributionSimilarity) {
            makeFragmentAnnotationFeatures();
        }
    }

    private HashMap<String, Float> makeBaseMap() {
        HashMap<String, Float> map = new HashMap<>();
        for (String s : Constants.fragmentIonHierarchy) {
            map.put(s, 0f);
        }
        return map;
    }

    //how to deal with this if ignored fragment ions types, so matchedIntensities and fragmentIonTypes not same length?
    //save masscalculator and annotateMZs
    private void makeFragmentAnnotationFeatures() {
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
        float totalPeaks = expIntensitiesList.size();
        float[] expIntensities = new float[expIntensitiesList.size()];
        float[] expMZs = new float[expMZsList.size()];
        for (int i = 0; i < expIntensities.length; i++) {
            expIntensities[i] = expIntensitiesList.get(i) / expTotalIntensity;
            expMZs[i] = expMZsList.get(i);
        }

        //pred
        ArrayList<Float> predIntensitiesList = new ArrayList<>();
        ArrayList<Float> predMZsList = new ArrayList<>();
//        maxInt = 0f;
//        for (float f : predInts) {
//            if (f > maxInt) {
//                maxInt = f;
//            }
//        }
//        minInt = maxInt / 100f * Constants.percentBasePeak;
        minInt = 10f * Constants.percentBasePeak; //change if we use anything besides predfull
        float predTotalIntensity = 0f;
        for (int i = 0; i < predInts.length; i++) {
            float f = predInts[i];
            if (f > minInt) {
                predIntensitiesList.add(f);
                predMZsList.add(predMZs[i]);
                predTotalIntensity += f;
            }
        }

        float[] predIntensities1 = new float[predIntensitiesList.size()];
        float[] predMZs1 = new float[predMZsList.size()];
        for (int i = 0; i < predIntensities1.length; i++) {
            predIntensities1[i] = predIntensitiesList.get(i) / predTotalIntensity;
            predMZs1[i] = predMZsList.get(i);
        }

        MassCalculator mc = new MassCalculator(name.split("\\|")[0], charge);
        String[] expFragmentIonTypes = mc.annotateMZs(expMZs)[1];
        String[] predFragmentIonTypes = mc.annotateMZs(predMZs1)[1];

        List<String> fragmentIonHierarchyList = Arrays.asList(Constants.fragmentIonHierarchy);
        for (int i = 0; i < expFragmentIonTypes.length; i++) {
            String presentType = expFragmentIonTypes[i];
            if (fragmentIonHierarchyList.contains(presentType)) {
                if (Constants.useMatchedIntensities || Constants.useIntensityDistributionSimilarity
                        || Constants.useIntensitiesDifference) {
                    matchedIntensities.put(presentType,
                            matchedIntensities.get(presentType) + (expIntensities[i]));
                }

                if (Constants.usePeakCounts) {
                    peakCounts.put(presentType,
                            peakCounts.get(presentType) + (1f / totalPeaks));
                }
            }
        }

        for (int i = 0; i < predFragmentIonTypes.length; i++) {
            String presentType = predFragmentIonTypes[i];
            if (fragmentIonHierarchyList.contains(presentType)) {
                if (Constants.useIntensityDistributionSimilarity || Constants.useIntensitiesDifference) {
                    predIntensities.put(presentType,
                            predIntensities.get(presentType) + predIntensities1[i]);
                }
            }
        }
        //intensitiesDifference.replaceAll((k, v) -> Math.abs(v));

//        if (Constants.usePredIntensities) {
//            float totalPredIntensity = 0f;
//            for (float f : spectralSimObj.predIntensities) {
//                totalPredIntensity += f;
//            }
//
//            fragmentIonTypes = mc.annotateMZs(spectralSimObj.predMZs)[1];
//            for (int i = 0; i < fragmentIonTypes.length; i++) {
//                String presentType = fragmentIonTypes[i];
//                predIntensities.put(presentType,
//                        predIntensities.get(presentType) + (spectralSimObj.predIntensities[i] / totalPredIntensity));
//            }
//        }

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
                        (float) new spectrumComparison(expMZs, expIntensities,
                                subsetPredMZsArray, subsetPredIntsArray,
                                Constants.useTopFragments, false).brayCurtis()); //already did base peak intensity filtering
            }
        }
    }
}
