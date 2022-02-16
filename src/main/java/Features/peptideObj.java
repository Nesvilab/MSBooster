package Features;

import java.util.ArrayList;
import java.util.HashMap;

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

    HashMap<String, Float> matchedIntensities = makeBaseMap();
    HashMap<String, Float> peakCounts = makeBaseMap();

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
        this.spectralSimObj = new spectrumComparison(scanNumObj.getExpMZs(), scanNumObj.getExpIntensities(),
                predMZs, predIntensities, Constants.useTopFragments, Constants.useBasePeak);
        this.RT = predRT;
        this.IM = predIM;
        if (Constants.useMatchedIntensities || Constants.usePeakCounts) {
            makeFragmentAnnotationFeatures();
        }
    }

    private HashMap<String, Float> makeBaseMap() {
        HashMap<String, Float> map = new HashMap<>();
        for (String s : MassCalculator.allowedFragmentIonTypes) {
            map.put(s, 0f);
        }
        return map;
    }

    //how to deal with this if ignored fragment ions types, so matchedIntensities and fragmentIonTypes not same length?
    //save masscalculator and annotateMZs
    private void makeFragmentAnnotationFeatures() {
        ArrayList<Float> expIntensitiesList = new ArrayList<>();
        ArrayList<Float> expMZsList = new ArrayList<>();
        float maxInt = 0f;
        for (float f : scanNumObj.getExpIntensities()) {
            if (f > maxInt) {
                maxInt = f;
            }
        }
        float minInt = maxInt / 100f * Constants.percentBasePeakExperimental;
        float totalIntensity = 0f;
        for (int i = 0; i < scanNumObj.getExpIntensities().length; i++) {
            float f = scanNumObj.getExpIntensities()[i];
            if (f > minInt) {
                expIntensitiesList.add(f);
                expMZsList.add(scanNumObj.getExpMZs()[i]);
                totalIntensity += f;
            }
        }
        float totalPeaks = expIntensitiesList.size();
        float[] expIntensities = new float[expIntensitiesList.size()];
        float[] expMZs = new float[expMZsList.size()];
        for (int i = 0; i < expIntensities.length; i++) {
            expIntensities[i] = expIntensitiesList.get(i);
            expMZs[i] = expMZsList.get(i);
        }

        MassCalculator mc = new MassCalculator(name.split("\\|")[0], charge);
        String[] fragmentIonTypes = mc.annotateMZs(expMZs)[1];

        for (int i = 0; i < fragmentIonTypes.length; i++) {
            String presentType = fragmentIonTypes[i];
            if (Constants.useMatchedIntensities) {
                matchedIntensities.put(presentType,
                        matchedIntensities.get(presentType) + (expIntensities[i] / totalIntensity));
            }
            if (Constants.usePeakCounts) {
                peakCounts.put(presentType,
                        peakCounts.get(presentType) + (1f / totalPeaks));
            }
        }
    }
}
