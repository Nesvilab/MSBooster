package Features;

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
    HashMap<String, Float> predIntensities = makeBaseMap();
    HashMap<String, Integer> peakCounts = makeBaseMapInt();

    spectrumComparison spectralSimObj;

    public peptideObj(mzmlScanNumber scanNumObj, String name, int rank, int targetORdecoy, String escore,
                      float[] predMZs, float[] predIntensities, float predRT, Float predIM,
                      String[] fragmentIonTypes) {
        this.name = name;
        this.charge = Integer.parseInt(name.split("\\|")[1]);
        this.rank = rank;
        this.scanNumObj = scanNumObj;
        this.scanNum = scanNumObj.scanNum;
        this.targetORdecoy = targetORdecoy;
        this.escore = escore;
        //TODO: provide with fragmentIonTypes
        this.spectralSimObj = new spectrumComparison(scanNumObj.getExpMZs(), scanNumObj.getExpIntensities(),
                predMZs, predIntensities, Constants.useTopFragments, Constants.useBasePeak);
        this.RT = predRT;
        this.IM = predIM;
        if (Constants.useMatchedIntensities || Constants.usePredIntensities || Constants.usePeakCounts) {
            if (fragmentIonTypes != null) {
                makeFragmentAnnotationFeatures(fragmentIonTypes, predMZs);
            }
        }
    }

    private HashMap<String, Float> makeBaseMap() {
        HashMap<String, Float> map = new HashMap<>();
        for (String s : MassCalculator.allowedFragmentIonTypes) {
            map.put(s, 0f);
        }
        return map;
    }
    private HashMap<String, Integer> makeBaseMapInt() {
        HashMap<String, Integer> map = new HashMap<>();
        for (String s : MassCalculator.allowedFragmentIonTypes) {
            map.put(s, 0);
        }
        return map;
    }

    //how to deal with this if ignored fragment ions types, so matchedIntensities and fragmentIonTypes not same length?
    //save masscalculator and annotateMZs
    private void makeFragmentAnnotationFeatures(String[] fragmentIonTypes, float[] predMZs) {
        float totalmatchedIntensity = 0f;
        if (Constants.useMatchedIntensities) {
            for (float f : spectralSimObj.matchedIntensities) {
                totalmatchedIntensity += f;
            }
        }
        float totalpredIntensity = 0f;
        if (Constants.usePredIntensities) {
            for (float f : spectralSimObj.predIntensities) {
                totalpredIntensity += f;
            }
        }

        //get intensity for each fragment type
        if (fragmentIonTypes.length != spectralSimObj.predMZs.length) {
            //for mz in spectralSimObj.predMZs, if equal to mz in predMZs, get that index of fragment ion types
            String[] newFragmentIonTypes = new String[spectralSimObj.predMZs.length];
            int index = 0;
            for (float mz : spectralSimObj.predMZs) {
                for (int i = 0; i < predMZs.length; i++) {
                    if (predMZs[i] == mz) {
                        newFragmentIonTypes[index] = fragmentIonTypes[i];
                        index += 1;
                        break;
                    }
                }
            }
            fragmentIonTypes = newFragmentIonTypes;
        }

        for (int i = 0; i < fragmentIonTypes.length; i++) {
            String presentType = fragmentIonTypes[i];
            if (Constants.useMatchedIntensities) {
                if (totalmatchedIntensity == 0f) {
                    matchedIntensities.put(presentType, 0f);
                } else {
                    matchedIntensities.put(presentType,
                            matchedIntensities.get(presentType) +
                                    (spectralSimObj.matchedIntensities[i] / totalmatchedIntensity));
                }
            }
            if (Constants.usePredIntensities) {
                predIntensities.put(presentType,
                        predIntensities.get(presentType) +
                                (spectralSimObj.predIntensities[i] / totalpredIntensity));
            }
            if (Constants.usePeakCounts) {
                peakCounts.put(presentType, peakCounts.get(presentType) + 1);
            }
        }
    }
}
