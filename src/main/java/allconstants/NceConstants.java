package allconstants;

import utils.CaseInsensitiveHashSet;

import java.util.HashMap;

public class NceConstants implements ConstantsInterface {
    public static CaseInsensitiveHashSet nceModels = new CaseInsensitiveHashSet(
            new String[] {
                    "PredFull", "Prosit", "PrositTMT", "alphapeptdeep", //externally supported
                    "AlphaPept_ms2_generic",
                    "Prosit_2019_intensity", "Prosit_2023_intensity_timsTOF",
                    "Prosit_2020_intensity_TMT", "Prosit_2020_intensity_HCD",
                    "UniSpec", "PredFull"});

    //these are default NCEs
    //set value as string because we often write it as a string in files
    //TODO: revisit when using mzmls with multiple fragmentation types
    public static HashMap<String, String> mzmlNCEs = new HashMap<>(); //key: fragmentation type, value: NCE
    public static String getNCE() {
        return mzmlNCEs.get(Constants.FragmentationType); //this is the default
    }

    public static Integer minNCE = 20;
    public static Integer maxNCE = 40;
    public static Boolean calibrateNCE = true;

    //these are optimized NCEs for each model
    //once models find their optimal NCE, place them here so they only have to be calibrated once
    public static HashMap<String, String> calibratedModels = new HashMap<>(); //key: model, value: NCE
    public static String getCalibratedNCE(String model) {
        String NCE = calibratedModels.get(model);
        if (NCE == null) {
            return getNCE();
        } else {
            return NCE;
        }
    }
}
