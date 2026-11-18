package allconstants;

import koinaclasses.KoinaException;
import koinaclasses.KoinaModelRegistry;
import koinaclasses.KoinaModelSpec;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;
import utils.CaseInsensitiveHashSet;

import java.util.HashMap;

public class NceConstants implements ConstantsInterface {
    private static final Logger log = LoggerFactory.getLogger(NceConstants.class);
    //models whose NCE is calibrated: the externally supported ones and the Koina models whose spec says so
    public static CaseInsensitiveHashSet nceModels = makeNceModels();
    private static CaseInsensitiveHashSet makeNceModels() {
        CaseInsensitiveHashSet models = new CaseInsensitiveHashSet(
                new String[] {"Prosit", "PrositTMT", "alphapeptdeep"}); //externally supported
        for (KoinaModelSpec spec : KoinaModelRegistry.all()) {
            if (spec.calibrateCe()) {
                models.add(spec.name());
            }
        }
        return models;
    }

    //these are default NCEs
    //set value as string because we often write it as a string in files
    //TODO: revisit when using mzmls with multiple fragmentation types
    public static HashMap<String, String> mzmlNCEs = new HashMap<>(); //key: fragmentation type, value: NCE
    public static String getNCE() {
        String NCE = mzmlNCEs.get(Constants.FragmentationType);
        if (NCE == null) { //Fragmentation type had to be switched
            String frag = mzmlNCEs.keySet().iterator().next();
            return mzmlNCEs.get(frag);
        }
        return NCE; //this is the default
    }

    public static Integer NCE = 25; //generic value used
    //the NCE range the calibration tries is the model's own (KoinaModelSpec minCe and maxCe). minNCE and maxNCE
    //are null unless set in the parameters; set, they replace that end of the range for every model
    public static Integer minNCE = null;
    public static Integer maxNCE = null;
    public static Boolean calibrateNCE = true;

    /** The lowest NCE the calibration tries for this model: minNCE if it is set, else the model's own. */
    public static int getMinNCE(String model) {
        if (minNCE != null) {
            return minNCE;
        }
        KoinaModelSpec spec = KoinaModelRegistry.find(model);
        return spec != null ? spec.minCe() : KoinaModelSpec.DEFAULT_MIN_CE;
    }

    /** The highest NCE the calibration tries for this model: maxNCE if it is set, else the model's own. */
    public static int getMaxNCE(String model) {
        if (maxNCE != null) {
            return maxNCE;
        }
        KoinaModelSpec spec = KoinaModelRegistry.find(model);
        return spec != null ? spec.maxCe() : KoinaModelSpec.DEFAULT_MAX_CE;
    }

    /**
     * Ends the run if this model's NCE is to be calibrated over an empty range: a minNCE above the model's highest
     * NCE, a maxNCE below its lowest, or a minNCE above the maxNCE. Only a Koina model whose spec calibrates its NCE
     * is calibrated; any other model, and no model (an empty name), pass.
     *
     * @throws KoinaException saying where each end of the range comes from
     */
    public static void requireCalibrationRange(String model) {
        KoinaModelSpec spec = model == null ? null : KoinaModelRegistry.find(model);
        if (!calibrateNCE || spec == null || !spec.calibrateCe()) {
            return;
        }
        int min = getMinNCE(model);
        int max = getMaxNCE(model);
        if (min <= max) {
            return;
        }
        throw new KoinaException("The NCE range to calibrate " + model + " over is empty: from " + min
                + (minNCE != null ? " (minNCE)" : " (the lowest NCE of " + model + ")") + " to " + max
                + (maxNCE != null ? " (maxNCE)" : " (the highest NCE of " + model + ")")
                + ". Set minNCE and maxNCE so that minNCE is not above maxNCE, or leave out the one that is not"
                + " needed");
    }

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
