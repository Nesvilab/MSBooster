package allconstants;

import koinaclasses.KoinaModelRegistry;
import koinaclasses.KoinaModelSpec;
import utils.CaseInsensitiveHashSet;
import utils.Print;

import java.util.*;
import java.util.function.Predicate;

public class ModelCollections implements ConstantsInterface {
    ////////////all Koina models, from koinaclasses.KoinaModelRegistry/////////////
    public static CaseInsensitiveHashSet KoinaRTmodels = koinaModels(spec -> spec.kind() == KoinaModelSpec.Kind.RT);
    public static CaseInsensitiveHashSet KoinaMS2models = koinaModels(spec -> spec.kind() == KoinaModelSpec.Kind.MS2);
    public static CaseInsensitiveHashSet KoinaIMmodels = koinaModels(spec -> spec.kind() == KoinaModelSpec.Kind.IM);
    public static CaseInsensitiveHashSet KoinaCCSmodels = koinaModels(spec -> spec.kind() == KoinaModelSpec.Kind.IM
            && spec.imConversion().unit() == KoinaModelSpec.ImUnit.CCS); // convert from ccs to 1/K0
    public static CaseInsensitiveHashSet KoinaModels = koinaModels(spec -> true);
    public static CaseInsensitiveHashSet KoinaIsoLabelmodels = koinaModels(KoinaModelSpec::isobaricLabel);

    private static CaseInsensitiveHashSet koinaModels(Predicate<KoinaModelSpec> filter) {
        CaseInsensitiveHashSet models = new CaseInsensitiveHashSet();
        for (KoinaModelSpec spec : KoinaModelRegistry.all()) {
            if (filter.test(spec)) {
                models.add(spec.name());
            }
        }
        return models;
    }

    ////////////individual collections/////////////
    //FragCast is added alongside DIA-NN as a findBest candidate (it does not replace DIA-NN).
    public static ArrayList<String> generalRTmodels = new ArrayList<>(
            List.of("DIA-NN", "FragCast", "AlphaPept_rt_generic", "Prosit_2019_irt", "Deeplc_hela_hf"));
    public static ArrayList<String> isolabelRTmodels = new ArrayList<>(
            List.of("DIA-NN", "FragCast", "Prosit_2020_irt_TMT", "AlphaPept_rt_generic", "Deeplc_hela_hf"));

    public static ArrayList<String> generalMS2models = new ArrayList<>(
            List.of("DIA-NN", "FragCast", "ms2pip_2021_HCD", "ms2pip_TTOF5600", "AlphaPept_ms2_generic",
                    "Prosit_2020_intensity_CID", "Prosit_2020_intensity_HCD"));
    public static ArrayList<String> isolabelMS2models = new ArrayList<>(
            List.of("Prosit_2020_intensity_TMT", "ms2pip_CID_TMT", "ms2pip_iTRAQphospho"));
    public static ArrayList<String> hlaMS2models = new ArrayList<>(
            List.of("DIA-NN", "FragCast", "ms2pip_Immuno_HCD", "AlphaPept_ms2_generic",
                    "Prosit_2020_intensity_CID", "Prosit_2020_intensity_HCD"));
    public static ArrayList<String> timstofMS2models = new ArrayList<>(
            List.of("DIA-NN", "FragCast", "ms2pip_timsTOF2024", "AlphaPept_ms2_generic", "Prosit_2023_intensity_timsTOF"));

    public static ArrayList<String> generalIMmodels = new ArrayList<>(
            List.of("DIA-NN", "FragCast", "AlphaPept_ccs_generic", "IM2Deep"));

    public static String rtCollection = "general";
    public static String ms2Collection = "general";
    public static String imCollection = "general";

    public static ArrayList<String> getRtCollection() {
        switch (rtCollection) {
            case "general":
                return generalRTmodels;
            case "isolabel":
                return isolabelRTmodels;
            default:
                Print.printError(rtCollection + " is not a supported collection. Exiting");
                System.exit(1);
                return new ArrayList<>();
        }
    }
    public static ArrayList<String> getMs2Collection() {
        switch (ms2Collection) {
            case "general":
                return generalMS2models;
            case "isolabel":
                return isolabelMS2models;
            case "hla":
                return hlaMS2models;
            case "timstof":
                return timstofMS2models;
            default:
                Print.printError(ms2Collection + " is not a supported collection. Exiting");
                System.exit(1);
                return new ArrayList<>();
        }
    }
    public static ArrayList<String> getImCollection() {
        switch (imCollection) {
            case "general":
                return generalIMmodels;
            default:
                Print.printError(imCollection + " is not a supported collection. Exiting");
                System.exit(1);
                return new ArrayList<>();
        }
    }
    ////////////user-specified, will override collection/////////////
    public static String rtSearchModelsString = "";
    public static String ms2SearchModelsString = "";
    public static String imSearchModelsString = "";

    ////////////allowed fragmentation types for models that require it, from the registry/////////////
    public static HashMap<String, HashSet<String>> allowedFragmentationTypes = new HashMap<>();
    static {
        for (KoinaModelSpec spec : KoinaModelRegistry.all()) {
            if (!spec.allowedFragmentationTypes().isEmpty()) {
                allowedFragmentationTypes.put(spec.name(), new HashSet<>(spec.allowedFragmentationTypes()));
            }
        }
    }
}
