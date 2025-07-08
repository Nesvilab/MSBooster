package allconstants;

import utils.CaseInsensitiveHashSet;
import utils.Print;

import java.util.*;

//TODO ms2pip tmt and itraq phospho models need to be corrected on koina
public class ModelCollections implements ConstantsInterface {
    ////////////all Koina models/////////////
    public static CaseInsensitiveHashSet KoinaRTmodels = new CaseInsensitiveHashSet(
            new String[] {
                    "AlphaPept_rt_generic",
                    "Prosit_2019_irt", "Prosit_2020_irt_TMT", "Prosit_2024_irt_cit",
                    "Deeplc_hela_hf"});
    public static CaseInsensitiveHashSet KoinaMS2models = new CaseInsensitiveHashSet(
            new String[] {
                    "ms2pip_2021_HCD", "ms2pip_timsTOF2024", "ms2pip_CID_TMT", "ms2pip_TTOF5600",
                    "ms2pip_Immuno_HCD", "ms2pip_iTRAQphospho",
                    "AlphaPept_ms2_generic",
                    "Prosit_2019_intensity", "Prosit_2020_intensity_CID", "Prosit_2020_intensity_TMT",
                    "Prosit_2020_intensity_HCD", "Prosit_2023_intensity_timsTOF", "Prosit_2024_intensity_cit",
                    "Prosit_2025_intensity_MultiFrag",
                    "UniSpec", "PredFull"});
    public static CaseInsensitiveHashSet KoinaIMmodels = new CaseInsensitiveHashSet(
            new String[] {"AlphaPept_ccs_generic", "IM2Deep"});
    public static CaseInsensitiveHashSet KoinaCCSmodels = new CaseInsensitiveHashSet(
            new String[] {"AlphaPept_ccs_generic", "IM2Deep"}); // convert from ccs to 1/K0
    public static CaseInsensitiveHashSet KoinaModels = new CaseInsensitiveHashSet();
    static {
        KoinaModels.addAll(KoinaRTmodels);
        KoinaModels.addAll(KoinaMS2models);
        KoinaModels.addAll(KoinaIMmodels);
    }
    public static CaseInsensitiveHashSet KoinaIsoLabelmodels = new CaseInsensitiveHashSet(
            new String[] {"Prosit_2020_irt_TMT", "Prosit_2020_intensity_TMT"});

    ////////////individual collections/////////////
    public static ArrayList<String> generalRTmodels = new ArrayList<>(
            List.of("DIA-NN", "AlphaPept_rt_generic", "Prosit_2019_irt", "Deeplc_hela_hf"));
    public static ArrayList<String> isolabelRTmodels = new ArrayList<>(
            List.of("DIA-NN", "Prosit_2020_irt_TMT", "AlphaPept_rt_generic", "Deeplc_hela_hf"));

    public static ArrayList<String> generalMS2models = new ArrayList<>(
            List.of("DIA-NN", "ms2pip_2021_HCD", "ms2pip_TTOF5600", "AlphaPept_ms2_generic",
                    "Prosit_2020_intensity_CID", "Prosit_2020_intensity_HCD"));
    public static ArrayList<String> isolabelMS2models = new ArrayList<>(
            List.of("Prosit_2020_intensity_TMT", "ms2pip_CID_TMT", "ms2pip_iTRAQphospho"));
    public static ArrayList<String> hlaMS2models = new ArrayList<>(
            List.of("DIA-NN", "ms2pip_Immuno_HCD", "AlphaPept_ms2_generic",
                    "Prosit_2020_intensity_CID", "Prosit_2020_intensity_HCD"));
    public static ArrayList<String> timstofMS2models = new ArrayList<>(
            List.of("DIA-NN", "ms2pip_timsTOF2024", "AlphaPept_ms2_generic", "Prosit_2023_intensity_timsTOF"));

    public static ArrayList<String> generalIMmodels = new ArrayList<>(
            List.of("DIA-NN", "AlphaPept_ccs_generic", "IM2Deep"));

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

    ////////////allowed fragmentation types for models that require it/////////////
    public static HashMap<String, HashSet<String>> allowedFragmentationTypes = new HashMap<>();
    static {
        allowedFragmentationTypes.put("Prosit_2020_intensity_TMT", new HashSet<>(Arrays.asList("HCD", "CID")));
        allowedFragmentationTypes.put("Prosit_2024_intensity_cit", new HashSet<>(Arrays.asList("HCD", "CID")));
        allowedFragmentationTypes.put("Prosit_2025_intensity_MultiFrag",
                new HashSet<>(Arrays.asList("HCD", "ECD", "EID", "UVPD", "ETciD")));
    }
}
