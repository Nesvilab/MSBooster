package allconstants;

import java.util.HashMap;

public class LowercaseModelMapper {
    public static HashMap<String, String> lowercaseToModel = new HashMap<>();
    static {
        lowercaseToModel.put("", "");
        lowercaseToModel.put("dia-nn", "DIA-NN");

        lowercaseToModel.put("alphapept_rt_generic", "AlphaPept_rt_generic");
        lowercaseToModel.put("prosit_2019_irt", "Prosit_2019_irt");
        lowercaseToModel.put("prosit_2020_irt_tmt", "Prosit_2020_irt_TMT");
        lowercaseToModel.put("deeplc_hela_hf", "Deeplc_hela_hf");

        lowercaseToModel.put("ms2pip_2021_hcd", "ms2pip_2021_HCD");
        lowercaseToModel.put("alphapept_ms2_generic", "AlphaPept_ms2_generic");
        lowercaseToModel.put("prosit_2019_intensity", "Prosit_2019_intensity");
        lowercaseToModel.put("prosit_2020_intensity_cid", "Prosit_2020_intensity_CID");
        lowercaseToModel.put("prosit_2020_intensity_tmt", "Prosit_2020_intensity_TMT");
        lowercaseToModel.put("prosit_2020_intensity_hcd", "Prosit_2020_intensity_HCD");
        lowercaseToModel.put("prosit_2023_intensity_timstof", "Prosit_2023_intensity_timsTOF");
        lowercaseToModel.put("unispec", "UniSpec");

        lowercaseToModel.put("alphapept_ccs_generic", "AlphaPept_ccs_generic");
    }
    public LowercaseModelMapper() {}
}
