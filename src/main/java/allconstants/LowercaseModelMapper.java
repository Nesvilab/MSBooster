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
        lowercaseToModel.put("ms2pip_timstof2024", "ms2pip_timsTOF2024");
        lowercaseToModel.put("ms2pip_cid_tmt", "ms2pip_CID_TMT");
        lowercaseToModel.put("ms2pip_ttof5600", "ms2pip_TTOF5600");
        lowercaseToModel.put("ms2pip_immuno_hcd", "ms2pip_Immuno_HCD");
        lowercaseToModel.put("ms2pip_itraqphospho", "ms2pip_iTRAQphospho");
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
