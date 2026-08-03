/*
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or
 * implied. See the License for the specific language governing
 * permissions and limitations under the License.
 */

package peptideptmformatting;

import allconstants.Constants;

public class PeptideSkipper {
    //provide peptide and see if it may be problematic
    //TODO: also require string argument fullpeptide in their language, to check if deeplc fullpeptide is longer than 60AA
    public static boolean skipPeptide(PeptideFormatter pf, String model) {
        model = model.toLowerCase();

        String stripped = pf.getStripped();
        String charge = pf.getCharge();
        String unimodName = "";
        if (model.contains("deeplc")) {
            unimodName = pf.getModel(model);
        }

        //letters
        if (model.contains("prosit") || model.contains("ms2pip") || model.contains("deeplc") ||
                model.contains("unispec") || model.contains("predfull") || model.contains("im2deep")) {
            for (char c : "OUBZJX".toCharArray()) {
                if (stripped.indexOf(c) != -1) {
                    return true;
                }
            }
        }
        //length
        if ((model.contains("ms2pip") || model.contains("prosit")) && stripped.length() > 30) {
            return true;
        }
        if (model.contains("unispec") && stripped.length() > 40) {
            return true;
        }
        if ((model.contains("deeplc") || model.contains("im2deep")) && stripped.length() > 60) {
            return true;
        }
        //charge
        int chargeInt = Integer.parseInt(charge);
        if (model.contains("unispec") && chargeInt > 5) {
            return true;
        }
        if (model.contains("prosit") && chargeInt > 6) { //predfull can handle charge up to 30
            return true;
        }
        //string length
        if (model.contains("deeplc")) {
            if (unimodName.length() > 60) {
                return true;
            }
        }
        return false;
    }
}
