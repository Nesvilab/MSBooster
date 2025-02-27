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

package peptideptmformatting;

import allconstants.Constants;

public class PeptideSkipper {
    //provide peptide and see if it may be problematic
    public static boolean skipPeptide(String stripped, String charge, String model) {
        model = model.toLowerCase();
        //letters
        if (model.contains("prosit") || model.contains("ms2pip") || model.contains("deeplc") ||
                model.contains("unispec") || model.contains("predfull")) {
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
        if (model.contains("deeplc") && stripped.length() > 60) {
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

        //library tsv is more lenient, let anything be missing
        if ((Constants.spectraPredFile != null && Constants.spectraPredFile.endsWith(".tsv")) ||
            (Constants.RTPredFile != null && Constants.RTPredFile.endsWith(".tsv")) ||
            (Constants.IMPredFile != null && Constants.IMPredFile.endsWith(".tsv")) ||
            (Constants.auxSpectraPredFile != null && Constants.auxSpectraPredFile.endsWith(".tsv"))) {
            return true;
        }
        return false;
    }
}
