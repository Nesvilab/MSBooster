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

package testutils;

import allconstants.Constants;
import allconstants.FragmentIonConstants;
import predictions.PredictionEntry;

/** What the in-source fragment tests share: a PSM built in a unit test, without a spectra model. */
public final class IsfFixtures {
    private IsfFixtures() {}

    /** A two-peak prediction with this RT. */
    public static PredictionEntry prediction(float rt) {
        PredictionEntry entry = new PredictionEntry(new float[]{150f, 250f}, new float[]{1f, 0.5f},
                new int[]{1, 2}, new int[]{1, 1}, new String[]{"y", "y"});
        entry.setRT(rt);
        return entry;
    }

    /**
     * Lets a PeptideObj be built without a SpectrumComparison, and hands back the useSpectra it
     * replaced for {@link #restoreSpectra}.
     */
    public static Boolean disableSpectra() {
        // PeptideObj's static baseMap iterates this; populate before PeptideObj loads.
        FragmentIonConstants.fragmentIonHierarchy = new String[]{"y", "b"};
        Constants.useIM = false;
        Constants.spectraModel = "";
        Constants.rtModel = "";
        Constants.imModel = "";
        Boolean originalUseSpectra = Constants.useSpectra;
        Constants.useSpectra = false;
        return originalUseSpectra;
    }

    public static void restoreSpectra(Boolean originalUseSpectra) {
        Constants.useSpectra = originalUseSpectra;
    }
}
