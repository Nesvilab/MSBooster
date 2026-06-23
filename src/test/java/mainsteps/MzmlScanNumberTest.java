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

package mainsteps;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertNotNull;

import allconstants.Constants;
import allconstants.FragmentIonConstants;
import org.junit.jupiter.api.Test;
import peptideptmformatting.PeptideFormatter;
import predictions.PredictionEntryHashMap;

class MzmlScanNumberTest {

    // A PSM whose precursor has no prediction in the (provided) library must NOT crash
    // MSBooster. Historically setPeptideObject dereferenced the null PredictionEntry and,
    // for a non-skippable peptide, called System.exit(1). It should instead degrade
    // gracefully: tally the PSM as skipped and return a zero-vector PeptideObj.
    @Test
    void missingPredictionIsSkippedNotFatal() throws Exception {
        // PeptideObj's static baseMap iterates this; populate before PeptideObj loads.
        FragmentIonConstants.fragmentIonHierarchy = new String[]{"y", "b"};
        Constants.useIM = false;
        Constants.spectraModel = "";  // skipPeptide() returns false -> the historical crash path
        Constants.rtModel = "";
        Constants.imModel = "";

        MzmlScanNumber scan = new MzmlScanNumber(
                1, new float[]{100f, 200f}, new float[]{10f, 20f}, 1.0f, 0f);

        // unmodified peptide; "base" format -> baseCharge "EEQEPETLAAWWNR|3"
        PeptideFormatter name = new PeptideFormatter("EEQEPETLAAWWNR", 3, "base");

        // empty prediction map -> allPreds.get("EEQEPETLAAWWNR|3") returns null (prediction missing)
        PredictionEntryHashMap allPreds = new PredictionEntryHashMap();

        // this fork is shared across tests; don't leak the useSpectra override
        Boolean originalUseSpectra = Constants.useSpectra;
        PeptideObj obj;
        try {
            Constants.useSpectra = false; // avoid building a SpectrumComparison in this unit test
            obj = scan.setPeptideObject(name, 1, 1, "0.01", allPreds, true);
        } finally {
            Constants.useSpectra = originalUseSpectra;
        }

        assertNotNull(obj, "a missing prediction should yield a zero-vector PeptideObj, not a crash");
        assertEquals(1, scan.skippedPSMs.get(),
                "a missing prediction should be counted as a skipped PSM");
        assertEquals(1, scan.peptideObjects.size(),
                "the placeholder PeptideObj should be registered on the scan");
    }
}
