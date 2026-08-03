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

package features.spectra;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertTrue;

import allconstants.Constants;
import allconstants.FragmentIonConstants;
import mainsteps.MzmlScanNumber;
import mainsteps.PeptideObj;
import org.junit.jupiter.api.BeforeEach;
import org.junit.jupiter.api.Test;

class SpectrumComparisonTest {

    private static final float[] EXP_MZS = {100f, 200f, 300f, 400f};
    private static final float[] PRED_INTENSITIES = {1f, 0.8f, 0.6f, 0.4f};

    @BeforeEach
    void setUp() {
        //PeptideObj's static baseMap iterates this; populate before PeptideObj loads
        FragmentIonConstants.fragmentIonHierarchy = new String[]{"y", "b"};
        FragmentIonConstants.primaryFragmentIonTypes.clear(); //-> matchWithDaltonsDefault decides
        Constants.matchWithDaltonsDefault = false;            //match on ppm
        Constants.ppmTolerance = 20f / 1000000f;
        Constants.minFragments = 2;
    }

    private SpectrumComparison compare(float[] expIntensities, float[] predMZs) {
        MzmlScanNumber scan = new MzmlScanNumber(1, EXP_MZS.clone(), expIntensities.clone(), 1.0f, 0f);
        PeptideObj pepObj = new PeptideObj();
        pepObj.scanNumObj = scan;
        return new SpectrumComparison(pepObj, scan.getExpMZs(), scan.getExpIntensities(),
                predMZs, PRED_INTENSITIES.clone(), new String[]{"y", "y", "y", "y"}, false);
    }

    // top6matchedIntensity() normalizes by the summed log intensity of every peak STRICTLY above
    // the spectrum mean. When all peaks share one intensity (a SCIEX TOF quantization artifact)
    // no peak clears the mean, so that denominator is 0 and the division escaped as Infinity/NaN
    // into the pin, breaking Percolator and crashing the xchart score histograms.
    @Test
    void flatSpectrumWithMatchesIsNotInfinity() {
        //every predicted fragment matches -> score > 0 -> score/0 was Infinity
        SpectrumComparison sc = compare(new float[]{48f, 48f, 48f, 48f}, new float[]{100f, 200f, 300f, 400f});

        double score = sc.top6matchedIntensity();

        assertTrue(Double.isFinite(score), "flat spectrum must not yield a non-finite score, was " + score);
        assertEquals(0.0, score, "an unnormalizable spectrum should fall back to 0, like the minFragments guard");
    }

    @Test
    void flatSpectrumWithoutMatchesIsNotNaN() {
        //no predicted fragment matches -> score == 0 -> 0.0/0.0 was NaN
        SpectrumComparison sc = compare(new float[]{48f, 48f, 48f, 48f}, new float[]{111f, 222f, 333f, 444f});

        double score = sc.top6matchedIntensity();

        assertTrue(Double.isFinite(score), "flat spectrum must not yield a non-finite score, was " + score);
        assertEquals(0.0, score, "an unnormalizable spectrum should fall back to 0, like the minFragments guard");
    }

    @Test
    void emptySpectrumIsNotNaN() {
        //mean() of an empty array is 0/0 = NaN, so no peak clears the threshold either
        MzmlScanNumber scan = new MzmlScanNumber(1, new float[0], new float[0], 1.0f, 0f);
        PeptideObj pepObj = new PeptideObj();
        pepObj.scanNumObj = scan;
        SpectrumComparison sc = new SpectrumComparison(pepObj, new float[0], new float[0],
                EXP_MZS.clone(), PRED_INTENSITIES.clone(), new String[]{"y", "y", "y", "y"}, false);

        double score = sc.top6matchedIntensity();

        assertTrue(Double.isFinite(score), "empty spectrum must not yield a non-finite score, was " + score);
        assertEquals(0.0, score);
    }

    // guard against the fix flattening ordinary spectra: one peak clears the mean here, so the
    // feature must still normalize as before rather than short-circuit to 0
    @Test
    void ordinarySpectrumStillNormalizes() {
        SpectrumComparison sc = compare(new float[]{10f, 20f, 30f, 100f}, new float[]{100f, 200f, 300f, 400f});

        double score = sc.top6matchedIntensity();

        //denominator = log(101); score = log(11)+log(21)+log(31)+log(101)
        double expected = (Math.log(11) + Math.log(21) + Math.log(31) + Math.log(101))
                / (float) Math.log(101);
        assertEquals(expected, score, 1e-6, "ordinary spectra must keep their existing normalization");
    }
}
