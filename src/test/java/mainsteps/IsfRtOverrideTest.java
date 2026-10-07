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

package mainsteps;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertFalse;
import static org.junit.jupiter.api.Assertions.assertTrue;

import org.junit.jupiter.api.AfterEach;
import org.junit.jupiter.api.BeforeEach;
import org.junit.jupiter.api.Test;
import peptideptmformatting.PeptideFormatter;
import predictions.PredictionEntryHashMap;
import readers.datareaders.IsfAnnotation;
import testutils.IsfFixtures;

/**
 * An in-source fragment (ISF) ion is made in the source from a co-eluting parent ion, so it elutes
 * with the parent, not where its own sequence would. Its expected RT must be the parent's
 * prediction; everything else about it (spectra, IM) stays its own.
 */
class IsfRtOverrideTest {
    private static final float ISF_OWN_RT = 11f;
    private static final float PARENT_RT = 42f;

    private Boolean originalUseSpectra;
    private PredictionEntryHashMap allPreds;
    private MzmlScanNumber scan;
    private IsfRtOverride override;

    @BeforeEach
    void setUp() throws Exception {
        originalUseSpectra = IsfFixtures.disableSpectra();

        scan = new MzmlScanNumber(7, new float[]{150f, 250f}, new float[]{10f, 20f}, 30f, 0f);
        scan.ms2LowerLimit = 100d;
        scan.ms2UpperLimit = 2000d;

        allPreds = new PredictionEntryHashMap();
        allPreds.put("EPTIDEK|2", IsfFixtures.prediction(ISF_OWN_RT));
        //the parent's own pin row is K.n[42.0106]PEPTM[15.9949]IDEK3.A
        allPreds.put(new PeptideFormatter("K.n[42.0106]PEPTM[15.9949]IDEK3.A", "3", "pin").getBaseCharge(),
                IsfFixtures.prediction(PARENT_RT));
        allPreds.put("PEPTIDEKR|2", IsfFixtures.prediction(0f));

        override = new IsfRtOverride();
    }

    @AfterEach
    void tearDown() {
        IsfFixtures.restoreSpectra(originalUseSpectra);
    }

    private PeptideObj psm(IsfAnnotation isf) throws Exception {
        return scan.setPeptideObject(new PeptideFormatter("EPTIDEK", 2, "base"), 1, 1, "0.01",
                allPreds, true, isf, override);
    }

    @Test
    void isfPsmTakesTheParentsPredictedRt() throws Exception {
        PeptideObj pep = psm(new IsfAnnotation("n[42.0106]PEPTM[15.9949]IDEK", 3));

        assertEquals(PARENT_RT, pep.RT);
        assertTrue(pep.isISF);
        assertEquals("EPTIDEK|2", pep.name, "everything but RT stays the ISF's own");
        assertEquals(new IsfAnnotation("n[42.0106]PEPTM[15.9949]IDEK", 3).parentBaseCharge, pep.rtPeptide,
                "RT calibration follows the peptide the RT was predicted for");
        assertEquals(1, override.matched());
        assertEquals(1, override.overridden());
        assertEquals(0, override.fallbacks());
    }

    @Test
    void nonIsfPsmKeepsItsOwnRt() throws Exception {
        PeptideObj pep = psm(null);

        assertEquals(ISF_OWN_RT, pep.RT);
        assertFalse(pep.isISF);
        assertEquals("EPTIDEK|2", pep.rtPeptide);
        assertEquals(0, override.matched());
    }

    @Test
    void missingParentPredictionFallsBackToOwnRt() throws Exception {
        PeptideObj pep = psm(new IsfAnnotation("PEPTIDEKAAA", 2));

        assertEquals(ISF_OWN_RT, pep.RT);
        assertTrue(pep.isISF, "still an ISF, so still kept out of the RT fits");
        assertEquals("EPTIDEK|2", pep.rtPeptide, "its own RT, so its own calibration curve");
        assertEquals(1, override.matched());
        assertEquals(0, override.overridden());
        assertEquals(1, override.fallbacks());
    }

    @Test
    void zeroParentRtFallsBackToOwnRt() throws Exception {
        PeptideObj pep = psm(new IsfAnnotation("PEPTIDEKR", 2));

        assertEquals(ISF_OWN_RT, pep.RT);
        assertEquals(1, override.fallbacks());
    }

    @Test
    void oldSignatureIsUnchanged() throws Exception {
        PeptideObj pep = scan.setPeptideObject(new PeptideFormatter("EPTIDEK", 2, "base"), 1, 1, "0.01",
                allPreds, true);

        assertEquals(ISF_OWN_RT, pep.RT);
        assertFalse(pep.isISF);
    }
}
