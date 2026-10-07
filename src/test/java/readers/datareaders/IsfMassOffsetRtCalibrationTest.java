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

package readers.datareaders;

import static org.junit.jupiter.api.Assertions.assertEquals;

import allconstants.Constants;
import features.rtandim.MassOffsetGroup;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.concurrent.ConcurrentSkipListMap;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import mainsteps.IsfRtOverride;
import mainsteps.MzmlScanNumber;
import mainsteps.PeptideObj;
import org.junit.jupiter.api.AfterEach;
import org.junit.jupiter.api.BeforeEach;
import org.junit.jupiter.api.Test;
import org.junit.jupiter.api.io.TempDir;
import peptideptmformatting.PeptideFormatter;
import predictions.PredictionEntryHashMap;
import readers.MgfFileReader;
import testutils.IsfFixtures;

/**
 * In a mass-offset search MSBooster fits one RT calibration per offset group, and a PSM is calibrated
 * with the curve of its group. An in-source fragment (ISF) that lost a labile offset (a glycan, say)
 * carries the parent's predicted RT, which lives on the parent's group curve, so it must be
 * calibrated with that curve and not with the one of its own (offset-free) name.
 */
class IsfMassOffsetRtCalibrationTest {
    private static final String GLYCAN = "203.0794";
    private static final String PARENT_PIN_PEPTIDE = "PEPTN[" + GLYCAN + "]IDEK";
    private static final String ISF_PEPTIDE = "PEPTNIDEK";
    private static final int CHARGE = 2;

    private static final float SCAN_RT = 10f;
    private static final float ISF_OWN_PRED_RT = 10f;
    private static final float PARENT_PRED_RT = 50f;
    //at the scan's RT the glycan curve expects the parent's prediction, the others curve the ISF's own
    private static final double GLYCAN_CURVE_AT_SCAN = PARENT_PRED_RT;
    private static final double OTHERS_CURVE_AT_SCAN = ISF_OWN_PRED_RT;
    private static final float BIN_IQR = 2f;
    private static final double EPS = 1e-9;

    @TempDir
    Path dir;

    private Boolean originalUseSpectra;
    private Integer originalRtBinMultiplier;
    private ExecutorService executor;

    @BeforeEach
    void setUp() {
        originalUseSpectra = IsfFixtures.disableSpectra();
        originalRtBinMultiplier = Constants.RTbinMultiplier;
        Constants.RTbinMultiplier = 1;
        executor = Executors.newSingleThreadExecutor();
    }

    @AfterEach
    void tearDown() {
        IsfFixtures.restoreSpectra(originalUseSpectra);
        Constants.RTbinMultiplier = originalRtBinMultiplier;
        executor.shutdownNow();
    }

    private static String parentBaseCharge() {
        return new IsfAnnotation(PARENT_PIN_PEPTIDE, CHARGE).parentBaseCharge;
    }

    /** One scan at SCAN_RT holding one ISF PSM, calibrated with a glycan curve and an "others" curve. */
    private MzmlReader runWithOneIsf(IsfAnnotation isf, boolean parentPredicted) throws Exception {
        Path mgfPath = dir.resolve("run.mgf");
        Files.write(mgfPath, new byte[0]);
        MgfFileReader mgf = new MgfFileReader(mgfPath.toString(), true, executor, "");
        mgf.scanNumberObjects.clear();

        PredictionEntryHashMap allPreds = new PredictionEntryHashMap();
        allPreds.put(ISF_PEPTIDE + "|" + CHARGE, IsfFixtures.prediction(ISF_OWN_PRED_RT));
        if (parentPredicted) {
            allPreds.put(parentBaseCharge(), IsfFixtures.prediction(PARENT_PRED_RT));
        }

        MzmlScanNumber msn = new MzmlScanNumber(1, new float[]{150f, 250f}, new float[]{10f, 20f}, SCAN_RT, 0f);
        msn.ms2LowerLimit = 100d;
        msn.ms2UpperLimit = 2000d;
        msn.setPeptideObject(new PeptideFormatter(ISF_PEPTIDE, CHARGE, "base"), 1, 1, "0.01", allPreds, true,
                isf, new IsfRtOverride());
        mgf.scanNumberObjects.put(1, msn);

        MzmlReader mzml = new MzmlReader(mgf);
        mzml.RTLOESS.put(GLYCAN, x -> GLYCAN_CURVE_AT_SCAN);
        mzml.RTLOESS.put("others", x -> OTHERS_CURVE_AT_SCAN);
        //predicted RT -> minutes: the glycan group puts the parent's prediction at the scan's RT
        mzml.irtToMinutes.put(GLYCAN, line(PARENT_PRED_RT, SCAN_RT));
        mzml.irtToMinutes.put("others", line(ISF_OWN_PRED_RT, SCAN_RT));
        mzml.RTbinStats = new float[(int) SCAN_RT + 1][3];
        mzml.RTbinStats[(int) SCAN_RT][2] = BIN_IQR;
        return mzml;
    }

    /** A predicted-RT-to-minutes map through the origin and (pred, minutes). */
    private static ConcurrentSkipListMap<Double, Double> line(double pred, double minutes) {
        ConcurrentSkipListMap<Double, Double> map = new ConcurrentSkipListMap<>();
        map.put(0d, 0d);
        map.put(pred, minutes);
        map.put(pred * 10, minutes * 10);
        return map;
    }

    private static PeptideObj onlyPsm(MzmlReader mzml) throws Exception {
        return mzml.getScanNumObject(1).peptideObjects.get(0);
    }

    @Test
    void fixtureParentKeyCarriesTheGlycanOffset() {
        assertEquals(GLYCAN, String.valueOf(
                MassOffsetGroup.deltaMasses(parentBaseCharge())[0]));
    }

    @Test
    void isfIsCalibratedWithItsParentsOffsetGroup() throws Exception {
        MzmlReader mzml = runWithOneIsf(new IsfAnnotation(PARENT_PIN_PEPTIDE, CHARGE), true);

        mzml.predictRTLOESS(executor);
        PeptideObj pep = onlyPsm(mzml);

        assertEquals(PARENT_PRED_RT, pep.RT, EPS);
        assertEquals(GLYCAN_CURVE_AT_SCAN, pep.calibratedRT, EPS);
        assertEquals(0d, pep.deltaRTLOESS, EPS);
        assertEquals(SCAN_RT, pep.predRTrealUnits, EPS);
        assertEquals(0d, pep.deltaRTLOESS_real, EPS);
    }

    @Test
    void isfNormalizedDeltaUsesItsParentsOffsetGroup() throws Exception {
        MzmlReader mzml = runWithOneIsf(new IsfAnnotation(PARENT_PIN_PEPTIDE, CHARGE), true);

        mzml.calculateDeltaRTLOESSnormalized(executor);

        assertEquals(0d, onlyPsm(mzml).deltaRTLOESSnormalized, EPS);
    }

    @Test
    void isfThatFellBackToItsOwnRtKeepsItsOwnGroup() throws Exception {
        MzmlReader mzml = runWithOneIsf(new IsfAnnotation(PARENT_PIN_PEPTIDE, CHARGE), false);

        mzml.predictRTLOESS(executor);
        PeptideObj pep = onlyPsm(mzml);

        assertEquals(ISF_OWN_PRED_RT, pep.RT, EPS);
        assertEquals(OTHERS_CURVE_AT_SCAN, pep.calibratedRT, EPS);
        assertEquals(0d, pep.deltaRTLOESS, EPS);
    }

    @Test
    void nonIsfKeepsItsOwnGroup() throws Exception {
        MzmlReader mzml = runWithOneIsf(null, true);

        mzml.predictRTLOESS(executor);
        PeptideObj pep = onlyPsm(mzml);

        assertEquals(ISF_OWN_PRED_RT, pep.RT, EPS);
        assertEquals(OTHERS_CURVE_AT_SCAN, pep.calibratedRT, EPS);
        assertEquals(0d, pep.deltaRTLOESS, EPS);
    }
}
