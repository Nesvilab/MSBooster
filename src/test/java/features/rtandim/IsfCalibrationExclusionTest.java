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

package features.rtandim;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertFalse;
import static org.junit.jupiter.api.Assertions.assertTrue;

import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.LinkedList;
import java.util.List;
import java.util.Set;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import mainsteps.IsfRtOverride;
import mainsteps.MzmlScanNumber;
import org.junit.jupiter.api.AfterEach;
import org.junit.jupiter.api.BeforeEach;
import org.junit.jupiter.api.Test;
import org.junit.jupiter.api.io.TempDir;
import peptideptmformatting.PeptideFormatter;
import predictions.PredictionEntryHashMap;
import readers.MgfFileReader;
import readers.datareaders.IsfAnnotation;
import readers.datareaders.MzmlReader;
import readers.datareaders.PinReader;
import testutils.IsfFixtures;

/**
 * In-source fragment (ISF) PSMs are scored against their parent's predicted RT, which by
 * construction agrees with where they elute. Letting them anchor an RT calibration or pick the RT
 * model would count the parent twice, so everything that builds or fits something from PSMs must
 * leave them out, while they are still scored.
 */
class IsfCalibrationExclusionTest {
    private static final int REGULAR_PSMS = 12;
    private static final int ISF_PSMS = 3;
    private static final float PARENT_RT = 999f;
    private static final String GOOD_ESCORE = "1.0E-5";
    private static final int USE_ALL = 0;
    private static final String HEADER =
            "SpecId\tLabel\tScanNr\tExpMass\tretentiontime\trank\tlog10_evalue\tPeptide\tProteins";

    @TempDir
    Path dir;

    private Boolean originalUseSpectra;
    private ExecutorService executor;
    private MzmlReader mzml;

    @BeforeEach
    void setUp() throws Exception {
        originalUseSpectra = IsfFixtures.disableSpectra();
        executor = Executors.newSingleThreadExecutor();
        mzml = buildRun();
    }

    @AfterEach
    void tearDown() {
        IsfFixtures.restoreSpectra(originalUseSpectra);
        executor.shutdownNow();
    }

    private static String regularPeptide(int i) {
        return "PEPTIDE" + new String(new char[i]).replace('\0', 'A') + "K";
    }

    private static String isfPeptide(int i) {
        return "ISF" + new String(new char[i]).replace('\0', 'G') + "R";
    }

    /** A run of REGULAR_PSMS ordinary PSMs (scans 1..12) and ISF_PSMS ISF PSMs (scans 13..15). */
    private MzmlReader buildRun() throws Exception {
        //an empty mgf gives a reader whose scans this test supplies itself
        Path mgfPath = dir.resolve("run.mgf");
        Files.write(mgfPath, new byte[0]);
        MgfFileReader mgf = new MgfFileReader(mgfPath.toString(), true, executor, "");
        mgf.scanNumberObjects.clear(); //an empty file still leaves a placeholder scan 0

        PredictionEntryHashMap allPreds = new PredictionEntryHashMap();
        allPreds.put("PARENTPEPTIDEK|2", IsfFixtures.prediction(PARENT_RT));
        IsfAnnotation isf = new IsfAnnotation("PARENTPEPTIDEK", 2);
        IsfRtOverride override = new IsfRtOverride();

        for (int scanNum = 1; scanNum <= REGULAR_PSMS + ISF_PSMS; scanNum++) {
            boolean isIsf = scanNum > REGULAR_PSMS;
            String peptide = isIsf ? isfPeptide(scanNum) : regularPeptide(scanNum);
            allPreds.put(peptide + "|2", IsfFixtures.prediction(scanNum * 10f));

            MzmlScanNumber msn = new MzmlScanNumber(scanNum, new float[]{150f, 250f}, new float[]{10f, 20f},
                    scanNum, 0f);
            msn.ms2LowerLimit = 100d;
            msn.ms2UpperLimit = 2000d;
            msn.setPeptideObject(new PeptideFormatter(peptide, 2, "base"), 1, 1, GOOD_ESCORE, allPreds, true,
                    isIsf ? isf : null, override);
            mgf.scanNumberObjects.put(scanNum, msn);
        }
        assertEquals(ISF_PSMS, override.overridden(), "fixture: every ISF PSM took the parent's RT");
        return new MzmlReader(mgf);
    }

    @Test
    @SuppressWarnings("unchecked")
    void rtCalibrationAnchorsLeaveIsfPsmsOut() throws Exception {
        Object[] arraysAndPeptides = LoessUtilities.getArrays(mzml, USE_ALL, "RT", 0);
        double[][] expAndPred = ((HashMap<String, double[][]>) arraysAndPeptides[0]).get("");
        ArrayList<String> peptides = ((HashMap<String, ArrayList<String>>) arraysAndPeptides[1]).get("");

        assertEquals(REGULAR_PSMS, expAndPred[0].length);
        assertEquals(REGULAR_PSMS, peptides.size());
        for (String peptide : peptides) {
            assertFalse(peptide.startsWith("ISF"), peptide + " must not anchor the RT calibration");
        }
        for (double pred : expAndPred[1]) {
            assertTrue(pred != PARENT_RT, "no anchor may carry an ISF's borrowed RT");
        }
    }

    @Test
    void imCalibrationStillUsesIsfPsms() throws Exception {
        //the ISF ion is made before the mobility cell, so its own IM is genuine
        for (int scanNum : mzml.getScanNums()) {
            MzmlScanNumber msn = mzml.getScanNumObject(scanNum);
            msn.IM = 1f + scanNum / 100f;
            msn.peptideObjects.get(0).IM = 1f + scanNum / 100f;
        }

        @SuppressWarnings("unchecked")
        HashMap<String, double[][]> expAndPred = (HashMap<String, double[][]>)
                LoessUtilities.getArrays(mzml, USE_ALL, "IM", 2)[0];

        assertEquals(REGULAR_PSMS + ISF_PSMS, expAndPred.get("")[0].length);
    }

    @Test
    void rtBinsLeaveIsfPsmsOut() throws Exception {
        ArrayList<Float>[] bins = RTFunctions.RTbins(mzml);

        int binned = 0;
        for (ArrayList<Float> bin : bins) {
            assertFalse(bin.contains(PARENT_RT), "RT bins (and the kernel densities built from them) "
                    + "must not hold an ISF's borrowed RT");
            binned += bin.isEmpty() ? 0 : 1;
        }
        assertEquals(REGULAR_PSMS, binned);
    }

    /** The top PSMs of a pin of these rows, written to the test folder and read back. */
    private LinkedList[] topPsms(String file, int num, Set<String> isfKeys, String... rows) throws Exception {
        Path pin = dir.resolve(file);
        List<String> lines = new ArrayList<>(List.of(HEADER));
        lines.addAll(List.of(rows));
        Files.write(pin, lines, StandardCharsets.UTF_8);
        PinReader reader = new PinReader(pin.toString());
        try {
            return reader.getTopPSMs(num, false, isfKeys);
        } finally {
            reader.close();
        }
    }

    @Test
    void topPsmsForModelSelectionSkipIsfPsms() throws Exception {
        LinkedList[] top = topPsms("run.pin", 2, Set.of(IsfAnnotation.key(1, 1), IsfAnnotation.key(3, 2)),
                "run.1.1.2_1\t1\t1\t800\t1.0\t1\t-10\tK.ISFPEPTIDE2.A\tsp|P1|P1",
                "run.2.2.2_1\t1\t2\t800\t2.0\t1\t-8\tK.PEPTIDEA2.A\tsp|P1|P1",
                "run.3.3.2_1\t1\t3\t800\t3.0\t1\t-6\tK.PEPTIDEAA2.A\tsp|P1|P1",
                "run.3.3.2_2\t1\t3\t800\t3.0\t2\t-9\tK.ISFPEPTIDEB2.A\tsp|P1|P1",
                "run.4.4.2_1\t1\t4\t800\t4.0\t1\t-2\tK.PEPTIDEAAA2.A\tsp|P1|P1");

        List<String> peptides = new ArrayList<>();
        for (Object pf : top[0]) {
            peptides.add(((PeptideFormatter) pf).getBase());
        }
        assertEquals(List.of("PEPTIDEA", "PEPTIDEAA"), peptides);
        assertEquals(List.of(2, 3), new ArrayList<>(top[1]));
    }

    @Test
    void topPsmsOfAPinWithOnlyIsfPsmsIsEmpty() throws Exception {
        LinkedList[] top = topPsms("only.pin", 5, Set.of(IsfAnnotation.key(1, 1)),
                "run.1.1.2_1\t1\t1\t800\t1.0\t1\t-10\tK.ISFPEPTIDE2.A\tsp|P1|P1");

        assertTrue(top[0].isEmpty());
        assertTrue(top[1].isEmpty());
    }
}
