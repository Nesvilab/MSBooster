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
import static org.junit.jupiter.api.Assertions.assertTrue;

import allconstants.Constants;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.Map;
import java.util.concurrent.ConcurrentSkipListMap;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import mainsteps.MzmlScanNumber;
import org.junit.jupiter.api.AfterEach;
import org.junit.jupiter.api.BeforeEach;
import org.junit.jupiter.api.Test;
import org.junit.jupiter.api.io.TempDir;
import readers.MgfFileReader;

/**
 * setInverseLoess tabulates the RT calibration curve (minutes to predicted RT) on a grid of minutes,
 * split across Constants.numThreads tasks, and stores it inverted (predicted RT to minutes). Where the
 * curve is flat, many minutes share one predicted RT, and the minute stored for it must be the one a
 * single-threaded loop keeps (the last, i.e. the largest), not the one of whichever task wrote last. The
 * grid itself must also be the one a single-threaded loop steps through, whatever Constants.numThreads is.
 */
class InverseRtCalibrationTest {
    private static final int THREADS = 8;
    private static final float MIN_RT = 0f;
    private static final float MAX_RT = 100f;
    //the curve is flat until this minute, so the first three of the eight tasks write the same key
    private static final double FLAT_UNTIL = 30d;
    private static final int REPEATS = 200;

    @TempDir
    Path dir;

    private Integer originalNumThreads;
    private String originalOutputDirectory;
    private ExecutorService sequential;
    private ExecutorService parallel;

    @BeforeEach
    void setUp() {
        originalNumThreads = Constants.numThreads;
        originalOutputDirectory = Constants.outputDirectory;
        Constants.numThreads = THREADS;
        Constants.outputDirectory = dir.toString();
        sequential = Executors.newSingleThreadExecutor();
        parallel = Executors.newFixedThreadPool(THREADS);
    }

    @AfterEach
    void tearDown() {
        Constants.numThreads = originalNumThreads;
        Constants.outputDirectory = originalOutputDirectory;
        sequential.shutdownNow();
        parallel.shutdownNow();
    }

    /** Scans from MIN_RT to MAX_RT, calibrated with a curve that is flat until FLAT_UNTIL. */
    private MzmlReader flatStartReader() throws Exception {
        Path mgfPath = dir.resolve("run.mgf");
        Files.write(mgfPath, new byte[0]);
        MgfFileReader mgf = new MgfFileReader(mgfPath.toString(), true, sequential, "");
        mgf.scanNumberObjects.clear();
        mgf.scanNumberObjects.put(1, new MzmlScanNumber(1, new float[0], new float[0], MIN_RT, 0f));
        mgf.scanNumberObjects.put(2, new MzmlScanNumber(2, new float[0], new float[0], MAX_RT, 0f));

        MzmlReader mzml = new MzmlReader(mgf);
        mzml.RTLOESS.put("", x -> Math.max(x, FLAT_UNTIL));
        return mzml;
    }

    @Test
    void flatCurveKeepsTheLastMinuteWhicheverTaskWritesLast() throws Exception {
        //one worker runs the tasks in the order they were submitted, as a single-threaded loop would
        MzmlReader sequentialMzml = flatStartReader();
        sequentialMzml.setInverseLoess(sequential);
        ConcurrentSkipListMap<Double, Double> expected = sequentialMzml.irtToMinutes.get("");

        Map.Entry<Double, Double> flat = expected.firstEntry();
        float increment = (MAX_RT - MIN_RT) / 10000f;
        assertEquals(FLAT_UNTIL, flat.getKey(), 0d);
        assertTrue(flat.getValue() <= FLAT_UNTIL && flat.getValue() > FLAT_UNTIL - 2 * increment,
                "flat part should map to its last minute, not " + flat.getValue());

        for (int i = 0; i < REPEATS; i++) {
            MzmlReader mzml = flatStartReader();
            mzml.setInverseLoess(parallel);
            ConcurrentSkipListMap<Double, Double> actual = mzml.irtToMinutes.get("");
            assertEquals(flat.getValue(), actual.get(FLAT_UNTIL), "repeat " + i);
            assertTrue(expected.equals(actual), "repeat " + i + ": the inverse curve differs from the sequential one");
        }
    }

    @Test
    void curveIsTheSameForAnyNumberOfThreads() throws Exception {
        Constants.numThreads = 1;
        MzmlReader singleTaskMzml = flatStartReader();
        singleTaskMzml.setInverseLoess(sequential);
        ConcurrentSkipListMap<Double, Double> expected = singleTaskMzml.irtToMinutes.get("");

        for (int threads : new int[] {2, 3, 7, THREADS, 11}) {
            Constants.numThreads = threads;
            MzmlReader mzml = flatStartReader();
            mzml.setInverseLoess(parallel);
            assertTrue(expected.equals(mzml.irtToMinutes.get("")),
                    threads + " threads: the inverse curve differs from the single-task one");
        }
    }
}
