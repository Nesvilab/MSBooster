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

package koinaclasses;

import static org.junit.jupiter.api.Assertions.assertDoesNotThrow;
import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertFalse;
import static org.junit.jupiter.api.Assertions.assertThrows;
import static org.junit.jupiter.api.Assertions.assertTrue;

import allconstants.Constants;
import allconstants.NceConstants;
import java.io.ByteArrayOutputStream;
import java.io.PrintStream;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.TreeMap;
import java.util.stream.Collectors;
import java.util.stream.Stream;
import org.junit.jupiter.api.AfterEach;
import org.junit.jupiter.api.BeforeEach;
import org.junit.jupiter.api.Test;
import org.junit.jupiter.api.io.TempDir;
import peptideptmformatting.PeptideFormatter;

/**
 * The NCE of a model is calibrated over its own range, or over minNCE and maxNCE where they are set (D27). A range
 * that is empty ends the run with a KoinaException that names the model and where each end comes from, before any
 * model is calibrated. Calibration PSMs the model cannot take leave nothing to compare: the model keeps the default
 * NCE, and whether the run has anything to send it is decided when its requests are written (D28).
 */
class NCEcalibratorTest {
    private static final String MODEL = "Prosit_2020_intensity_HCD"; //20 to 40

    @TempDir
    Path dir;

    private Integer minNCE;
    private Integer maxNCE;
    private Boolean calibrateNCE;
    private Map<String, String> mzmlNCEs;
    private String fragmentationType;
    private String auxSpectraModel;

    @BeforeEach
    void setUp() {
        minNCE = NceConstants.minNCE;
        maxNCE = NceConstants.maxNCE;
        calibrateNCE = NceConstants.calibrateNCE;
        mzmlNCEs = new HashMap<>(NceConstants.mzmlNCEs);
        fragmentationType = Constants.FragmentationType;
        auxSpectraModel = Constants.auxSpectraModel;
        NceConstants.minNCE = null;
        NceConstants.maxNCE = null;
        NceConstants.calibrateNCE = true;
        Constants.auxSpectraModel = "";
    }

    @AfterEach
    void tearDown() {
        NceConstants.minNCE = minNCE;
        NceConstants.maxNCE = maxNCE;
        NceConstants.calibrateNCE = calibrateNCE;
        NceConstants.mzmlNCEs.clear();
        NceConstants.mzmlNCEs.putAll(mzmlNCEs);
        Constants.FragmentationType = fragmentationType;
        Constants.auxSpectraModel = auxSpectraModel;
    }

    @Test
    void anEmptyRangeEndsTheRunSayingWhereEachEndComesFrom() {
        NceConstants.minNCE = 45; //to widen the scan of the timsTOF model, whose own range goes to 70
        assertDoesNotThrow(() -> NceConstants.requireCalibrationRange("Prosit_2023_intensity_timsTOF"));
        KoinaException e = assertThrows(KoinaException.class, () -> NceConstants.requireCalibrationRange(MODEL));
        assertEquals("The NCE range to calibrate " + MODEL + " over is empty: from 45 (minNCE) to 40 (the highest NCE"
                + " of " + MODEL + "). Set minNCE and maxNCE so that minNCE is not above maxNCE, or leave out the one"
                + " that is not needed", e.getMessage());

        NceConstants.minNCE = null;
        NceConstants.maxNCE = 15;
        e = assertThrows(KoinaException.class, () -> NceConstants.requireCalibrationRange(MODEL));
        assertTrue(e.getMessage().contains("from 20 (the lowest NCE of " + MODEL + ") to 15 (maxNCE)"),
                e.getMessage());

        NceConstants.minNCE = 30;
        NceConstants.maxNCE = 25;
        e = assertThrows(KoinaException.class, () -> NceConstants.requireCalibrationRange(MODEL));
        assertTrue(e.getMessage().contains("from 30 (minNCE) to 25 (maxNCE)"), e.getMessage());

        //the calibration itself says so too, as an exception MainClass ends the run with, not an exit
        e = assertThrows(KoinaException.class,
                () -> NCEcalibrator.calibrateNCE(MODEL, new KoinaMethods(null), dir.toString(), false));
        assertTrue(e.getMessage().startsWith("The NCE range to calibrate " + MODEL + " over is empty"),
                e.getMessage());
    }

    @Test
    void aModelWhoseNceIsNotCalibratedPasses() {
        NceConstants.minNCE = 45;
        assertDoesNotThrow(() -> NceConstants.requireCalibrationRange("Prosit_2019_irt"), "an RT model");
        assertDoesNotThrow(() -> NceConstants.requireCalibrationRange(""), "no model");
        assertDoesNotThrow(() -> NceConstants.requireCalibrationRange("DIA-NN"));
        assertDoesNotThrow(() -> NceConstants.requireCalibrationRange("alphapeptdeep"),
                "a model predicted outside Koina, whose NCE calibrateNCE does not calibrate");
        NceConstants.calibrateNCE = false;
        assertDoesNotThrow(() -> NceConstants.requireCalibrationRange(MODEL), "calibrateNCE = false");
    }

    @Test
    void calibrationPsmsTheModelCannotTakeLeaveItsDefaultNce() throws Exception {
        NceConstants.mzmlNCEs.clear();
        NceConstants.mzmlNCEs.put("HCD", "27.5");
        Constants.FragmentationType = "HCD";
        KoinaMethods km = new KoinaMethods(null);
        km.peptideArraylist.add(new PeptideFormatter("GDUVLPAEK", "2", "base")); //residue U: never sent
        km.peptideArraylist.add(new PeptideFormatter("ADUGLR", "3", "base"));

        Object[] result = NCEcalibrator.calibrateNCE(MODEL, km, dir.resolve("nce").toString(), true);
        assertTrue(((TreeMap<?, ?>) result[0]).isEmpty(), "nothing compared");
        assertEquals(0d, result[1]);
        assertEquals(27, result[2], "the default NCE, as the best-model search sends a model it does not calibrate");

        //no chart for a model that was not calibrated, and no exception
        String figures = Constants.figureDirectory;
        Constants.figureDirectory = dir.resolve("figures").toString();
        try {
            NCEcalibrator.plotNCEchart(MODEL, new TreeMap<>());
            TreeMap<Integer, ArrayList<Double>> empty = new TreeMap<>();
            empty.put(30, new ArrayList<>());
            NCEcalibrator.plotNCEchart(MODEL, empty);
            assertFalse(Files.exists(dir.resolve("figures")));
        } finally {
            Constants.figureDirectory = figures;
        }
    }

    //D27: the calibration tries the model's own NCE range, and minNCE or maxNCE only where one is set
    @Test
    void theCalibrationTriesTheModelsOwnRangeUnlessMinOrMaxNceIsSet() throws Exception {
        String timsTof = "Prosit_2023_intensity_timsTOF";
        assertEquals("Calibrating NCE for " + timsTof + " from 20 to 70", calibrationRange(timsTof));
        assertEquals("Calibrating NCE for " + MODEL + " from 20 to 40", calibrationRange(MODEL));
        NceConstants.maxNCE = 50;
        assertEquals("Calibrating NCE for " + timsTof + " from 20 to 50", calibrationRange(timsTof));
        assertEquals("Calibrating NCE for " + MODEL + " from 20 to 50", calibrationRange(MODEL));
        NceConstants.maxNCE = null;
        NceConstants.minNCE = 30;
        assertEquals("Calibrating NCE for " + timsTof + " from 30 to 70", calibrationRange(timsTof));
        assertEquals("Calibrating NCE for " + MODEL + " from 30 to 40", calibrationRange(MODEL));
    }

    //the line calibrateNCE prints for the range it tries, for PSMs it cannot send, so that nothing is predicted
    private String calibrationRange(String model) throws Exception {
        NceConstants.mzmlNCEs.clear();
        NceConstants.mzmlNCEs.put("HCD", "27.5");
        Constants.FragmentationType = "HCD";
        KoinaMethods km = new KoinaMethods(null);
        km.peptideArraylist.add(new PeptideFormatter("GDUVLPAEK", "2", "base"));

        PrintStream out = System.out;
        ByteArrayOutputStream printed = new ByteArrayOutputStream();
        try {
            System.setOut(new PrintStream(printed, true, StandardCharsets.UTF_8));
            NCEcalibrator.calibrateNCE(model, km, dir.resolve("nce-" + model).toString(), false);
        } finally {
            System.setOut(out);
        }
        String range = null;
        for (String line : printed.toString(StandardCharsets.UTF_8).split("\\R")) {
            int at = line.indexOf("Calibrating NCE for ");
            if (at >= 0) {
                range = line.substring(at);
            }
        }
        return range;
    }

    //the charts of a calibration cover the range it tried, six NCEs to a chart
    @Test
    void theChartsCoverTheModelsOwnRange() throws Exception {
        String timsTof = "Prosit_2023_intensity_timsTOF";
        TreeMap<Integer, ArrayList<Double>> similarities = new TreeMap<>();
        for (int nce : new int[] {21, 65, 70}) {
            similarities.put(nce, new ArrayList<>(List.of(0.4, 0.5, 0.6)));
        }
        assertEquals(List.of("NCE_calibration_" + timsTof + "20to25.png", "NCE_calibration_" + timsTof + "62to67.png",
                "NCE_calibration_" + timsTof + "68to70.png"), charts(timsTof, similarities, "own"));
        //a model of 20 to 40 has no chart above 40
        assertEquals(List.of("NCE_calibration_" + MODEL + "20to25.png"), charts(MODEL, similarities, "hcd"));
        NceConstants.maxNCE = 66;
        assertEquals(List.of("NCE_calibration_" + timsTof + "20to25.png", "NCE_calibration_" + timsTof + "62to66.png"),
                charts(timsTof, similarities, "max"));
    }

    private List<String> charts(String model, TreeMap<Integer, ArrayList<Double>> similarities, String folder)
            throws Exception {
        String figures = Constants.figureDirectory;
        String extension = Constants.plotExtension;
        Path figureDir = dir.resolve("figures-" + folder);
        Constants.figureDirectory = figureDir.toString();
        Constants.plotExtension = "png";
        try {
            NCEcalibrator.plotNCEchart(model, similarities);
        } finally {
            Constants.figureDirectory = figures;
            Constants.plotExtension = extension;
        }
        try (Stream<Path> files = Files.list(figureDir)) {
            return files.map(f -> f.getFileName().toString()).sorted().collect(Collectors.toList());
        }
    }
}
