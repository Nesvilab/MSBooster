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

package utils;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertFalse;
import static org.junit.jupiter.api.Assertions.assertTrue;

import allconstants.Constants;
import java.io.ByteArrayOutputStream;
import java.io.PrintStream;
import java.lang.reflect.Proxy;
import java.nio.charset.StandardCharsets;
import java.util.Collections;
import org.junit.jupiter.api.AfterEach;
import org.junit.jupiter.api.BeforeEach;
import org.junit.jupiter.api.Test;
import umich.ms.datatypes.lcmsrun.LCMSRunInfo;
import umich.ms.datatypes.scan.props.Instrument;
import umich.ms.datatypes.scancollection.impl.ScanCollectionDefault;
import umich.ms.fileio.filetypes.LCMSDataSource;

// getInstrument reads the instrument model and analyzer from the run header of the first mzML. The header
// MSFragger writes for an Orbitrap Astral Zoom .raw reads back as model "Model N/A" with an empty analyzer,
// and the log said only "Could not detect instrument type", which hid what was read.
class InstrumentUtilsTest {

    private String originalInstrument;
    private String printed;

    @BeforeEach
    void forgetTheDetectedInstrument() {
        originalInstrument = Constants.instrument;
        Constants.instrument = "";
    }

    @AfterEach
    void restoreTheDetectedInstrument() {
        Constants.instrument = originalInstrument;
    }

    // The analyzer is read for every keyword before the Astral one, so a null analyzer used to throw
    // before the model got there.
    @Test
    void anAstralModelIsAThermoTofEvenWithoutAnAnalyzer() {
        Instrument instrument = new Instrument("Thermo", "Orbitrap Astral Zoom", null);

        assertEquals("THERMOTOF", getInstrument(scansWith(instrument)));
        assertEquals("THERMOTOF", Constants.instrument);
        assertTrue(printed.contains("Instrument detected: Astral"), printed);
    }

    // What batmass-io 1.38.x reads from MSFragger's mzML of an Orbitrap Astral Zoom run.
    @Test
    void theLogNamesTheModelAndAnalyzerThatMatchedNothing() {
        Instrument instrument = new Instrument("Vendor N/A", "Model N/A", "");

        assertEquals("QE", getInstrument(scansWith(instrument)));
        assertEquals("QE", Constants.instrument);
        assertTrue(printed.contains("Could not detect instrument type from model [Model N/A], analyzer []. "
                + "Setting to QE. If a different instrument was used, specify using '--instrument' via the "
                + "command line or 'instrument=' in the param file."), printed);
    }

    @Test
    void aNullAnalyzerFallsBackToQeAndTheLogNamesTheModel() {
        Instrument instrument = new Instrument();
        instrument.setAnalyzer(null);

        assertEquals("QE", getInstrument(scansWith(instrument)));
        assertEquals("QE", Constants.instrument);
        assertTrue(printed.contains("from model [Model N/A], analyzer []. Setting to QE."), printed);
    }

    // A run header without a componentList leaves the analyzer at its placeholder.
    @Test
    void theLogNamesTheAnalyzerPlaceholder() {
        Instrument instrument = new Instrument("Vendor N/A", "Model N/A", "Analyzer N/A");

        assertEquals("QE", getInstrument(scansWith(instrument)));
        assertTrue(printed.contains("from model [Model N/A], analyzer [Analyzer N/A]. Setting to QE."), printed);
    }

    // A keyword in the analyzer alone still decides, as it did before.
    @Test
    void anAnalyzerKeywordIsEnough() {
        Instrument instrument = new Instrument("Bruker", "Model N/A", "timsTOF");

        assertEquals("TIMSTOF", getInstrument(scansWith(instrument)));
        assertTrue(printed.contains("Instrument detected: TIMSTOF"), printed);
    }

    // A null model used to throw on the first keyword check and fall back to QE.
    @Test
    void aNullModelLeavesTheAnalyzerToDecide() {
        Instrument instrument = new Instrument("Bruker", null, "timsTOF");

        assertEquals("TIMSTOF", getInstrument(scansWith(instrument)));
        assertTrue(printed.contains("Instrument detected: TIMSTOF"), printed);
    }

    @Test
    void aRunWithoutAnInstrumentFallsBackToQe() {
        assertEquals("QE", getInstrument(scansWith(null)));
        assertTrue(printed.contains("from model [], analyzer []. Setting to QE."), printed);
    }

    // The scans of an MGF input (MzmlReader(MgfFileReader)) are null.
    @Test
    void noScansFallBackToQe() {
        assertEquals("QE", getInstrument(null));
        assertTrue(printed.contains("from model [], analyzer []. Setting to QE."), printed);
    }

    @Test
    void anInstrumentAlreadySetIsKept() {
        Constants.instrument = "LUMOS";

        assertEquals("LUMOS", getInstrument(scansWith(new Instrument("Thermo", "Orbitrap Astral", null))));
        assertEquals("LUMOS", Constants.instrument);
        assertFalse(printed.contains("Instrument detected"), printed);
    }

    private String getInstrument(ScanCollectionDefault scans) {
        String[] instrument = new String[1];
        printed = stdout(() -> instrument[0] = InstrumentUtils.getInstrument(scans));
        return instrument[0];
    }

    // scans whose run header holds this instrument, or no instrument when it is null. Parsing the run
    // yields no scans; getInstrument only needs the header.
    private static ScanCollectionDefault scansWith(Instrument instrument) {
        LCMSRunInfo runInfo = new LCMSRunInfo();
        if (instrument != null) {
            runInfo.addInstrument(instrument, "IC1");
        }
        LCMSDataSource<?> source = (LCMSDataSource<?>) Proxy.newProxyInstance(
                LCMSDataSource.class.getClassLoader(), new Class<?>[]{LCMSDataSource.class},
                (proxy, method, args) -> switch (method.getName()) {
                    case "parse" -> Collections.emptyList();
                    case "getRunInfo", "fetchRunInfo", "parseRunInfo" -> runInfo;
                    case "getName", "toString" -> "test.mzML";
                    case "hashCode" -> System.identityHashCode(proxy);
                    case "equals" -> proxy == args[0];
                    default -> throw new UnsupportedOperationException(method.getName());
                });
        ScanCollectionDefault scans = new ScanCollectionDefault(true);
        scans.setDataSource(source);
        return scans;
    }

    private static String stdout(Runnable code) {
        PrintStream out = System.out;
        ByteArrayOutputStream bytes = new ByteArrayOutputStream();
        try {
            System.setOut(new PrintStream(bytes, true, StandardCharsets.UTF_8));
            code.run();
        } finally {
            System.setOut(out);
        }
        return bytes.toString(StandardCharsets.UTF_8);
    }
}
