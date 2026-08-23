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

package modelcallers;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertTrue;

import java.util.Arrays;
import java.util.Collections;
import java.util.List;

import org.junit.jupiter.api.Test;

// The number a native process dies by tells the person reading the log nothing. The DIA-NN caller
// turns its own three into sentences; these are FragCast's, and the one that matters is the code a
// library too large to hold in memory dies with after hours of work.
public class FragCastProcessTest {
    private static FragCastProcess.Result result(int exitCode, String... output) {
        return new FragCastProcess.Result(exitCode, Arrays.asList(output));
    }

    @Test
    public void aRunThatWorkedHasNothingToExplain() {
        assertEquals("", result(0, "Wrote 2306447 rows to lib.parquet (parquet)").diagnosis());
    }

    // The real failure: 15.8 million precursors, four hours, then this. The exit code alone
    // (0xC0000409) is Rust aborting, which any un-unwindable panic produces - what names it as memory
    // is the line FragCast printed, so that is what is read.
    @Test
    public void anAllocationFailureIsNamedAsMemoryAndSaysWhatToCut() {
        String d = result(-1073740791,
                "[INFO] Building library: 15863768 peptides, 15863768 precursors",
                "memory allocation of 1610612736 bytes failed").diagnosis();

        assertTrue(d.contains("out of memory"), d);
        assertTrue(d.contains("charge range") && d.contains("fragCastTopN"),
                "a user who just lost four hours needs the two knobs, not just the diagnosis: " + d);
    }

    // Same code without that line is some other abort, and must not be reported as memory.
    @Test
    public void anAbortWithNoAllocationLineIsNotCalledAMemoryFailure() {
        String d = result(-1073740791, "[INFO] Spec backend: conformer").diagnosis();

        assertTrue(d.contains("aborted"), d);
        assertTrue(!d.contains("memory"), "guessed at memory with no evidence for it: " + d);
    }

    @Test
    public void theCodesAWindowsProcessCannotStartWithAreTranslated() {
        assertTrue(result(-1073741515).diagnosis().contains("vc-redist"),
                "a missing runtime DLL should point at the redistributable, as the DIA-NN caller does");
        assertTrue(result(-1073741819).diagnosis().contains("access violation"));
        assertTrue(result(137).diagnosis().contains("memory"));
    }

    // Inventing an explanation for a code nobody has seen would be worse than saying nothing: the
    // caller prints the raw code either way.
    @Test
    public void anUnrecognizedFailureExplainsNothing() {
        assertEquals("", result(42, "something went wrong").diagnosis());
    }

    // --- how a line reaches the console -------------------------------------------------

    // FragCast stamps its own timestamp and level. Echoed whole, every line carried two of each and
    // the second was MSBooster's, so a FragCast warning arrived labelled [INFO].
    @Test
    public void fragCastsOwnPrefixIsTakenOffBeforePrinting() {
        assertEquals("Building library: 15863768 peptides",
                FragCastProcess.message("[2026-08-23 12:19:31 INFO] Building library: 15863768 peptides"));
    }

    @Test
    public void aWarningPrintsAsAProblemAndInformationDoesNot() {
        assertTrue(FragCastProcess.isProblem("[2026-08-23 10:17:37 WARN] Nothing is written until the whole library is built"));
        assertTrue(FragCastProcess.isProblem("[2026-08-23 10:17:37 ERROR] build-library: read ONNX weights"));
        assertTrue(!FragCastProcess.isProblem("[2026-08-23 12:19:31 INFO] Predicting: 10%"));
        assertTrue(!FragCastProcess.isProblem("[2026-08-23 12:19:31 DEBUG] bucket lv=12"));
    }

    // A line FragCast did not stamp - a panic, an allocation failure, anything on plain stderr - is
    // printed exactly as it came. Guessing at a level for it would either hide it or cry wolf.
    @Test
    public void anUnstampedLineIsPassedThroughUntouched() {
        String raw = "memory allocation of 1610612736 bytes failed";
        assertEquals(raw, FragCastProcess.message(raw));
        assertTrue(!FragCastProcess.isProblem(raw));

        assertEquals("", FragCastProcess.message(""));
        assertEquals("[not a log line]", FragCastProcess.message("[not a log line]"));
    }

    @Test
    public void aFailureThatPrintedNothingStillAnswers() {
        List<String> nothing = Collections.emptyList();
        assertEquals("", new FragCastProcess.Result(42, nothing).diagnosis());
        assertTrue(new FragCastProcess.Result(137, nothing).diagnosis().contains("memory"));
    }
}
