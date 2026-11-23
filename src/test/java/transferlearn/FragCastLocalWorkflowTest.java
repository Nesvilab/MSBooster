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

package transferlearn;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertFalse;
import static org.junit.jupiter.api.Assertions.assertNotEquals;
import static org.junit.jupiter.api.Assertions.assertNotNull;
import static org.junit.jupiter.api.Assertions.assertThrows;

import java.io.File;
import java.lang.reflect.Field;
import java.lang.reflect.Modifier;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.HashMap;
import java.util.Map;

import allconstants.Constants;
import org.junit.jupiter.api.AfterEach;
import org.junit.jupiter.api.BeforeEach;
import org.junit.jupiter.api.Test;
import org.junit.jupiter.api.io.TempDir;
import testutils.ExitTrap;

// This entry point initialises the parameters with requirePinMzml=false, which leaves
// Constants.outputDirectory empty - and new File("", name) is not a relative path, it is the ROOT of
// the filesystem. Everything this workflow writes (the split libraries, the fine-tuned ONNX models)
// goes through the directory chosen here, so getting it wrong scatters files at C:\ or fails to
// create them at all.
public class FragCastLocalWorkflowTest {
    /** The directory a bare {@code new File("", name)} would land in: the root of the filesystem. */
    private static File filesystemRoot() {
        return new File("", "fragcast_transfer").getAbsoluteFile();
    }

    @Test
    public void withoutAnOutputDirectoryTheOutputLandsBesideTheLibrary(@TempDir Path dir) {
        File library = dir.resolve("library.tsv").toFile();
        File resolved = FragCastLocalWorkflow.resolveOutputDirectory("", library.getPath());

        assertNotNull(resolved.getParentFile());
        assertEquals(dir.toFile().getAbsoluteFile(), resolved.getParentFile().getAbsoluteFile(),
                "the output was not written beside the library");
        assertNotEquals(filesystemRoot(), resolved.getAbsoluteFile(),
                "the output was written to the root of the filesystem");
    }

    // A bare filename has no parent directory of its own, which is exactly the case that used to
    // resolve to the root.
    @Test
    public void arelativeLibraryPathStillResolvesSomewhereWritable() {
        File resolved = FragCastLocalWorkflow.resolveOutputDirectory("", "library.tsv");
        assertNotEquals(filesystemRoot(), resolved.getAbsoluteFile(),
                "a bare library filename sent the output to the root of the filesystem");
        assertEquals(new File(System.getProperty("user.dir")).getAbsoluteFile(),
                resolved.getParentFile().getAbsoluteFile());
    }

    // --output-dir is the only thing that relocates the output, so it must beat the default rather
    // than be merged with it.
    @Test
    public void anexplicitOutputDirectoryWins(@TempDir Path dir) {
        assertEquals(dir.resolve("typed").toFile(),
                FragCastLocalWorkflow.resolveOutputDirectory(dir.resolve("typed").toString(),
                        dir.resolve("library.tsv").toString()));
    }

    // --- the exit code FragPipe sees --------------------------------------------------------------
    //
    // FragPipe reads a step's outcome from its exit code alone. A FragCast that cannot be started
    // stands in for any failed task here: every task fails, and the exit code is what the options
    // decide. These drive the real parameter machinery, which writes into public statics across the
    // whole of Constants, so every one of them is put back afterwards.

    private final Map<Field, Object> constantsBefore = new HashMap<>();

    @BeforeEach
    public void saveConstants() throws Exception {
        constantsBefore.clear();
        for (Field f : Constants.class.getFields()) {
            if (!Modifier.isFinal(f.getModifiers()) && Modifier.isStatic(f.getModifiers())) {
                constantsBefore.put(f, f.get(null));
            }
        }
    }

    @AfterEach
    public void restoreConstants() throws Exception {
        for (Map.Entry<Field, Object> e : constantsBefore.entrySet()) {
            e.getKey().set(null, e.getValue());
        }
    }

    /** The status the workflow exits with on a library it cannot train on, given these options. */
    private static int exitStatus(Path dir, String cmdOpts) throws Exception {
        File library = dir.resolve("library.tsv").toFile();
        Files.write(library.toPath(), "ModifiedPeptideSequence\tPrecursorCharge\n"
                .getBytes(StandardCharsets.UTF_8));
        File params = dir.resolve("params.txt").toFile();
        Files.write(params.toPath(), String.join("\n",
                "FragCast = " + dir.resolve("no-such-fragcast"),
                "FragCastModelDir = " + dir,
                "numThreads = 1",
                "FragCastFineTuneCmdOpts = " + cmdOpts).getBytes(StandardCharsets.UTF_8));
        try (ExitTrap trap = ExitTrap.install()) {
            return assertThrows(ExitTrap.ExitCalled.class,
                    () -> FragCastLocalWorkflow.fineTuneOnLibrary(params.getPath(), library.getPath(),
                            dir.resolve("out").toString(), "")).status;
        }
    }

    // Without options, the run that trained nothing still ends in success, as it always has: the
    // pretrained models are unchanged, and the run says so.
    @Test
    public void withoutOptionsARunThatTrainedNothingExitsZero(@TempDir Path dir) throws Exception {
        assertEquals(0, exitStatus(dir, ""));
    }

    // With options, the same failure is a failed step: the user asked for training those options
    // describe, and exiting 0 would tell FragPipe they got it.
    @Test
    public void withOptionsAFailedTaskExitsOne(@TempDir Path dir) throws Exception {
        assertEquals(1, exitStatus(dir, "--epochs 5"));
    }

    // Options that cannot be split stop the run before any task starts or any output is made.
    @Test
    public void optionsWithAnUnclosedQuoteStopTheRunBeforeAnyTask(@TempDir Path dir) throws Exception {
        assertEquals(1, exitStatus(dir, "--eval \"C:\\My data\\eval.tsv"));
        assertFalse(dir.resolve("out").toFile().exists(), "the run went on to make its output directory");
    }
}
