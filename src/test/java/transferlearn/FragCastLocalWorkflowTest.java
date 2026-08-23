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
import static org.junit.jupiter.api.Assertions.assertNotEquals;
import static org.junit.jupiter.api.Assertions.assertNotNull;

import java.io.File;
import java.nio.file.Path;

import org.junit.jupiter.api.Test;
import org.junit.jupiter.api.io.TempDir;

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
}
