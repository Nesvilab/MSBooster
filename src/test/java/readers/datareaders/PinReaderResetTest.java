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
import static org.junit.jupiter.api.Assertions.assertFalse;
import static org.junit.jupiter.api.Assertions.assertTrue;

import java.io.BufferedReader;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.Arrays;
import org.junit.jupiter.api.Test;
import org.junit.jupiter.api.io.TempDir;

/**
 * {@link PinReader#reset()} reopens the pin. It must close the reader it replaces: otherwise every
 * reset leaks a file handle, which on Windows also keeps the pin from being deleted or renamed.
 */
class PinReaderResetTest {

    @TempDir
    Path dir;

    private static boolean isClosed(BufferedReader reader) {
        try {
            reader.ready();
            return false;
        } catch (java.io.IOException e) {
            return true; // "Stream closed"
        }
    }

    @Test
    void resetClosesTheReaderItReplaces() throws Exception {
        Path file = dir.resolve("run.pin");
        Files.write(file, Arrays.asList("SpecId\tLabel\tScanNr", "run.1.1.2_1\t1\t1", "run.2.2.2_1\t1\t2"),
                StandardCharsets.UTF_8);
        PinReader pin = new PinReader(file.toString());
        BufferedReader first = pin.in;

        pin.reset();

        assertTrue(isClosed(first), "the replaced reader must be closed");
        assertFalse(isClosed(pin.in));
        assertTrue(pin.next(true), "reading starts over after the header");
        assertEquals(1, pin.getScanNum());
        pin.close();
    }
}
