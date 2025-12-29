/*
 * This file is part of MSBooster.
 *
 * MSBooster is free software: you can redistribute it and/or modify it under
 * the terms of the GNU Lesser General Public License as published by the Free
 * Software Foundation, either version 3 of the License, or (at your
 * option) any later version.
 *
 * MSBooster is distributed in the hope that it will be useful, but WITHOUT
 * ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
 * FITNESS FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public License
 * for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with MSBooster. If not, see <https://www.gnu.org/licenses/>.
 */

package readers;

import static org.junit.jupiter.api.Assertions.assertNotNull;

import allconstants.Constants;
import java.nio.file.Paths;
import java.util.Objects;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.ScheduledThreadPoolExecutor;

import org.junit.jupiter.api.Test;
import readers.predictionreaders.LibraryTsvReader;


class LibraryTsvReaderTest {

      @Test
      void testLibraryReader() throws Exception {
          Constants.unimodObo = Paths.get(Objects.requireNonNull(LibraryTsvReaderTest.class.getResource("/unimod.obo")).toURI()).toString();
          Runtime run = Runtime.getRuntime();
          ExecutorService executorService = new ScheduledThreadPoolExecutor(
                  run.availableProcessors() - 1);
          LibraryTsvReader libraryTsvReader = new LibraryTsvReader(Paths.get(Objects.requireNonNull(
                  LibraryTsvReader.class.getResource("/library_1.tsv")).toURI()).toString(),
                  executorService, "unimod.obo");
          assertNotNull(libraryTsvReader);
      }
}