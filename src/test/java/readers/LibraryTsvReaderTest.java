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