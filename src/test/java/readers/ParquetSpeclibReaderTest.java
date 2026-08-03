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

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertFalse;
import static org.junit.jupiter.api.Assertions.assertTrue;

import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashSet;
import java.util.List;
import java.util.Objects;
import java.util.concurrent.ScheduledThreadPoolExecutor;

import allconstants.Constants;
import org.junit.jupiter.api.Test;
import predictions.PredictionEntry;
import predictions.PredictionEntryHashMap;
import readers.predictionreaders.LibraryTsvReader;
import readers.predictionreaders.ParquetSpeclibReader;

// Proves ParquetSpeclibReader produces the same predictions as LibraryTsvReader for the same
// FragCast build-library run (fragcast_lib.tsv and fragcast_lib.parquet hold identical content).
class ParquetSpeclibReaderTest {

    private static String res(String name) throws Exception {
        return Paths.get(Objects.requireNonNull(
                ParquetSpeclibReaderTest.class.getResource("/" + name)).toURI()).toString();
    }

    private static List<String> fragTuples(PredictionEntry pe) {
        float[] mz = pe.getMzs();
        float[] in = pe.getIntensities();
        int[] num = pe.getFragNums();
        int[] ch = pe.getCharges();
        String[] type = pe.getFragmentIonTypes();
        List<String> out = new ArrayList<>();
        for (int i = 0; i < mz.length; i++) {
            out.add(String.format("%s|%d|%d|%.4f|%.4f", type[i], num[i], ch[i], mz[i], in[i]));
        }
        Collections.sort(out);
        return out;
    }

    @Test
    void parquetMatchesTsv() throws Exception {
        Constants.unimodObo = res("unimod.obo");
        ScheduledThreadPoolExecutor es = new ScheduledThreadPoolExecutor(
                Math.max(1, Runtime.getRuntime().availableProcessors() - 1));

        PredictionEntryHashMap tsv =
                new LibraryTsvReader(res("fragcast_lib.tsv"), es, "unimod.obo").getPreds();
        PredictionEntryHashMap pq =
                new ParquetSpeclibReader(res("fragcast_lib.parquet"), es, new HashSet<>()).getPreds();

        assertFalse(tsv.isEmpty(), "TSV reader produced no predictions");
        assertEquals(tsv.keySet(), pq.keySet(), "precursor key sets differ");

        for (String key : tsv.keySet()) {
            PredictionEntry a = tsv.get(key);
            PredictionEntry b = pq.get(key);
            assertEquals(a.RT, b.RT, 1e-4f, "RT differs for " + key);
            assertEquals(a.IM, b.IM, 1e-4f, "IM differs for " + key);
            assertEquals(fragTuples(a), fragTuples(b), "fragments differ for " + key);
        }
        System.out.println("Parquet == TSV for " + tsv.keySet().size() + " precursors");

        // allowedPrecursors filtering: restrict to a single precursor
        String one = tsv.keySet().iterator().next();
        HashSet<String> allowed = new HashSet<>();
        allowed.add(one);
        PredictionEntryHashMap filtered =
                new ParquetSpeclibReader(res("fragcast_lib.parquet"), es, allowed).getPreds();
        assertEquals(1, filtered.keySet().size(), "allowedPrecursors filter not applied");
        assertTrue(filtered.containsKey(one));
        es.shutdown();
    }
}
