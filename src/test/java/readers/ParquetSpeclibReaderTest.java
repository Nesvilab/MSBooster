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
