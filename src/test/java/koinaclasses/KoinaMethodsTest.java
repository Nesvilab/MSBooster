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

import static org.junit.jupiter.api.Assertions.assertEquals;

import allconstants.Constants;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.List;
import java.util.Set;
import org.junit.jupiter.api.AfterEach;
import org.junit.jupiter.api.BeforeEach;
import org.junit.jupiter.api.Test;
import org.junit.jupiter.api.io.TempDir;
import peptideptmformatting.PeptideFormatter;
import readers.datareaders.IsfAnnotation;
import readers.datareaders.PinReader;

/**
 * The top PSMs of a pin pick and calibrate the models. An in-source fragment (ISF) PSM borrows its
 * parent's RT, so it must not help pick the RT model; its spectrum is a genuine fragmentation of its
 * own peptide, so it stays in the list the MS2 model and the NCE are picked with.
 */
class KoinaMethodsTest {
    private static final String MZML = "run.mzML";
    private static final String HEADER =
            "SpecId\tLabel\tScanNr\tExpMass\tretentiontime\trank\tlog10_evalue\tPeptide\tProteins";

    @TempDir
    Path dir;

    private Boolean originalUseRT;
    private Boolean originalUseIM;

    @BeforeEach
    void setUp() {
        originalUseRT = Constants.useRT;
        originalUseIM = Constants.useIM;
        Constants.useRT = true;
        Constants.useIM = false;
    }

    @AfterEach
    void tearDown() {
        Constants.useRT = originalUseRT;
        Constants.useIM = originalUseIM;
    }

    private static String row(int scan, String peptide, int log10Evalue) {
        return "run." + scan + "." + scan + ".2_1\t1\t" + scan + "\t800\t" + scan + ".0\t1\t" + log10Evalue
                + "\tK." + peptide + "2.A\tsp|P1|P1";
    }

    private static List<String> bases(List<PeptideFormatter> peptides) {
        List<String> bases = new ArrayList<>();
        for (PeptideFormatter pf : peptides) {
            bases.add(pf.getBase());
        }
        return bases;
    }

    @Test
    void isfPsmsStayInTheMs2ListAndLeaveTheRtList() throws Exception {
        Path pin = dir.resolve("run.pin");
        Files.write(pin, List.of(HEADER,
                row(1, "ISFPEPTIDE", -10),
                row(2, "PEPTIDEA", -8),
                row(3, "PEPTIDEAA", -6)), StandardCharsets.UTF_8);
        KoinaMethods km = new KoinaMethods(null);

        PinReader reader = new PinReader(pin.toString());
        try {
            km.addTopPSMs(MZML, reader, 2, Set.of(IsfAnnotation.key(1, 1)));
        } finally {
            reader.close();
        }

        assertEquals(List.of(1, 2), new ArrayList<>(km.scanNums.get(MZML)), "MS2/NCE keep the ISF PSM");
        assertEquals(List.of("ISFPEPTIDE", "PEPTIDEA"), bases(km.peptideArraylist));
        assertEquals(List.of(2, 3), new ArrayList<>(km.scanNumsRT.get(MZML)), "RT leaves the ISF PSM out");
        assertEquals(List.of("PEPTIDEA", "PEPTIDEAA"), bases(km.peptidesRT.get(MZML)));
        assertEquals(List.of("PEPTIDEA", "PEPTIDEAA"), bases(km.peptideArrayListRT));
    }
}
