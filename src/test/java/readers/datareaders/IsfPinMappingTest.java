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

import java.io.File;
import java.io.IOException;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.List;
import java.util.Map;
import org.junit.jupiter.api.Test;
import org.junit.jupiter.api.io.TempDir;

/**
 * pepXML annotations reach pin rows through (ScanNr, rank): the pin's ScanNr column and its rank
 * column, or the SpecId's "_rank" suffix when a pin has no rank column.
 */
class IsfPinMappingTest {

    @TempDir
    Path dir;

    private static final String PEPXML_HEAD = "<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n"
            + "<msms_pipeline_analysis xmlns=\"http://regis-web.systemsbiology.net/pepXML\">\n"
            + "<msms_run_summary base_name=\"run\">\n";
    private static final String PEPXML_TAIL = "</msms_run_summary>\n</msms_pipeline_analysis>\n";

    private static String isfQuery(int scan, int rank, String parent, int charge) {
        return "<spectrum_query start_scan=\"" + scan + "\" end_scan=\"" + scan + "\">\n<search_result>\n"
                + "<search_hit peptide=\"X\" hit_rank=\"" + rank + "\">\n"
                + "<search_score name=\"isf_parent_peptide\" value=\"" + parent + "\"/>\n"
                + "<search_score name=\"isf_parent_charge\" value=\"" + charge + "\"/>\n"
                + "</search_hit>\n</search_result>\n</spectrum_query>\n";
    }

    private File writePin(String header, String... rows) throws IOException {
        StringBuilder sb = new StringBuilder(header).append('\n');
        for (String row : rows) {
            sb.append(row).append('\n');
        }
        Path pin = dir.resolve("run.pin");
        Files.write(pin, sb.toString().getBytes(StandardCharsets.UTF_8));
        return pin.toFile();
    }

    private static List<String> annotatedRows(File pin, Map<String, IsfAnnotation> annotations) throws IOException {
        List<String> matched = new ArrayList<>();
        PinReader reader = new PinReader(pin.getAbsolutePath());
        try {
            while (reader.next(true)) {
                IsfAnnotation annotation = annotations.get(IsfAnnotation.key(reader.getScanNum(), reader.getRank()));
                if (annotation != null) {
                    matched.add(reader.getRow()[reader.specIdx] + "->" + annotation.parentBaseCharge);
                }
            }
        } finally {
            reader.close();
        }
        return matched;
    }

    @Test
    void ddaRowsMatchByScanNrAndRankColumn() throws Exception {
        Files.write(dir.resolve("run.pepXML"), (PEPXML_HEAD
                + isfQuery(10, 2, "PEPTIDEKR", 3)
                + isfQuery(11, 1, "AAAPEPTIDEK", 2)
                + PEPXML_TAIL).getBytes(StandardCharsets.UTF_8));
        File pin = writePin("SpecId\tLabel\tScanNr\tExpMass\tretentiontime\trank\tPeptide\tProteins",
                "run.10.10.2_1\t1\t10\t800.0\t30.0\t1\tK.EPTIDEKR2.A\tsp|P1|P1",
                "run.10.10.2_2\t1\t10\t800.0\t30.0\t2\tK.PEPTIDE2.K\tsp|P1|P1",
                "run.11.11.2_1\t-1\t11\t700.0\t30.1\t1\tK.PEPTIDEK2.A\tsp|P1|P1",
                "run.12.12.2_1\t1\t12\t700.0\t30.1\t1\tK.AAAPEPTIDEK2.A\tsp|P1|P1");

        List<String> matched = annotatedRows(pin, IsfPepXmlReader.readForPin(pin));

        assertEquals(List.of("run.10.10.2_2->PEPTIDEKR|3", "run.11.11.2_1->AAAPEPTIDEK|2"), matched,
                "decoys are matched too: the ISF verdict never looks at the target/decoy label");
    }

    @Test
    void diaRowsMatchTheRankOfTheirRankFile() throws Exception {
        Files.write(dir.resolve("run_rank1.pepXML"), (PEPXML_HEAD
                + PEPXML_TAIL).getBytes(StandardCharsets.UTF_8));
        Files.write(dir.resolve("run_rank3.pepXML"), (PEPXML_HEAD
                + isfQuery(500, 1, "PEPTIDEKR", 3)
                + PEPXML_TAIL).getBytes(StandardCharsets.UTF_8));
        //no rank column: the rank comes from the SpecId suffix
        File pin = writePin("SpecId\tLabel\tScanNr\tExpMass\tretentiontime\tPeptide\tProteins",
                "run.500.500.2_1\t1\t500\t800.0\t30.0\tK.PEPTIDEKR3.A\tsp|P1|P1",
                "run.500.500.2_3\t1\t500\t800.0\t30.0\tK.PEPTIDE2.K\tsp|P1|P1");

        List<String> matched = annotatedRows(pin, IsfPepXmlReader.readForPin(pin));

        assertEquals(List.of("run.500.500.2_3->PEPTIDEKR|3"), matched);
    }

    @Test
    void rankFallsBackToTheSpecIdSuffix() {
        assertEquals(4, PinReader.rankOf(new String[]{"run.1.1.2_4", "1", "1"}, -1, 0));
        assertEquals(2, PinReader.rankOf(new String[]{"run.1.1.2_4", "2"}, 1, 0));
    }
}
