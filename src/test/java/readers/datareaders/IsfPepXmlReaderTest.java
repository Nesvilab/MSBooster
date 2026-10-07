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
import static org.junit.jupiter.api.Assertions.assertNotNull;
import static org.junit.jupiter.api.Assertions.assertNull;
import static org.junit.jupiter.api.Assertions.assertTrue;

import java.io.File;
import java.io.IOException;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.List;
import java.util.Map;
import java.util.stream.Collectors;
import org.junit.jupiter.api.Test;
import org.junit.jupiter.api.io.TempDir;
import peptideptmformatting.PeptideFormatter;

/**
 * MSFragger writes in-source fragment (ISF) annotations as extra {@code isf_*} search scores in the
 * pepXML it writes next to the pin. These tests pin down how MSBooster finds those files for a pin
 * and turns them into (scan, rank) keyed annotations.
 */
class IsfPepXmlReaderTest {

    @TempDir
    Path dir;

    private static String header() {
        return "<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n"
                + "<?xml-stylesheet type=\"text/xsl\" href=\"pepXML_std.xsl\"?>\n"
                + "<msms_pipeline_analysis date=\"2026-10-01\" xmlns=\"http://regis-web.systemsbiology.net/pepXML\">\n"
                + "<msms_run_summary base_name=\"run\" raw_data_type=\"mzML\" raw_data=\"mzML\">\n";
    }

    private static String footer() {
        return "</msms_run_summary>\n</msms_pipeline_analysis>\n";
    }

    private static String query(int scan, String... hits) {
        StringBuilder sb = new StringBuilder();
        sb.append("<spectrum_query start_scan=\"").append(scan).append("\" end_scan=\"").append(scan)
                .append("\" spectrum=\"run.").append(scan).append('.').append(scan)
                .append(".2\" assumed_charge=\"2\" index=\"1\">\n<search_result>\n");
        for (String hit : hits) {
            sb.append(hit);
        }
        sb.append("</search_result>\n</spectrum_query>\n");
        return sb.toString();
    }

    private static String plainHit(int rank, String peptide) {
        return "<search_hit peptide=\"" + peptide + "\" hit_rank=\"" + rank + "\" protein=\"sp|P1|P1\">\n"
                + "<search_score name=\"hyperscore\" value=\"20.0\"/>\n"
                + "<search_score name=\"fI_nterm\" value=\"0.0\"/>\n"
                + "</search_hit>\n";
    }

    private static String isfHit(int rank, String peptide, String parentModified, String parentCharge) {
        return "<search_hit peptide=\"" + peptide + "\" hit_rank=\"" + rank + "\" protein=\"sp|P1|P1\">\n"
                + "<modification_info modified_peptide=\"" + peptide + "\"/>\n"
                + "<search_score name=\"hyperscore\" value=\"20.0\"/>\n"
                + "<search_score name=\"fI_nterm\" value=\"0.0\"/>\n"
                + "<search_score name=\"isf_parent_peptide\" value=\"" + parentModified + "\"/>\n"
                + "<search_score name=\"isf_parent_charge\" value=\"" + parentCharge + "\"/>\n"
                + "<search_score name=\"isf_apex_rt_delta\" value=\"0.0123\"/>\n"
                + "<ptm_result ptm=\"x\"/>\n"
                + "</search_hit>\n";
    }

    private void write(String name, String content) throws IOException {
        Files.write(dir.resolve(name), content.getBytes(StandardCharsets.UTF_8));
    }

    private File pin(String base) throws IOException {
        Path pin = dir.resolve(base + ".pin");
        Files.write(pin, "SpecId\tLabel\tScanNr\n".getBytes(StandardCharsets.UTF_8));
        return pin.toFile();
    }

    @Test
    void ddaPepXmlUsesHitRank() throws Exception {
        write("run.pepXML", header()
                + query(100, plainHit(1, "PEPTIDEKR"), isfHit(2, "PEPTIDE", "PEPTIDEKR", "3"))
                + query(200, isfHit(1, "EPTIDEK", "n[42.0106]PEPTM[15.9949]IDEK", "2"))
                + footer());

        Map<String, IsfAnnotation> annotations = IsfPepXmlReader.readForPin(pin("run"));

        assertEquals(2, annotations.size());
        assertNull(annotations.get(IsfAnnotation.key(100, 1)), "a hit without isf_* scores is not an ISF");
        IsfAnnotation rank2 = annotations.get(IsfAnnotation.key(100, 2));
        assertNotNull(rank2);
        assertEquals("PEPTIDEKR", rank2.parentModifiedPeptide);
        assertEquals(3, rank2.parentCharge);
        IsfAnnotation modified = annotations.get(IsfAnnotation.key(200, 1));
        assertNotNull(modified);
        assertEquals("n[42.0106]PEPTM[15.9949]IDEK", modified.parentModifiedPeptide);
        assertEquals(2, modified.parentCharge);
    }

    @Test
    void diaRankFilesTakeTheRankFromTheFileName() throws Exception {
        // hit_rank is forced to 1 in every _rank<N> file, so N must come from the name
        write("run_rank1.pepXML", header()
                + query(300, isfHit(1, "PEPTIDE", "PEPTIDEKR", "2"))
                + footer());
        write("run_rank2.pepXML", header()
                + query(300, isfHit(1, "EPTIDEKR", "PEPTIDEKR", "3"))
                + query(301, plainHit(1, "AAAAK"))
                + footer());

        Map<String, IsfAnnotation> annotations = IsfPepXmlReader.readForPin(pin("run"));

        assertEquals(2, annotations.size());
        assertEquals(2, annotations.get(IsfAnnotation.key(300, 1)).parentCharge);
        assertEquals(3, annotations.get(IsfAnnotation.key(300, 2)).parentCharge);
    }

    @Test
    void locatesPlainAndRankPepXmlsOfThePinOnly() throws Exception {
        write("run.pepXML", header() + footer());
        write("run_rank1.pepXML", header() + footer());
        write("run_rank12.pepXML", header() + footer());
        write("run_rankX.pepXML", header() + footer());
        // the base is a prefix, a suffix or an infix of these names, never the whole name
        write("run2.pepXML", header() + footer());
        write("run2_rank1.pepXML", header() + footer());
        write("other_run.pepXML", header() + footer());
        write("a_run_b.pepXML", header() + footer());
        write("run.pin.pepXML", header() + footer());

        List<String> names = IsfPepXmlReader.locate(pin("run")).entrySet().stream()
                .map(e -> e.getKey().getFileName() + "@" + e.getValue())
                .sorted()
                .collect(Collectors.toList());

        assertEquals(List.of("run.pepXML@0", "run_rank1.pepXML@1", "run_rank12.pepXML@12"), names);
    }

    @Test
    void noPepXmlMeansNoAnnotations() throws Exception {
        assertTrue(IsfPepXmlReader.readForPin(pin("run")).isEmpty());
    }

    @Test
    void pepXmlWithoutIsfScoresMeansNoAnnotations() throws Exception {
        write("run.pepXML", header() + query(1, plainHit(1, "PEPTIDEK")) + footer());

        assertTrue(IsfPepXmlReader.readForPin(pin("run")).isEmpty());
    }

    @Test
    void malformedPepXmlIsIgnoredNotFatal() throws Exception {
        write("run.pepXML", header()
                + query(1, isfHit(1, "PEPTIDE", "PEPTIDEK", "2"))
                + "<spectrum_query start_scan=\"2\"><search_result>"); // truncated file

        assertTrue(IsfPepXmlReader.readForPin(pin("run")).isEmpty(),
                "a file that cannot be parsed to the end contributes nothing, not a partial map");
    }

    @Test
    void hitsWithIncompleteOrInvalidIsfScoresAreSkipped() throws Exception {
        String missingCharge = "<search_hit peptide=\"PEPTIDE\" hit_rank=\"1\">\n"
                + "<search_score name=\"isf_parent_peptide\" value=\"PEPTIDEK\"/>\n"
                + "</search_hit>\n";
        write("run.pepXML", header()
                + query(1, missingCharge)
                + query(2, isfHit(1, "PEPTIDE", "PEPTIDEK", "two"))
                + query(3, isfHit(1, "PEPTIDE", "", "2"))
                + query(4, isfHit(1, "PEPTIDE", "PEPTIDEK", "0"))
                + query(5, isfHit(1, "PEPTIDE", "PEPTIDEK", "2"))
                + footer());

        Map<String, IsfAnnotation> annotations = IsfPepXmlReader.readForPin(pin("run"));

        assertEquals(1, annotations.size());
        assertNotNull(annotations.get(IsfAnnotation.key(5, 1)));
    }

    @Test
    void parentKeyMatchesTheKeyOfTheParentPinRow() {
        IsfAnnotation annotation = new IsfAnnotation("n[42.0106]PEPTM[15.9949]IDEKc[0.9840]", 2);
        // the parent hit is itself a pin row; its Peptide column is prev.peptide+charge.next
        PeptideFormatter parentRow = new PeptideFormatter("K.n[42.0106]PEPTM[15.9949]IDEKc[0.9840]2.A", "2", "pin");

        assertEquals(parentRow.getBaseCharge(), annotation.parentBaseCharge);
        assertEquals("[42.0106]PEPTM[15.9949]IDEK[0.9840]|2", annotation.parentBaseCharge);
    }

    @Test
    void unmodifiedParentKeyIsThePlainSequence() {
        assertEquals("PEPTIDEK|3", new IsfAnnotation("PEPTIDEK", 3).parentBaseCharge);
    }
}
