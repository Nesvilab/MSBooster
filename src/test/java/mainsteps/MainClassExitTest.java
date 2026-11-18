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

package mainsteps;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertFalse;
import static org.junit.jupiter.api.Assertions.assertTrue;
import static org.junit.jupiter.api.Assertions.fail;
import static org.junit.jupiter.api.Assumptions.assumeFalse;

import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.List;
import java.util.concurrent.TimeUnit;
import org.junit.jupiter.api.Test;
import org.junit.jupiter.api.io.TempDir;

/**
 * A run that cannot go on ends with exit status 1 and says why, without a stack trace: a model name MSBooster does
 * not know (D6), IM2Deep for a search with an isobaric label (D25), an empty NCE range to calibrate over (D27), and a
 * model that can take none of the PSMs (D28). The checks throw, and MainClass turns what they throw into exit
 * status 1. Each run is a JVM of its own, on a
 * synthetic mzML and pin file, and ends before any model is called; no Koina URL is set, so no run could reach Koina.
 */
class MainClassExitTest {
    //one MS2 scan of a Lumos, with an ion mobility so that useIM stays on
    private static final String MZML = String.join("\n",
            "<?xml version=\"1.0\" encoding=\"utf-8\"?>",
            "<mzML xmlns=\"http://psi.hupo.org/ms/mzml\" version=\"1.1.0\" id=\"run\">",
            "  <cvList count=\"1\">",
            "    <cv id=\"MS\" fullName=\"Proteomics Standards Initiative Mass Spectrometry Ontology\" version=\"4.1.0\""
                    + " URI=\"https://raw.githubusercontent.com/HUPO-PSI/psi-ms-CV/master/psi-ms.obo\"/>",
            "  </cvList>",
            "  <fileDescription>",
            "    <fileContent>",
            "      <cvParam cvRef=\"MS\" accession=\"MS:1000580\" name=\"MSn spectrum\" value=\"\"/>",
            "    </fileContent>",
            "  </fileDescription>",
            "  <softwareList count=\"1\">",
            "    <software id=\"synthetic\" version=\"1\"/>",
            "  </softwareList>",
            "  <instrumentConfigurationList count=\"1\">",
            "    <instrumentConfiguration id=\"IC1\">",
            "      <cvParam cvRef=\"MS\" accession=\"MS:1002416\" name=\"Orbitrap Fusion Lumos\" value=\"\"/>",
            "    </instrumentConfiguration>",
            "  </instrumentConfigurationList>",
            "  <dataProcessingList count=\"1\">",
            "    <dataProcessing id=\"dp\">",
            "      <processingMethod order=\"0\" softwareRef=\"synthetic\">",
            "        <cvParam cvRef=\"MS\" accession=\"MS:1000544\" name=\"Conversion to mzML\" value=\"\"/>",
            "      </processingMethod>",
            "    </dataProcessing>",
            "  </dataProcessingList>",
            "  <run id=\"run\" defaultInstrumentConfigurationRef=\"IC1\">",
            "    <spectrumList count=\"1\" defaultDataProcessingRef=\"dp\">",
            "      <spectrum index=\"0\" id=\"controllerType=0 controllerNumber=1 scan=1\" defaultArrayLength=\"2\">",
            "        <cvParam cvRef=\"MS\" accession=\"MS:1000511\" name=\"ms level\" value=\"2\"/>",
            "        <cvParam cvRef=\"MS\" accession=\"MS:1000580\" name=\"MSn spectrum\" value=\"\"/>",
            "        <cvParam cvRef=\"MS\" accession=\"MS:1000127\" name=\"centroid spectrum\" value=\"\"/>",
            "        <scanList count=\"1\">",
            "          <cvParam cvRef=\"MS\" accession=\"MS:1000795\" name=\"no combination\" value=\"\"/>",
            "          <scan>",
            "            <cvParam cvRef=\"MS\" accession=\"MS:1000016\" name=\"scan start time\" value=\"0.2\""
                    + " unitCvRef=\"UO\" unitAccession=\"UO:0000031\" unitName=\"minute\"/>",
            "            <cvParam cvRef=\"MS\" accession=\"MS:1002815\" name=\"inverse reduced ion mobility\""
                    + " value=\"0.95\" unitCvRef=\"MS\" unitAccession=\"MS:1002814\""
                    + " unitName=\"volt-second per square centimeter\"/>",
            "            <cvParam cvRef=\"MS\" accession=\"MS:1000512\" name=\"filter string\""
                    + " value=\"FTMS + p NSI d Full ms2 500.00@hcd30.00 [100.0000-2000.0000]\"/>",
            "          </scan>",
            "        </scanList>",
            "        <precursorList count=\"1\">",
            "          <precursor>",
            "            <selectedIonList count=\"1\">",
            "              <selectedIon>",
            "                <cvParam cvRef=\"MS\" accession=\"MS:1000744\" name=\"selected ion m/z\" value=\"500.0\""
                    + " unitCvRef=\"MS\" unitAccession=\"MS:1000040\" unitName=\"m/z\"/>",
            "                <cvParam cvRef=\"MS\" accession=\"MS:1000041\" name=\"charge state\" value=\"2\"/>",
            "              </selectedIon>",
            "            </selectedIonList>",
            "            <activation>",
            "              <cvParam cvRef=\"MS\" accession=\"MS:1000422\""
                    + " name=\"beam-type collision-induced dissociation\" value=\"\"/>",
            "              <cvParam cvRef=\"MS\" accession=\"MS:1000045\" name=\"collision energy\" value=\"30.0\"/>",
            "            </activation>",
            "          </precursor>",
            "        </precursorList>",
            "        <binaryDataArrayList count=\"2\">",
            "          <binaryDataArray encodedLength=\"24\">",
            "            <cvParam cvRef=\"MS\" accession=\"MS:1000523\" name=\"64-bit float\" value=\"\"/>",
            "            <cvParam cvRef=\"MS\" accession=\"MS:1000576\" name=\"no compression\" value=\"\"/>",
            "            <cvParam cvRef=\"MS\" accession=\"MS:1000514\" name=\"m/z array\" value=\"\""
                    + " unitCvRef=\"MS\" unitAccession=\"MS:1000040\" unitName=\"m/z\"/>",
            "            <binary>AAAAAAAAWUAAAAAAAABpQA==</binary>", //100, 200
            "          </binaryDataArray>",
            "          <binaryDataArray encodedLength=\"24\">",
            "            <cvParam cvRef=\"MS\" accession=\"MS:1000523\" name=\"64-bit float\" value=\"\"/>",
            "            <cvParam cvRef=\"MS\" accession=\"MS:1000576\" name=\"no compression\" value=\"\"/>",
            "            <cvParam cvRef=\"MS\" accession=\"MS:1000515\" name=\"intensity array\" value=\"\""
                    + " unitCvRef=\"MS\" unitAccession=\"MS:1000131\" unitName=\"number of detector counts\"/>",
            "            <binary>AAAAAAAAJEAAAAAAAAA0QA==</binary>", //10, 20
            "          </binaryDataArray>",
            "        </binaryDataArrayList>",
            "      </spectrum>",
            "    </spectrumList>",
            "  </run>",
            "</mzML>",
            "");

    private static final String PIN = pin("K.LGEPTIDEK.A");

    //a pin file of one PSM, of the mzML's one scan
    private static String pin(String peptide) {
        return "SpecId\tLabel\tScanNr\tExpMass\tretentiontime\trank\thyperscore\tlog10_evalue\tPeptide\tProteins\n"
                + "run.1.1.2_1\t1\t1\t998.0\t12.0\t1\t30.0\t-8.0\t" + peptide + "\tsp|SYNTH|SYNTH_PEPTIDE\n";
    }

    private static final String NO_LABEL = "fragment_mass_tolerance = 20\nfragment_mass_units = 1\n"
            + "add_C_cysteine = 57.021464\n";
    private static final String TMT = "fragment_mass_tolerance = 20\nfragment_mass_units = 1\n"
            + "add_K_lysine = 229.162932\nadd_Nterm_peptide = 229.162932\n";

    @TempDir
    Path dir;

    @Test
    void anUnknownModelEndsTheRunWithStatus1() throws Exception {
        Run run = run(NO_LABEL, "useSpectra = true", "spectraModel = PredFull");
        run.assertEndedWith("No spectra model called PredFull (spectraModel). The Koina spectra models MSBooster "
                + "can call are ");
    }

    @Test
    void im2DeepForASearchWithTmtEndsTheRunWithStatus1() throws Exception {
        Run run = run(TMT, "useIM = true", "imModel = IM2Deep");
        run.assertEndedWith("IM2Deep cannot encode isobaric labels (TMT/iTRAQ), and the search in fragger.params "
                + "has one (229.1629). Use AlphaPept_ccs_generic instead");
    }

    @Test
    void anEmptyNceRangeEndsTheRunWithStatus1() throws Exception {
        Run run = run(NO_LABEL, "useSpectra = true", "spectraModel = Prosit_2020_intensity_HCD", "minNCE = 45");
        run.assertEndedWith("The NCE range to calibrate Prosit_2020_intensity_HCD over is empty: from 45 (minNCE) to "
                + "40 (the highest NCE of Prosit_2020_intensity_HCD)");
    }

    //the check runs where the request file is written (PeptideFileCreator), before Koina is asked anything
    @Test
    void aModelThatCanTakeNoneOfThePsmsEndsTheRunWithStatus1() throws Exception {
        String residues51 = "ELVSAPLGTNEIYDHFAQWTGLPVEDAKLGRSTEPVNGHAKYLDTEGAIPK";
        assertEquals(51, residues51.length());
        Run run = runOn(pin("K." + residues51 + ".A"), NO_LABEL, "useRT = true", "rtModel = Chronologer_RT");
        run.assertEndedWith("Koina model Chronologer_RT cannot take any of the 1 PSMs: longer than 50 residues (1), "
                + "so it has nothing to predict. Choose another model");
    }

    private record Run(int status, String out) {
        //exit status 1, the reason as an error line, then "Exiting", and no stack trace
        void assertEndedWith(String reason) {
            assertEquals(1, status, out);
            List<String> errors = new ArrayList<>();
            for (String line : out.split("\\R")) {
                int at = line.indexOf("[ERROR] - ");
                if (at >= 0) {
                    errors.add(line.substring(at + "[ERROR] - ".length()));
                }
            }
            assertEquals(2, errors.size(), out);
            assertTrue(errors.get(0).startsWith(reason), out);
            assertEquals("Exiting", errors.get(1), out);
            assertFalse(out.contains("\tat "), out);
        }
    }

    //runs MSBooster in a JVM of its own on the synthetic files, with these parameters besides the files
    private Run run(String fraggerParams, String... params) throws Exception {
        return runOn(PIN, fraggerParams, params);
    }

    private Run runOn(String pin, String fraggerParams, String... params) throws Exception {
        assumeFalse(dir.toString().contains(" "), "MSBooster reads a list of paths separated by spaces");
        Files.writeString(dir.resolve("run.mzML"), MZML, StandardCharsets.UTF_8);
        Files.writeString(dir.resolve("run.pin"), pin, StandardCharsets.UTF_8);
        Files.createDirectories(dir.resolve("out")); //FragPipe hands MSBooster an output directory that exists
        Files.writeString(dir.resolve("fragger.params"), fraggerParams, StandardCharsets.UTF_8);
        List<String> lines = new ArrayList<>(List.of(
                "pinPepXMLDirectory = " + dir.resolve("run.pin"),
                "mzmlDirectory = " + dir.resolve("run.mzML"),
                "outputDirectory = " + dir.resolve("out"),
                "fragger = " + dir.resolve("fragger.params"),
                "instrument = LUMOS",
                "FragmentationType = HCD",
                "numThreads = 1",
                "useSpectra = false",
                "useRT = false",
                "useIM = false"));
        lines.addAll(List.of(params)); //a later line sets a parameter again
        Path paramsFile = dir.resolve("msbooster_params.txt");
        Files.write(paramsFile, lines, StandardCharsets.UTF_8);

        String classpath = System.getProperty("surefire.test.class.path", System.getProperty("java.class.path"));
        Path java = Path.of(System.getProperty("java.home"), "bin", "java");
        Path out = dir.resolve("stdout.txt");
        Process process = new ProcessBuilder(java.toString(), "-Djava.awt.headless=true", "-cp", classpath,
                MainClass.class.getName(), "--paramsList", paramsFile.toString())
                .directory(dir.toFile())
                .redirectErrorStream(true)
                .redirectOutput(out.toFile())
                .start();
        if (!process.waitFor(2, TimeUnit.MINUTES)) {
            process.destroyForcibly();
            fail("MSBooster did not end within 2 minutes:\n" + Files.readString(out, StandardCharsets.UTF_8));
        }
        return new Run(process.exitValue(), Files.readString(out, StandardCharsets.UTF_8));
    }
}
