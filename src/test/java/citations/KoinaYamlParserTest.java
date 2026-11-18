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

package citations;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertFalse;
import static org.junit.jupiter.api.Assertions.assertNotNull;
import static org.junit.jupiter.api.Assertions.assertNull;
import static org.junit.jupiter.api.Assertions.assertTrue;

import allconstants.Constants;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;
import java.util.Set;
import koinaclasses.KoinaModelRegistry;
import koinaclasses.KoinaModelSpec;
import org.junit.jupiter.api.Test;
import org.junit.jupiter.api.io.TempDir;
import utils.Model;

class KoinaYamlParserTest {
    //the citation is looked up by the spec's key as it is, so every key must be a path of koina_citations.yml with a
    //citation to print
    @Test
    void everyKoinaModelHasACitation() {
        Map<String, String> keys = new LinkedHashMap<>();
        for (KoinaModelSpec spec : KoinaModelRegistry.all()) {
            keys.put(spec.name(), spec.citationKey());
        }

        KoinaYamlParser parser = new KoinaYamlParser();
        List<String> missing = new ArrayList<>();
        List<String> blank = new ArrayList<>();
        for (Map.Entry<String, String> e : keys.entrySet()) {
            String citation = parser.getCitation(e.getValue());
            if (citation == null) {
                missing.add(e.getKey() + " (" + e.getValue() + ")");
            } else if (citation.isBlank()) {
                blank.add(e.getKey() + " (" + e.getValue() + ")");
            }
        }
        assertTrue(missing.isEmpty(), "no path in koina_citations.yml for " + missing);
        assertTrue(blank.isEmpty(), "an empty citation in koina_citations.yml for " + blank);
    }

    @Test
    void theCitationIsTheTextOfItsSection() {
        KoinaYamlParser parser = new KoinaYamlParser();
        String multiFrag = parser.getCitation("Prosit_2025_intensity_MultiFrag");
        assertTrue(multiFrag.contains("doi.org/10.1101/2025.05.28.656555"), multiFrag);
        assertFalse(multiFrag.contains("<br>") || multiFrag.contains("\n"), multiFrag);
        assertTrue(parser.getCitation("Prosit_2025_intensity_lac").contains("Prosit_2025_intensity_lac"),
                "a model with no publication yet is named in the line printed for it");
    }

    @Test
    void anUnknownKeyHasNoCitation() {
        assertNull(new KoinaYamlParser().getCitation("No_such_model"));
        assertNull(new KoinaYamlParser().getCitation("AlphaPept_ms2_generic"),
                "the AlphaPeptDeep models are cited under their AlphaPeptDeep_* path, which their specs name");
    }

    //D6: the run's Koina models are cited in the params file under their specs' citation keys, each citation once;
    //a model that is not a Koina model is cited by no line, and a run without a Koina model writes none
    @Test
    void theCitationsOfTheRunsKoinaModelsAreAddedToTheParamsFile(@TempDir Path dir) throws Exception {
        Path params = dir.resolve("msbooster_params.txt");
        String saved = Constants.paramsList;
        try {
            Constants.paramsList = params.toString();
            KoinaYamlParser parser = new KoinaYamlParser();

            Files.writeString(params, "spectraModel=DIA-NN\n", StandardCharsets.UTF_8);
            parser.writeCitations(new ArrayList<>(List.of(new Model("DIA-NN", "spectra"),
                    new Model("FragCast", "RT"))));
            assertEquals(List.of("spectraModel=DIA-NN"), Files.readAllLines(params, StandardCharsets.UTF_8));

            //ms2pip_2021_HCD and the AlphaPept models have no path of their own name in koina_citations.yml
            Files.writeString(params, "spectraModel=ms2pip_2021_HCD\n", StandardCharsets.UTF_8);
            List<String> models = List.of("ms2pip_2021_HCD", "AlphaPept_rt_generic", "AlphaPept_ccs_generic",
                    "Altimeter_2024_intensities");
            ArrayList<Model> run = new ArrayList<>();
            run.add(new Model("ms2pip_2021_HCD", "spectra"));
            run.add(new Model("AlphaPept_rt_generic", "RT"));
            run.add(new Model("AlphaPept_ccs_generic", "IM"));
            run.add(new Model("Altimeter_2024_intensities", "spectra")); //as the auxiliary spectra model
            run.add(new Model("DIA-NN", "RT"));
            parser.writeCitations(run);

            List<String> lines = Files.readAllLines(params, StandardCharsets.UTF_8);
            assertEquals("spectraModel=ms2pip_2021_HCD", lines.get(0));
            assertEquals("#CITATIONS: If using models from Koina, please cite the following:", lines.get(1));
            assertTrue(lines.get(2).startsWith("#CITATIONS: Lautenbacher L, Yang KL, Kockmann T, et al. Koina:"),
                    lines.get(2));
            Set<String> expected = new HashSet<>();
            for (String model : models) {
                String citation = parser.getCitation(KoinaModelRegistry.require(model).citationKey());
                assertNotNull(citation, model);
                expected.add("#CITATIONS: " + citation);
            }
            //the AlphaPeptDeep models share one citation
            assertEquals(models.size() - 1, expected.size(), expected.toString());
            List<String> cited = lines.subList(3, lines.size());
            assertEquals(expected.size(), cited.size(), "each citation once: " + cited);
            assertEquals(expected, new HashSet<>(cited));
            assertTrue(cited.stream().anyMatch(line -> line.contains("doi:10.1093/nar/gkad335")), //MS2PIP 2023
                    cited.toString());
        } finally {
            Constants.paramsList = saved;
        }
    }
}
