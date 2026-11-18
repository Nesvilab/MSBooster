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

import static org.junit.jupiter.api.Assertions.assertDoesNotThrow;
import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertFalse;
import static org.junit.jupiter.api.Assertions.assertThrows;
import static org.junit.jupiter.api.Assertions.assertTrue;

import allconstants.Constants;
import allconstants.ModelCollections;
import bestmodelsearch.ModelCollectionDecider;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.HashMap;
import java.util.List;
import koinaclasses.KoinaException;
import koinaclasses.KoinaModelRegistry;
import koinaclasses.KoinaPreCheck;
import org.junit.jupiter.api.AfterEach;
import org.junit.jupiter.api.BeforeEach;
import org.junit.jupiter.api.Test;
import org.junit.jupiter.api.io.TempDir;
import peptideptmformatting.PTMhandler;

/**
 * IM2Deep cannot encode isobaric labels, so a search with TMT or iTRAQ in fragger.params cannot use it for IM (D25):
 * the run ends, pointing to AlphaPept_ccs_generic, and the search for the best IM model leaves it out.
 */
class IsobaricLabelCheckTest {
    @TempDir
    Path tmp;

    private String rtCollection;
    private String ms2Collection;
    private String imSearchModels;
    private Boolean matchWithDaltons;
    private Boolean matchWithDaltonsAux;
    private Boolean matchWithDaltonsDefault;
    private Float daTolerance;

    @BeforeEach
    void save() {
        rtCollection = ModelCollections.rtCollection;
        ms2Collection = ModelCollections.ms2Collection;
        imSearchModels = ModelCollections.imSearchModelsString;
        matchWithDaltons = Constants.matchWithDaltons;
        matchWithDaltonsAux = Constants.matchWithDaltonsAux;
        matchWithDaltonsDefault = Constants.matchWithDaltonsDefault;
        daTolerance = Constants.DaTolerance;
        PTMhandler.resetIsobaricLabel();
    }

    @AfterEach
    void restore() {
        ModelCollections.rtCollection = rtCollection;
        ModelCollections.ms2Collection = ms2Collection;
        ModelCollections.imSearchModelsString = imSearchModels;
        Constants.matchWithDaltons = matchWithDaltons;
        Constants.matchWithDaltonsAux = matchWithDaltonsAux;
        Constants.matchWithDaltonsDefault = matchWithDaltonsDefault;
        Constants.DaTolerance = daTolerance;
        PTMhandler.resetIsobaricLabel();
    }

    private void readFraggerParams(String... lines) throws Exception {
        Path fragger = tmp.resolve("fragger.params");
        Files.writeString(fragger, String.join("\n", lines) + "\n", StandardCharsets.UTF_8);
        HashMap<String, String> params = new HashMap<>();
        params.put("fragger", fragger.toString());
        ParameterUtils.readFraggerParams(params);
    }

    @Test
    void im2DeepWithATmtLabelEndsTheRunAndNamesTheModelToUse() throws Exception {
        readFraggerParams("fragment_mass_tolerance = 20", "fragment_mass_units = 1", "add_K_lysine = 229.162932",
                "add_Nterm_peptide = 229.162932");
        assertTrue(PTMhandler.hasIsobaricLabel());

        KoinaException e = assertThrows(KoinaException.class, () -> KoinaPreCheck.checkIsobaricLabel("IM2Deep"));
        assertEquals("IM2Deep cannot encode isobaric labels (TMT/iTRAQ), and the search in fragger.params has one "
                + "(229.1629). Use AlphaPept_ccs_generic instead", e.getMessage());
        assertDoesNotThrow(() -> KoinaPreCheck.checkIsobaricLabel("AlphaPept_ccs_generic"));
        assertDoesNotThrow(() -> KoinaPreCheck.checkIsobaricLabel("DIA-NN"), "not a Koina model");
        assertDoesNotThrow(() -> KoinaPreCheck.checkIsobaricLabel("Prosit_2020_intensity_TMT"));
    }

    @Test
    void aVariableItraqLabelEndsTheRunToo() throws Exception {
        readFraggerParams("fragment_mass_tolerance = 20", "fragment_mass_units = 1",
                "variable_mod_01 = 144.10207 n^ 1");
        KoinaException e = assertThrows(KoinaException.class, () -> KoinaPreCheck.checkIsobaricLabel("IM2Deep"));
        assertTrue(e.getMessage().contains("(144.1021)"), e.getMessage());
    }

    @Test
    void im2DeepWithoutALabelIsUsed() throws Exception {
        readFraggerParams("fragment_mass_tolerance = 20", "fragment_mass_units = 1", "add_C_cysteine = 57.021464",
                "variable_mod_01 = 15.9949 M 3");
        assertFalse(PTMhandler.hasIsobaricLabel());
        assertDoesNotThrow(() -> KoinaPreCheck.checkIsobaricLabel("IM2Deep"));
        assertTrue(ModelCollectionDecider.decideCollection("IM").contains("IM2Deep"));
    }

    @Test
    void theSearchForTheBestImModelLeavesIm2DeepOut() throws Exception {
        readFraggerParams("fragment_mass_tolerance = 20", "fragment_mass_units = 1", "add_K_lysine = 304.207146");
        assertTrue(PTMhandler.hasIsobaricLabel(), "TMTpro");

        List<String> considered = ModelCollectionDecider.decideCollection("IM");
        assertFalse(considered.contains("IM2Deep"), considered.toString());
        assertTrue(considered.contains("AlphaPept_ccs_generic"), considered.toString());
        assertTrue(ModelCollections.getImCollection().contains("IM2Deep"), "the collection itself is left as it is");

        ModelCollections.imSearchModelsString = "IM2Deep";
        KoinaException e = assertThrows(KoinaException.class, () -> ModelCollectionDecider.decideCollection("IM"));
        assertTrue(e.getMessage().startsWith("No IM model is left to search: IM2Deep cannot encode isobaric labels"),
                e.getMessage());
    }

    @Test
    void onlyIm2DeepCannotEncodeIsobaricLabels() {
        for (var spec : KoinaModelRegistry.all()) {
            assertEquals(!spec.name().equals("IM2Deep"), spec.takesIsobaricLabels(), spec.name());
        }
    }
}
