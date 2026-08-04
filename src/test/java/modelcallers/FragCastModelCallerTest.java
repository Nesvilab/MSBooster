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

package modelcallers;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertFalse;
import static org.junit.jupiter.api.Assertions.assertNotEquals;
import static org.junit.jupiter.api.Assertions.assertTrue;

import java.util.HashMap;
import java.util.List;

import allconstants.Constants;
import allconstants.FragCastModels;
import allconstants.LowercaseModelMapper;
import allconstants.ModelCollections;
import bestmodelsearch.ModelCollectionDecider;
import org.junit.jupiter.api.BeforeEach;
import org.junit.jupiter.api.Test;
import utils.Model;

// Pins the two contracts that let MSBooster reach FragCast's small/fast Spec model: the model names
// a user may select, and the command line handed to the FragCast executable. Neither test starts a
// process - buildCommand is the seam.
public class FragCastModelCallerTest {

    @BeforeEach
    public void setExecutablePath() {
        Constants.FragCast = "fragcast";
    }

    @Test
    public void bothSpecVariantsAreFragCastAndOnlyOneIsFast() {
        assertTrue(FragCastModels.isFragCast(FragCastModels.CONFORMER));
        assertTrue(FragCastModels.isFragCast(FragCastModels.FAST));
        assertFalse(FragCastModels.isFragCast("DIA-NN"));
        assertFalse(FragCastModels.isFragCast("Prosit_2019_intensity"));
        assertFalse(FragCastModels.isFragCast(null));

        assertTrue(FragCastModels.usesFastSpecModel(FragCastModels.FAST));
        assertFalse(FragCastModels.usesFastSpecModel(FragCastModels.CONFORMER));
    }

    // MainClass resolves user-supplied model names through this map and exits when one is missing,
    // so an unmapped name makes the model unselectable however well the rest is wired.
    @Test
    public void theFastModelNameIsSelectableInAnyCasing() {
        HashMap<String, String> mapper = new LowercaseModelMapper().getLowercaseToModel();
        for (String typed : new String[] {"FragCast-Fast", "fragcast-fast", "FRAGCAST-FAST"}) {
            assertEquals(FragCastModels.FAST, mapper.get(typed.toLowerCase()),
                    "not selectable as " + typed);
        }
        assertEquals(FragCastModels.CONFORMER, mapper.get("fragcast"));
    }

    // The fast model is opt-in: it costs a second FragCast run to compare against the Conformer, and
    // RT/IM would be identical, so it is deliberately not a default findBest candidate.
    @Test
    public void theFastModelIsNotADefaultFindBestCandidate() {
        for (List<String> collection : List.of(
                ModelCollections.generalMS2models, ModelCollections.hlaMS2models,
                ModelCollections.timstofMS2models, ModelCollections.generalRTmodels,
                ModelCollections.generalIMmodels)) {
            assertTrue(collection.contains(FragCastModels.CONFORMER));
            assertFalse(collection.contains(FragCastModels.FAST),
                    "the fast model became a default candidate: " + collection);
        }
    }

    // Being absent from the default collections does not put the fast model out of findBest's reach:
    // the *SearchModelsString params resolve names through LowercaseModelMapper, so listing it there
    // makes findBest score it against the Conformer.
    @Test
    public void theFastModelCanBeNamedInTheFindBestSearchList() {
        String saved = ModelCollections.ms2SearchModelsString;
        try {
            ModelCollections.ms2SearchModelsString = "FragCast,fragcast-fast";
            assertEquals(List.of(FragCastModels.CONFORMER, FragCastModels.FAST),
                    ModelCollectionDecider.decideCollection("MS2"));
        } finally {
            ModelCollections.ms2SearchModelsString = saved;
        }
    }

    // One MSBooster job runs FragCast at most once, so both Spec variants cache their library under
    // the same key. Other predictors keep keying on themselves.
    @Test
    public void bothSpecVariantsShareOnePredictionKey() {
        assertEquals(FragCastModels.predictionKey(FragCastModels.CONFORMER),
                FragCastModels.predictionKey(FragCastModels.FAST));
        assertEquals("DIA-NN", FragCastModels.predictionKey("DIA-NN"));
        assertNotEquals(FragCastModels.predictionKey("DIA-NN"),
                FragCastModels.predictionKey(FragCastModels.CONFORMER));
    }

    // Only the MS2 intensities depend on which Spec weights load, so the spectra model chooses them
    // for the whole job and RT/IM follow: the numbers they read are identical either way.
    @Test
    public void theSpectraModelDecidesTheJobsSpecVariant() {
        assertTrue(FragCastModels.jobUsesFastSpecModel(List.of(
                new Model(FragCastModels.FAST, "spectra"),
                new Model(FragCastModels.CONFORMER, "RT"),
                new Model(FragCastModels.CONFORMER, "IM"))));
        assertFalse(FragCastModels.jobUsesFastSpecModel(List.of(
                new Model(FragCastModels.CONFORMER, "spectra"),
                new Model(FragCastModels.FAST, "RT"),
                new Model(FragCastModels.FAST, "IM"))));
    }

    // With no FragCast spectra model the choice cannot change a number, so the first FragCast model
    // named wins - the job never loads weights (FragCast-Spec-Fast.onnx) the user did not ask for.
    @Test
    public void withoutAFragCastSpectraModelTheFirstFragCastModelDecides() {
        assertTrue(FragCastModels.jobUsesFastSpecModel(List.of(
                new Model("Prosit_2020_intensity_HCD", "spectra"),
                new Model(FragCastModels.FAST, "RT"))));
        assertFalse(FragCastModels.jobUsesFastSpecModel(List.of(
                new Model("Prosit_2020_intensity_HCD", "spectra"),
                new Model(FragCastModels.CONFORMER, "RT"),
                new Model(FragCastModels.FAST, "IM"))));
        assertFalse(FragCastModels.jobUsesFastSpecModel(List.of(new Model("DIA-NN", "RT"))));
    }

    @Test
    public void fastAddsTheFlagToTheBuildLibraryCall() {
        List<String> fast = FragCastModelCaller.buildCommand("spectraRT.tsv", "out.parquet", true);
        assertEquals(List.of("fragcast", "--task", "build-library", "--fast", "--in", "spectraRT.tsv"),
                fast.subList(0, 6));
        assertTrue(fast.containsAll(List.of("--out", "out.parquet", "--format", "parquet")));
    }

    @Test
    public void theDefaultCallIsUnchangedAndPassesNoFastFlag() {
        List<String> conformer = FragCastModelCaller.buildCommand("spectraRT.tsv", "out.parquet", false);
        assertFalse(conformer.contains("--fast"));
        assertEquals(List.of("fragcast", "--task", "build-library", "--in", "spectraRT.tsv"),
                conformer.subList(0, 5));
        // the fragment-filter flags FragCast's build-library task reads, and their param sources
        assertTrue(conformer.containsAll(List.of("--top-n", "--min-frag-mz", "--min-rel-intensity",
                "--min-frag-size", "--threads")));
        assertEquals(String.valueOf(Constants.fragCastTopN),
                conformer.get(conformer.indexOf("--top-n") + 1));
    }

    // A job runs one variant, but findBest scores both against each other and a rerun may switch
    // variants while keeping its predictions; sharing an output path would mix the two libraries up.
    @Test
    public void theTwoVariantsWriteToDifferentFiles() {
        String conformer = FragCastModelCaller.predictionFile("spectraRT.tsv", false);
        String fast = FragCastModelCaller.predictionFile("spectraRT.tsv", true);
        assertEquals("spectraRT.predicted.parquet", conformer);
        assertEquals("spectraRT.fast.predicted.parquet", fast);
        assertNotEquals(conformer, fast);
        // ParquetSpeclibReader is chosen by extension, so both must stay .parquet
        assertTrue(fast.endsWith(".parquet"));
    }
}
