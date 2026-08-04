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

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;

import allconstants.FragCastModels;
import org.junit.jupiter.api.Test;
import utils.Model;

// getPredictionFiles runs each local executable predictor once per job and reuses its library for
// every property that model serves. MainUtils.localPredFile is that caching, with the run itself
// injected so these tests need no FragCast or DIA-NN executable.
public class MainUtilsTest {

    //stands in for a real predictor call: records what was run, returns the file it would write.
    //Mirrors the real runner in that only FragCast reads the fast flag.
    private static class RecordingPredictor implements MainUtils.LocalPredictorRunner {
        final ArrayList<String> runs = new ArrayList<>();

        @Override
        public String run(String model, boolean fast) {
            if (!FragCastModels.isFragCast(model)) { //DIA-NN
                runs.add(model);
                return "spectraRT.predicted.bin";
            }
            runs.add(model + (fast ? " --fast" : ""));
            return "spectraRT" + (fast ? ".fast" : "") + ".predicted.parquet";
        }
    }

    private static ArrayList<String> resolveAll(List<Model> models, boolean fragCastFast,
                                                MainUtils.LocalPredictorRunner predictor)
            throws Exception {
        HashMap<String, String> cache = new HashMap<>();
        ArrayList<String> predFiles = new ArrayList<>();
        for (Model model : models) {
            predFiles.add(MainUtils.localPredFile(model.name, fragCastFast, cache, predictor));
        }
        return predFiles;
    }

    // The case this caching exists for: fast spectra with Conformer RT and IM. RT and IM are
    // identical under either set of Spec weights, so a second FragCast run would only reproduce
    // numbers the first already wrote - all three properties come from the one library.
    @Test
    public void mixedFragCastSpecVariantsRunFragCastOnce() throws Exception {
        List<Model> models = List.of(
                new Model(FragCastModels.FAST, "spectra"),
                new Model(FragCastModels.CONFORMER, "RT"),
                new Model(FragCastModels.CONFORMER, "IM"));
        RecordingPredictor predictor = new RecordingPredictor();

        ArrayList<String> predFiles = resolveAll(models,
                FragCastModels.jobUsesFastSpecModel(models), predictor);

        assertEquals(List.of("FragCast-Fast --fast"), predictor.runs, "FragCast ran more than once");
        assertEquals(List.of("spectraRT.fast.predicted.parquet", "spectraRT.fast.predicted.parquet",
                "spectraRT.fast.predicted.parquet"), predFiles,
                "spectra, RT and IM must read the same library");
    }

    // The same holds for a single variant used throughout, which is how FragCast already behaved.
    @Test
    public void oneFragCastVariantThroughoutStillRunsOnce() throws Exception {
        List<Model> models = List.of(
                new Model(FragCastModels.CONFORMER, "spectra"),
                new Model(FragCastModels.CONFORMER, "RT"),
                new Model(FragCastModels.CONFORMER, "IM"));
        RecordingPredictor predictor = new RecordingPredictor();

        ArrayList<String> predFiles = resolveAll(models,
                FragCastModels.jobUsesFastSpecModel(models), predictor);

        assertEquals(List.of("FragCast"), predictor.runs);
        assertEquals(3, predFiles.size());
        assertEquals("spectraRT.predicted.parquet", predFiles.get(0));
    }

    // FragCast and DIA-NN are separate executables writing separate libraries, so they must not
    // collapse into one cache entry the way the two FragCast Spec variants do.
    @Test
    public void differentLocalPredictorsEachGetTheirOwnRun() throws Exception {
        List<Model> models = List.of(
                new Model(FragCastModels.FAST, "spectra"),
                new Model("DIA-NN", "RT"),
                new Model("DIA-NN", "IM"));
        RecordingPredictor predictor = new RecordingPredictor();

        ArrayList<String> predFiles = resolveAll(models,
                FragCastModels.jobUsesFastSpecModel(models), predictor);

        assertEquals(List.of("FragCast-Fast --fast", "DIA-NN"), predictor.runs);
        assertEquals(List.of("spectraRT.fast.predicted.parquet", "spectraRT.predicted.bin",
                "spectraRT.predicted.bin"), predFiles);
    }

    private static ArrayList<String> inputFileModelNames(List<Model> models) {
        ArrayList<String> names = new ArrayList<>();
        for (Model model : MainUtils.inputFileModels(models)) {
            names.add(model.name);
        }
        return names;
    }

    // makeInputFiles dedupes on the same prediction key getPredictionFiles caches on. Both FragCast
    // Spec variants write spectraRT.tsv from the same peptides, so a job naming both must build it
    // once; keying on the model name instead built the identical file twice.
    @Test
    public void mixedFragCastSpecVariantsBuildInputFileOnce() {
        List<Model> models = List.of(
                new Model(FragCastModels.FAST, "spectra"),
                new Model(FragCastModels.CONFORMER, "RT"),
                new Model(FragCastModels.CONFORMER, "IM"));

        assertEquals(List.of(FragCastModels.FAST), inputFileModelNames(models),
                "FragCast input file was built more than once");
    }

    // The RT model naming FragCast first must not drag the job onto a name it did not ask for: the
    // file is the same either way, and jobUsesFastSpecModel - not this list - picks the Spec weights.
    @Test
    public void firstFragCastModelNamedStandsInForBothVariants() {
        List<Model> models = List.of(
                new Model(FragCastModels.CONFORMER, "RT"),
                new Model(FragCastModels.FAST, "spectra"));

        assertEquals(List.of(FragCastModels.CONFORMER), inputFileModelNames(models));
    }

    // Every other model keys on its own name, so deduping is unchanged for them: repeats collapse,
    // distinct models each keep their own input file.
    @Test
    public void nonFragCastModelsDedupeByNameAsBefore() {
        List<Model> models = List.of(
                new Model("Prosit_2020_intensity_HCD", "spectra"),
                new Model("Prosit_2020_intensity_HCD", "RT"),
                new Model("DIA-NN", "IM"));

        assertEquals(List.of("Prosit_2020_intensity_HCD", "DIA-NN"), inputFileModelNames(models));
    }

    // A FragCast job alongside another predictor keeps both input files.
    @Test
    public void fragCastAndDiannBothBuildInputFiles() {
        List<Model> models = List.of(
                new Model(FragCastModels.FAST, "spectra"),
                new Model("DIA-NN", "RT"),
                new Model("DIA-NN", "IM"));

        assertEquals(List.of(FragCastModels.FAST, "DIA-NN"), inputFileModelNames(models));
    }
}
