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
import static org.junit.jupiter.api.Assertions.assertThrows;
import static org.junit.jupiter.api.Assertions.assertTrue;

import java.io.File;
import java.io.IOException;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.List;

import allconstants.Constants;
import allconstants.FragCastWeights;
import modelcallers.FragCastFineTuneJob;
import modelcallers.FragCastFineTuneOutcome;
import org.junit.jupiter.api.AfterEach;
import org.junit.jupiter.api.BeforeEach;
import org.junit.jupiter.api.Test;
import org.junit.jupiter.api.io.TempDir;

// Two things this class decides. A fine-tune produces loose ONNX files and a bundle holding copies
// of them, and what may be deleted afterwards is the first: getting it wrong removes weights the
// installation - or a model the user pointed a parameter at - still needs. The second is which
// failed task ends the run, and that depends on whether the user handed FragCast options of their
// own; FragCast never runs here, a fake stands in for each task.
public class FragCastTransferLearnerTest {
    private String cmdOptsBefore;
    private String modelDirBefore;

    @TempDir
    Path tmp;

    @BeforeEach
    public void saveConstants() {
        cmdOptsBefore = Constants.FragCastFineTuneCmdOpts;
        modelDirBefore = Constants.FragCastModelDir;
        //named outright, so the starting weights do not depend on where a FragCast happens to sit
        Constants.FragCastModelDir = tmp.resolve("models").toString();
    }

    @AfterEach
    public void restoreConstants() {
        Constants.FragCastFineTuneCmdOpts = cmdOptsBefore;
        Constants.FragCastModelDir = modelDirBefore;
    }

    // The bundle is the product; the loose ONNX files are what went into it. Deleting a file the
    // fine-tune did not write - the shipped pretrained weights, or a model a parameter names -
    // would take away weights that are not this run's to remove.
    @Test
    public void onlyTheModelsThisRunWroteAreRemoved(@TempDir Path dir) throws Exception {
        Path pretrained = dir.resolve("FragCast-IM.onnx");
        Path finetunedRt = dir.resolve("FragCast-RT.finetuned.onnx");
        Path finetunedSpec = dir.resolve("FragCast-SPEC.finetuned.onnx");
        for (Path p : new Path[] {pretrained, finetunedRt, finetunedSpec}) {
            Files.write(p, "onnx".getBytes(StandardCharsets.UTF_8));
        }

        FragCastWeights before = FragCastWeights.base().withTask("im", pretrained.toString());
        FragCastWeights accepted = before
                .withTask("rt", finetunedRt.toString())
                .withTask("spec", finetunedSpec.toString());

        assertEquals(2, FragCastTransferLearner.deleteLooseModels(accepted, before));
        assertFalse(Files.exists(finetunedRt), "the fine-tuned RT model was kept");
        assertFalse(Files.exists(finetunedSpec), "the fine-tuned MS2 model was kept");
        assertTrue(Files.exists(pretrained),
                "a model this run only started from was deleted along with its own output");
    }

    // A task that produced nothing still names whatever it started from, so a run where every task
    // was skipped names only weights it did not write. None of those may be removed.
    @Test
    public void aRunThatFineTunedNothingDeletesNothing(@TempDir Path dir) throws Exception {
        Path weights = dir.resolve("FragCast-RT.onnx");
        Files.write(weights, "onnx".getBytes(StandardCharsets.UTF_8));
        FragCastWeights before = FragCastWeights.base().withTask("rt", weights.toString());

        assertEquals(0, FragCastTransferLearner.deleteLooseModels(before, before));
        assertTrue(Files.exists(weights));
    }

    // A path naming a file that is not there is passed over rather than reported: FragCast can
    // finish a task without exporting anything, and the caller has already said so.
    @Test
    public void aPathThatNamesNoFileIsPassedOver(@TempDir Path dir) {
        FragCastWeights before = FragCastWeights.base();
        FragCastWeights accepted = before.withTask("rt", dir.resolve("absent.onnx").toString());
        assertEquals(0, FragCastTransferLearner.deleteLooseModels(accepted, before));
    }

    // --- which failed task ends the run ---------------------------------------------------------

    /** What FragCast does for each task, in place of running it; the tasks asked for are recorded. */
    private static final class FakeFragCast implements FragCastTransferLearner.Trainer {
        final List<String> ran = new ArrayList<>();
        private final String failing;
        private final FragCastTransferLearner.Trainer failure;

        /** Every task exports, except {@code failing}, which does what {@code failure} does. */
        FakeFragCast(String failing, FragCastTransferLearner.Trainer failure) {
            this.failing = failing;
            this.failure = failure;
        }

        @Override
        public FragCastFineTuneOutcome fineTune(FragCastFineTuneJob job) throws Exception {
            ran.add(job.task);
            return job.task.equals(failing) ? failure.fineTune(job) : exported(job);
        }
    }

    private static FragCastFineTuneOutcome exported(FragCastFineTuneJob job) throws IOException {
        Files.write(job.saveOnnx.toPath(), "onnx".getBytes(StandardCharsets.UTF_8));
        return FragCastFineTuneOutcome.parse(job.task, 0, job.saveOnnx,
                List.of("[t INFO] Wrote fine-tuned model -> " + job.saveOnnx));
    }

    private static FragCastTransferLearner.Trainer exitsWith(int code) {
        return job -> FragCastFineTuneOutcome.parse(job.task, code, job.saveOnnx,
                List.of("[t ERROR] unknown arg: --bogus"));
    }

    /** Exit 0 and no export line, as an rt --grid sweep ends. */
    private static FragCastFineTuneOutcome wroteNothing(FragCastFineTuneJob job) {
        return FragCastFineTuneOutcome.parse(job.task, 0, job.saveOnnx, List.of("[t INFO] grid done"));
    }

    private static FragCastFineTuneOutcome cannotStart(FragCastFineTuneJob job) throws IOException {
        throw new IOException("Cannot run program \"fragcast\"");
    }

    private FragCastWeights fineTuneAll(FragCastWeights before, FakeFragCast fake)
            throws FragCastTransferLearner.CmdOptsFailure {
        return FragCastTransferLearner.fineTuneAll(before, tmp.resolve("library.tsv").toFile(),
                tmp.toFile(), false, fake);
    }

    private File saved(String task) {
        return tmp.resolve("FragCast-" + task.toUpperCase() + ".finetuned.onnx").toFile();
    }

    // A fast MS2 fine-tune announces its export in words of its own. Read as no export, it would
    // stop a run with options set, blaming them for a model that was written.
    @Test
    public void withOptionsAFastMs2FineTuneThatExportedIsKept() throws Exception {
        Constants.FragCastFineTuneCmdOpts = "--epochs 5";
        FakeFragCast fake = new FakeFragCast("spec", job -> {
            assertTrue(job.baseOnnx.endsWith("FragCast-Spec-Fast.onnx"), job.baseOnnx);
            Files.write(job.saveOnnx.toPath(), "onnx".getBytes(StandardCharsets.UTF_8));
            return FragCastFineTuneOutcome.parse(job.task, 0, job.saveOnnx,
                    List.of("[t INFO] Wrote fine-tuned FAST model -> " + job.saveOnnx));
        });
        FragCastWeights accepted = FragCastTransferLearner.fineTuneAll(FragCastWeights.base(),
                tmp.resolve("library.tsv").toFile(), tmp.toFile(), true, fake);

        assertEquals(List.of(FragCastWeights.TASKS), fake.ran);
        assertEquals(saved("spec").getAbsolutePath(), accepted.forTask("spec"));
    }

    // With no options of the user's own, a failed task costs only that task, however it failed: the
    // others still run and their models are kept, and the pretrained weights stand in for the rest.
    @Test
    public void withoutOptionsAFailedTaskCostsOnlyThatTask() throws Exception {
        for (String blank : new String[] {"", "   "}) {
            Constants.FragCastFineTuneCmdOpts = blank;
            for (FragCastTransferLearner.Trainer failure : List.<FragCastTransferLearner.Trainer>of(
                    FragCastTransferLearnerTest::cannotStart, exitsWith(2),
                    FragCastTransferLearnerTest::wroteNothing)) {
                FakeFragCast fake = new FakeFragCast("im", failure);
                FragCastWeights accepted = fineTuneAll(FragCastWeights.base(), fake);

                assertEquals(List.of(FragCastWeights.TASKS), fake.ran, "a failed task stopped the others");
                assertEquals(saved("rt").getAbsolutePath(), accepted.forTask("rt"));
                assertEquals("", accepted.forTask("im"), "the failed task named a model");
                assertEquals(saved("spec").getAbsolutePath(), accepted.forTask("spec"));
            }
        }
    }

    // With options, they are the likeliest cause of a failure and every task gets the same ones, so
    // the first failure ends the run, naming them, before another task spends hours failing the same
    // way - whether FragCast could not be run, exited non-zero, or wrote no usable weights.
    @Test
    public void withOptionsAFailedTaskStopsTheRunAndNamesThem() {
        Constants.FragCastFineTuneCmdOpts = " --epochs 5 --fast ";
        for (FragCastTransferLearner.Trainer failure : List.<FragCastTransferLearner.Trainer>of(
                FragCastTransferLearnerTest::cannotStart, exitsWith(2), exitsWith(101),
                FragCastTransferLearnerTest::wroteNothing)) {
            FakeFragCast fake = new FakeFragCast("im", failure);
            FragCastTransferLearner.CmdOptsFailure e = assertThrows(
                    FragCastTransferLearner.CmdOptsFailure.class,
                    () -> fineTuneAll(FragCastWeights.base(), fake));

            assertEquals(List.of("rt", "im"), fake.ran, "a task ran after the run should have stopped");
            String message = e.getMessage();
            assertTrue(message.contains("FragCast im fine-tune"), message);
            assertTrue(message.contains("FragCastFineTuneCmdOpts"), message);
            assertTrue(message.contains("likely cause"), message);
            assertTrue(message.endsWith("--epochs 5 --fast"), "the options were not named: " + message);
        }
    }

    // The reason the task failed is part of the message, not lost behind it.
    @Test
    public void theStopSaysHowTheTaskFailed() {
        Constants.FragCastFineTuneCmdOpts = "--bogus 1";
        assertTrue(assertThrows(FragCastTransferLearner.CmdOptsFailure.class,
                () -> fineTuneAll(FragCastWeights.base(), new FakeFragCast("rt", exitsWith(2))))
                .getMessage().contains("FragCast exited with code 2"));
        assertTrue(assertThrows(FragCastTransferLearner.CmdOptsFailure.class,
                () -> fineTuneAll(FragCastWeights.base(),
                        new FakeFragCast("rt", FragCastTransferLearnerTest::cannotStart)))
                .getMessage().contains("Cannot run program"));
    }

    // A library with nothing to train a task on - no ion mobility, most often - says nothing about
    // the options, so it is skipped as it always was.
    @Test
    public void withOptionsALibraryWithNothingToTrainOnIsStillSkipped() throws Exception {
        Constants.FragCastFineTuneCmdOpts = "--epochs 5";
        FakeFragCast fake = new FakeFragCast("im",
                exitsWith(FragCastFineTuneOutcome.EXIT_NO_TRAINING_DATA));
        FragCastWeights accepted = fineTuneAll(FragCastWeights.base(), fake);

        assertEquals(List.of(FragCastWeights.TASKS), fake.ran);
        assertEquals("", accepted.forTask("im"));
        assertEquals(saved("spec").getAbsolutePath(), accepted.forTask("spec"));
    }

    // Nor does a task that would overwrite the weights it starts from; FragCast is never run for it.
    @Test
    public void withOptionsATaskThatWouldOverwriteItsOwnWeightsIsStillSkipped() throws Exception {
        Constants.FragCastFineTuneCmdOpts = "--epochs 5";
        FragCastWeights before = FragCastWeights.base().withTask("rt", saved("rt").getPath());
        FakeFragCast fake = new FakeFragCast("rt", job -> {
            throw new AssertionError("the rt task was run over its own weights");
        });
        FragCastWeights accepted = fineTuneAll(before, fake);

        assertEquals(List.of("im", "spec"), fake.ran);
        assertEquals(saved("rt").getPath(), accepted.forTask("rt"));
    }
}
