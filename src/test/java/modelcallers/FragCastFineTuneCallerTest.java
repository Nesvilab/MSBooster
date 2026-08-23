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
import static org.junit.jupiter.api.Assertions.assertNotNull;
import static org.junit.jupiter.api.Assertions.assertNull;
import static org.junit.jupiter.api.Assertions.assertTrue;

import java.io.File;
import java.util.Arrays;
import java.util.List;

import allconstants.Constants;
import allconstants.FragCastWeights;
import org.junit.jupiter.api.BeforeEach;
import org.junit.jupiter.api.Test;

// Two contracts, neither of which starts a process: the command line FragCast is handed for a
// fine-tune, and the reading of what it prints back. The log fixtures below are verbatim output
// from a real FragCast run, so a change to its message format fails here rather than in the field,
// where an unparsed export line would silently mean "this run produced nothing to predict with".
public class FragCastFineTuneCallerTest {
    private static final File TRAIN = new File("train.tsv");
    private static final File SAVE = new File("out.onnx");

    @BeforeEach
    public void setUp() {
        Constants.FragCast = "fragcast";
        Constants.numThreads = 8;
    }

    private static List<String> argv(String task, String baseOnnx) {
        return FragCastFineTuneCaller.buildCommand(
                FragCastFineTuneJob.forTask(task, TRAIN, baseOnnx, SAVE));
    }

    @Test
    public void everyTaskTrainsAndExports() {
        for (String task : FragCastWeights.TASKS) {
            List<String> command = argv(task, "/models/x.onnx");
            assertEquals(Arrays.asList("fragcast", "--task", task), command.subList(0, 3), task);
            assertTrue(command.containsAll(Arrays.asList("--train", TRAIN.getAbsolutePath(),
                    "--save-onnx", SAVE.getAbsolutePath())), task);
        }
    }

    // Each of these finishes successfully having written no usable model, so none may ever appear.
    @Test
    public void flagsThatSilentlyDiscardTheExportAreNeverPassed() {
        for (String task : FragCastWeights.TASKS) {
            List<String> command = argv(task, "/models/x.onnx");
            for (String flag : new String[] {"--fast", "--lora", "--grid", "--extra-blocks", "--ce"}) {
                assertFalse(command.contains(flag), task + " passed " + flag + ": " + command);
            }
        }
    }

    // FragCast derives the ion-mobility default from the retention-time path by string
    // substitution, so the weights are always named outright.
    @Test
    public void weightsAreNamedExplicitlyWhenTheyCanBeLocated() {
        List<String> command = argv("im", "/models/FragCast-IM.onnx");
        int i = command.indexOf("--onnx");
        assertTrue(i >= 0, "--onnx was omitted: " + command);
        assertEquals("/models/FragCast-IM.onnx", command.get(i + 1));

        //with nowhere to point, FragCast resolves its own bundled file rather than being handed a guess
        assertFalse(argv("im", "").contains("--onnx"));
    }

    // One library is all FragCast needs: it holds out its own --val-frac slice of --train to select
    // the best epoch against. Naming --val-frac here - and this once named it as 0, to reclaim the
    // rows for training - would leave it nothing to select on.
    @Test
    public void noevalLibraryAndNoHoldoutFractionArePassed() {
        for (String task : FragCastWeights.TASKS) {
            List<String> command = argv(task, "/models/x.onnx");
            assertFalse(command.contains("--eval"), task + " passed --eval: " + command);
            assertFalse(command.contains("--val-frac"), task + " pinned --val-frac: " + command);
        }
    }

    // Every hyperparameter belongs to FragCast. Pinning one here would freeze it at whatever was
    // true when this integration was written, and each model's tuned value moves with the executable.
    @Test
    public void nohyperparameterIsNamedOnTheCommandLine() {
        for (String task : FragCastWeights.TASKS) {
            List<String> command = argv(task, "/models/x.onnx");
            for (String flag : new String[] {"--topk", "--epochs", "--batch", "--lr", "--wd",
                    "--warmup", "--seed", "--train-blocks", "--optimizer", "--patience"}) {
                assertFalse(command.contains(flag),
                        task + " named the hyperparameter " + flag + ": " + command);
            }
        }
    }

    // Both MS2 models are fine-tunable; --fast is never passed, because FragCast refuses it
    // alongside --onnx and reads the architecture out of the weights file instead. Naming the fast
    // weights is therefore the whole of what makes a fine-tune a fast one.
    @Test
    public void thefastMs2ModelIsFineTunedByNamingItsWeights() {
        assertEquals("FragCast-Spec-Fast.onnx", FragCastFineTuneJob.weightsFileName("spec", true));
        assertEquals("FragCast-Spec.onnx", FragCastFineTuneJob.weightsFileName("spec", false));
        //the other two tasks have one model each, whichever Spec model the job predicts with
        for (boolean fast : new boolean[] {true, false}) {
            assertEquals("FragCast-RT.onnx", FragCastFineTuneJob.weightsFileName("rt", fast));
            assertEquals("FragCast-IM.onnx", FragCastFineTuneJob.weightsFileName("im", fast));
        }
        List<String> command = argv("spec", "/models/FragCast-Spec-Fast.onnx");
        assertFalse(command.contains("--fast"), "--fast was passed alongside --onnx: " + command);
        assertEquals("/models/FragCast-Spec-Fast.onnx", command.get(command.indexOf("--onnx") + 1));
    }

    // What MSBooster does still name is what only it knows.
    @Test
    public void whatMsboosterKnowsIsStillNamed() {
        List<String> command = argv("rt", "/models/x.onnx");
        assertEquals("8", command.get(command.indexOf("--threads") + 1),
                "the thread budget was left to FragCast: " + command);
        assertTrue(command.contains("--save-onnx"), command.toString());
    }

    private static final List<String> RT_LOG = Arrays.asList(
            "[2026-08-18 14:39:49 INFO] Fine-tune: task=rt params=334337 base=FragCast-RT.onnx",
            "[2026-08-18 14:39:49 INFO] No --eval given: holding out 10% of --train (303 precursors) for model selection",
            "[2026-08-18 14:39:49 INFO] Loaded train n=2697 (0.16s), eval n=303",
            "[2026-08-18 14:39:49 INFO] Baseline: Pearson=0.7800 MAE=34.688s",
            "[2026-08-18 14:39:58 INFO] Fine-tuned: Pearson=0.9200 Spearman=0.9100 MAE=12.400s RMSE=18.100s",
            "[2026-08-18 14:39:58 INFO] Wrote fine-tuned model (de-normalized to physical units, "
                    + "mu=49.4343 sd=40.0534) -> out.onnx");

    private static final List<String> SPEC_LOG = Arrays.asList(
            "[2026-08-18 14:40:26 INFO] Spec: train precursors=2684 eval precursors=302 (read 0.3s)",
            "[2026-08-18 14:40:26 INFO] Spec baseline: cosine=0.3381 spectral_angle=0.2235 (eval 0.7s, 451 prec/s)",
            "[2026-08-18 14:41:14 INFO] Spec fine-tuned: cosine=0.4083 spectral_angle=0.2724",
            "[2026-08-18 14:41:14 INFO] Wrote fine-tuned model -> out.onnx");

    @Test
    public void aRegressionRunIsReadOutOfItsLog() {
        FragCastFineTuneOutcome outcome = FragCastFineTuneOutcome.parse("rt", 0, TRAIN, RT_LOG);
        assertEquals(2697, outcome.trainCount);
        assertEquals(303, outcome.evalCount);
        assertEquals(0.78, outcome.baselineMetric, 1e-9);
        assertEquals(0.92, outcome.finalMetric, 1e-9);
        assertEquals(34.688, outcome.baselineMae, 1e-9);
        assertEquals(12.4, outcome.finalMae, 1e-9);
        assertEquals("s", outcome.maeUnit);
        assertEquals(49.4343, outcome.labelMean, 1e-9);
        assertEquals(40.0534, outcome.labelSd, 1e-9);
        assertTrue(outcome.exported);
        assertEquals("Pearson", outcome.metricName());
    }

    // The ion-mobility unit is appended to the error with no separator ("MAE=0.1491/K0"), which a
    // less careful pattern reads as 0.14910.
    @Test
    public void theIonMobilityUnitDoesNotBleedIntoTheError() {
        FragCastFineTuneOutcome outcome = FragCastFineTuneOutcome.parse("im", 0, TRAIN, Arrays.asList(
                "[2026-08-18 14:40:04 INFO] Loaded train n=2697 (0.19s), eval n=303",
                "[2026-08-18 14:40:04 INFO] Baseline: Pearson=0.8272 MAE=0.1491/K0",
                "[2026-08-18 14:40:13 INFO] Fine-tuned: Pearson=0.8694 Spearman=0.8600 "
                        + "MAE=0.1201/K0 RMSE=0.1731/K0"));
        assertEquals(0.8272, outcome.baselineMetric, 1e-9);
        assertEquals(0.8694, outcome.finalMetric, 1e-9);
        //FragCast prints the error to three decimals and then appends "1/K0" with no separator,
        //so the true values are 0.149 and 0.120 - reading 0.1491 would be reading the unit
        assertEquals(0.149, outcome.baselineMae, 1e-9);
        assertEquals(0.120, outcome.finalMae, 1e-9);
        assertEquals("1/K0", outcome.maeUnit, "the unit was read as part of the number");
        assertEquals("im: train=2697 eval=303 Pearson 0.8272 -> 0.8694 MAE 0.149 -> 0.120 1/K0",
                outcome.summary());
    }

    // The MS2 task labels its lines differently and reports cosine rather than correlation.
    @Test
    public void anMs2RunIsReadOutOfItsDifferentlyLabelledLog() {
        FragCastFineTuneOutcome outcome = FragCastFineTuneOutcome.parse("spec", 0, TRAIN, SPEC_LOG);
        assertEquals(2684, outcome.trainCount);
        assertEquals(302, outcome.evalCount);
        assertEquals(0.3381, outcome.baselineMetric, 1e-9);
        assertEquals(0.4083, outcome.finalMetric, 1e-9);
        assertTrue(outcome.exported);
        assertEquals("cosine", outcome.metricName());
        assertEquals("spec: train=2684 eval=302 cosine 0.3381 -> 0.4083", outcome.summary(),
                "MS2 reports no error, so the summary must not claim one");
    }

    @Test
    public void aGenuineImprovementIsUsable() {
        FragCastFineTuneOutcome rt = FragCastFineTuneOutcome.parse("rt", 0, existingFile(), RT_LOG);
        assertNull(rt.unusableReason(), "a real improvement was refused");
        FragCastFineTuneOutcome spec = FragCastFineTuneOutcome.parse("spec", 0, existingFile(), SPEC_LOG);
        assertNull(spec.unusableReason());
    }

    // Correlation is scale-invariant, so it can improve while the predictions drift away from the
    // measurements. The error is still read and still printed - that is the whole reason to print
    // it - but a model the user asked for and FragCast wrote is a model MSBooster predicts with;
    // how good it is, is not MSBooster's call.
    @Test
    public void abetterCorrelationWithAworseErrorIsStillUsed() {
        FragCastFineTuneOutcome outcome = FragCastFineTuneOutcome.parse("rt", 0, existingFile(), Arrays.asList(
                "[t INFO] Loaded train n=2697 (0.16s), eval n=303",
                "[t INFO] Baseline: Pearson=-0.1391 MAE=29.450s",
                "[t INFO] Fine-tuned: Pearson=-0.1008 Spearman=-0.1000 MAE=29.865s RMSE=30.000s",
                "[t INFO] Wrote fine-tuned model (de-normalized to physical units, mu=1.0 sd=1.0) -> out.onnx"));
        assertNull(outcome.unusableReason(), "a written model was refused over its error");
        //the numbers are still parsed and still reported: they are the run log, not a verdict
        assertEquals(29.450, outcome.baselineMae, 1e-9);
        assertEquals(29.865, outcome.finalMae, 1e-9);
        assertTrue(outcome.summary().endsWith("MAE 29.450 -> 29.865 s"),
                "the worsening error was left out of the summary: " + outcome.summary());
    }

    // FragCast selects the best epoch against its own holdout. If what it exports still scores below
    // the baseline, that is its judgement and the user's to read, not a reason to discard the file.
    @Test
    public void amodelWhoseMetricFellIsStillUsed() {
        FragCastFineTuneOutcome outcome = FragCastFineTuneOutcome.parse("spec", 0, existingFile(),
                Arrays.asList(
                        "[t INFO] Spec: train precursors=2684 eval precursors=302 (read 0.3s)",
                        "[t INFO] Spec baseline: cosine=0.4083 spectral_angle=0.2724 (eval 0.7s, 451 prec/s)",
                        "[t INFO] Spec fine-tuned: cosine=0.3381 spectral_angle=0.2235",
                        "[t INFO] Wrote fine-tuned model -> out.onnx"));
        assertTrue(outcome.finalMetric < outcome.baselineMetric, "the fixture no longer falls");
        assertNull(outcome.unusableReason(), "a model was refused for scoring below its baseline");
    }

    // The test above covers spec, the one task whose export states no scale and so never reaches the
    // mu/sd check. Retention time and ion mobility do reach it, so they need their own fixture: a
    // fallen metric alongside a sound standardization must still be used, or the check that survives
    // to guard the scale quietly becomes the merit gate that was removed.
    @Test
    public void anrtModelWhoseMetricFellIsStillUsedWhenItsScaleIsSound() {
        FragCastFineTuneOutcome outcome = FragCastFineTuneOutcome.parse("rt", 0, existingFile(),
                Arrays.asList(
                        "[t INFO] Loaded train n=2697 (0.16s), eval n=303",
                        "[t INFO] Baseline: Pearson=0.9200 MAE=12.400s",
                        "[t INFO] Fine-tuned: Pearson=0.7800 Spearman=0.7700 MAE=34.688s RMSE=40.100s",
                        "[t INFO] Wrote fine-tuned model (de-normalized to physical units, "
                                + "mu=49.4343 sd=40.0534) -> out.onnx"));
        assertTrue(outcome.denormalizedExport, "the fixture no longer reaches the mu/sd check");
        assertTrue(outcome.finalMetric < outcome.baselineMetric, "the fixture no longer falls");
        assertNull(outcome.unusableReason(), "an rt model was refused for scoring below its baseline");
    }

    // The metrics were once a gate of their own: a run whose numbers would not parse was refused
    // outright, which meant any drift in how FragCast labels its metric lines threw away a model it
    // had genuinely written. Only the export decides now, so unreadable numbers make for a quiet
    // summary and nothing more.
    @Test
    public void aRunWhoseMetricsCouldNotBeReadIsStillUsable() {
        FragCastFineTuneOutcome outcome = FragCastFineTuneOutcome.parse("rt", 0, existingFile(),
                Arrays.asList(
                        "[t INFO] Loaded train n=2697 (0.16s), eval n=303",
                        "[t INFO] Baseline correlation 0.7800, MAE 34.688s",
                        "[t INFO] After fine-tuning: correlation 0.9200, MAE 12.400s",
                        "[t INFO] Wrote fine-tuned model (de-normalized to physical units, "
                                + "mu=49.4343 sd=40.0534) -> out.onnx"));
        assertTrue(Double.isNaN(outcome.baselineMetric), "the fixture's metrics still parse");
        assertTrue(Double.isNaN(outcome.finalMetric), "the fixture's metrics still parse");
        assertTrue(outcome.exported);
        assertNull(outcome.unusableReason(), "a written model was refused over an unreadable metric");
    }

    // Only the de-normalized export line names the standardization, so a run that took the plain
    // export path reports mu and sd as NaN. That is a quiet log, not a corrupt model: refusing it
    // would throw away a perfectly good fine-tune over a scale FragCast never stated.
    @Test
    public void aplainExportIsUsedEvenThoughItNamesNoScale() {
        FragCastFineTuneOutcome outcome = FragCastFineTuneOutcome.parse("rt", 0, existingFile(),
                Arrays.asList(
                        "[t INFO] Loaded train n=2697 (0.16s), eval n=303",
                        "[t INFO] Baseline: Pearson=0.7800 MAE=34.688s",
                        "[t INFO] Fine-tuned: Pearson=0.9200 Spearman=0.9100 MAE=12.400s RMSE=18.100s",
                        "[t INFO] Wrote fine-tuned model -> out.onnx"));
        assertTrue(outcome.exported);
        assertFalse(outcome.denormalizedExport, "a plain export was read as a de-normalized one");
        assertTrue(Double.isNaN(outcome.labelMean));
        assertNull(outcome.unusableReason(), "a plain export was refused over a scale it never named");
    }

    // Standardizing an empty target set gives a NaN mean, which the exporter folds into the weights
    // and then reports success. Every prediction that file makes is NaN, so there is nothing to use.
    @Test
    public void aDegenerateStandardizationIsUnusable() {
        FragCastFineTuneOutcome outcome = FragCastFineTuneOutcome.parse("rt", 0, existingFile(), Arrays.asList(
                "[t INFO] Loaded train n=0 (0.01s), eval n=303",
                "[t INFO] Baseline: Pearson=0.7800 MAE=34.688s",
                "[t INFO] Fine-tuned: Pearson=0.9900 Spearman=0.9900 MAE=1.000s RMSE=1.000s",
                "[t INFO] Wrote fine-tuned model (de-normalized to physical units, mu=NaN sd=0.0000) -> out.onnx"));
        assertNotNull(outcome.unusableReason(), "a NaN-standardized model was used");
    }

    @Test
    public void aRunThatWroteNothingIsUnusable() {
        List<String> noExport = RT_LOG.subList(0, RT_LOG.size() - 1);
        FragCastFineTuneOutcome outcome = FragCastFineTuneOutcome.parse("rt", 0, existingFile(), noExport);
        assertFalse(outcome.exported);
        assertNotNull(outcome.unusableReason());

        assertNotNull(FragCastFineTuneOutcome.parse("rt", 1, existingFile(), RT_LOG)
                .unusableReason(), "a failing exit code was used");
        assertNotNull(FragCastFineTuneOutcome.parse("rt", 0, new File("nope.onnx"), RT_LOG)
                .unusableReason(), "a missing weights file was used");
    }

    // A library with nothing to train a task on is an ordinary thing - no ion mobility in the data
    // means no ion mobility in the library - and FragCast says so with its own exit code. Telling
    // that apart from a broken run by reading the log text would make the wording an interface.
    @Test
    public void alibraryWithNothingToTrainOnIsSkippedNotFailed() {
        FragCastFineTuneOutcome outcome = FragCastFineTuneOutcome.parse("im",
                FragCastFineTuneOutcome.EXIT_NO_TRAINING_DATA, existingFile(), Arrays.asList(
                        "[t WARN] lib.tsv: skipped 26468 of 32023 precursors that FragCast cannot train on",
                        "[t ERROR] --train has 0 precursor(s) FragCast can train on; at least 2 are needed"));
        assertTrue(outcome.noTrainingData(), "the no-data exit code was read as something else");
    }

    // Any other non-zero exit is a failure, and must not be quietly skipped as though the library
    // simply had nothing in it.
    @Test
    public void anordinaryFailureIsNotMistakenForAnEmptyLibrary() {
        for (int code : new int[] {1, 2, 101, -1}) {
            assertFalse(FragCastFineTuneOutcome.parse("rt", code, existingFile(), RT_LOG).noTrainingData(),
                    "exit code " + code + " was read as an empty library");
        }
    }

    // FragCast decides whether a training set is big enough to learn from, by selecting the best
    // epoch against its own held-out slice. A size floor here would override that with a number
    // MSBooster invented, so even a handful of precursors produces a model like any other.
    @Test
    public void asmallTrainingSetIsUsedWhateverItsSize() {
        FragCastFineTuneOutcome outcome = FragCastFineTuneOutcome.parse("rt", 0, existingFile(),
                Arrays.asList(
                        "[t INFO] Loaded train n=7 (0.01s), eval n=1",
                        "[t INFO] Baseline: Pearson=0.7800 MAE=34.688s",
                        "[t INFO] Fine-tuned: Pearson=0.7830 Spearman=0.7800 MAE=34.100s RMSE=40.100s",
                        "[t INFO] Wrote fine-tuned model (de-normalized to physical units, "
                                + "mu=49.4343 sd=40.0534) -> out.onnx"));
        assertEquals(7, outcome.trainCount);
        assertNull(outcome.unusableReason(), "a 7-precursor run was refused on size");
    }

    // The held-out slice is FragCast's own, and an empty one means it selected the best epoch
    // against nothing. The count is worth reading in the log, but the weights it wrote are still
    // the weights the user asked for.
    @Test
    public void amodelScoredAgainstAnEmptyHoldoutIsStillUsed() {
        FragCastFineTuneOutcome outcome = FragCastFineTuneOutcome.parse("rt", 0, existingFile(),
                Arrays.asList(
                        "[t INFO] Loaded train n=2697 (0.16s), eval n=0",
                        "[t INFO] Baseline: Pearson=0.7800 MAE=34.688s",
                        "[t INFO] Fine-tuned: Pearson=0.9200 Spearman=0.9100 MAE=12.400s RMSE=18.100s",
                        "[t INFO] Wrote fine-tuned model (de-normalized to physical units, "
                                + "mu=49.4343 sd=40.0534) -> out.onnx"));
        assertEquals(0, outcome.evalCount);
        assertNull(outcome.unusableReason(), "a model was refused for the size of FragCast's holdout");
    }

    /** A real, non-empty file standing in for exported weights. */
    private static File existingFile() {
        return new File("pom.xml");
    }
}
