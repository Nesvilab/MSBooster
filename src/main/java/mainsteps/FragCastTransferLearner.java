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

import allconstants.FragCastWeights;
import modelcallers.FragCastFineTuneCaller;
import modelcallers.FragCastFineTuneJob;
import modelcallers.FragCastFineTuneOutcome;
import modelcallers.FragCastProcess;
import transferlearn.fragcast.FragCastModelBundle;

import java.io.File;
import java.io.IOException;

import static utils.Print.printError;
import static utils.Print.printInfo;

/**
 * Fine-tunes FragCast on a spectral library the user supplies, locally.
 *
 * <p>This is the second transfer-learning option MSBooster offers. The first, in
 * {@code transferlearn}, uploads a spectral library to a training server, waits in its queue, and
 * downloads a model. This one needs no server, no account and no network: FragCast trains on the
 * CPU of the machine it is running on.
 *
 * <p>The work of one fine-tune is here; who asks for it is {@link
 * transferlearn.FragCastLocalWorkflow}. Each task FragCast can serve is trained separately, and what
 * it exports is gathered into one zip. Predicting with those models is a separate run, told to load
 * that zip - so this class produces a model and nothing else.
 *
 * <p>Every model FragCast exports is used. Whether a fine-tune was worth doing is FragCast's
 * judgement and the user's, not MSBooster's: FragCast selects the best epoch against its own
 * held-out slice, so nothing here sets a training-set floor, an improvement threshold or a holdout
 * fraction - each of those would override a decision already made, with a number MSBooster invented,
 * and would leave the user asking for a fine-tuned model and silently getting the pretrained one.
 * The metrics FragCast reports are printed for the run log; the only thing checked is that there is
 * a weights file to name.
 *
 * <p>The supplied library's targets are taken as they stand, which puts one thing in the user's
 * hands: {@code NormalizedRetentionTime} is trained on directly, so a library holding gradient
 * minutes teaches the model to predict that gradient's minutes rather than an index meaning the same
 * thing in every run. Nothing downstream will catch that; the library is where it has to be right.
 *
 * <p>Which of a library's entries can be trained on at all is FragCast's to decide and to report:
 * it refuses a peptide it cannot represent as written rather than quietly mangling it, and says how
 * many it refused.
 */
public class FragCastTransferLearner {
    /** Stem of the one-file model a fine-tune writes, without the {@code .zip}. */
    public static final String DEFAULT_BUNDLE_BASENAME = "FragCast-finetuned";

    /**
     * Run every task FragCast can be fine-tuned for, and take every model it exported. There is no
     * subset to choose: a property the library carries no measurement for is FragCast's to skip and
     * to report, which it does per task, so naming one here would only be a second opinion on a
     * question already answered.
     *
     * @param trainLibrary the library itself, passed to FragCast unchanged
     * @param fast start the MS2 task from FragCast's small/fast Spec weights rather than the
     *             Conformer ones. Both are fine-tunable; this only decides which file it patches.
     * @return {@code before} with each fine-tuned model substituted in
     */
    public static FragCastWeights fineTuneAll(FragCastWeights before, File trainLibrary, File dir,
                                              boolean fast) {
        String modelDir = FragCastProcess.resolveModelDir();
        FragCastWeights accepted = before;
        for (String task : FragCastWeights.TASKS) {
            File saveOnnx = new File(dir, "FragCast-" + task.toUpperCase() + ".finetuned.onnx");
            String startFrom = startingWeights(task, before, modelDir, fast);
            if (sameFile(startFrom, saveOnnx)) {
                //FragCast writes the export as a patched copy of the model it started from, so a run
                //pointed at its own output would overwrite the user's weights with an unvalidated
                //model - and the path would be unchanged, so nothing downstream could tell
                printError("Skipping the FragCast " + task + " model: it would be fine-tuned from " +
                        startFrom + " and written back over it. Point the fine-tune somewhere else " +
                        "with --output-dir, or unset the parameter naming that model.");
                continue;
            }
            FragCastFineTuneJob job = FragCastFineTuneJob.forTask(task, trainLibrary,
                    startFrom, saveOnnx);
            try {
                //echoed rather than captured: FragCast logs a line per epoch, and fine-tuning is
                //long enough that a silent one cannot be told from a stalled one
                FragCastFineTuneOutcome outcome = FragCastFineTuneCaller.fineTune(job, true);
                if (outcome.noTrainingData()) {
                    //the library holds nothing this model can learn from - no ion mobility in it,
                    //most often. Not a failure: there is simply no fine-tuned model to be had.
                    printInfo("  skipped: the library carries nothing the " + task +
                            " model can be trained on");
                    continue;
                }
                printInfo("  " + outcome.summary());
                String reason = outcome.unusableReason();
                if (reason != null) {
                    //not a verdict on the fine-tune: there is no weights file to point the next
                    //prediction pass at, so this run has nothing new to predict this property with
                    printInfo("  no usable fine-tuned " + task + " model was written: " + reason);
                    continue;
                }
                accepted = accepted.withTask(task, saveOnnx.getAbsolutePath());
            } catch (Exception e) {
                printError("  fine-tuning the FragCast " + task + " model failed (" + e +
                        "); this run has no fine-tuned " + task + " model to predict with");
            }
        }
        return accepted;
    }

    /**
     * The weights a task starts from: whatever this job predicts with, which is the pretrained model
     * unless the user pointed a parameter at one of their own. Naming the file explicitly matters -
     * FragCast's own default for the ion-mobility task is derived from the retention-time path by
     * string substitution.
     */
    private static String startingWeights(String task, FragCastWeights before, String modelDir,
                                          boolean fast) {
        String configured = before.forTask(task);
        if (!configured.isEmpty()) {
            return configured;
        }
        if (modelDir == null || modelDir.isEmpty()) {
            return ""; //let FragCast resolve its own bundled weights
        }
        return new File(modelDir, FragCastFineTuneJob.weightsFileName(task, fast)).getPath();
    }

    /** Do these two paths name the same file on disk, whatever they look like as strings? */
    private static boolean sameFile(String path, File file) {
        if (path == null || path.isEmpty()) {
            return false;
        }
        try {
            return new File(path).getCanonicalFile().equals(file.getCanonicalFile());
        } catch (IOException e) {
            return new File(path).getAbsoluteFile().equals(file.getAbsoluteFile());
        }
    }

    /**
     * Gather the accepted models into one zip beside them, and say where it is.
     *
     * <p>What this buys is that a fine-tune produces one thing rather than three. What a user
     * copies to another machine, or hands to the MSBooster panel in FragPipe, is one file - and the
     * loose ONNX files that went into it are removed once it exists ({@link #deleteLooseModels}).
     *
     * <p>Failing to write it is reported and not raised: the loose models are still on disk, and the
     * caller names those instead.
     *
     * @return the bundle, or {@code null} if there was nothing to bundle or it could not be written
     */
    public static File writeBundle(FragCastWeights accepted, File dir, String basename) {
        try {
            return FragCastModelBundle.write(accepted, new File(dir, basename + ".zip"));
        } catch (IOException e) {
            printError("Could not write the FragCast model bundle (" + e.getMessage() +
                    "); the fine-tuned models are still there individually");
            return null;
        }
    }

    /**
     * Remove the loose ONNX files a fine-tune wrote, now that a bundle holds them.
     *
     * <p>Only ever called once a bundle has actually been written. Without one the loose files are
     * the only copy of what the fine-tune produced, so the caller names those instead.
     *
     * <p>Only files this run produced are removed. A task that was skipped still names whatever it
     * started from - the shipped pretrained weights, or a model the user pointed a parameter at -
     * and those are not this run's to delete, which is what comparing against {@code before} decides.
     *
     * @return how many were removed
     */
    public static int deleteLooseModels(FragCastWeights accepted, FragCastWeights before) {
        int removed = 0;
        for (String task : FragCastWeights.TASKS) {
            String path = accepted.forTask(task);
            if (path.isEmpty() || path.equals(before.forTask(task))) {
                continue;
            }
            File file = new File(path);
            if (!file.isFile()) {
                continue;
            }
            if (file.delete()) {
                removed++;
            } else {
                printError("Could not remove " + file.getPath() + "; it is a copy of what the " +
                        "bundle already holds, and can be deleted by hand");
            }
        }
        return removed;
    }
}
