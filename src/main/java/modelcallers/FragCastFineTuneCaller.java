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

import allconstants.Constants;

import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

import static utils.Print.printInfo;

/**
 * Runs FragCast's fine-tuning tasks - the local, CPU-only counterpart to the server-based transfer
 * learning in {@code transferlearn}.
 *
 * <p>The same executable that predicts a spectral library also trains on one: {@code --task rt},
 * {@code --task im} and {@code --task spec} each read one training library, hold out a slice of it
 * to select the best epoch against, report the model's accuracy before and after, and write the
 * fine-tuned weights to {@code --save-onnx} as a patched copy of the model they started from. Nothing leaves the machine, and the result drops
 * straight back into {@code --task build-library} through {@code --rt-onnx}, {@code --im-onnx} or
 * {@code --spec-onnx}.
 *
 * <p>{@link #buildCommand} is the seam, in the same spirit as
 * {@link FragCastModelCaller#buildCommand}: the whole command-line contract is exercisable without
 * an installed executable. The command it builds is deliberately short: MSBooster names the task,
 * the library, the weights and the thread count, and leaves every hyperparameter - including how
 * much of the library to hold out - to FragCast's own per-task defaults.
 *
 * <p>Four flags are not merely omitted but must never be added, each for a reason FragCast documents
 * or demonstrates - {@code --fast} names an inference-only model and is refused alongside
 * {@code --save-onnx}; {@code --lora} and {@code --grid} train and then return before the export ever
 * runs, so they finish successfully having written no weights; {@code --extra-blocks} appends layers
 * that have no representation in the ONNX file and skips the export; and {@code --ce} is a no-op in
 * the production MS2 model.
 */
public class FragCastFineTuneCaller {
    /** The FragCast command line for one fine-tuning run. */
    static List<String> buildCommand(FragCastFineTuneJob job) {
        List<String> command = new ArrayList<>();
        command.add(Constants.FragCast);
        command.add("--task");
        command.add(job.task);
        command.add("--train");
        command.add(job.train.getAbsolutePath());
        if (!job.baseOnnx.isEmpty()) {
            //naming the weights explicitly also sidesteps FragCast's default for the ion-mobility
            //task, which is the retention-time path with its file name rewritten by substring
            command.add("--onnx");
            command.add(job.baseOnnx);
        }
        command.add("--save-onnx");
        command.add(job.saveOnnx.getAbsolutePath());
        command.add("--threads");
        command.add(String.valueOf(job.threads));
        //--val-frac is deliberately absent along with the rest. Without --eval it is what FragCast
        //holds out of --train to select the best epoch against, so pinning it - and pinning it to
        //zero, as this once did to reclaim training rows - would leave nothing to select on
        return command;
    }

    /**
     * Run one fine-tuning job to completion and parse what it reported.
     *
     * <p>Nothing here exits the JVM. A failed fine-tune is a reason to keep predicting with the
     * pretrained model, not a reason to abandon a rescoring run that is otherwise fine.
     */
    public static FragCastFineTuneOutcome fineTune(FragCastFineTuneJob job, boolean verbose)
            throws IOException, InterruptedException {
        List<String> command = buildCommand(job);
        printInfo("Fine-tuning the FragCast " + job.task + " model on " + job.train.getName() +
                ", with FragCast's own hyperparameters for the task. FragCast reports how many of " +
                "its precursors it could train on.");
        FragCastProcess.Result result = FragCastProcess.run(command, verbose);
        return FragCastFineTuneOutcome.parse(job.task, result.exitCode, job.saveOnnx, result.output);
    }

}
