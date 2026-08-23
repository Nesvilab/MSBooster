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

import java.io.File;

/**
 * One FragCast fine-tuning run, fully described before anything is executed.
 *
 * <p>The hyperparameters are FragCast's own. MSBooster names no epoch count, learning rate, batch
 * size, weight decay, warmup, seed, training depth or {@code --topk}: it omits every one of those
 * flags, so each run trains on whatever the executable considers right for the task. FragCast is the
 * component that knows what its models need, and the tuned values move with it rather than being
 * pinned to whatever was true when this integration was written.
 *
 * <p>What MSBooster still names is what only it knows: which task to run, which library to read,
 * which weights to start from, where to write the result, and how many threads it may use.
 *
 * <p>There is no eval library either. FragCast holds out its own {@code --val-frac} slice of
 * {@code --train} and selects the best epoch against it, so supplying a second library would only
 * duplicate a decision it already makes - and would require choosing how much to hold out.
 *
 * <p>The one hyperparameter this used to set deliberately was {@code --warmup}, because FragCast's
 * default of 100 optimizer steps can exceed the entire step count of a small fine-tune, leaving the
 * learning rate on its ramp for the whole run. FragCast now selects the best epoch against its own
 * held-out slice, so which epoch a short run ends on is a decision it makes with the ramp it chose
 * rather than one MSBooster forces by pinning a warmup it cannot tune per task.
 */
public class FragCastFineTuneJob {
    public final String task;
    public final File train;
    /** The weights to start from, or empty to let FragCast resolve its own pretrained file. */
    public final String baseOnnx;
    public final File saveOnnx;
    public final int threads;

    FragCastFineTuneJob(String task, File train, String baseOnnx, File saveOnnx, int threads) {
        this.task = task;
        this.train = train;
        this.baseOnnx = baseOnnx == null ? "" : baseOnnx;
        this.saveOnnx = saveOnnx;
        this.threads = threads;
    }

    /**
     * The bundled weights file name for a task, used to build an explicit {@code --onnx} path.
     *
     * <p>The MS2 task has two: naming the fast file is how a fast job is fine-tuned, since FragCast
     * refuses {@code --fast} alongside {@code --onnx} and reads the architecture out of the weights
     * instead.
     */
    public static String weightsFileName(String task, boolean fast) {
        switch (task) {
            case "im":
                return "FragCast-IM.onnx";
            case "spec":
                return fast ? "FragCast-Spec-Fast.onnx" : "FragCast-Spec.onnx";
            default:
                return "FragCast-RT.onnx";
        }
    }

    /**
     * Build the job for one task.
     *
     * @param task            {@code rt}, {@code im} or {@code spec}
     * @param train           the training library
     * @param baseOnnx        weights to start from, or empty for FragCast's own default
     * @param saveOnnx        where the fine-tuned weights are written
     */
    public static FragCastFineTuneJob forTask(String task, File train, String baseOnnx,
                                              File saveOnnx) {
        return new FragCastFineTuneJob(task, train, baseOnnx, saveOnnx,
                Math.max(1, Constants.numThreads));
    }
}
