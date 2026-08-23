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

package allconstants;

import java.util.List;

import utils.Model;

/**
 * The FragCast model names MSBooster exposes, and what separates them.
 *
 * <p>One FragCast executable ships two MS2 (Spec) predictors: the default Conformer
 * ({@code FragCast-Spec.onnx}) and a smaller, ~5x cheaper cnn_bigru model
 * ({@code FragCast-Spec-Fast.onnx}) selected with the {@code --fast} flag. FragCast has no
 * per-property invocation to pick between them — a single {@code --task build-library} run predicts
 * RT, IM and MS2 together — so the fast predictor is exposed here as its own model name, and the
 * {@code --fast} flag is derived from whichever name the user asked for.
 *
 * <p>Only the Spec weights differ: RT and IM predictions are identical either way. {@link #FAST} is
 * therefore a spectra choice, and is deliberately absent from the RT/IM findBest candidate lists in
 * {@link ModelCollections}, where it would cost a second FragCast run to reproduce numbers
 * {@link #CONFORMER} already produced. For the same reason a job that names both variants
 * (fast spectra, Conformer RT/IM) still runs FragCast once: {@link #jobUsesFastSpecModel} picks the
 * job's single set of Spec weights and {@link #predictionKey} makes both names share its library.
 */
public class FragCastModels {
    /** Default FragCast: Conformer Spec weights, the most accurate of the two. */
    public static final String CONFORMER = "FragCast";
    /**
     * Small/fast Spec weights — ~5x cheaper MS2 forward, ~0.01-0.02 lower cosine. The name is the
     * model a user selects, not the weights file it loads ({@code FragCast-Spec-Fast.onnx}).
     */
    public static final String FAST = "FragCast-Fast";

    /** Is this model name served by the FragCast executable, under either set of Spec weights? */
    public static boolean isFragCast(String model) {
        return CONFORMER.equals(model) || FAST.equals(model);
    }

    /** Does this model name require FragCast's {@code --fast} flag? */
    public static boolean usesFastSpecModel(String model) {
        return FAST.equals(model);
    }

    /** The cache key both Spec variants share, so one job asks the executable for one library. */
    private static final String EXECUTABLE_KEY = "FragCast-executable";

    /**
     * The key a job caches a model's prediction file under. Both Spec variants collapse to one key,
     * so a job that names both (fast spectra with Conformer RT/IM) runs FragCast once and reads
     * every property out of that run's library — legitimate because RT and IM are identical under
     * either set of Spec weights, so the second run could only reproduce them. Every other model,
     * DIA-NN included, keys on its own name and keeps its own run.
     */
    public static String predictionKey(String model) {
        return isFragCast(model) ? EXECUTABLE_KEY : model;
    }

    /**
     * The key a job caches a model's prediction <i>file</i> under once weights can differ.
     *
     * <p>{@link #predictionKey(String)} answers "which executable produces this model's library",
     * which is the right question for the peptide input file — that file is the same whatever
     * weights are loaded. It is the wrong question for the library itself, whose <i>name</i> carries
     * {@link FragCastWeights#fileTag()}: a key that left the weights out would stand for two
     * different files, so a model loaded under one weight set could be served a library predicted
     * under another. Folding the same tag in keeps the key and the file name saying the same thing,
     * and leaves the pretrained key byte-identical to what it always was.
     */
    public static String predictionKey(String model, FragCastWeights weights) {
        String key = predictionKey(model);
        return isFragCast(model) ? key + weights.fileTag() : key;
    }

    /**
     * Which Spec weights the job's single FragCast run should load, given the models it will run.
     *
     * <p>Only the MS2 intensities depend on the choice, so the spectra model makes it. When no
     * FragCast model serves spectra the choice cannot change a number, and the first FragCast model
     * named decides — never loading {@code FragCast-Spec-Fast.onnx} for a job that did not ask for
     * it. Returns {@code false} when no model here is a FragCast one.
     */
    public static boolean jobUsesFastSpecModel(List<Model> models) {
        String firstFragCast = null;
        for (Model model : models) {
            if (!isFragCast(model.name)) {
                continue;
            }
            if ("spectra".equals(model.modelType)) {
                return usesFastSpecModel(model.name);
            }
            if (firstFragCast == null) {
                firstFragCast = model.name;
            }
        }
        return usesFastSpecModel(firstFragCast);
    }
}
