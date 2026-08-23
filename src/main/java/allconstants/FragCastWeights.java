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
import java.util.Objects;

/**
 * Which FragCast weights a prediction run loads: one ONNX path per property, each empty when that
 * property should use the pretrained model in {@link Constants#FragCastModelDir}.
 *
 * <p>FragCast selects weights per property on its {@code build-library} command line
 * ({@code --rt-onnx}, {@code --im-onnx}, {@code --spec-onnx}), so a job may mix a fine-tuned RT
 * model with pretrained IM and MS2 models. That is exactly what a fine-tune produces when only some
 * tasks exported a model to load, so the three paths travel together as one immutable value rather
 * than as three loose strings.
 *
 * <p>The value also names the library a run produces. Fine-tuning is its own run and rescoring
 * with the result is another, so one output directory holds both, and {@link #fileTag()} is what
 * keeps them apart: a library predicted from a fine-tuned model never shares a path with the
 * pretrained one, and neither can be mistaken for or overwritten by the other.
 * {@link FragCastModels#predictionKey(String, FragCastWeights)} folds the same tag into the key a
 * job caches that file under, so the two always agree on what identifies a library. The pretrained
 * set tags nothing, so a job naming no custom weights writes exactly the names it always did.
 */
public final class FragCastWeights {
    /**
     * The properties a weight set covers, in the order a run works through them. Every task-keyed
     * loop in the codebase iterates this rather than writing the three names out again.
     */
    public static final String[] TASKS = {"rt", "im", "spec"};

    private static final FragCastWeights BASE = new FragCastWeights("", "", "");

    public final String rtOnnx;
    public final String imOnnx;
    public final String specOnnx;

    private FragCastWeights(String rtOnnx, String imOnnx, String specOnnx) {
        this.rtOnnx = normalize(rtOnnx);
        this.imOnnx = normalize(imOnnx);
        this.specOnnx = normalize(specOnnx);
    }

    //a param left at "null" in a params file arrives as the literal 4-character string, and an
    //unset one arrives as null; both mean "use the pretrained weights"
    private static String normalize(String path) {
        if (path == null) {
            return "";
        }
        String trimmed = path.trim();
        return trimmed.equals("null") ? "" : trimmed;
    }

    /** The pretrained weights: no {@code --rt-onnx}/{@code --im-onnx}/{@code --spec-onnx} at all. */
    public static FragCastWeights base() {
        return BASE;
    }

    /** The weights the params file asks for. */
    public static FragCastWeights fromConstants() {
        return of(Constants.FragCastRtOnnx, Constants.FragCastImOnnx, Constants.FragCastSpecOnnx);
    }

    public static FragCastWeights of(String rtOnnx, String imOnnx, String specOnnx) {
        FragCastWeights weights = new FragCastWeights(rtOnnx, imOnnx, specOnnx);
        return weights.isBase() ? BASE : weights;
    }

    public FragCastWeights withRt(String path) {
        return of(path, imOnnx, specOnnx);
    }

    public FragCastWeights withIm(String path) {
        return of(rtOnnx, path, specOnnx);
    }

    public FragCastWeights withSpec(String path) {
        return of(rtOnnx, imOnnx, path);
    }

    /** The path this set names for a task, or empty when that property uses the pretrained model. */
    public String forTask(String task) {
        switch (task) {
            case "rt":
                return rtOnnx;
            case "im":
                return imOnnx;
            case "spec":
                return specOnnx;
            default:
                return "";
        }
    }

    /** This set with one task substituted. Unknown tasks change nothing. */
    public FragCastWeights withTask(String task, String path) {
        switch (task) {
            case "rt":
                return withRt(path);
            case "im":
                return withIm(path);
            case "spec":
                return withSpec(path);
            default:
                return this;
        }
    }

    /** Does this name any custom weights at all? */
    public boolean isBase() {
        return rtOnnx.isEmpty() && imOnnx.isEmpty() && specOnnx.isEmpty();
    }

    /**
     * The marker a prediction file carries so libraries written under different weights never share
     * a path. Empty for {@link #base()}, so the pretrained run keeps writing the file name it
     * always has; otherwise {@code .ft-} plus eight stable hex digits of the three paths.
     */
    public String fileTag() {
        if (isBase()) {
            return "";
        }
        //String.hashCode is specified exactly, so the same weights give the same tag on any JVM
        int hash = (rtOnnx + '\u0000' + imOnnx + '\u0000' + specOnnx).hashCode();
        return ".ft-" + String.format("%08x", hash);
    }

    /**
     * Reject a weight set FragCast itself would reject, so the failure names the MSBooster parameter
     * rather than arriving as a Rust error. FragCast treats {@code --fast} and {@code --spec-onnx}
     * as a hard conflict: both name the Spec weights, and letting either win would silently ignore
     * one of them.
     *
     * @param fast whether the job runs FragCast's small/fast Spec model
     *             ({@link FragCastModels#jobUsesFastSpecModel})
     */
    public boolean supersedesFastFlag(boolean fast) {
        return fast && !specOnnx.isEmpty();
    }

    /**
     * Append the weight flags this set names to a FragCast {@code build-library} command line. A
     * property with no custom weights contributes nothing, leaving FragCast to resolve its own
     * pretrained file from {@code FRAGCAST_MODEL_DIR}.
     */
    public void appendFlags(List<String> command) {
        if (!rtOnnx.isEmpty()) {
            command.add("--rt-onnx");
            command.add(rtOnnx);
        }
        if (!imOnnx.isEmpty()) {
            command.add("--im-onnx");
            command.add(imOnnx);
        }
        if (!specOnnx.isEmpty()) {
            command.add("--spec-onnx");
            command.add(specOnnx);
        }
    }

    @Override
    public boolean equals(Object o) {
        if (this == o) {
            return true;
        }
        if (!(o instanceof FragCastWeights)) {
            return false;
        }
        FragCastWeights other = (FragCastWeights) o;
        return rtOnnx.equals(other.rtOnnx) && imOnnx.equals(other.imOnnx) && specOnnx.equals(other.specOnnx);
    }

    @Override
    public int hashCode() {
        return Objects.hash(rtOnnx, imOnnx, specOnnx);
    }

    @Override
    public String toString() {
        return isBase() ? "pretrained FragCast weights"
                : "FragCast weights [rt=" + describe(rtOnnx) + ", im=" + describe(imOnnx) +
                        ", spec=" + describe(specOnnx) + "]";
    }

    private static String describe(String path) {
        return path.isEmpty() ? "pretrained" : path;
    }
}
