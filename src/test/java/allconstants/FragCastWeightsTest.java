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

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertFalse;
import static org.junit.jupiter.api.Assertions.assertNotEquals;
import static org.junit.jupiter.api.Assertions.assertThrows;
import static org.junit.jupiter.api.Assertions.assertTrue;

import java.util.ArrayList;
import java.util.List;

import org.junit.jupiter.api.Test;

// FragCastWeights is what keeps a fine-tuned FragCast run from being confused with the pretrained
// one, so the two properties that matter are that pretrained weights change nothing at all, and
// that any other weight set is distinguishable from them and from each other.
public class FragCastWeightsTest {
    // The backward-compatibility pin: a job that fine-tunes nothing must still write
    // spectraRT.predicted.parquet and pass FragCast no weight flags.
    @Test
    public void pretrainedWeightsAreInvisible() {
        FragCastWeights base = FragCastWeights.base();
        assertTrue(base.isBase());
        assertEquals("", base.fileTag(), "the pretrained run's file name changed");

        List<String> command = new ArrayList<>();
        base.appendFlags(command);
        assertTrue(command.isEmpty(), "pretrained weights added flags: " + command);
    }

    // A parameter file may leave a path unset, blank, or at the literal "null" the parser uses to
    // mean "keep the default"; all three mean the pretrained model.
    @Test
    public void unsetPathsAllMeanPretrained() {
        assertTrue(FragCastWeights.of(null, "", "null").isBase());
        assertTrue(FragCastWeights.of("  ", null, null).isBase());
        assertEquals(FragCastWeights.base(), FragCastWeights.of("", "", ""));
    }

    @Test
    public void customWeightsGetADistinctStableTag() {
        FragCastWeights rt = FragCastWeights.of("/models/rt.onnx", "", "");
        assertFalse(rt.isBase());
        assertEquals(rt.fileTag(), FragCastWeights.of("/models/rt.onnx", "", "").fileTag(),
                "the same weights produced two different file names");
        assertTrue(rt.fileTag().startsWith(".ft-"), rt.fileTag());

        FragCastWeights im = FragCastWeights.of("", "/models/rt.onnx", "");
        assertNotEquals(rt.fileTag(), im.fileTag(),
                "the same path on a different property produced the same file name");
    }

    // Each property is swapped independently, because a fine-tune may only have models worth
    // keeping for some of them.
    @Test
    public void propertiesAreSwappedOneAtATime() {
        FragCastWeights weights = FragCastWeights.base().withRt("/models/rt.onnx").withSpec("/models/spec.onnx");
        assertEquals("/models/rt.onnx", weights.rtOnnx);
        assertEquals("", weights.imOnnx, "the ion-mobility model was changed unasked");
        assertEquals("/models/spec.onnx", weights.specOnnx);

        List<String> command = new ArrayList<>();
        weights.appendFlags(command);
        assertEquals(List.of("--rt-onnx", "/models/rt.onnx", "--spec-onnx", "/models/spec.onnx"), command);
    }

    // FragCast refuses --fast alongside --spec-onnx, because both name the Spec weights. The
    // weights win: the ONNX declares its own architecture, so a fine-tuned fast model still runs
    // fast when named this way - which is the only way to predict with one.
    @Test
    public void namedMs2WeightsSupersedeTheFastFlag() {
        assertTrue(FragCastWeights.base().withSpec("/models/spec.onnx").supersedesFastFlag(true));

        //without the fast model there is no flag to supersede
        assertFalse(FragCastWeights.base().withSpec("/models/spec.onnx").supersedesFastFlag(false));
        //and RT weights say nothing about which Spec model runs
        assertFalse(FragCastWeights.base().withRt("/models/rt.onnx").supersedesFastFlag(true));
        assertFalse(FragCastWeights.base().supersedesFastFlag(true));
    }

    // Both FragCast Spec variants share one prediction key by name, because they read the same
    // library; different weights write different libraries, so those must not share a key.
    @Test
    public void predictionKeySeparatesWeightsButNotSpecVariants() {
        FragCastWeights base = FragCastWeights.base();
        FragCastWeights tuned = base.withRt("/models/rt.onnx");

        assertEquals(FragCastModels.predictionKey(FragCastModels.CONFORMER, base),
                FragCastModels.predictionKey(FragCastModels.FAST, base));
        assertNotEquals(FragCastModels.predictionKey(FragCastModels.CONFORMER, base),
                FragCastModels.predictionKey(FragCastModels.CONFORMER, tuned),
                "a fine-tuned pass would have been served the pretrained library");
        assertEquals(FragCastModels.predictionKey(FragCastModels.CONFORMER),
                FragCastModels.predictionKey(FragCastModels.CONFORMER, base),
                "the pretrained key stopped matching the one-argument form");
        assertEquals("DIA-NN", FragCastModels.predictionKey("DIA-NN", tuned),
                "FragCast weights leaked into another predictor's key");
    }
}
