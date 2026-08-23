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
import static org.junit.jupiter.api.Assertions.assertTrue;

import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;

import allconstants.FragCastWeights;
import org.junit.jupiter.api.Test;
import org.junit.jupiter.api.io.TempDir;

// A fine-tune produces loose ONNX files and a bundle holding copies of them. What may be deleted
// afterwards is the whole of what this class decides, and getting it wrong removes weights the
// installation - or a model the user pointed a parameter at - still needs.
public class FragCastTransferLearnerTest {

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
}
