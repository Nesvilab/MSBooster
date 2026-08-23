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

import static org.junit.jupiter.api.Assertions.assertFalse;
import static org.junit.jupiter.api.Assertions.assertTrue;

import java.io.File;
import java.lang.reflect.Field;
import java.lang.reflect.Modifier;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashSet;
import java.util.List;

import allconstants.Constants;
import allconstants.FragmentIonConstants;
import allconstants.ModelCollections;
import allconstants.NceConstants;
import org.junit.jupiter.api.Test;

// An unrecognized key in a parameter file is a hard error that ends the run, so a typo in the
// shipped template is not a documentation bug - it is a file that cannot be used. This checks the
// template against the fields the parser actually resolves, without mutating any of them.
public class ParameterTemplateTest {
    private static final File TEMPLATE = new File("msbooster_params.txt");

    /** Every name updateConstants will accept: public, non-final fields of the four holders. */
    private static HashSet<String> settableParameters() {
        HashSet<String> names = new HashSet<>();
        for (Class<?> holder : Arrays.asList(Constants.class, FragmentIonConstants.class,
                NceConstants.class, ModelCollections.class)) {
            for (Field field : holder.getFields()) {
                if (!Modifier.isFinal(field.getModifiers())) {
                    names.add(field.getName());
                }
            }
        }
        return names;
    }

    /** The keys the template sets, read the way processParameterList reads them. */
    private static List<String> templateKeys() throws Exception {
        ArrayList<String> keys = new ArrayList<>();
        for (String line : Files.readAllLines(TEMPLATE.toPath(), StandardCharsets.UTF_8)) {
            String stripped = line.split("#", 2)[0];
            if (!stripped.contains("=")) {
                continue;
            }
            keys.add(stripped.split("=", 2)[0].trim());
        }
        return keys;
    }

    @Test
    public void everyKeyInTheShippedTemplateIsARealParameter() throws Exception {
        assertTrue(TEMPLATE.isFile(), "the parameter template moved: " + TEMPLATE.getAbsolutePath());
        HashSet<String> settable = settableParameters();
        ArrayList<String> unknown = new ArrayList<>();
        for (String key : templateKeys()) {
            if (!settable.contains(key)) {
                unknown.add(key);
            }
        }
        assertTrue(unknown.isEmpty(), "msbooster_params.txt sets parameters that do not exist, so " +
                "MSBooster would exit on it: " + unknown);
    }

    // The template is how a user discovers a feature; a parameter absent from it is one nobody finds.
    @Test
    public void theTemplateDocumentsTheLocalTransferLearningParameters() throws Exception {
        List<String> keys = templateKeys();
        for (String key : Arrays.asList("FragCastRtOnnx", "FragCastImOnnx", "FragCastSpecOnnx",
                "FragCastModelZip")) {
            assertTrue(keys.contains(key), key + " is missing from msbooster_params.txt");
        }
    }

    // Nothing about how the fine-tune is run belongs to MSBooster. Pinning a hyperparameter here
    // would freeze it at whatever was true when this integration was written, and choosing a
    // threshold or a holdout fraction would override a decision FragCast already makes against its
    // own held-out slice. A parameter reintroduced by accident would be silently authoritative.
    @Test
    public void nothingAboutHowTheFineTuneRunsIsExposed() throws Exception {
        HashSet<String> settable = settableParameters();
        ArrayList<String> exposed = new ArrayList<>();
        for (String key : Arrays.asList("fragCastTransferRtEpochs", "fragCastTransferImEpochs",
                "fragCastTransferSpecEpochs", "fragCastTransferRtLr", "fragCastTransferImLr",
                "fragCastTransferSpecLr", "fragCastTransferBatch", "fragCastTransferWeightDecay",
                "fragCastTransferWarmup", "fragCastTransferSeed", "fragCastTransferMaxTrain",
                "fragCastTransferRtTrainBlocks", "fragCastTransferImTrainBlocks",
                "fragCastTransferMinTrain", "fragCastTransferEvalFraction",
                "fragCastTransferEvalLibrary", "fragCastTransferMinImprovement",
                "fragCastTransferProbeSize", "fragCastTransferMaxScaleDrift")) {
            if (settable.contains(key)) {
                exposed.add(key);
            }
        }
        assertTrue(exposed.isEmpty(), "these belong to FragCast, not to MSBooster: " + exposed);
    }

}
