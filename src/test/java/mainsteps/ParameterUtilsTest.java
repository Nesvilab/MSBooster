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

import java.io.File;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.HashMap;

import org.junit.jupiter.api.Test;
import org.junit.jupiter.api.io.TempDir;

// Pins how the decoy tag reaches Constants.decoyPrefix. Every consumer compares it against bare
// protein labels (pin proteins, library ProteinId, FragCast's --decoy-tag, the prediction server's
// decoy_tag), so any path that stores something other than the bare tag silently unmarks every
// decoy downstream.
public class ParameterUtilsTest {

    @TempDir
    Path tmp;

    private File fileOf(String name, String... lines) throws Exception {
        File f = tmp.resolve(name).toFile();
        Files.write(f.toPath(), String.join("\n", lines).getBytes(StandardCharsets.UTF_8));
        return f;
    }

    // decoy_prefix used to be stored as ">" + tag, a FASTA-header-shaped value only ever compared
    // by dead code; ">rev_" matches no protein label, so every run naming a fragger.params had its
    // decoys silently unmarked.
    @Test
    public void fraggerDecoyPrefixIsStoredAsTheBareTag() throws Exception {
        //the tolerance lines every real fragger.params carries; readFraggerParams parses the
        //tolerance unconditionally at the end
        File fragger = fileOf("fragger.params", "decoy_prefix = rev_",
                "fragment_mass_tolerance = 20", "fragment_mass_units = 1");
        HashMap<String, String> params = new HashMap<>();
        params.put("fragger", fragger.getAbsolutePath());

        ParameterUtils.readFraggerParams(params);

        assertEquals("rev_", params.get("decoyPrefix"));
    }

    // The ordering Predictor and FragCastPredictor rely on: the params file is expanded the moment
    // --paramsList is read, so a --decoyPrefix placed after it outranks the file's own decoyPrefix
    // line. Predictor used to assign --decoy-tag to Constants before the file was read, which let
    // the file overwrite it.
    @Test
    public void aDecoyPrefixAfterParamsListOutranksTheFile() throws Exception {
        File params = fileOf("msbooster_params.txt", "decoyPrefix = DECOY_");
        HashMap<String, String> map = new HashMap<>();

        ParameterUtils.processCommandLineInputs(new String[]{
                "--paramsList", params.getAbsolutePath(),
                "--requirePinMzml", "false",
                "--decoyPrefix", "rev_"}, map);

        assertEquals("rev_", map.get("decoyPrefix"));
    }
}
