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

package peptideptmformatting;

import java.util.HashMap;
//TODO: delete later, as only pepxml uses. pDeep3 can use new PeptideFormatter
public class ModFormatter {
    int[] pos;
    double[] mass;
    final HashMap <Double, String> modifications= new HashMap<Double, String>() {{
        put(160.03065, "Carbamidomethyl[C]"); //msfragger 3.1
        //modifications.put(160.0307, "Carbamidomethyl[C]"); //msfragger 3.0
        put(147.0354, "Oxidation[M]");
    }};

    public ModFormatter(int[] positions, double[] masses) {
        pos = positions;
        mass = masses;
    }

    public String format() {
        int len = pos.length;
        String formatted = "";
        for (int i = 0; i < len; i++) {
            String toAdd = pos[i] + "," + modifications.get(mass[i]) + ";";
            formatted += toAdd;
        }
        return formatted;
    }
}
