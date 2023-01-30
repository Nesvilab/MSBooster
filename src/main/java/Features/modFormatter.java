/*
 * This file is part of MSBooster.
 *
 * MSBooster is free software: you can redistribute it and/or modify it under
 * the terms of the GNU Lesser General Public License as published by the Free
 * Software Foundation, either version 3 of the License, or (at your
 * option) any later version.
 *
 * MSBooster is distributed in the hope that it will be useful, but WITHOUT
 * ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
 * FITNESS FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public License
 * for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with MSBooster. If not, see <https://www.gnu.org/licenses/>.
 */

package Features;

import java.util.HashMap;
//TODO: delete later, as only pepxml uses. pDeep3 can use new PeptideFormatter
public class modFormatter {
    int[] pos;
    double[] mass;
    final HashMap <Double, String> modifications= new HashMap<Double, String>() {{
        put(160.03065, "Carbamidomethyl[C]"); //msfragger 3.1
        //modifications.put(160.0307, "Carbamidomethyl[C]"); //msfragger 3.0
        put(147.0354, "Oxidation[M]");
    }};

    public modFormatter(int[] positions, double[] masses) {
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
