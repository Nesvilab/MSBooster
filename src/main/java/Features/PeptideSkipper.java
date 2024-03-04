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

public class PeptideSkipper {
    //provide peptide and see if it may be problematic
    public static boolean skipPeptide(String stripped, String charge) {
        //letters
        for (char c : "OUBZJX".toCharArray()) {
            if (stripped.indexOf(c) != -1) {
                return true;
            }
        }
        //length
        if (stripped.length() < 7 || stripped.length() > 20) {
            return true;
        }
        //charge
        return Integer.parseInt(charge) > 6;
    }
}
