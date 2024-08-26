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

package utils;

import static utils.Print.printInfo;

import java.util.Arrays;

public class RandomThings {

    public RandomThings(){
    }

    public static void main(String[] args) {
        float[] a = new float[1];
        a[0] = 1f;
        float[] b = a;
        b[0] = 2f;
        printInfo(Arrays.toString(a));
        printInfo(Arrays.toString(b));
    }
}

