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

import java.util.ArrayList;
import java.util.List;

public class NumericUtils {
    //move to another class if used more often
    public static float[] doubleToFloat(double[] arr){
        float[] newArr = new float[arr.length];
        for (int i = 0; i < arr.length; i++){
            newArr[i] = (float) arr[i];
        }
        return newArr;
    }

    public static double[] floatToDouble(float[] arr){
        double[] newArr = new double[arr.length];
        for (int i = 0; i < arr.length; i++){
            newArr[i] = arr[i];
        }
        return newArr;
    }

    public static double[] floatToDouble(ArrayList<Float> arr){
        double[] newArr = new double[arr.size()];
        for (int i = 0; i < arr.size(); i++){
            newArr[i] = (double) arr.get(i);
        }
        return newArr;
    }

    public static double sum(List<Double> arr) {
        double sum = 0f;
        for (double f : arr) {
            sum += f;
        }
        return sum;
    }

    public static boolean isUppercaseOrDigit(char c) {
        return Character.isUpperCase(c) || Character.isDigit(c);
    }

    public static boolean massesCloseEnough(double mass1, double mass2) { //TODO: use this method more
        double errorTol = 0.001;
        return Math.abs(mass2 - mass1) < errorTol;
    }
}
