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
import java.util.Comparator;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

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

    public static int intSum(int[] arr) {
        int sum = 0;
        for (int f : arr) {
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
        
    //method from chat gpt to return ranked indices of two arrays
    public static void getRanks(float[] arr1, float[] arr2, List<Integer> rankList1, List<Integer> rankList2) {
        List<float[]> combined = new ArrayList<>();

        // Add arr1 elements with indices
        for (int i = 0; i < arr1.length; i++) {
            combined.add(new float[]{arr1[i], i, 1}); // 1 means from arr1
        }
        // Add arr2 elements with indices
        for (int i = 0; i < arr2.length; i++) {
            combined.add(new float[]{arr2[i], i, 2}); // 2 means from arr2
        }

        // Sort the list by float values
        combined.sort(Comparator.comparingDouble(a -> a[0]));

        // Assign ranks
        Map<Integer, List<Integer>> rankMap1 = new HashMap<>();
        Map<Integer, List<Integer>> rankMap2 = new HashMap<>();
        for (int rank = 0; rank < combined.size(); rank++) {
            float[] item = combined.get(rank);
            int index = (int) item[1];
            int arrayType = (int) item[2];

            if (arrayType == 1) {
                rankMap1.computeIfAbsent(index, k -> new ArrayList<>()).add(rank);
            } else {
                rankMap2.computeIfAbsent(index, k -> new ArrayList<>()).add(rank);
            }
        }

        // Populate rank lists
        for (int i = 0; i < arr1.length; i++) {
            rankList1.add(rankMap1.get(i).get(0));
        }
        for (int i = 0; i < arr2.length; i++) {
            rankList2.add(rankMap2.get(i).get(0));
        }
    }
}
