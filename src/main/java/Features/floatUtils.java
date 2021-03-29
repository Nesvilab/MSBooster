package Features;

import java.util.ArrayList;

public class floatUtils {
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
            newArr[i] = (double) arr[i];
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
}
