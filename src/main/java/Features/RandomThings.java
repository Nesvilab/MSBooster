package Features;

import umich.ms.fileio.exceptions.FileParsingException;

import java.io.IOException;
import java.lang.reflect.InvocationTargetException;
import java.util.Arrays;
import java.util.HashMap;
import java.util.concurrent.ExecutionException;

public class RandomThings {

    public RandomThings(){
    }

    public static void main(String[] args) {
        HashMap<String, String[]> a = new HashMap<>();
        //a.put("a", new String[1]);
        String[] b = a.get("a");
        System.out.println(b != null);
    }
}

