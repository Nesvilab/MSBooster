package Features;

import Exceptions.UnsupportedInputException;

import java.lang.reflect.Field;
import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Map;

//this is what I use in the java jar file
public class MainClass {
    public static void main(String[] args) throws Exception {
        //to do: take constants as input file

        //accept command line inputs
        HashSet<String> fields = new HashSet<>();
        for (Field f : Constants.class.getDeclaredFields()) {
            fields.add(f.getName());
        }

        HashMap<String, String> params = new HashMap<String, String>();

        //setting new values
        for (int i = 0; i < args.length; i++) {
            String key = args[i].substring(2); //remove --
            i++;
            params.put(key, args[i]);
        }

        //adding to constants class
        Constants c = new Constants();
        for (Map.Entry<String, String> entry : params.entrySet()) {
            String key = entry.getKey();
            if (! fields.contains(key)) {
                throw new Exception(entry.getKey() + " is not a valid parameter");
            } else {
                //get class of field
                Field field = Constants.class.getField(key);
                Class<?> myClass = field.getType();

                //parse to appropriate type
                field.set(c, myClass.getConstructor(String.class).newInstance(entry.getValue()));
            }
        }

        //check that features are allowed
        boolean allFeatures = false;
        boolean autoFeatures = false;
        String[] featuresArray = Constants.features.split(",");
        for (String f : featuresArray) {
            if (! Constants.allowedFeatures.contains(f)) {
                throw new UnsupportedInputException("UnsupportedInputException", f + " is not an allowed feature. " +
                        "Please choose from the following: " + Constants.allowedFeatures);
            }
            if (f.equals("all")) {
                allFeatures = true;
                break;
            }
            if (f.equals("auto")) {
                autoFeatures = true;
                break;
            }
        }

        //create file for spectral and RT prediction
        boolean createSpectraRTPredFile = false;
        boolean createDetectPredFile = false;

        //check which ones we need
        if (allFeatures || autoFeatures) {
            createSpectraRTPredFile = true;
            createDetectPredFile = true;
        } else {
            HashSet<String> featureSet = new HashSet<>(Arrays.asList(featuresArray));
            featureSet.retainAll(Constants.spectraRTFeatures);
            if (featureSet.size() > 0) {
                createSpectraRTPredFile = true;
            }

            featureSet = new HashSet<>(Arrays.asList(featuresArray));
            featureSet.retainAll(Constants.detectFeatures);
            if (featureSet.size() > 0) {
                createDetectPredFile = true;
            }
        }

        //call prediction models
    }
}
