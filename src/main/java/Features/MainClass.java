package Features;

import java.lang.reflect.Field;
import java.util.HashMap;
import java.util.Map;
import java.util.Set;

//this is what I use in the java jar file
public class MainClass {
    public static void main(String[] args) throws Exception {
        //take constants as input file
        HashMap<String, String> hash1 = new HashMap<String, String>() {{
            put("ppmTolerance", "20");
            put("useTopFragments", "true");
            put("topFragments", "12");
            put("removeRankPeaks", "true");
            put("RTregressionSize", "5000");
            put("uniformPriorPercentile", "10");
            put("predFileFormat", "bin");
            put("predFile", null);
            put("usedFeatures", "brayCurtis euclideanDistance cosineSimilarity spectralContrastAngle pearsonCorr " +
                    "dotProduct deltaRTlinear deltaRTbins RTzscore RTprobability RTprobabilityUnifPrior");
            put("outputFolder", null);
            put("params", null);
        }};

        HashMap<String, String> hash2 = new HashMap<String, String>();

            //setting new values
        for (int i = 0; i < args.length; i++) {
            if (args[i].startsWith("--")) {
                String key = args[i].substring(2);

                //in case of usedFeatures, could span multiple elements
                StringBuilder sb = new StringBuilder(" ");
                while (! args[i + 1].substring(0, 2).equals("--")) {
                    i++;
                    sb.append(args[i]);

                    if (i + 1 >= args.length) { //last argument
                        break;
                    }
                }

                hash2.put(key, sb.toString());
            }
        }

            //adding to constants class
        Set<String> hash1Keys = hash1.keySet();
        Constants c = new Constants();
        for (Map.Entry<String, String> entry : hash2.entrySet()) {
            String key = entry.getKey();
            if (! hash1Keys.contains(key)) {
                throw new Exception(entry.getKey() + " is not a valid parameter");
            } else {
                //get class of field
                Field field = Constants.class.getField(key);
                Class<?> myClass = field.getType();

                //parse to appropriate type

                //field.set(c, entry.getValue());
            }
            System.out.println(Constants.ppmTolerance);
        }

        //create file for spectral and RT prediction
    }
}
