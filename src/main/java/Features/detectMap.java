package Features;

import java.io.*;
import java.util.HashMap;

//only using DeepMSPeptide so far for prediction
public class detectMap {
    HashMap<String, Float> detectabilities = new HashMap<>();

    public detectMap(String detectFile) throws IOException {
        BufferedReader br = new BufferedReader(new FileReader(new File(detectFile)));
        br.readLine(); //header
        String line;
        while ((line = br.readLine()) != null) {
            String[] info = line.split("\t");
            detectabilities.put(info[0], Float.valueOf(info[1]));
        }
    }

    public float getDetectability(String pep) {
        //try to intelligently reformat peptide to one the Hashmap recognizes
        try {
            float d = detectabilities.get(pep);
            return d;
        } catch (Exception e) {
            return detectabilities.get(pep.split("\\|")[0]);
        }
    }
}