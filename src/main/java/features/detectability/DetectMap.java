/*
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or
 * implied. See the License for the specific language governing
 * permissions and limitations under the License.
 */

package features.detectability;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.HashMap;

//only using DeepMSPeptide so far for prediction
public class DetectMap {
    HashMap<String, Float> detectabilities = new HashMap<>();

    public DetectMap(String detectFile) throws IOException {
//        long startTime = System.nanoTime();
        BufferedReader br = new BufferedReader(new FileReader(new File(detectFile)));
        br.readLine(); //header
        String line;
        while ((line = br.readLine()) != null) {
            String[] info = line.split("\t");
            detectabilities.put(info[0], Float.valueOf(info[1]));
        }
//        long endTime = System.nanoTime();
//        long duration = (endTime - startTime);
//        printInfo("Detectability map loading took " + duration / 1000000 +" milliseconds");
    }

    public float getDetectability(String pep) {
        //try to intelligently reformat peptide to one the Hashmap recognizes
        try {
            return detectabilities.get(pep.split("\\|")[0]);
        } catch (Exception e) {
            return detectabilities.get(pep);
        }
    }

    public void clear() {
        detectabilities.clear();
    }
}
