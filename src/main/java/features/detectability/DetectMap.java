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
