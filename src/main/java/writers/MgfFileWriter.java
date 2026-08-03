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

package writers;

import predictions.PredictionEntry;
import predictions.PredictionEntryHashMap;
import utils.NumericUtils;

import java.io.BufferedWriter;
import java.io.FileWriter;
import java.io.IOException;
import java.util.Map;

import static utils.Print.printInfo;

public class MgfFileWriter {
    PredictionEntryHashMap allPreds;
    public MgfFileWriter(PredictionEntryHashMap allPreds) {
        this.allPreds = allPreds;
    }

    public void write(String outfile) throws IOException {
        printInfo("Writing " + outfile);
        BufferedWriter bw = new BufferedWriter(new FileWriter(outfile));
        for (Map.Entry<String, PredictionEntry> entry : allPreds.entrySet()) {
            String[] peptide = entry.getKey().split("\\|");
            PredictionEntry pe = entry.getValue();
            if (pe.numFragments() == 0) { //was not able to have its mzs predicted
                continue;
            }

            //isotopic information
            int isotopeSum = NumericUtils.intSum(pe.isotopes);

            bw.write("BEGIN IONS" + "\n");
            bw.write("TITLE=" + peptide[0] + "\n");
            bw.write("CHARGE=" + peptide[1] + "\n");
            bw.write("RT=" + pe.RT + "\n");
            bw.write("1/K0=" + pe.IM + "\n");
            for (int i = 0; i < pe.numFragments(); i++) {
                //no need to filter by intensity since that's already done
                if (pe.getIntensity(i) != 0) {
                    if (isotopeSum > 0) {
                        bw.write(pe.getMz(i) + "\t" + pe.getIntensity(i) + " " + pe.getIonTypeString(i) + " " + pe.isotopes[i] + "\n");
                    } else {
                        bw.write(pe.getMz(i) + "\t" + pe.getIntensity(i) + " " + pe.getIonTypeString(i) + "\n");
                    }
                }
            }
            bw.write("END IONS" + "\n");
        }
        bw.close();
    }
}
