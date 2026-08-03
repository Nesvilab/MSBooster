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

package features.rtandim;

import allconstants.Constants;
import mainsteps.MzmlScanNumber;
import mainsteps.PeptideObj;
import readers.datareaders.MzmlReader;
import umich.ms.fileio.exceptions.FileParsingException;

import java.io.IOException;
import java.util.ArrayList;

public class IMFunctions {
    public static int numCharges = 7;

    public static ArrayList<Float>[][] IMbins(MzmlReader mzml) throws IOException, FileParsingException {
        //hard coded as 2, but if there are higher IM values, this can change
        int numBins = 2 * Constants.IMbinMultiplier;

        ArrayList<Float>[][] predIMround = new ArrayList[numCharges][numBins + 1];
        for (int c = 0; c < numCharges; c++) {
            for (int col = 0; col < numBins + 1; col++) {
                predIMround[c][col] = new ArrayList<Float>();
            }
        }

        //iterate through scanNumbers
        for (int scanNum : mzml.getScanNums()) {
            MzmlScanNumber scanNumObj = mzml.getScanNumObject(scanNum);
            int round = (int) (scanNumObj.IM * Constants.IMbinMultiplier); //experimental RT for this scan, assume in minutes

            //iterate through PSMs
            for (int i = 1; i < scanNumObj.peptideObjects.size() + 1; i++) {
                PeptideObj pep = scanNumObj.getPeptideObject(i);
                if (pep == null) {
                    break;
                }
                int charge = pep.charge - 1;

                int instances = Math.max(1, -1 * (int) Math.ceil(Math.log10(Double.parseDouble(pep.escore)))); //this version avoids empty bins
                for (int j = 0; j < instances; j++) {
                    predIMround[charge][round].add(pep.IM);
                }
            }
        }
        return predIMround;
    }
}
