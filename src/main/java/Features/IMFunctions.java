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

package Features;

import umich.ms.fileio.exceptions.FileParsingException;

import java.io.IOException;
import java.util.ArrayList;
import java.util.Comparator;
import java.util.HashMap;
import java.util.stream.IntStream;

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
