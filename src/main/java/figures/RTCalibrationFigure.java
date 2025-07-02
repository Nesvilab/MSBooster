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

package figures;

import readers.datareaders.MzmlReader;
import kotlin.jvm.functions.Function1;

import java.io.IOException;
import java.util.HashMap;

public class RTCalibrationFigure extends CalibrationFigure {

    private static final String FOLDER_STRING = "RT_calibration_curves";
    private static final String MODE = "RT";

    public RTCalibrationFigure(MzmlReader mzml, String outFile, float opacity,
                               HashMap<String, double[][]> massToData,
                               HashMap<String, Function1<Double, Double>> loessFunctions) throws IOException {
        super();
        super.folderString = FOLDER_STRING;
        super.mode = MODE;
        curves = plotFigure(mzml, outFile, opacity, massToData, loessFunctions);
    }

    //this version is calibrated
    public RTCalibrationFigure(MzmlReader mzml, String outFile, float opacity,
                               HashMap<String, double[][]> massToData) throws IOException {
        super();
        super.folderString = FOLDER_STRING;
        super.mode = MODE;
        super.yaxislabel = "calibrated";
        super.regressionLabel = "y = x";

        //return y = x
        HashMap<String, Function1<Double, Double>> identityFunctions = new HashMap<>();
        for (String key : massToData.keySet()) {
            identityFunctions.put(key, identityFunction());
        }
        plotFigure(mzml, outFile.substring(0, outFile.length() - 4) + "_calibrated.pin", opacity,
                massToData, identityFunctions);
    }

    private static Function1<Double, Double> identityFunction() {
        return x -> {
            return x;  // y = x
        };
    }
}