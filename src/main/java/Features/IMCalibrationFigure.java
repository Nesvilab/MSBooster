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

import kotlin.jvm.functions.Function1;

import java.io.IOException;
import java.util.HashMap;

public class IMCalibrationFigure extends CalibrationFigure {

    private static final String FOLDER_STRING = "IM_calibration_curves";
    private static final String MODE = "IM";

    public IMCalibrationFigure(MzmlReader mzml, String outFile, float opacity,
                               HashMap<String, double[][]> massToData,
                               HashMap<String, Function1<Double, Double>> loessFunctions,
                               int charge) throws IOException {
        super();
        super.folderString = FOLDER_STRING;
        super.mode = MODE;
        super.charge = String.valueOf(charge);
        plotFigure(mzml, outFile, opacity, massToData, loessFunctions);
    }
}
