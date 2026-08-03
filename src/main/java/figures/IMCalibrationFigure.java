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

package figures;

import readers.datareaders.MzmlReader;
import java.util.function.DoubleUnaryOperator;

import java.io.IOException;
import java.util.HashMap;

public class IMCalibrationFigure extends CalibrationFigure {

    private static final String FOLDER_STRING = "IM_calibration_curves";
    private static final String MODE = "IM";

    public IMCalibrationFigure(MzmlReader mzml, String outFile, float opacity,
                               HashMap<String, double[][]> massToData,
                               HashMap<String, DoubleUnaryOperator> loessFunctions,
                               int charge) throws IOException {
        super();
        super.folderString = FOLDER_STRING;
        super.mode = MODE;
        super.charge = String.valueOf(charge);
        plotFigure(mzml, outFile, opacity, massToData, loessFunctions);
    }
}
