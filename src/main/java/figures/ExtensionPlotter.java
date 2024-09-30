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

import allconstants.Constants;

import static utils.Print.printError;

import org.knowm.xchart.BitmapEncoder;
import org.knowm.xchart.VectorGraphicsEncoder;
import org.knowm.xchart.internal.chartpart.Chart;

import java.io.IOException;

public class ExtensionPlotter {
    public static void plot(Chart<?, ?> chart, String basename) throws IOException {
        Constants.plotExtension = Constants.plotExtension.toLowerCase();
        switch (Constants.plotExtension) {
            case "png":
                BitmapEncoder.saveBitmap(chart, basename,
                        BitmapEncoder.BitmapFormat.PNG);
                break;
            case "pdf":
                VectorGraphicsEncoder.saveVectorGraphic(chart, basename,
                        VectorGraphicsEncoder.VectorGraphicsFormat.PDF);
                break;
            default:
                printError(Constants.plotExtension + " not supported for plotting. Exiting");
                System.exit(1);
        }
    }
}
