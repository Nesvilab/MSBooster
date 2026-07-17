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

import java.io.File;
import java.io.IOException;

public class ExtensionPlotter {
    @FunctionalInterface
    private interface ChartWriter {
        void write() throws IOException;
    }

    //Both xchart encoders open their output file before they render the chart, so a chart that
    //throws while painting truncates the file and leaves it behind at 0 bytes. That empty figure
    //is indistinguishable from a real one on disk and hides which feature actually failed, so
    //drop it and let the caller see the original failure.
    private static void writeOrDeletePartial(String path, ChartWriter writer) throws IOException {
        try {
            writer.write();
        } catch (Throwable t) { //Errors too, a partial figure must not outlive any failed render
            new File(path).delete();
            throw t;
        }
    }

    public static void plot(Chart<?, ?> chart, String basename) throws IOException {
        Constants.plotExtension = Constants.plotExtension.toLowerCase();
        switch (Constants.plotExtension) {
            case "png":
                writeOrDeletePartial(
                        BitmapEncoder.addFileExtension(basename, BitmapEncoder.BitmapFormat.PNG),
                        () -> BitmapEncoder.saveBitmap(chart, basename,
                                BitmapEncoder.BitmapFormat.PNG));
                break;
            case "pdf":
                writeOrDeletePartial(
                        VectorGraphicsEncoder.addFileExtension(basename,
                                VectorGraphicsEncoder.VectorGraphicsFormat.PDF),
                        () -> VectorGraphicsEncoder.saveVectorGraphic(chart, basename,
                                VectorGraphicsEncoder.VectorGraphicsFormat.PDF));
                break;
            default:
                printError(Constants.plotExtension + " not supported for plotting. Exiting");
                System.exit(1);
        }
    }
}
