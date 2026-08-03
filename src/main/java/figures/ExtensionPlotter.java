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
