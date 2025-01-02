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

import allconstants.Constants;

import java.io.BufferedReader;
import java.io.File;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;

import static utils.Print.printInfo;

public class DeepMSPeptideModelCaller {
    public static void callModel() {
        long startTime = System.nanoTime();
        try {
            printInfo("Generating DeepMSPeptide predictions");
            ProcessBuilder builder = new ProcessBuilder("." + File.separator + "DeepMSPeptideRevised.exe",
                    Constants.detectPredInput);
            printInfo(String.join(" ", builder.command()));
            builder.redirectErrorStream(true);
            Process process = builder.start();
            InputStream is = process.getInputStream();
            BufferedReader reader = new BufferedReader(new InputStreamReader(is));

            String line = null;
            while ((line = reader.readLine()) != null) {
                printInfo(line);
            }

            Constants.detectPredFile =
                    Constants.detectPredInput.substring(0, Constants.detectPredInput.length() - 4) +
                            "_Predictions.txt";
            printInfo("Done generating DeepMSPeptide predictions");
        } catch (IOException e) {
            e.printStackTrace();
            System.exit(1);
        }
        long endTime = System.nanoTime();
        long duration = (endTime - startTime);
        printInfo("Model running took " + duration / 1000000 +" milliseconds");
    }
}
