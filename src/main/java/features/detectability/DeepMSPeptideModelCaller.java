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
