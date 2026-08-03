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

package modelcallers;

import static utils.Print.printInfo;

import allconstants.Constants;

import java.io.BufferedReader;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;

public class AlphaPeptDeepModelCaller {
    public static void callModel(String mode) {
        long startTime = System.nanoTime();
        try {
            if (mode.equals("transfer")) {
                printInfo("Generating AlphaPeptDeep transfer learned model");
                ProcessBuilder builder = new ProcessBuilder(Constants.AlphaPeptDeep,
                        "transfer", Constants.paramsList);
                printInfo(String.join(" ", builder.command()));
                builder.redirectErrorStream(true);
                Process process = builder.start();
                InputStream is = process.getInputStream();
                BufferedReader reader = new BufferedReader(new InputStreamReader(is));

                String line = null;
                while ((line = reader.readLine()) != null) {
                    printInfo(line);
                }
            } else if (mode.equals("predict")) {
                printInfo("Generating AlphaPeptDeep predictions");
                ProcessBuilder builder = new ProcessBuilder(Constants.AlphaPeptDeep,
                        "predict", Constants.paramsList);
                printInfo(String.join(" ", builder.command()));
                builder.redirectErrorStream(true);
                Process process = builder.start();
                InputStream is = process.getInputStream();
                BufferedReader reader = new BufferedReader(new InputStreamReader(is));

                String line = null;
                while ((line = reader.readLine()) != null) {
                    printInfo(line);
                }
            }
        } catch (IOException e) {
            e.printStackTrace();
        }
        long endTime = System.nanoTime();
        long duration = (endTime - startTime);
        printInfo("Model running took " + duration / 1000000 +" milliseconds");
    }
}
