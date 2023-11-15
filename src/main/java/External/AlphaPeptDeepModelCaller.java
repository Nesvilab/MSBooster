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

package External;

import Features.Constants;

import java.io.BufferedReader;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;

public class AlphaPeptDeepModelCaller {
    public static void callModel(String mode) {
        long startTime = System.nanoTime();
        try {
            if (mode.equals("transfer")) {
                System.out.println("Generating AlphaPeptDeep transfer learned model");
                ProcessBuilder builder = new ProcessBuilder(Constants.AlphaPeptDeep,
                        "transfer", Constants.paramsList);
                System.out.println(String.join(" ", builder.command()));
                builder.redirectErrorStream(true);
                Process process = builder.start();
                InputStream is = process.getInputStream();
                BufferedReader reader = new BufferedReader(new InputStreamReader(is));

                String line = null;
                while ((line = reader.readLine()) != null) {
                    System.out.println(line);
                }
            } else if (mode.equals("predict")) {
                System.out.println("Generating AlphaPeptDeep predictions");
                ProcessBuilder builder = new ProcessBuilder(Constants.AlphaPeptDeep,
                        "predict", Constants.paramsList);
                System.out.println(String.join(" ", builder.command()));
                builder.redirectErrorStream(true);
                Process process = builder.start();
                InputStream is = process.getInputStream();
                BufferedReader reader = new BufferedReader(new InputStreamReader(is));

                String line = null;
                while ((line = reader.readLine()) != null) {
                    System.out.println(line);
                }
            }
        } catch (IOException e) {
            e.printStackTrace();
        }
        long endTime = System.nanoTime();
        long duration = (endTime - startTime);
        System.out.println("Model running took " + duration / 1000000 +" milliseconds");
    }
}
