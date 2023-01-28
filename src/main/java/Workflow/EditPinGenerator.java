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

package Workflow;

import Features.Constants;
import org.apache.commons.lang3.ArrayUtils;

import java.io.*;
import java.nio.file.Paths;

public class EditPinGenerator {
    public static void leaveOneOut(String editedPin) throws IOException {
        //what columns have features?
        int colMin = 28;
        int colMax = 41;
        String[] features = Constants.features.split(",");

        //set file paths
        File oldPin = new File(editedPin);
        String parentDirectory = oldPin.getParent();

        for (int i = colMin; i < colMax + 1; i++) {
            //what feature is left out?
            String excludedFeature = features[i - colMin];
            System.out.println(excludedFeature);

            //make folder for excluded feature
            File newDir = new File(parentDirectory + File.separator + "exclude" + excludedFeature);
            if (!newDir.exists()){
                newDir.mkdirs();
            }

            //read in file, but leave out one column
            System.out.println("creating temporary pin");
            BufferedReader oldPinReader = new BufferedReader(new FileReader(editedPin));
            BufferedWriter newPin = new BufferedWriter(new FileWriter(parentDirectory + File.separator + "tmp.pin"));
            String line;
            while( (line = oldPinReader.readLine()) != null) {
                String[] columns = line.split("\t");
                columns = ArrayUtils.remove(columns, i);
                newPin.write(String.join("\t", columns));
                newPin.write("\n");
            }
            oldPinReader.close();
            newPin.close();

            System.out.println("running through percolator");
            ProcessBuilder builder = new ProcessBuilder("C:/Users/kevin/Downloads/proteomics/prediction_tools/percolator-v3-05/bin/percolator.exe",
                    "--num-threads", "9",
                    "--only-psms",
                    "--no-terminate",
                    "--post-processing-tdc",
                    "--results-psms", newDir.getCanonicalPath() + File.separator + excludedFeature + "_exclude.tsv",
                    "--decoy-results-psms", newDir.getCanonicalPath() + File.separator + excludedFeature + "_excludeD.tsv",
                    parentDirectory + File.separator + "tmp.pin");
            System.out.println(String.join(" ", builder.command()));
            builder.redirectErrorStream(true);
            Process process = builder.start();
            InputStream is = process.getInputStream();
            BufferedReader reader = new BufferedReader(new InputStreamReader(is));

            line = null;
            while ((line = reader.readLine()) != null) {
                System.out.println(line);
            }

            //create interact files
            System.out.println("creating interact.pep.xml");
            PercolatorOutputToPepXML.percolatorToPepXML(Paths.get(editedPin),
                    parentDirectory + File.separator + oldPin.getName().substring(7, oldPin.getName().length() - 4),
                    Paths.get(newDir.getCanonicalPath() + File.separator + excludedFeature + "_exclude.tsv"),
                    Paths.get(newDir.getCanonicalPath() + File.separator + excludedFeature + "_excludeD.tsv"),
                    Paths.get(newDir.getCanonicalPath() + File.separator + "interact-" +
                            oldPin.getName().substring(7, oldPin.getName().length() - 4)),
                    "DIA");
        }
    }

    public static void leaveOneIn(String editedPin) throws IOException {
        //what columns have features?
        int colMin = 28;
        int colMax = 41;
        String[] features = Constants.features.split(",");

        //set file paths
        File oldPin = new File(editedPin);
        String parentDirectory = oldPin.getParent();

        for (int i = colMin; i < colMax + 1; i++) {
            //what feature is left out?
            String includedFeature = features[i - colMin];
            System.out.println(includedFeature);

            //make folder for excluded feature
            File newDir = new File(parentDirectory + File.separator + includedFeature);
            if (!newDir.exists()){
                newDir.mkdirs();
            }

            //read in file, but leave out one column
            System.out.println("creating temporary pin");
            BufferedReader oldPinReader = new BufferedReader(new FileReader(editedPin));
            BufferedWriter newPin = new BufferedWriter(new FileWriter(parentDirectory + File.separator + "tmp.pin"));
            String line;
            while( (line = oldPinReader.readLine()) != null) {
                String[] columns = line.split("\t");
                for (int j = colMax; j > colMin - 1; j--) {
                    if (j != i) {
                        columns = ArrayUtils.remove(columns, j);
                    }
                }
                newPin.write(String.join("\t", columns));
                newPin.write("\n");
            }
            oldPinReader.close();
            newPin.close();

            System.out.println("running through percolator");
            ProcessBuilder builder = new ProcessBuilder("C:/Users/kevin/Downloads/proteomics/prediction_tools/percolator-v3-05/bin/percolator.exe",
                    "--num-threads", "9",
                    "--only-psms",
                    "--no-terminate",
                    "--post-processing-tdc",
                    "--results-psms", newDir.getCanonicalPath() + File.separator + includedFeature + ".tsv",
                    "--decoy-results-psms", newDir.getCanonicalPath() + File.separator + includedFeature + "D.tsv",
                    parentDirectory + File.separator + "tmp.pin");
            System.out.println(String.join(" ", builder.command()));
            builder.redirectErrorStream(true);
            Process process = builder.start();
            InputStream is = process.getInputStream();
            BufferedReader reader = new BufferedReader(new InputStreamReader(is));

            line = null;
            while ((line = reader.readLine()) != null) {
                System.out.println(line);
            }

            //create interact files
            System.out.println("creating interact.pep.xml");
            PercolatorOutputToPepXML.percolatorToPepXML(Paths.get(editedPin),
                    parentDirectory + File.separator + oldPin.getName().substring(7, oldPin.getName().length() - 4),
                    Paths.get(newDir.getCanonicalPath() + File.separator + includedFeature + ".tsv"),
                    Paths.get(newDir.getCanonicalPath() + File.separator + includedFeature + "D.tsv"),
                    Paths.get(newDir.getCanonicalPath() + File.separator + "interact-" +
                            oldPin.getName().substring(7, oldPin.getName().length() - 4)),
                    "DIA");
        }

        //remove tmp.pin
        File tmp = new File(parentDirectory + File.separator + "tmp.pin");
        tmp.delete();
    }

    public static void main(String[] args) throws IOException {
        //leaveOneOut("C:/Users/kevin/Downloads/proteomics/cptac/2021-2-21/pep1XML1tmp/percToPep/edited_CPTAC_CCRCC_W_JHU_LUMOS_C3L-01665_T.pin");
        leaveOneIn("C:/Users/kevin/Downloads/proteomics/cptac/2021-2-21/pep1XML1tmp/percToPep/edited_CPTAC_CCRCC_W_JHU_LUMOS_C3L-01665_T.pin");
    }
}
