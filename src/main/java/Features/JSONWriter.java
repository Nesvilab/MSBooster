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

import static utils.Print.printError;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.*;
import java.util.concurrent.ExecutionException;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Future;
import org.json.JSONArray;
import org.json.JSONObject;

public class JSONWriter {
    String url;
    String model;
    String property;
    String[] peptides;
    int[] charges;
    float[] nces;
    String[] instruments;
    String[] fragmentations;
    boolean TMT = false;
    int numFiles;
    int iteration = -1;
    public List<Future> futureList = new ArrayList<>(Constants.numThreads);

    public static final int maxJsonLength = 1000;

    public JSONWriter(String url, HashSet<String> entries, boolean fullSet) {
        this.url = url;
        this.model = url.toLowerCase().split("_")[0];
        if (Constants.KoinaRTmodels.contains(url)) {
            property = "rt";
        } else if (Constants.KoinaMS2models.contains(url)) {
            property = "ms2";
        }
        if (Constants.KoinaTMTmodels.contains(url)) {
            TMT = true;
        }

        //set entries
        this.peptides = new String[entries.size()];
        this.charges = new int[entries.size()];
        this.nces = new float[entries.size()];
        this.instruments = new String[entries.size()];
        this.fragmentations = new String[entries.size()];

        //fullset is to indicate if the stripped peptide length is there
        //this is true if using pinreader.createjson
        if ((model.contains("alphapept")) && (fullSet)) {
            List<String> sortedList = new ArrayList<>(entries);

            // Custom comparator to sort strings based on their stripped length
            Comparator<String> customComparator = (str1, str2) -> {
                // Extract the numerical part from both strings.
                int num1 = Integer.parseInt(str1.substring(str1.lastIndexOf(',') + 1));
                int num2 = Integer.parseInt(str2.substring(str2.lastIndexOf(',') + 1));
                // Compare these numerical parts
                return Integer.compare(num1, num2);
            };

            Collections.sort(sortedList, customComparator);

            int i = 0;
            for (String entry : sortedList) {
                String[] info = entry.split(",");
                this.peptides[i] = info[0];
                this.charges[i] = Integer.parseInt(info[1]);
                this.nces[i] = Float.parseFloat(info[2]);
                this.instruments[i] = info[3].toUpperCase();
                this.fragmentations[i] = info[4];
                i++;
            }
        } else {
            int i = 0;
            for (String entry : entries) {
                String[] info = entry.split(",");
                this.peptides[i] = info[0];
                this.charges[i] = Integer.parseInt(info[1]);
                this.nces[i] = Float.parseFloat(info[2]);
                this.instruments[i] = info[3].toUpperCase();
                this.fragmentations[i] = info[4];
                i++;
            }
        }

        numFiles = (int) Math.ceil((float) this.peptides.length / (float) maxJsonLength);
    }

    //this version if we need to spawn smaller ones
    public JSONWriter(JSONWriter parent, int iteration) {
        this.iteration = iteration;
        int start = iteration * maxJsonLength;
        this.url = parent.url;
        this.model = url.toLowerCase().split("_")[0];
        if (Constants.KoinaRTmodels.contains(url)) {
            property = "rt";
        } else if (Constants.KoinaMS2models.contains(url)) {
            property = "ms2";
        }
        if (Constants.KoinaTMTmodels.contains(url)) {
            TMT = true;
        }

        //take subarrays
        this.peptides = Arrays.copyOfRange(parent.peptides, start,
                Math.min(start + maxJsonLength, parent.peptides.length));
        this.charges = Arrays.copyOfRange(parent.charges, start,
                Math.min(start + maxJsonLength, parent.peptides.length));
        this.nces = Arrays.copyOfRange(parent.nces, start,
                Math.min(start + maxJsonLength, parent.peptides.length));
        this.instruments = Arrays.copyOfRange(parent.instruments, start,
                Math.min(start + maxJsonLength, parent.peptides.length));
        this.fragmentations = Arrays.copyOfRange(parent.fragmentations, start,
                Math.min(start + maxJsonLength, parent.peptides.length));
    }

    public String write(boolean createDir, String jsonOutFolder, ExecutorService executorService)
            throws IOException, ExecutionException, InterruptedException {
        if (createDir) {
            MyFileUtils.createWholeDirectory(jsonOutFolder);
        } else {
            if (! Files.exists(Paths.get(jsonOutFolder))) {
                Files.createDirectories(Paths.get(jsonOutFolder));
            }
        }

        String fileName = "";
        if (numFiles > 1) {
            futureList.clear();
            for (int i = 0; i < Constants.numThreads; i++) {
                int start = (int) (numFiles * (long) i) / Constants.numThreads;
                int end = (int) (numFiles * (long) (i + 1)) / Constants.numThreads;
                futureList.add(executorService.submit(() -> {
                    for (int rep = start; rep < end; rep ++) {
                        JSONWriter jw = new JSONWriter(this, rep);
                        try {
                            jw.write(false, jsonOutFolder, executorService);
                        } catch (IOException | ExecutionException | InterruptedException e) {
                            e.printStackTrace();
                        }
                    }
                }));
            }
            for (Future future : futureList) {
                future.get();
            }
        } else {
            switch (model) {
                case "alphapept":
                case "prosit":
                case "ms2pip":
                case "deeplc":
                    // Create the JSON data structure
                    JSONObject jsonData = new JSONObject();
                    jsonData.put("id", "0");

                    JSONArray inputsArray = new JSONArray();
                    JSONObject peptideObject = new JSONObject();
                    peptideObject.put("name", "peptide_sequences");
                    peptideObject.put("shape", new JSONArray("[" + peptides.length + ",1]"));
                    peptideObject.put("datatype", "BYTES");

                    JSONArray peptideData = new JSONArray();
                    JSONArray innerArray;

//                    if (TMT) {
//                        for (int i = 0; i < peptides.length; i++) {
//                            peptides[i] = "[UNIMOD:737]-" + peptides[i];
//                        }
//                    }
                    for (String pep : peptides) {
                        innerArray = new JSONArray();
                        innerArray.put(pep);
                        peptideData.put(innerArray);
                    }

                    peptideObject.put("data", peptideData);
                    inputsArray.put(peptideObject);

                    if (property.equals("ms2")) {
                        //precursor charges
                        JSONObject chargeObject = new JSONObject();
                        chargeObject.put("name", "precursor_charges");
                        chargeObject.put("shape", new JSONArray("[" + peptides.length + ",1]"));
                        chargeObject.put("datatype", "INT32");

                        JSONArray chargeData = new JSONArray();

                        for (int charge : charges) {
                            innerArray = new JSONArray();
                            innerArray.put(charge);
                            chargeData.put(innerArray);
                        }

                        chargeObject.put("data", chargeData);
                        inputsArray.put(chargeObject);

                        //collision energies
                        if (!url.equals("Prosit_2020_intensity_CID") && !url.equals("ms2pip_2021_HCD")) {
                            JSONObject nceObject = new JSONObject();
                            nceObject.put("name", "collision_energies");
                            nceObject.put("shape", new JSONArray("[" + peptides.length + ",1]"));
                            nceObject.put("datatype", "FP32");

                            JSONArray nceData = new JSONArray();

                            for (float nce : nces) {
                                innerArray = new JSONArray();
                                innerArray.put(nce);
                                nceData.put(innerArray);
                            }

                            nceObject.put("data", nceData);
                            inputsArray.put(nceObject);
                        }

                        //instrument types
                        if (model.equals("alphapept")) {
                            JSONObject instrumentObject = new JSONObject();
                            instrumentObject.put("name", "instrument_types");
                            instrumentObject.put("shape", new JSONArray("[" + peptides.length + ",1]"));
                            instrumentObject.put("datatype", "BYTES");

                            JSONArray instrumentData = new JSONArray();

                            for (String instrument : instruments) {
                                innerArray = new JSONArray();
                                innerArray.put(instrument);
                                instrumentData.put(innerArray);
                            }

                            instrumentObject.put("data", instrumentData);
                            inputsArray.put(instrumentObject);
                        }
                        if (TMT) {
                            JSONObject fragmentationObject = new JSONObject();
                            fragmentationObject.put("name", "fragmentation_types");
                            fragmentationObject.put("shape", new JSONArray("[" + peptides.length + ",1]"));
                            fragmentationObject.put("datatype", "BYTES");

                            JSONArray fragmentationData = new JSONArray();

                            for (String fragmentation : fragmentations) {
                                innerArray = new JSONArray();
                                innerArray.put(fragmentation);
                                fragmentationData.put(innerArray);
                            }

                            fragmentationObject.put("data", fragmentationData);
                            inputsArray.put(fragmentationObject);
                        }
                    }

                    jsonData.put("inputs", inputsArray);

                    // Write JSON data to a file
                    if (iteration != -1) {
                        fileName = jsonOutFolder + File.separator +
                                "spectraRT_" + url + iteration + property + ".json";
                    } else {
                        fileName = jsonOutFolder + File.separator +
                                "spectraRT_" + url + property + ".json";
                    }
                    try (FileWriter fileWriter = new FileWriter(fileName)) {
                        fileWriter.write(jsonData.toString(4)); // Using 4 for indentation
                    } catch (IOException e) {
                        e.printStackTrace();
                    }
                    break;
                default:
                    printError(model + " is not supported by Koina. Exiting now.");
                    System.exit(0);
            }
        }
        return jsonOutFolder;
    }
}
