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

package writers;

import allconstants.Constants;
import koinaclasses.KoinaException;
import koinaclasses.KoinaModelRegistry;
import koinaclasses.KoinaModelSpec;
import org.json.JSONArray;
import org.json.JSONObject;
import utils.Multithreader;
import utils.MyFileUtils;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.*;
import java.util.concurrent.ExecutionException;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Future;

public class JSONWriter {
    String url;
    KoinaModelSpec spec;
    String property;
    String[] peptides;
    int[] charges;
    float[] nces;
    String[] instruments;
    String[] fragmentations;
    int numFiles;
    int iteration = -1;
    //most precursors in one request file
    int batchSize;
    public List<Future> futureList = new ArrayList<>(Constants.numThreads);

    /**
     * @param batchSize the most precursors one request file holds: KoinaModelCheck.batchSize, which is the
     *                  model's spec batch size unless the Koina server takes smaller batches
     * @throws KoinaException if {@code url} is not a Koina model MSBooster can call
     */
    public JSONWriter(String url, Set<String> entries, boolean fullSet, int batchSize) {
        this.url = url;
        this.spec = KoinaModelRegistry.find(url);
        if (spec == null) {
            throw new KoinaException(url + " is not a Koina model MSBooster can call");
        }
        if (batchSize < 1) {
            throw new IllegalArgumentException("A Koina request holds at least one precursor, not " + batchSize);
        }
        this.property = spec.property(Constants.auxSpectraModel.equals(url));
        this.batchSize = batchSize;

        //set entries
        this.peptides = new String[entries.size()];
        this.charges = new int[entries.size()];
        this.nces = new float[entries.size()];
        this.instruments = new String[entries.size()];
        this.fragmentations = new String[entries.size()];

        //fullset is to indicate if the stripped peptide length is there
        //this is true if using pinreader.createjson
        if (spec.sortBatchByLength() && fullSet) {
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

        numFiles = (int) Math.ceil((float) this.peptides.length / (float) batchSize);
    }

    //this version if we need to spawn smaller ones
    public JSONWriter(JSONWriter parent, int iteration) {
        this.iteration = iteration;
        this.url = parent.url;
        this.spec = parent.spec;
        this.property = parent.property;
        this.batchSize = parent.batchSize;
        int start = iteration * batchSize;

        //take subarrays
        this.peptides = Arrays.copyOfRange(parent.peptides, start,
                Math.min(start + batchSize, parent.peptides.length));
        this.charges = Arrays.copyOfRange(parent.charges, start,
                Math.min(start + batchSize, parent.peptides.length));
        this.nces = Arrays.copyOfRange(parent.nces, start,
                Math.min(start + batchSize, parent.peptides.length));
        this.instruments = Arrays.copyOfRange(parent.instruments, start,
                Math.min(start + batchSize, parent.peptides.length));
        this.fragmentations = Arrays.copyOfRange(parent.fragmentations, start,
                Math.min(start + batchSize, parent.peptides.length));
    }

    public String write(boolean createDir, String jsonOutFolder, ExecutorService executorService)
            throws IOException, ExecutionException, InterruptedException {
        if (createDir) {
            MyFileUtils.deleteWholeDirectory(jsonOutFolder);
            MyFileUtils.createWholeDirectory(jsonOutFolder);
        } else {
            if (! Files.exists(Paths.get(jsonOutFolder))) {
                Files.createDirectories(Paths.get(jsonOutFolder));
            }
        }

        String fileName = "";
        if (iteration == -1 && peptides.length == 0) {
            //no precursor to send, so no request: Koina refuses a request without one
        } else if (numFiles > 1) {
            futureList.clear();
            Multithreader mt = new Multithreader(numFiles, Constants.numThreads);
            for (int i = 0; i < Constants.numThreads; i++) {
                int finalI = i;
                futureList.add(executorService.submit(() -> {
                    for (int rep = mt.indices[finalI]; rep < mt.indices[finalI + 1]; rep++) {
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
            // Create the JSON data structure
            JSONObject jsonData = new JSONObject();
            jsonData.put("id", "0");

            JSONArray inputsArray = new JSONArray();
            JSONArray innerArray;

            JSONArray peptideData = new JSONArray();
            for (String pep : peptides) {
                innerArray = new JSONArray();
                innerArray.put(pep);
                peptideData.put(innerArray);
            }
            inputsArray.put(inputObject("peptide_sequences", "BYTES", peptideData));

            //the other inputs the model takes, in the order of the spec
            for (KoinaModelSpec.Input input : spec.inputs()) {
                JSONArray data = new JSONArray();
                switch (input) {
                    case PRECURSOR_CHARGES:
                        for (int charge : charges) {
                            innerArray = new JSONArray();
                            innerArray.put(charge);
                            data.put(innerArray);
                        }
                        break;
                    case COLLISION_ENERGIES:
                        for (float nce : nces) {
                            innerArray = new JSONArray();
                            innerArray.put(nce);
                            data.put(innerArray);
                        }
                        break;
                    case INSTRUMENT_TYPES:
                        for (String instrument : instruments) {
                            innerArray = new JSONArray();
                            innerArray.put(instrument);
                            data.put(innerArray);
                        }
                        break;
                    case FRAGMENTATION_TYPES:
                        for (String fragmentation : fragmentations) {
                            innerArray = new JSONArray();
                            innerArray.put(fragmentation);
                            data.put(innerArray);
                        }
                        break;
                }
                inputsArray.put(inputObject(input.koinaName(), input.datatype(), data));
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
        }
        return jsonOutFolder;
    }

    //one input of the request: a column of one value per precursor
    private JSONObject inputObject(String name, String datatype, JSONArray data) {
        JSONObject inputObject = new JSONObject();
        inputObject.put("name", name);
        inputObject.put("shape", new JSONArray("[" + peptides.length + ",1]"));
        inputObject.put("datatype", datatype);
        inputObject.put("data", data);
        return inputObject;
    }
}
