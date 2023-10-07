package Features;

import org.apache.commons.io.FileUtils;
import org.json.JSONArray;
import org.json.JSONObject;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.lang.reflect.Array;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashSet;
import java.util.List;
import java.util.concurrent.*;

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

    private HashSet<String> RTmodels = new HashSet<>(Arrays.asList("AlphaPept_rt_generic",
            "Prosit_2019_irt", "Prosit_2020_irt_TMT"));
    private HashSet<String> MS2models = new HashSet<>(Arrays.asList("ms2pip_2021_HCD",
            "AlphaPept_ms2_generic", "Prosit_2019_intensity", "Prosit_2020_intensity_CID",
            "Prosit_2020_intensity_TMT", "Prosit_2020_intensity_HCD", "Prosit_2023_intensity_TOF"));
    private HashSet<String> TMTmodels = new HashSet<>(Arrays.asList("Prosit_2020_irt_TMT",
            "Prosit_2020_intensity_TMT"));
    private final int maxJsonLength = 1000;

    public JSONWriter(String url, HashSet<String> entries) {
        this.url = url;
        this.model = url.toLowerCase().split("_")[0];
        if (RTmodels.contains(url)) {
            property = "rt";
        } else if (MS2models.contains(url)) {
            property = "ms2";
        }
        if (TMTmodels.contains(url)) {
            TMT = true;
        }

        //set entries
        this.peptides = new String[entries.size()];
        this.charges = new int[entries.size()];
        this.nces = new float[entries.size()];
        this.instruments = new String[entries.size()];
        this.fragmentations = new String[entries.size()];
        int i = 0;
        for (String entry : entries) {
            String[] info = entry.split(",");
            this.peptides[i] = info[0];
            this.charges[i] = Integer.parseInt(info[1]);
            this.nces[i] = Float.parseFloat(info[2]);
            this.instruments[i] = info[3];
            this.fragmentations[i] = info[4];
            i++;
        }

        numFiles = (int) Math.ceil((float) this.peptides.length / (float) maxJsonLength);
    }

    //this version if we need to spawn smaller ones
    public JSONWriter(JSONWriter parent, int iteration) {
        this.iteration = iteration;
        int start = iteration * maxJsonLength;
        this.url = parent.url;
        this.model = url.toLowerCase().split("_")[0];
        if (RTmodels.contains(url)) {
            property = "rt";
        } else if (MS2models.contains(url)) {
            property = "ms2";
        }
        if (TMTmodels.contains(url)) {
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

    public String write(boolean createDir) throws IOException, ExecutionException, InterruptedException {
        String jsonOutFolder = Constants.outputDirectory + File.separator + "jsonFiles";
        if (createDir) {
            if (Files.exists(Paths.get(jsonOutFolder))) {
                FileUtils.deleteDirectory(new File(jsonOutFolder));
            }
            Files.createDirectories(Paths.get(jsonOutFolder));
        }

        String fileName = "";
        if (numFiles > 1) {
            futureList.clear();
            ExecutorService executorService = Executors.newFixedThreadPool(Constants.numThreads);
            for (int i = 0; i < Constants.numThreads; i++) {
                int start = (int) (numFiles * (long) i) / Constants.numThreads;
                int end = (int) (numFiles * (long) (i + 1)) / Constants.numThreads;
                futureList.add(executorService.submit(() -> {
                    for (int rep = start; rep < end; rep ++) {
                        JSONWriter jw = new JSONWriter(this, rep);
                        try {
                            jw.write(false);
                        } catch (IOException | ExecutionException | InterruptedException e) {
                            e.printStackTrace();
                        }
                    }
                }));
            }
            for (Future future : futureList) {
                future.get();
            }
            executorService.shutdown();
        } else {
            switch (model) {
                case "alphapept":
                case "prosit":
                case "ms2pip":
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

                    if (TMT) {
                        for (int i = 0; i < peptides.length; i++) {
                            peptides[i] = "[UNIMOD:737]-" + peptides[i];
                        }
                    }
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
                                "spectraRT_" + url + iteration + ".json";
                    } else {
                        fileName = jsonOutFolder + File.separator +
                                "spectraRT_" + url + ".json";
                    }
                    try (FileWriter fileWriter = new FileWriter(fileName)) {
                        fileWriter.write(jsonData.toString(4)); // Using 4 for indentation
                    } catch (IOException e) {
                        e.printStackTrace();
                    }
                    break;
                default:
                    System.out.println(model + " is not supported by Koina. Exiting now.");
                    System.exit(0);
            }
        }
        return jsonOutFolder;
    }
}
