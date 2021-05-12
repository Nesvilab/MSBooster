package Features;

import External.ExternalModelCaller;

import java.io.*;
import java.lang.reflect.Field;
import java.lang.reflect.Modifier;
import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Map;

//this is what I use in the java jar file
public class MainClass {
    public static void main(String[] args) throws Exception {
        //to do: take constants as input file

        //accept command line inputs
        HashSet<String> fields = new HashSet<>();
        for (Field f : Constants.class.getDeclaredFields()) {
            fields.add(f.getName());
        }

        HashMap<String, String> params = new HashMap<String, String>();

        //setting new values
        //TODO: description of all flags. Make pinDirectory and mzmlDirectory positional arguments?
        for (int i = 0; i < args.length; i++) {
            String key = args[i].substring(2); //remove --
            if (key.equals("help")) { //help message
                System.out.println("Usage: java -jar MSFraggerDIA_postprocess-1.0-SNAPSHOT-jar-with-dependencies.jar [flags]");
                System.exit(0);
            }
            i++;
            params.put(key, args[i]);
        }

        //adding to constants class
        if (params.containsKey("paramsList")) { //override previous params input
                                                //params with nulls are left out
            String line;
            BufferedReader reader = new BufferedReader(new FileReader(params.get("paramsList")));
            while ((line = reader.readLine()) != null) {
                String[] lineSplit = line.split("=", 2);
                params.put(lineSplit[0].trim(), lineSplit[1].trim());
            }
            reader.close();
        }

        Constants c = new Constants();
        for (Map.Entry<String, String> entry : params.entrySet()) {
            String key = entry.getKey();
            if (! fields.contains(key)) {
                throw new Exception(entry.getKey() + " is not a valid parameter");
            } else {
                //get class of field
                Field field = Constants.class.getField(key);
                Class<?> myClass = field.getType();

                //parse to appropriate type
                field.set(c, myClass.getConstructor(String.class).newInstance(entry.getValue()));
            }
        }

        //defining num threads
        Runtime run  = Runtime.getRuntime();
        if (Constants.numThreads <= 0) {
            Constants.numThreads = run.availableProcessors();
        } //otherwise use user-defined
        System.out.println("Using " + Constants.numThreads + " threads");

        //check that at least pinPepXMLDirectory and mzmlDirectory are provided
        if (Constants.pinPepXMLDirectory == null) {
            throw new IllegalArgumentException("pinPepXMLDirectory must be provided");
        }
        if (Constants.mzmlDirectory == null) {
            throw new IllegalArgumentException("mzmlDirectory must be provided");
        }

        //check that features are allowed
        boolean allFeatures = false;
        boolean autoFeatures = false;
        String[] featuresArray = Constants.features.split(",");
        for (String f : featuresArray) {
            if (! Constants.allowedFeatures.contains(f)) {
                throw new IllegalArgumentException(f + " is not an allowed feature. " +
                        "Please choose from the following: " + Constants.allowedFeatures);
            }
            if (f.equals("all")) {
                allFeatures = true;
                featuresArray = (String[]) Constants.allowedFeatures.toArray();
                break;
            }
            if (f.equals("auto")) {
                autoFeatures = true;
                break;
            }
        }
        HashSet<String> featureSet = new HashSet<>(Arrays.asList(featuresArray));

        //if detectFractionGreater, need fasta
        boolean createDetectAllPredFile = false;
        if (featureSet.contains("detectFractionGreater") || autoFeatures) {
            if (Constants.fasta == null) {
                throw new IllegalArgumentException("Using current combination of features, " +
                        "detectFractionGreater is calculated and needs a fasta provided using " +
                        "--fasta <fasta file location>");
            }
            //need to have another boolean that supersedes regular createDetectPredFile
            createDetectAllPredFile = true;
        }

        //create file for spectral and RT prediction
        //ignore if files already created
        boolean createSpectraRTPredFile = false;
        boolean createDetectPredFile = false;
        boolean createSpectraRTPredFile2 = false;
        boolean createDetectPredFile2 = false;

        //check which ones we need
        if (allFeatures || autoFeatures) {
            createSpectraRTPredFile = true;
            createDetectPredFile = true;
            createSpectraRTPredFile2 = true;
            createDetectPredFile2 = true;
        } else {
            featureSet.retainAll(Constants.spectraRTFeatures);
            if (featureSet.size() > 0) {
                createSpectraRTPredFile = true;
                createSpectraRTPredFile2 = true;
            }

            featureSet = new HashSet<>(Arrays.asList(featuresArray));
            featureSet.retainAll(Constants.detectFeatures);
            if (featureSet.size() > 0) {
                createDetectPredFile = true;
                createDetectPredFile2 = true;
            }
        }
        //overriding if intermediate files already made
        if (Constants.spectraRTPredInput != null || Constants.spectraRTPredFile != null) {
            createSpectraRTPredFile = false;
        }
        if (Constants.detectPredInput != null || Constants.detectPredFile != null) {
            createDetectPredFile = false;
            createDetectAllPredFile = false;
        }

        c.updatePaths(); //setting null paths

        //generate files for prediction models
        if (createSpectraRTPredFile) {
            if (Constants.spectraRTPredModel.equals("DIA-NN")) {
                if (Constants.DiaNN == null) {
                    throw new IllegalArgumentException("path to DIA-NN executable must be provided");
                }
                System.out.println("Generating input file for DIA-NN");
                peptideFileCreator.createPeptideFile(Constants.pinPepXMLDirectory, Constants.spectraRTPredInput, "Diann", "pin");
            } else if (Constants.spectraRTPredModel.equals("pDeep3")) {
                System.out.println("Generating input file for pDeep3");
                peptideFileCreator.createPeptideFile(Constants.pinPepXMLDirectory, Constants.spectraRTPredInput, "pDeep3", "pin");
            }
        }
        if (createDetectAllPredFile) {
            System.out.println("Generating input file for DeepMSPeptide");
            peptideFileCreator.createPeptideFile(Constants.pinPepXMLDirectory, Constants.detectPredInput, "DeepMSPeptideAll", "pin");
        } else if (createDetectPredFile) {
            System.out.println("Generating input file for DeepMSPeptide");
            peptideFileCreator.createPeptideFile(Constants.pinPepXMLDirectory, Constants.detectPredInput, "DeepMSPeptide", "pin");
        }

        //generate predictions
        if ((Constants.spectraRTPredFile == null) && (createSpectraRTPredFile2)) {
            ExternalModelCaller.callModel(run, "DIA-NN");
        }
        if ((Constants.detectPredFile == null) && (createDetectPredFile2)) {
            ExternalModelCaller.callModel(run, "DeepMSPeptide");
        }

        //print new params file
        try {
            BufferedWriter buffer = new BufferedWriter(new FileWriter(Constants.outputDirectory +
                    File.separator + "finalParams.txt"));

            Field[] f = Constants.class.getDeclaredFields();
            for (Field field : f) {
                if ((field.getModifiers() & Modifier.FINAL) != Modifier.FINAL) {
                    if (!field.getName().equals("paramsList")) {
                        buffer.write(field.getName() + " = " + field.get(c) + "\n");
                    }
                }
            }
            buffer.close();
        } catch (Exception e) {
            e.getStackTrace();
        }

        //create new pin file with features
        System.out.println("Generating edited pin with following features: " + Arrays.toString(featuresArray));
        percolatorFormatter.editPin(Constants.pinPepXMLDirectory, Constants.mzmlDirectory, Constants.spectraRTPredFile,
                Constants.detectPredFile, featuresArray, Constants.editedPin);

        //TODO: how to deal with auto
    }
}
