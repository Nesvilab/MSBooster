package Features;

import External.ExternalModelCaller;

import java.lang.reflect.Field;
import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Map;

//this is what I use in the java jar file
public class MainClass {
    public static void main(String[] args) throws Exception {
        //to do: take constants as input file

        //defining num threads
        Runtime run  = Runtime.getRuntime();
        if (Constants.numThreads <= 0) {
            Constants.numThreads = run.availableProcessors();
        } //otherwise use user-defined
        System.out.println("Using " + Constants.numThreads + " threads");

        //accept command line inputs
        HashSet<String> fields = new HashSet<>();
        for (Field f : Constants.class.getDeclaredFields()) {
            fields.add(f.getName());
        }

        HashMap<String, String> params = new HashMap<String, String>();

        //setting new values
        for (int i = 0; i < args.length; i++) {
            String key = args[i].substring(2); //remove --
            i++;
            params.put(key, args[i]);
        }

        //adding to constants class
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

        //check that at least pinPepXMLDirectory and mzmlDirectory are provided
        if (Constants.pinPepXMLDirectory == null) {
            throw new IllegalArgumentException("pinPepXMLDirectory must be provided");
        }
        if (Constants.mzmlDirectory == null) {
            throw new IllegalArgumentException("mzmlDirectory must be provided");
        }
        c.updatePaths();

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

        //create file for spectral and RT prediction
        //ignore if files already created
        boolean createSpectraRTPredFile = false;
        boolean createDetectPredFile = false;

        //check which ones we need
        if (allFeatures || autoFeatures) {
            createSpectraRTPredFile = true;
            createDetectPredFile = true;
        } else {
            HashSet<String> featureSet = new HashSet<>(Arrays.asList(featuresArray));
            featureSet.retainAll(Constants.spectraRTFeatures);
            if (featureSet.size() > 0) {
                createSpectraRTPredFile = true;
            }

            featureSet = new HashSet<>(Arrays.asList(featuresArray));
            featureSet.retainAll(Constants.detectFeatures);
            if (featureSet.size() > 0) {
                createDetectPredFile = true;
            }
        }

        //generate files for prediction models
        if (createSpectraRTPredFile) {
            if (Constants.spectraRTPredModel.equals("DIA-NN")) {
                if (Constants.DiaNN == null) {
                    throw new IllegalArgumentException("path to DIA-NN executable must be provided");
                }
                System.out.println("Generating input file for DIA-NN");
                peptideFileCreator.createPeptideFile(Constants.pinPepXMLDirectory, Constants.spectraRTPredInput, "Diann");
                ExternalModelCaller.callModel(run, "DIA-NN");
            } else if (Constants.spectraRTPredModel.equals("pDeep3")) {
                System.out.println("Generating input file for pDeep3");
                peptideFileCreator.createPeptideFile(Constants.pinPepXMLDirectory, Constants.spectraRTPredInput, "pDeep3");
                //TODO: run pDeep3 and add to ExternalModelCaller
            }
        }
        if (createDetectPredFile) {
            System.out.println("Generating input file for DeepMSPeptide");
            peptideFileCreator.createPeptideFile(Constants.pinPepXMLDirectory, Constants.detectPredInput, "DeepMSPeptide");
            //TODO: run DeepMSPeptide and add to ExternalModelCaller
        }

        //create new pin file with features
        System.out.println("Generating edited pin with following features: " + Arrays.toString(featuresArray));
        percolatorFormatter.editPin(Constants.pinPepXMLDirectory, Constants.mzmlDirectory, Constants.spectraRTPredFile,
                Constants.detectPredFile, featuresArray, Constants.editedPin);

        //TODO: how to deal with auto
    }
}