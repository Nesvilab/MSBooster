package bestmodelsearch;

import allconstants.LowercaseModelMapper;
import allconstants.ModelCollections;
import koinaclasses.KoinaException;
import koinaclasses.KoinaModelRegistry;
import koinaclasses.KoinaPreCheck;

import java.util.ArrayList;
import java.util.HashMap;

import static allconstants.ModelCollections.*;
import static utils.Print.printError;
import static utils.Print.printInfo;

public class ModelCollectionDecider {
    public static ArrayList<String> decideCollection(String mode) {
        String searchString = "";
        ArrayList<String> collection = null;

        switch (mode) {
            case "RT":
                searchString = rtSearchModelsString;
                collection = getRtCollection();
                break;
            case "MS2":
                searchString = ms2SearchModelsString;
                collection = getMs2Collection();
                break;
            case "IM":
                searchString = imSearchModelsString;
                collection = getImCollection();
                break;
        }

        printInfo("Searching for best " + mode + " model for your data");
        ArrayList<String> consideredModels = new ArrayList<>();
        if (! searchString.isEmpty()) {
            String[] searchModels = searchString.split(",");

            LowercaseModelMapper lmm = new LowercaseModelMapper();
            HashMap<String, String> lowercaseMapper = lmm.getLowercaseToModel();
            for (String model : searchModels) {
                String lowerModel = model.toLowerCase();
                if (lowercaseMapper.containsKey(lowerModel)) {
                    consideredModels.add(lowercaseMapper.get(lowerModel));
                } else {
                    printError("No " + mode + " model called " + model + ". Exiting.");
                    System.exit(1);
                }
            }
        } else {
            consideredModels.addAll(collection);
        }

        //a model that cannot encode the search's isobaric label could not be used, so it is not tried
        ArrayList<String> labelledOut = new ArrayList<>();
        for (String model : consideredModels) {
            if (!KoinaPreCheck.takesSearchLabel(model)) {
                labelledOut.add(model);
            }
        }
        if (!labelledOut.isEmpty()) {
            for (String model : labelledOut) {
                printInfo(KoinaPreCheck.isobaricLabelMessage(KoinaModelRegistry.require(model))
                        + ". It is left out of the search for the best " + mode + " model");
            }
            consideredModels.removeAll(labelledOut);
            if (consideredModels.isEmpty()) {
                throw new KoinaException("No " + mode + " model is left to search: "
                        + KoinaPreCheck.isobaricLabelMessage(KoinaModelRegistry.require(labelledOut.get(0))));
            }
        }
        printInfo("Searching the following models: " + consideredModels);
        return consideredModels;
    }
}
