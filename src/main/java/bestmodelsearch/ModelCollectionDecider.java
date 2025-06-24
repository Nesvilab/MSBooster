package bestmodelsearch;

import allconstants.LowercaseModelMapper;
import allconstants.ModelCollections;

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
        if (! searchString.isEmpty()) {
            ArrayList<String> consideredModels = new ArrayList<>();
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
            printInfo("Searching the following models: " + consideredModels);
            return consideredModels;
        } else {
            printInfo("Searching the following models: " + collection);
            return collection;
        }
    }
}
