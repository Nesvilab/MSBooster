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

package readers.predictionreaders;

import allconstants.Constants;
import allconstants.ModelCollections;
import predictions.PredictionEntryHashMap;

import static utils.Print.printError;

public class KoinaLibReader implements LibraryPredictionMapper {
    public boolean failed = false;
    public PredictionEntryHashMap allPreds = new PredictionEntryHashMap();
    public String finalModel;
    public String modelType; //TODO: may need better way to handle this i.e. prosit and prosit cit
    public String property;
    public boolean useFullAnnotation = false;

    public KoinaLibReader(String model) {
        finalModel = model;
        modelType = model.toLowerCase().split("_")[0];
        if (modelType.equals("prosit") && model.contains("TMT")) {
            modelType = "prosittmt";
        } else if (modelType.equals("prosit") && model.contains("_cit")) {
            modelType = "prosit_cit";
        }

        if (modelType.contains("unispec") || modelType.contains("predfull")) {
            useFullAnnotation = true;
        }

        //decide if this is RT or MS2 model
        if (ModelCollections.KoinaRTmodels.contains(model)) {
            property = "rt";
        } else if (ModelCollections.KoinaMS2models.contains(model)) {
            if (Constants.auxSpectraModel.equals(model)) {
                property = "ms2_aux";
            } else {
                property = "ms2";
            }
        } else if (ModelCollections.KoinaIMmodels.contains(model)) {
            property = "im";
        } else {
            printError(model + " not in Koina models");
            System.exit(1);
        }
    }

    public PredictionEntryHashMap getPreds() {return allPreds;}
    public void setPreds(PredictionEntryHashMap preds) {
        allPreds = preds;
    }

    public void clear() {
        allPreds.clear();
    }
}
