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
