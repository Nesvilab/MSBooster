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
import koinaclasses.KoinaModelRegistry;
import koinaclasses.KoinaModelSpec;
import peptideptmformatting.PeptideFormatter;
import predictions.PredictionEntryHashMap;

import java.util.Map;
import java.util.concurrent.ConcurrentHashMap;

import static utils.Print.printError;

public class KoinaLibReader implements LibraryPredictionMapper {
    public boolean failed = false;
    public PredictionEntryHashMap allPreds = new PredictionEntryHashMap();
    public String finalModel;
    //the model's spec, which PeptideFormatter, PeptideSkipper and transferKoinaPreds take
    public final KoinaModelSpec spec;
    public String property;
    //the precursors Koina refused (HTTP 400) when each was sent alone, in mass notation, with what Koina said; see
    //refuse
    private final Map<String, String> refused = new ConcurrentHashMap<>();

    public KoinaLibReader(String model) {
        spec = KoinaModelRegistry.find(model);
        if (spec == null) {
            printError(model + " not in Koina models");
            System.exit(1);
        }
        finalModel = model;

        //decide if this is an RT, MS2 or IM model
        property = spec.property(Constants.auxSpectraModel.equals(model));
    }

    public PredictionEntryHashMap getPreds() {return allPreds;}
    public void setPreds(PredictionEntryHashMap preds) {
        allPreds = preds;
    }

    public void clear() {
        allPreds.clear();
    }

    /**
     * Records a precursor Koina refused (HTTP 400) when it was sent alone (KoinaBisection). It has no prediction from
     * the model, and transferKoinaPreds scores it as a precursor the model cannot take.
     * <p>
     * It is recorded in mass notation, as the model's predictions are, so it is found whichever UniMod code the
     * peptide is written with: a model that remembers the codes it finds (PeptideFormatter.foundUnimods) may write a
     * mass with two codes one way in the request and the other way when transferKoinaPreds writes it again.
     *
     * @param modelPeptide the peptide as the model was sent it
     * @param charge       its charge as the request wrote it; ignored, and may be null, for a model that is sent no
     *                     charge
     * @param message      what Koina answered
     * @return whether it was not recorded before
     */
    public boolean refuse(String modelPeptide, String charge, String message) {
        PeptideFormatter pf = new PeptideFormatter(modelPeptide, charge == null ? "0" : charge, spec);
        return refused.putIfAbsent(key(pf), message) == null;
    }

    /**
     * Whether Koina refused this precursor ({@link #refuse}).
     *
     * @param modelPf the peptide as the model was sent it, read with the model's spec
     *                ({@link PeptideFormatter#PeptideFormatter(String, Object, KoinaModelSpec)}), and its charge
     */
    public boolean isRefused(PeptideFormatter modelPf) {
        return !refused.isEmpty() && refused.containsKey(key(modelPf));
    }

    /** Whether Koina refused this precursor ({@link #refuse}), the peptide as the model was sent it. */
    public boolean isRefused(String modelPeptide, String charge) {
        return !refused.isEmpty() && isRefused(new PeptideFormatter(modelPeptide, charge, spec));
    }

    /** How many precursors Koina refused ({@link #refuse}). */
    public int refusedCount() {
        return refused.size();
    }

    //an RT model is sent no charge, so Koina refuses its peptide whatever the charge
    private String key(PeptideFormatter modelPf) {
        return spec.sends(KoinaModelSpec.Input.PRECURSOR_CHARGES) ? modelPf.getBaseCharge() : modelPf.getBase();
    }
}
