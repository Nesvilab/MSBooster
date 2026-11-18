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

package peptideptmformatting;

import allconstants.FragCastCharges;
import koinaclasses.KoinaModelRegistry;
import koinaclasses.KoinaModelSpec;
import koinaclasses.KoinaModelSpec.Limits;

//whether a model cannot predict a precursor, so it is not sent to the model and its missing prediction is expected
public class PeptideSkipper {
    //the local Prosit and PrositTMT predictions (MspReader) take what the Koina Prosit models take
    private static final Limits LOCAL_PROSIT_LIMITS = Limits.NONE
            .withForbiddenResidues(Limits.NONSTANDARD_RESIDUES).withLength(0, 30).withCharge(0, 6);

    /** Whether a Koina model cannot take this precursor ({@link KoinaModelSpec#limits()}). */
    public static boolean skipPeptide(PeptideFormatter pf, KoinaModelSpec spec) {
        return skipReason(pf, spec) != null;
    }

    /**
     * Why a Koina model cannot take this precursor ({@link KoinaModelSpec#limits()}), such as "residue U" or
     * "longer than 30 residues", or null if it can.
     */
    public static String skipReason(PeptideFormatter pf, KoinaModelSpec spec) {
        Limits limits = spec.limits();
        //the peptide is written for the model before the other checks, so a model that remembers the UniMod codes
        //it finds (PeptideFormatter.foundUnimods) remembers those of a skipped peptide too
        String modelPeptide = limits.needsModelPeptide() ? pf.getModel(spec) : null;
        return limits.skipReason(baseWithoutCtermMod(pf), pf.getStripped(), Integer.parseInt(pf.getCharge()),
                modelPeptide);
    }

    //the base peptide without its C-term mod, which is written after the last residue ("LGEPTIDEC[-0.9840]" with
    //cterm), so that it is not read as a mod of a last C
    private static String baseWithoutCtermMod(PeptideFormatter pf) {
        String base = pf.getBase();
        if (pf.cterm && base.endsWith("]")) {
            return base.substring(0, base.lastIndexOf('['));
        }
        return base;
    }

    /**
     * Whether this model cannot take this precursor: a Koina model by its spec, the local Prosit and PrositTMT
     * predictions by the Prosit limits, and FragCast by the charges it encodes. Other models take every precursor.
     */
    public static boolean skipPeptide(PeptideFormatter pf, String model) {
        KoinaModelSpec spec = KoinaModelRegistry.find(model);
        if (spec != null) {
            return skipPeptide(pf, spec);
        }

        int chargeInt = Integer.parseInt(pf.getCharge());
        if (model.equalsIgnoreCase("Prosit") || model.equalsIgnoreCase("PrositTMT")) {
            return LOCAL_PROSIT_LIMITS.skipReason(pf.getStripped(), chargeInt, null) != null;
        }
        //FragCast encodes charge as a six-way one-hot and its build-library refuses a range past it,
        //so these precursors are never predicted and must be recognized as unsupported rather than
        //looked up and found missing
        if (model.toLowerCase().contains("fragcast")
                && (chargeInt < FragCastCharges.LOWEST || chargeInt > FragCastCharges.HIGHEST)) {
            return true;
        }
        return false;
    }

    /**
     * Whether any of these models cannot take this precursor, for a PSM whose prediction comes from several
     * models (spectra, RT and IM). An empty name is a property no model predicts.
     */
    public static boolean skipPeptideForAnyModel(PeptideFormatter pf, String... models) {
        for (String model : models) {
            if (model != null && !model.isEmpty() && skipPeptide(pf, model)) {
                return true;
            }
        }
        return false;
    }
}
