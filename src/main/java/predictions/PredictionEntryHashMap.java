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

package predictions;

import allconstants.Constants;
import features.spectra.MassCalculator;
import peptideptmformatting.PeptideFormatter;
import peptideptmformatting.PeptideSkipper;

import static utils.Print.printError;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.Map;
import java.util.concurrent.ConcurrentHashMap;
import java.util.concurrent.ExecutionException;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Future;

public class PredictionEntryHashMap extends ConcurrentHashMap<String, PredictionEntry> {
    public String modelType;
    public String property;
    public void filterTopFragments(ExecutorService executorService)
            throws ExecutionException, InterruptedException {
        String[] peptides = new String[this.size()];
        PredictionEntry[] predictions = new PredictionEntry[this.size()];
        int entryIdx = 0;
        for (Entry<String, PredictionEntry> entry : this.entrySet()) {
            peptides[entryIdx] = entry.getKey();
            predictions[entryIdx] = entry.getValue();
            entryIdx++;
        }

        List<Future> futureList = new ArrayList<>(Constants.numThreads);
        for (int i = 0; i < Constants.numThreads; i++) {
            int start = (int) (this.size() * (long) i) / Constants.numThreads;
            int end = (int) (this.size() * (long) (i + 1)) / Constants.numThreads;
            futureList.add(executorService.submit(() -> {
                for (int j = start; j < end; j++) {
                    PredictionEntry pe = predictions[j];
                    pe.filterFragments();
                    this.put(peptides[j], pe);
                }
            }));
        }
        for (Future future : futureList) {
            future.get();
        }
    }

    public float getMaxPredRT() {
        float maxRT = 0f;
        for (PredictionEntry entry : this.values()) {
            if (entry.RT > maxRT) {
                maxRT = entry.RT;
            }
        }
        return maxRT;
    };

    //fullLib is allPreds, or the library all properties should be transfered onto
    //only need if reading in prepredicted libraries from different sources
    //example would be IM from DIANN, spectra and RT from Koina
    //for both files, they already have all the peptides needed, no transferring needed
    public void mergeIntoLibrary(PredictionEntryHashMap fullLib, String mode) {
        //switch case for spectra, RT, IM, aux
        for (Entry<String, PredictionEntry> entry : this.entrySet()) {
            String peptide = entry.getKey();
            PredictionEntry pe2 = entry.getValue();
            if (fullLib.containsKey(peptide)) {
                PredictionEntry pe1 = fullLib.get(peptide);
                switch (mode) {
                    //for RT and IM, just replace that
                    case "RT":
                        pe1.setRT(pe2.getRT());
                        break;
                    case "IM":
                        pe1.setIM(pe2.getIM());
                        break;
                    //for spectra and aux, probably need two separate vectors to be saved
                    case "spectra":
                        pe1.mzs = pe2.mzs;
                        pe1.intensities = pe2.intensities;
                        pe1.fragNums = pe2.fragNums;
                        pe1.charges = pe2.charges;
                        pe1.fragmentIonTypes = pe2.fragmentIonTypes;
                        pe1.flags = pe2.flags;
                        pe1.filtered = pe2.filtered;
                        break;
                    case "auxSpectra":
                        break;
                }

                fullLib.put(peptide, pe1); //entry in library1, but getting property for it from library2
            } else {
                fullLib.put(peptide, pe2); //if entry missing in library1, add it from library2
            }
        }

        //delete old one at end and save library
        clear();
    }

    //use this when we need don't already have predictions for all final peptides
    //models may have different PTM restrictions, so need to transfer their predictions onto the ones we actually need
    //for when msbooster is calling models for the first time, not when it is reading in prepredicted libraries
    //TODO: also need to consider model. Unispec and Predfull transfer differently
    //TODO: can save aby fragments separately from other ones
    public void transferKoinaPreds(ArrayList<PredictionEntryHashMap> predMaps, String fulltsv) throws IOException {
        //iterate through entries of full tsv
        BufferedReader TSVReader = new BufferedReader(new FileReader(fulltsv));
        String l;

        while ((l = TSVReader.readLine()) != null) {
            String[] line = l.split("\t");
            String baseCharge = line[0] + "|" + line[1];

            //make new prediction entry
            PredictionEntry newPred = new PredictionEntry();

            //for each peptide, get the translation from each pred model
            for (PredictionEntryHashMap predMap : predMaps) {
                PredictionEntry oldPred;
                float[] newMZs = new float[0];
                boolean skipPeptide = false;
                String modelSpecificBaseCharge = "";

                if (predMap.containsKey(baseCharge)) {
                    oldPred = predMap.get(baseCharge);
                } else { //need to translate
                    PeptideFormatter pf = null;
                    oldPred = null;
                    switch (predMap.modelType) {
                        case "alphapept":
                        case "ms2pip":
                        case "deeplc":
                        case "unispec":
                        case "prosit":
                        case "prosittmt":
                        case "predfull":
                            pf = new PeptideFormatter(
                                    new PeptideFormatter(line[0], line[1], "base").getModel(predMap.modelType),
                                    line[1], predMap.modelType);
                            modelSpecificBaseCharge = pf.getBaseCharge();
                            oldPred = predMap.get(modelSpecificBaseCharge);
                            break;
                        default:
                            printError(predMap.modelType + " not supported by Koina");
                            System.exit(1);
                    }

                    if (oldPred == null) {
                        //check if missing
                        if (PeptideSkipper.skipPeptide(pf.getStripped(), pf.getCharge(), predMap.modelType)) {
                            skipPeptide = true;
                        } else {
                            printError("Missing peptide to transfer prediction onto " + l + ": " + modelSpecificBaseCharge);
                            printError("Exiting now.");
                            System.exit(1);
                        }
                    } else { //don't want to worry about calculating all the different NLs
                        if (predMap.property.equals("ms2")) {
                            MassCalculator mc = new MassCalculator(line[0], line[1]);
                            newMZs = new float[oldPred.getMzs().length];

                            MassCalculator oldMc = new MassCalculator(pf.getBase(), pf.getCharge());

                            for (int i = 0; i < newMZs.length; i++) {
                                newMZs[i] = oldPred.getMzs()[i] +
                                        mc.compareModMasses(oldMc, oldPred.getFragNums()[i],
                                                oldPred.getFragmentIonTypes()[i],
                                                oldPred.getCharges()[i], oldPred.getFullAnnotations()[i]);
                            }
                        }
                    }
                }

                switch (predMap.property) {
                    case "ms2":
                        if (skipPeptide) {
                            newPred.mzs = new float[]{0};
                            newPred.intensities = new float[]{0};
                            newPred.fragNums = new int[]{0};
                            newPred.charges = new int[]{0};
                            newPred.fragmentIonTypes = new String[]{"y"};
                            newPred.flags = new int[]{1};
                        } else {
                            if (newMZs.length == 0) {
                                newMZs = oldPred.mzs;
                            }
                            newPred = new PredictionEntry(newMZs, oldPred.intensities, oldPred.fragNums,
                                    oldPred.charges, oldPred.fragmentIonTypes, oldPred.flags, oldPred.fullAnnotations);
                        }
                        break;
                    case "rt":
                        if (skipPeptide) {
                            newPred.setRT(0);
                        } else {
                            newPred.setRT(oldPred.RT);
                        }
                        break;
                    case "im":
                        if (skipPeptide) {
                            newPred.setIM(0);
                        } else {
                            newPred.setIM(oldPred.IM);
                        }
                        break;
                    case "auxms2":
                        break;
                    default:
                        printError(predMap.property + " is not supported. Exiting");
                        System.exit(1);
                }
            }

            baseCharge = baseCharge.replace("cterm", "");
            this.put(baseCharge, newPred);
        }
    }
}
