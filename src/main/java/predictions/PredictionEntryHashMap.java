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
import readers.predictionreaders.KoinaLibReader;
import utils.Multithreader;

import java.io.BufferedReader;
import java.io.FileReader;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;
import java.util.concurrent.ConcurrentHashMap;
import java.util.concurrent.ExecutionException;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Future;

import static utils.Print.printError;

public class PredictionEntryHashMap extends ConcurrentHashMap<String, PredictionEntry> {
    public PredictionEntryHashMap() { super(); }
    public PredictionEntryHashMap(int initialCapacity) { super(initialCapacity); }

    public void preprocessPredictedSpectra(ExecutorService executorService, HashSet<String> primaryTypes, HashSet<String> auxTypes)
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
        Multithreader mt = new Multithreader(this.size(), Constants.numThreads);
        for (int i = 0; i < Constants.numThreads; i++) {
            int finalI = i;
            futureList.add(executorService.submit(() -> {
                for (int j = mt.indices[finalI]; j < mt.indices[finalI + 1]; j++) {
                    PredictionEntry pe = predictions[j];
                    pe.preprocessFragments(primaryTypes, Constants.topFragments);
                    if (pe.auxSpectra != null) {
                        pe.auxSpectra.preprocessFragments(auxTypes, Constants.topFragments + 10);
                    }
                    this.put(peptides[j], pe);
                }
            }));
        }
        for (Future future : futureList) {
            future.get();
        }
    }

    public float getMaxPredRT() {
        float maxRT = Float.MIN_VALUE;
        for (PredictionEntry entry : this.values()) {
            if (entry.RT > maxRT) {
                maxRT = entry.RT;
            }
        }
        return maxRT;
    };

    //fullLib is allPreds, or the library all properties should be transfered onto
    //only need if reading in prepredicted libraries from different sources
    //example would be IM from FragPred, spectra and RT from Koina
    //for both files, they already have all the peptides needed, no transferring needed
    public void mergeIntoLibrary(PredictionEntryHashMap fullLib, String mode) {
        //switch case for spectra, RT, IM, aux
        for (Entry<String, PredictionEntry> entry : this.entrySet()) {
            String peptide = entry.getKey();
            PredictionEntry pe2 = entry.getValue(); //library being added
            if (fullLib.containsKey(peptide)) {
                PredictionEntry pe1 = fullLib.get(peptide); //allpreds
                switch (mode) {
                    //for RT and IM, just replace that
                    case "RT":
                        pe1.setRT(pe2.getRT());
                        break;
                    case "IM":
                        pe1.setIM(pe2.getIM());
                        break;
                    case "spectra":
                        pe1.fragments = pe2.fragments;
                        break;
                    case "auxSpectra":
                        if (pe1.numFragments() <= 1) { //populate missing main spectra
                            if (Constants.auxSpectraModel.equalsIgnoreCase("unispec")) {
                                //filter to use only 0 isotopes
                                ArrayList<Integer> zeroIsotopes = new ArrayList<>();
                                for (int i = 0; i < pe2.isotopes.length; i++) {
                                    if (pe2.isotopes[i] == 0) {
                                        zeroIsotopes.add(i);
                                    }
                                }

                                int sz = zeroIsotopes.size();
                                pe1.fragments = new int[PredictionEntry.STRIDE * sz];
                                pe1.fullAnnotations = pe2.fullAnnotations.length > 0 ? new String[sz] : pe2.fullAnnotations;
                                pe1.isotopes = pe2.isotopes.length > 0 ? new int[sz] : pe2.isotopes;

                                int i = 0;
                                for (int idx : zeroIsotopes) {
                                    pe1.fragments[PredictionEntry.STRIDE * i]     = pe2.fragments[PredictionEntry.STRIDE * idx];
                                    pe1.fragments[PredictionEntry.STRIDE * i + 1] = pe2.fragments[PredictionEntry.STRIDE * idx + 1];
                                    pe1.fragments[PredictionEntry.STRIDE * i + 2] = pe2.fragments[PredictionEntry.STRIDE * idx + 2];
                                    if (pe2.fullAnnotations.length > 0) {
                                        pe1.fullAnnotations[i] = pe2.fullAnnotations[idx];
                                    }
                                    if (pe2.isotopes.length > 0) {
                                        pe1.isotopes[i] = pe2.isotopes[idx];
                                    }
                                    i++;
                                }
                            } else { //PredFull
                                pe1.fragments = pe2.fragments;
                                pe1.fullAnnotations = pe2.fullAnnotations;
                                pe1.isotopes = pe2.isotopes;

                                //prediction entry completely from predfull, need to use dalton matching
                                pe1.daltonMatching = true;
                            }
                        }
                        pe1.auxSpectra = pe2;
                        break;
                }

                fullLib.put(peptide, pe1); //entry in library1, but getting property for it from library2
            } else { //first entry loaded
                //if running whole workflow, even if one model does not support all peptides, will have empty entry in mgf
                fullLib.put(peptide, pe2);
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
    public void transferKoinaPreds(ArrayList<KoinaLibReader> klrs, String fulltsv) throws Exception {
        //iterate through entries of full tsv
        BufferedReader TSVReader = new BufferedReader(new FileReader(fulltsv));
        String l;

        while ((l = TSVReader.readLine()) != null) {
            String[] line = l.split("\t");
            String baseCharge = line[0] + "|" + line[1];

            //make new prediction entry
            PredictionEntry newPred = new PredictionEntry();

            //for each peptide, get the translation from each pred model
            for (KoinaLibReader klr : klrs) {
                PredictionEntryHashMap predMap = klr.getPreds();
                PredictionEntry oldPred;
                float[] newMZs = new float[0];
                boolean skipPeptide = false;
                String modelSpecificBaseCharge = "";

                if (predMap.containsKey(baseCharge)) {
                    oldPred = predMap.get(baseCharge);
                } else { //need to translate
                    PeptideFormatter pf = null;
                    oldPred = null;
                    switch (klr.modelType) {
                        case "alphapept":
                        case "ms2pip":
                        case "deeplc":
                        case "im2deep":
                        case "unispec":
                        case "prosit":
                        case "prosittmt":
                        case "prosit_cit":
                        case "predfull":
                            pf = new PeptideFormatter(
                                    new PeptideFormatter(line[0], line[1], "base").getModel(klr.modelType),
                                    line[1], klr.modelType);
                            modelSpecificBaseCharge = pf.getBaseCharge();
                            oldPred = predMap.get(modelSpecificBaseCharge);
                            break;
                        default:
                            printError(klr.modelType + " not supported by Koina");
                            System.exit(1);
                    }

                    if (oldPred == null) {
                        //check if missing
                        if (PeptideSkipper.skipPeptide(pf, klr.modelType)) {
                            skipPeptide = true;
                        } else {
                            printError("Missing peptide to transfer prediction onto " + l + ": " + modelSpecificBaseCharge);
                            printError("Exiting now.");
                            System.exit(1);
                        }
                    } else { //don't want to worry about calculating all the different NLs
                        if (klr.property.equals("ms2") || klr.property.equals("ms2_aux")) {
                            MassCalculator mc = new MassCalculator(line[0], line[1]);
                            newMZs = new float[oldPred.numFragments()];

                            MassCalculator oldMc = new MassCalculator(pf.getBase(), pf.getCharge());

                            for (int i = 0; i < newMZs.length; i++) {
                                newMZs[i] = oldPred.getMz(i) +
                                        mc.compareModMasses(oldMc, oldPred.getFragNum(i),
                                                oldPred.getIonTypeString(i),
                                                oldPred.getCharge(i),
                                                oldPred.fullAnnotations.length > 0 ? oldPred.fullAnnotations[i] : "");
                            }
                        }
                    }
                }

                switch (klr.property) {
                    case "ms2":
                        if (skipPeptide) {
                            newPred.fragments = new int[PredictionEntry.STRIDE];
                            PredictionEntry.packSlot(newPred.fragments, 0, 0f, 0f, 0,
                                    allconstants.FragmentIonConstants.ION_INDEX.get("y"), 0);
                        } else {
                            if (newMZs.length == 0) { //this entry is already in old pred hashmap
                                newMZs = oldPred.getMzs();
                            }
                            if (klr.modelType.contains("unispec") || klr.modelType.contains("predfull")) {
                                newPred = new PredictionEntry(newMZs, oldPred.getIntensities(), oldPred.getFragNums(),
                                        oldPred.getCharges(), oldPred.getFragmentIonTypes(), oldPred.fullAnnotations);
                            } else {
                                newPred = new PredictionEntry(newMZs, oldPred.getIntensities(), oldPred.getFragNums(),
                                        oldPred.getCharges(), oldPred.getFragmentIonTypes());
                            }
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
                    case "ms2_aux":
                        PredictionEntry auxPe = new PredictionEntry();
                        if (skipPeptide) {
                            auxPe.fragments = new int[PredictionEntry.STRIDE];
                            PredictionEntry.packSlot(auxPe.fragments, 0, 0f, 0f, 0,
                                    allconstants.FragmentIonConstants.ION_INDEX.get("y"), 0);
                        } else {
                            if (newMZs.length == 0) {
                                newMZs = oldPred.getMzs();
                            }
                            if (klr.modelType.contains("unispec") || klr.modelType.contains("predfull")) {
                                auxPe = new PredictionEntry(newMZs, oldPred.getIntensities(), oldPred.getFragNums(),
                                        oldPred.getCharges(), oldPred.getFragmentIonTypes(), oldPred.fullAnnotations);
                            } else {
                                auxPe = new PredictionEntry(newMZs, oldPred.getIntensities(), oldPred.getFragNums(),
                                        oldPred.getCharges(), oldPred.getFragmentIonTypes());
                            }
                        }
                        newPred.auxSpectra = auxPe;
                        break;
                    default:
                        printError(klr.property + " is not supported. Exiting");
                        System.exit(1);
                }
            }

            baseCharge = baseCharge.replace("cterm", "");
            this.put(baseCharge, newPred);
        }
    }
}
