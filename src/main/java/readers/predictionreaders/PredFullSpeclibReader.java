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

import allconstants.FragmentIonConstants;
import features.spectra.MassCalculator;
import peptideptmformatting.PeptideFormatter;
import peptideptmformatting.PeptideSkipper;
import predictions.PredictionEntry;
import readers.MgfFileReader;
import readers.datareaders.PinReader;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.net.URISyntaxException;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.Set;
import java.util.concurrent.ExecutorService;

import static features.spectra.MassCalculator.allNeutralLossMasses;
import static utils.Print.printInfo;

public class PredFullSpeclibReader extends MgfFileReader {
    HashSet<String> peptides = new HashSet<>();

    public PredFullSpeclibReader(String file, boolean createScanNumObjects, File[] pinFiles,
                                 ExecutorService executorService)
            throws IOException, URISyntaxException {
        //initially start with mgf reading
        super(file, createScanNumObjects, executorService, "PredFull");

        //read in peptides from pinFiles
        for (File pinFile : pinFiles) {
            PinReader pin = new PinReader(pinFile.getCanonicalPath());

            //add to counter
            while (pin.next(true)) {
                peptides.add(pin.getPep().getBaseCharge());
            }
            pin.close();
        }

        printInfo("Shifting peaks for PredFull predictions");
        //also need to add to allPreds the peptides with PTMs not supported by PredFull
        //this requires that mgf file name has same prefix as full input file
        String fullFile = file.substring(0, file.length() - 4) + "_full.tsv";

        //for each peptide, check if base contain PTMs not supported by predfull
        //if so, need to start shifting
        //can do this by enumerating all possible fragment m/z (method in mass calculator)
        BufferedReader TSVReader = new BufferedReader(new FileReader(fullFile));
        String l;

        //filter here so don't need to annotate everything
        Set<String> ignoredFragmentIonTypesSet = FragmentIonConstants.makeIgnoredFragmentIonTypes();
        this.getPreds();
        while ((l = TSVReader.readLine()) != null) {
            //doing fragment annotation for everything, not just modified ones
            //start shifting and annotating
            //first, get unmodified peptide prediction from dictionary
            String[] lSplit = l.split("\t");

            //see if this peptide is in this dataset. If not, don't process it
            PeptideFormatter formattedPeptide = new PeptideFormatter(lSplit[0], lSplit[1], "base");
            if ((peptides.size() != 0) & (!peptides.contains(formattedPeptide.getBaseCharge()))) {
                continue;
            }

            PeptideFormatter peptideToSearch = new PeptideFormatter(formattedPeptide.getPredfull(),
                    lSplit[1], "predfull");

            PredictionEntry pe = this.allPredsHashMap.get(peptideToSearch.getBaseCharge());

            if (pe == null) { //valid reasons to be empty
                if (PeptideSkipper.skipPeptide(peptideToSearch, "predfull")) {
                    continue;
                }
            }
            //for each fragment ion, see if it can be annotated
            MassCalculator mc = new MassCalculator(peptideToSearch.getBase(), peptideToSearch.getCharge());
            //annotate ion
            //any calculated annotated fragment is considered
            String[][] info = mc.annotateMZs(pe.mzs, "default", true);
            String[] annotations = info[0];
            String[] fragmentIonTypes = info[1];

            //for annotated fragment ions, shift their predicted peak
            //shift only if including the amino acid
            //create new massCalculator object, and change m/z to new one
            //give String in annotation to calcMass and round to closest tenth of m/z
            //adjust in newMZs
            MassCalculator shiftedMC = new MassCalculator(lSplit[0], lSplit[1]);
            int index = 0;
            ArrayList<Float> finalMZs = new ArrayList<>();
            ArrayList<Float> finalIntensities = new ArrayList<>();
            ArrayList<String> finalFragmentIonTypes = new ArrayList<>();
            for (int j = 0; j < annotations.length; j++) {
                //skip fragment ion if of ignored type
                if (ignoredFragmentIonTypesSet.contains(fragmentIonTypes[j])) {
                    index += 1;
                    continue;
                }

                String fragment = annotations[j];
                Set<Float> newMZ = new HashSet<>();

                //possibly multiple identities, and only accept fragment if all or none are shifted
                String[] fragments = fragment.split(";");

                for (String frag : fragments) {
                    //split by + to get change. If no split, assume charge 1
                    String[] chargeSplit = frag.split("\\+");
                    int charge = 1;
                    if (chargeSplit.length == 2) {
                        charge = Integer.parseInt(chargeSplit[1]);
                    }

                    //split by - to get neutral loss. If no split, assume no neutral loss
                    String[] nlSplit = chargeSplit[0].split("-");
                    String nl = "";
                    if (nlSplit.length == 2) {
                        nl = nlSplit[1];
                    }

                    //switch case (MH, immonium, extract letters and numbers)
                    String[] nlSplitSplit = nlSplit[0].split(":");
                    String ionName = nlSplitSplit[0];
                    switch (ionName) {
                        case "p":
                            newMZ.add(shiftedMC.calcMass(shiftedMC.modMasses.size() - 2, "y", charge, allNeutralLossMasses.get(nl), 0));
                            break;
                        case "imm":
                            //need to check if mox of modified and unmodified, and exclude if so
                            //don't know if both mod and unmod m/z will be that intensity, safer to exclude
                            Set<Double> mods = new HashSet<>();
                            char aa = nlSplitSplit[1].charAt(0);
                            for (int i = 0; i < shiftedMC.peptide.length(); i++) {
                                if (shiftedMC.peptide.charAt(i) == aa) {
                                    mods.add(shiftedMC.modMasses.get(i + 1));
                                }
                            }
                            for (int i = 0; i < mods.size(); i++) {
                                newMZ.add((float) (mods.iterator().next() + shiftedMC.AAmap.get(aa) - 26.99));
                            }
                            break;
                        case "unknown":
                            newMZ.add(pe.mzs[index]);
                            break;
                        default: //backbone and internal
                            if (ionName.contains("dot")) {
                                newMZ.add(shiftedMC.calcMass(Integer.parseInt(ionName.substring(ionName.length() - 1)),
                                        ionName.substring(0, ionName.length() - 1), charge, 0));
                                break;
                            }
                            StringBuilder sb = new StringBuilder();
                            ArrayList<Integer> letterPositions = new ArrayList<>();
                            ArrayList<Integer> numbers = new ArrayList<>();
                            letterPositions.add(ionName.length());
                            //get letters first
                            for (int i = ionName.length() - 2; i > -1; i--) { //last char will always be number
                                char c = ionName.charAt(i);
                                if (Character.isLetter(c)) {
                                    sb.append(c);
                                    letterPositions.add(i);
                                }
                            }
                            //get numbers
                            for (int i = letterPositions.size() - 1; i > 0; i--) {
                                numbers.add(Integer.parseInt(ionName.substring(letterPositions.get(i) + 1, letterPositions.get(i - 1))));
                            }
                            if (numbers.size() == 1) {
                                newMZ.add(shiftedMC.calcMass(numbers.get(0), sb.toString(), charge, allNeutralLossMasses.get(nl), 0));
                            } else {
                                newMZ.add(shiftedMC.calcMass(numbers.get(0), numbers.get(1), sb.toString(), allNeutralLossMasses.get(nl), 0));
                            }
                            break;
                    }
                }
                if (newMZ.size() == 1) { //problem if one possibility is shifted and other isn't
                    finalMZs.add(newMZ.iterator().next());
                    finalIntensities.add(pe.intensities[index]);
                    finalFragmentIonTypes.add(fragmentIonTypes[index]);
                } else {
                    float totalFloat = 0f;
                    float compareFloat1 = newMZ.iterator().next();
                    totalFloat += compareFloat1;
                    boolean problem = false;
                    for (int i = 0; i < newMZ.size() - 1; i++) {
                        float compareFloat2 = newMZ.iterator().next();
                        totalFloat += compareFloat2;
                        if (Math.abs(compareFloat1 - compareFloat2) > 0.1f) {
                            problem = true;
                            break;
                        }
                    }
                    if (! problem) {
                        finalMZs.add(totalFloat / newMZ.size());
                        finalIntensities.add(pe.intensities[index]);
                        finalFragmentIonTypes.add(fragmentIonTypes[index]);
                    }
                }
                index += 1;
            }

            //now can assign new PredictionEntry
            float[] mzArray = new float[finalMZs.size()];
            float[] intArray = new float[finalMZs.size()];
            String[] fragmentArray = new String[finalMZs.size()];
            for (int i = 0; i < finalMZs.size(); i++) {
                mzArray[i] = finalMZs.get(i);
                intArray[i] = finalIntensities.get(i);
                fragmentArray[i] = finalFragmentIonTypes.get(i);
            }

            PredictionEntry newPred = new PredictionEntry(mzArray, intArray,
                    new int[0], new int[0], fragmentArray);
            newPred.setRT(pe.RT);
            newPred.setIM(pe.IM);
            this.allPredsHashMap.put(lSplit[0] + "|" + lSplit[1], newPred);
        }
    }
}
