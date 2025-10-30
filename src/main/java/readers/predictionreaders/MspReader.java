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
import features.spectra.MassCalculator;
import peptideptmformatting.PeptideFormatter;
import peptideptmformatting.PeptideSkipper;
import predictions.PredictionEntry;
import predictions.PredictionEntryHashMap;

import static peptideptmformatting.PTMhandler.prositAAMods;
import static utils.Print.printError;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;

public class MspReader implements LibraryPredictionMapper {
    final ArrayList<String> filenames;
    PredictionEntryHashMap allPreds = new PredictionEntryHashMap();

    public MspReader(String files) {
        File predsDirectory = new File(files);
        String[] predsFiles = predsDirectory.list();
        filenames = new ArrayList<String>();

        if (predsFiles == null) { //if user provided a file, not a directory
            filenames.add(files);
        } else { //user provides directory
            for (String predsFile : predsFiles) {
                if (predsFile.contains(".msp")) {
                    filenames.add(files + File.separator + predsFile);
                }
            }
        }

        String line;
        for (String f : filenames) {
            try (BufferedReader msp = new BufferedReader(new FileReader(f))) {

                while ((line = msp.readLine()) != null) { //Name
                    StringBuilder sb = new StringBuilder();
                    String[] lineSplit = line.split(": ")[1].split("/");
                    sb.append(lineSplit[0]).append("|");
                    String charge = lineSplit[1];

                    msp.readLine(); //MW
                    line = msp.readLine(); //comment
                    lineSplit = line.split(" ");

                    String pep;
                    float RT = 0f;

                    for (String s : lineSplit) {
                        String[] commentContents = s.split("=");
                        switch (commentContents[0]) {
                            case "Mods":
                                sb.append(commentContents[1]);
                                break;
                            case "iRT":
                                RT = Float.parseFloat(commentContents[1]);
                                break;
                        }
                    }
                    pep = new PeptideFormatter(sb.toString(), charge, "msp").getBaseCharge();

                    //fragments
                    //if it keeps reporting 174 peaks, set that as initial arraylist size
                    line = msp.readLine(); //Num peaks
                    int numPeaks = Integer.parseInt(line.split(" ")[2]);
                    ArrayList<Float> mzsList = new ArrayList<>(numPeaks);
                    ArrayList<Float> intsList = new ArrayList<>(numPeaks);
                    ArrayList<Integer> fragNumsList = new ArrayList<>(numPeaks);
                    ArrayList<String> fragmentIonTypesList = new ArrayList<>(numPeaks);
                    ArrayList<Integer> chargesList = new ArrayList<>(numPeaks);

                    for (int i = 0; i < numPeaks; i++) {
                        line = msp.readLine();
                        lineSplit = line.split("\t");
                        float tmpInt = Float.parseFloat(lineSplit[1]);
                        if (tmpInt == 0f) {
                            continue;
                        }
                        String fragment = lineSplit[2].substring(1, lineSplit[2].length() - 1).split("/")[0];
                        String[] chargeSplit = fragment.split("\\^");

                        //add fragment
                        int fragCharge = 1;
                        if (chargeSplit.length != 1) {
                            fragCharge = Integer.parseInt(chargeSplit[1]);
                        }
                        if (fragCharge > Constants.maxPredictedFragmentCharge) {
                            continue;
                        }
                        int fragNum = Integer.parseInt(chargeSplit[0].substring(1));
                        if (fragNum < Constants.minPredictedFragmentNum) {
                            continue;
                        }

                        mzsList.add(Float.parseFloat(lineSplit[0]));
                        intsList.add(tmpInt);
                        chargesList.add(fragCharge);
                        fragmentIonTypesList.add(fragment.substring(0, 1));
                        fragNumsList.add(fragNum);
                    }

                    float[] mzs = new float[mzsList.size()];
                    float[] ints = new float[intsList.size()];
                    int[] fragNums = new int[intsList.size()];
                    String[] fragmentIonTypes = new String[intsList.size()];
                    int[] charges = new int[intsList.size()];
                    for (int i = 0; i < mzs.length; i++) {
                        mzs[i] = mzsList.get(i);
                        ints[i] = intsList.get(i);
                        fragNums[i] = fragNumsList.get(i);
                        fragmentIonTypes[i] = fragmentIonTypesList.get(i);
                        charges[i] = chargesList.get(i);
                    }

                    PredictionEntry newPred = new PredictionEntry(mzs, ints,
                            fragNums, charges, fragmentIonTypes);
                    newPred.setRT(RT);
                    newPred.setIM(0f);
                    allPreds.put(pep, newPred);
                }

                //read in extra files
                BufferedReader TSVReader = new BufferedReader(new FileReader(
                        f.substring(0, f.length() - 4) + "_full.tsv"));
                while ((line = TSVReader.readLine()) != null) {
                    String[] lineSplit2 = line.split("\t");
                    PeptideFormatter pf = null;
                    if (Constants.spectraModel.equals("Prosit") || Constants.rtModel.equals("Prosit")) {
                        pf = new PeptideFormatter(
                                new PeptideFormatter(lineSplit2[0], lineSplit2[1], "base").getProsit(prositAAMods),
                                lineSplit2[1], "prosit");
                    } else if (Constants.spectraModel.equals("PrositTMT") || Constants.rtModel.equals("PrositTMT")) {
                        pf = new PeptideFormatter(
                                new PeptideFormatter(lineSplit2[0], lineSplit2[1], "base").getPrositTMT(),
                                lineSplit2[1], "prosit");
                    } else {
                        printError("spectraRTPredModel must either be Prosit or PrositTMT");
                        System.exit(1);
                    }

                    if (! pf.getBase().equals(lineSplit2[0])) {
                        //get predictionEntry
                        PredictionEntry tmp = allPreds.get(pf.getBaseCharge());

                        if (tmp == null) { //valid reasons to be empty
                            if (PeptideSkipper.skipPeptide(pf,
                                    Constants.spectraModel + Constants.rtModel + Constants.imModel)) {
                                continue;
                            }
                        }

                        MassCalculator mc = new MassCalculator(lineSplit2[0], lineSplit2[1]);
                        float[] newMZs = new float[tmp.mzs.length];
                        for (int i = 0; i < newMZs.length; i++) {
                            newMZs[i] = mc.calcMass(tmp.fragNums[i], tmp.fragmentIonTypes[i], tmp.charges[i],
                                    tmp.isotopes[i]);
                        }

                        //add to hashmap
                        PredictionEntry newPred = new PredictionEntry(newMZs, tmp.intensities,
                                tmp.fragNums, tmp.charges, tmp.fragmentIonTypes);
                        newPred.setRT(tmp.RT);
                        newPred.setIM(0f);
                        allPreds.put(mc.fullPeptide, newPred);
                    }
                }
            } catch (IOException e) {
                e.printStackTrace();
                System.exit(1);
            }
        }
    }

    public PredictionEntryHashMap getPreds() { return allPreds; }
    public void setPreds(PredictionEntryHashMap preds) {
        allPreds = preds;
    }

    public void clear() {
        allPreds.clear();
    }
}
