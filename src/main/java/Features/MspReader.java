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

package Features;

import static utils.Print.printError;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.concurrent.ConcurrentHashMap;

public class MspReader implements SpectralPredictionMapper {
    final ArrayList<String> filenames;
    PredictionEntryHashMap allPreds = new PredictionEntryHashMap();

    //convert int flag to fragment ion type
    private static HashMap<Integer, String> makeFlagTOion() {
        HashMap<Integer, String> map = new HashMap<>();
        map.put(0, "b");
        map.put(1, "y");
        return map;
    }
    HashMap<Integer, String> flagTOion = makeFlagTOion();

    public MspReader(String files) throws FileNotFoundException {
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
                    pep = new PeptideFormatter(sb.toString(), charge, "msp").baseCharge;

                    //fragments
                    //if it keeps reporting 174 peaks, set that as initial arraylist size
                    line = msp.readLine(); //Num peaks
                    int numPeaks = Integer.parseInt(line.split(" ")[2]);
                    ArrayList<Float> mzsList = new ArrayList<>(numPeaks);
                    ArrayList<Float> intsList = new ArrayList<>(numPeaks);
                    ArrayList<Integer> fragNumsList = new ArrayList<>(numPeaks);
                    ArrayList<Integer> flagsList = new ArrayList<>(numPeaks);
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
                        if (fragment.charAt(0) == 'y') {
                            flagsList.add(1);
                        } else {
                            flagsList.add(0);
                        }
                        fragNumsList.add(fragNum);
                    }

                    float[] mzs = new float[mzsList.size()];
                    float[] ints = new float[intsList.size()];
                    int[] fragNums = new int[intsList.size()];
                    int[] flags = new int[intsList.size()];
                    int[] charges = new int[intsList.size()];
                    for (int i = 0; i < mzs.length; i++) {
                        mzs[i] = mzsList.get(i);
                        ints[i] = intsList.get(i);
                        fragNums[i] = fragNumsList.get(i);
                        flags[i] = flagsList.get(i);
                        charges[i] = chargesList.get(i);
                    }

                    PredictionEntry newPred = new PredictionEntry(mzs, ints,
                            fragNums, charges, new String[0]);
                    newPred.setRT(RT);
                    newPred.setIM(0f);
                    newPred.setFlags(flags);
                    allPreds.put(pep, newPred);
                }

                //read in extra files
                BufferedReader TSVReader = new BufferedReader(new FileReader(
                        f.substring(0, f.length() - 4) + "_full.tsv"));
                while ((line = TSVReader.readLine()) != null) {
                    String[] lineSplit2 = line.split("\t");
                    PeptideFormatter pf = null;
                    if (Constants.spectraRTPredModel.equals("Prosit")) {
                        pf = new PeptideFormatter(
                                new PeptideFormatter(lineSplit2[0], lineSplit2[1], "base").prosit,
                                lineSplit2[1], "prosit");
                    } else if (Constants.spectraRTPredModel.equals("PrositTMT")) {
                        pf = new PeptideFormatter(
                                new PeptideFormatter(lineSplit2[0], lineSplit2[1], "base").prositTMT,
                                lineSplit2[1], "prosit");
                    } else {
                        printError("spectraRTPredModel must either be Prosit or PrositTMT");
                        System.exit(1);
                    }

                    if (! pf.base.equals(lineSplit2[0])) {
                        //get predictionEntry
                        PredictionEntry tmp = allPreds.get(pf.baseCharge);

                        if (tmp == null) { //valid reasons to be empty
                            if (PeptideSkipper.skipPeptide(pf.stripped, pf.charge)) {
                                continue;
                            }
                        }

                        MassCalculator mc = new MassCalculator(lineSplit2[0], lineSplit2[1]);
                        float[] newMZs = new float[tmp.mzs.length];
                        for (int i = 0; i < newMZs.length; i++) {
                            newMZs[i] = mc.calcMass(tmp.fragNums[i], flagTOion.get(tmp.flags[i]), tmp.charges[i]);
                        }

                        //add to hashmap
                        PredictionEntry newPred = new PredictionEntry(newMZs, tmp.intensities,
                                tmp.fragNums, tmp.charges, new String[0]);
                        newPred.setFlags(tmp.flags);
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

    public float getMaxPredRT() {
        float maxRT = 0f;
        for (PredictionEntry entry : allPreds.values()) {
            if (entry.RT > maxRT) {
                maxRT = entry.RT;
            }
        }
        return maxRT;
    }

    public void clear() {
        allPreds.clear();
    }
}
