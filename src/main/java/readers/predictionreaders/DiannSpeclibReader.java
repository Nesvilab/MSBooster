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

import features.spectra.MassCalculator;
import peptideptmformatting.PeptideFormatter;
import predictions.PredictionEntry;
import predictions.PredictionEntryHashMap;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.io.InputStream;
import java.nio.ByteBuffer;
import java.nio.ByteOrder;
import java.util.ArrayList;

import static utils.Print.printError;

public class DiannSpeclibReader implements LibraryPredictionMapper {
    final ArrayList<String> filenames;
    PredictionEntryHashMap allPreds = new PredictionEntryHashMap();

    //https://stackoverflow.com/questions/46163114/get-bit-values-from-byte-array
    //https://www.geeksforgeeks.org/bitwise-operators-in-java/
    public DiannSpeclibReader(String binFile) throws FileNotFoundException {
        File predsDirectory = new File(binFile);
        String[] predsFiles = predsDirectory.list();
        filenames = new ArrayList<String>();

        if (predsFiles == null) { //if user provided a file, not a directory
            filenames.add(binFile);
        } else { //user provides directory
            for (String predsFile : predsFiles) {
                if (predsFile.contains("predicted.bin")) {
                    filenames.add(binFile + File.separator + predsFile);
                }
            }
        }

        for (String bFile : filenames) {
            //try to infer binary file name from text file
            int splitDot = bFile.indexOf("predicted.bin");
            if (! new File(bFile).exists()) {
                printError("Error: no prediction file available at: " + bFile);
                System.exit(1);
            }
            String textFile = bFile.substring(0, splitDot) + "tsv"; //enforces tsv naming convention
            if (! new File(textFile).exists()) {
                printError("Error: no prediction file available at: " + textFile);
                System.exit(1);
            }

            try{
                InputStream is = new FileInputStream(bFile);
                BufferedReader TSVReader = new BufferedReader(new FileReader(textFile));

                int len; //holds length of bytes
                byte[] buffer1 = new byte[12];
                String[] line = TSVReader.readLine().split("\t"); //header

                while ((len = is.read(buffer1)) != -1) {
                    line = TSVReader.readLine().split("\t");
                    MassCalculator mc = new MassCalculator(new PeptideFormatter(line[0], line[1], "diann").getBase(), line[1]);

                    //get data for precursor
                    int numFrags = ByteBuffer.wrap(buffer1, 0, 4).order(ByteOrder.LITTLE_ENDIAN).getInt();
                    float iRT = ByteBuffer.wrap(buffer1, 4, 4).order(ByteOrder.LITTLE_ENDIAN).getFloat();
                    float IM = ByteBuffer.wrap(buffer1, 8, 4).order(ByteOrder.LITTLE_ENDIAN).getFloat();

                    //arrays for hashmap
                    float[] mzs = new float[numFrags];
                    float[] intensities = new float[numFrags];
                    int[] fragNums = new int[numFrags];
                    int[] flags = new int[numFrags];
                    String[] fragmentIonTypes = new String[numFrags];
                    int[] charges = new int[numFrags];

                    //load fragment info
                    byte[] buffer2 = new byte[4 * numFrags];
                    len = is.read(buffer2);

                    //iterate through fragments
                    for (int i = 0; i < numFrags; i++) {
                        int fragInt = ByteBuffer.wrap(buffer2, i * 4, 4).order(ByteOrder.LITTLE_ENDIAN).getInt();
                        int intensity = bits(fragInt, 0, 16);
                        int fragNum = bits(fragInt, 16, 8);
                        int flag = bits(fragInt, 29, 1); //y 1, b 0
                        int charge = bits(fragInt, 30, 2) + 1; //start from end

                        //get fragment m/z
                        String ionType = MassCalculator.flagTOion.get(flag);
                        float fragMZ = mc.calcMass(fragNum, ionType, charge, 0);

                        //add to arrays
                        mzs[i] = fragMZ;
                        intensities[i] = intensity;
                        fragNums[i] = fragNum;
                        flags[i] = flag;
                        fragmentIonTypes[i] = ionType;
                        charges[i] = charge;
                    }

                    //add to hashmap
                    PredictionEntry newPred = new PredictionEntry(mzs, intensities,
                            fragNums, charges, fragmentIonTypes, flags);
                    newPred.setRT(iRT);
                    newPred.setIM(IM);
                    allPreds.put(mc.fullPeptide, newPred);
                }
                is.close();

                if (TSVReader.readLine() != null) {
                    printError("Prediction file is missing some entries. Please rerun MSBooster");
                    System.exit(1);
                }
                TSVReader.close();

                textFile = bFile.substring(0, splitDot - 1) + "_full.tsv";
                TSVReader = new BufferedReader(new FileReader(textFile));
                String l;

                while ((l = TSVReader.readLine()) != null) {
                    line = l.split("\t");
                    //check if diann to base results in same base peptide
                    PeptideFormatter pf = new PeptideFormatter(
                            new PeptideFormatter(line[0], line[1], "base").getDiann(), line[1], "diann");

                    if (! pf.getBase().equals(line[0])) {
                        //get predictionEntry
                        PredictionEntry tmp = allPreds.get(pf.getBaseCharge());
                        MassCalculator mc = new MassCalculator(line[0], line[1]);
                        float[] newMZs = new float[tmp.mzs.length];
                        for (int i = 0; i < newMZs.length; i++) {
                            newMZs[i] = mc.calcMass(tmp.fragNums[i], MassCalculator.flagTOion.get(tmp.flags[i]),
                                    tmp.charges[i], tmp.isotopes[i]);
                        }

                        //add to hashmap
                        PredictionEntry newPred = new PredictionEntry(newMZs, tmp.intensities,
                                tmp.fragNums, tmp.charges, tmp.fragmentIonTypes, tmp.flags);
                        newPred.setRT(tmp.RT);
                        newPred.setIM(tmp.IM);
                        allPreds.put(mc.fullPeptide, newPred);
                    }
                }
            } catch (IOException e) {
                e.printStackTrace();
                System.exit(1);
            }
        }
    }

    //for extracting info from fragments
    private static int bits(int n, int offset, int length) {
        return n >> (32 - offset - length) & ~(-1 << length);
    }

    public PredictionEntryHashMap getPreds() { return allPreds; }
    public void setPreds(PredictionEntryHashMap preds) {
        allPreds = preds;
    }

    public void clear() {
        allPreds.clear();
    }
}
