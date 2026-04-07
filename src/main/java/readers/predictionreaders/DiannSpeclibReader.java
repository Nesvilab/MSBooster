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
import predictions.PredictionEntry;
import predictions.PredictionEntryHashMap;
import utils.Print;
import utils.ProgressReporter;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.io.InputStream;
import java.util.ArrayList;

import static utils.Print.printError;

//Note: DiannSpeclibReader class is the .predicted.bin file and accompanying files from DIA-NN
public class DiannSpeclibReader implements LibraryPredictionMapper {
    final ArrayList<String> filenames;
    PredictionEntryHashMap allPreds = new PredictionEntryHashMap();

    private static final String[] flagTOion = {"b", "y"}; // 0=b, 1=y

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

        // Pre-count total entries across all files to pre-size the HashMap (avoids costly rehashing)
        int totalAllEntries = 0;
        for (String bFile : filenames) {
            int splitDot = bFile.indexOf("predicted.bin");
            String textFile = bFile.substring(0, splitDot) + "tsv";
            try (BufferedReader counter = new BufferedReader(new FileReader(textFile))) {
                counter.readLine(); // skip header
                while (counter.readLine() != null) totalAllEntries++;
            } catch (IOException e) {
                // validation below will catch missing files
            }
        }
        allPreds = new PredictionEntryHashMap(totalAllEntries);

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
                // Count entries for progress reporting
                int totalEntries = -1; // -1 to subtract header
                try (BufferedReader counter = new BufferedReader(new FileReader(textFile))) {
                    while (counter.readLine() != null) totalEntries++;
                }
                Print.printInfo("Reading " + totalEntries + " entries from binary spectral library");
                ProgressReporter pr = new ProgressReporter(totalEntries);

                InputStream is = new FileInputStream(bFile);
                BufferedReader TSVReader = new BufferedReader(new FileReader(textFile));

                int len; //holds length of bytes
                byte[] buffer1 = new byte[12];
                byte[] buffer2 = new byte[4 * 100]; // pre-allocated; reused each peptide (max 100 frags)
                TSVReader.readLine(); //header

                while ((len = is.read(buffer1)) != -1) {
                    String rawLine = TSVReader.readLine();
                    int tab = rawLine.indexOf('\t');
                    String peptideSeq = rawLine.substring(0, tab);
                    String chargeStr = rawLine.substring(tab + 1);
                    pr.progress();
                    MassCalculator mc = new MassCalculator(new PeptideFormatter(peptideSeq, chargeStr, "diann").getBase(), chargeStr);

                    //get data for precursor — read directly from buffer1 bytes (little-endian)
                    int numFrags = (buffer1[0] & 0xFF) | ((buffer1[1] & 0xFF) << 8)
                            | ((buffer1[2] & 0xFF) << 16) | (buffer1[3] << 24);
                    int iRTbits = (buffer1[4] & 0xFF) | ((buffer1[5] & 0xFF) << 8)
                            | ((buffer1[6] & 0xFF) << 16) | (buffer1[7] << 24);
                    float iRT = Float.intBitsToFloat(iRTbits);
                    int IMbits = (buffer1[8] & 0xFF) | ((buffer1[9] & 0xFF) << 8)
                            | ((buffer1[10] & 0xFF) << 16) | (buffer1[11] << 24);
                    float IM = Float.intBitsToFloat(IMbits);

                    //arrays for hashmap
                    float[] mzs = new float[numFrags];
                    float[] intensities = new float[numFrags];
                    int[] fragNums = new int[numFrags];
                    String[] fragmentIonTypes = new String[numFrags];
                    int[] charges = new int[numFrags];

                    //load fragment info — reuse buffer2 (grow only if needed)
                    int needed = 4 * numFrags;
                    if (needed > buffer2.length) {
                        buffer2 = new byte[needed];
                    }
                    len = is.read(buffer2, 0, needed);

                    //iterate through fragments — read bytes directly, inline bit extraction
                    for (int i = 0; i < numFrags; i++) {
                        int off = i * 4;
                        int fragInt = (buffer2[off] & 0xFF) | ((buffer2[off + 1] & 0xFF) << 8)
                                | ((buffer2[off + 2] & 0xFF) << 16) | (buffer2[off + 3] << 24);
                        int intensity = (fragInt >>> 16) & 0xFFFF;   // bits[31..16]
                        int fragNum   = (fragInt >>> 8)  & 0xFF;     // bits[15..8]
                        int flag      = (fragInt >>> 2)  & 0x1;      // bit[2]  y=1, b=0
                        int charge    = (fragInt & 0x3) + 1;         // bits[1..0]

                        //get fragment m/z
                        String ionType = flagTOion[flag];
                        float fragMZ = mc.calcMass(fragNum, ionType, charge, 0);

                        //add to arrays
                        mzs[i] = fragMZ;
                        intensities[i] = intensity;
                        fragNums[i] = fragNum;
                        fragmentIonTypes[i] = ionType;
                        charges[i] = charge;
                    }

                    //add to hashmap
                    PredictionEntry newPred = new PredictionEntry(mzs, intensities,
                            fragNums, charges, fragmentIonTypes);
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
                int totalFullEntries = 0;
                try (BufferedReader counter = new BufferedReader(new FileReader(textFile))) {
                    while (counter.readLine() != null) totalFullEntries++;
                }
                Print.printInfo("Reading " + totalFullEntries + " entries from full TSV");
                ProgressReporter prFull = new ProgressReporter(totalFullEntries);

                TSVReader = new BufferedReader(new FileReader(textFile));
                String l;

                while ((l = TSVReader.readLine()) != null) {
                    int tab = l.indexOf('\t');
                    String peptideSeq = l.substring(0, tab);
                    String chargeStr = l.substring(tab + 1);
                    prFull.progress();
                    //check if diann to base results in same base peptide
                    PeptideFormatter pf = new PeptideFormatter(
                            new PeptideFormatter(peptideSeq, chargeStr, "base").getDiann(), chargeStr, "diann");

                    if (! pf.getBase().equals(peptideSeq)) {
                        //get predictionEntry
                        PredictionEntry tmp = allPreds.get(pf.getBaseCharge());
                        MassCalculator mc = new MassCalculator(peptideSeq, chargeStr);
                        float[] newMZs = new float[tmp.mzs.length];
                        for (int i = 0; i < newMZs.length; i++) {
                            newMZs[i] = mc.calcMass(tmp.fragNums[i], tmp.fragmentIonTypes[i],
                                    tmp.charges[i], tmp.isotopes.length > 0 ? tmp.isotopes[i] : 0);
                        }

                        //add to hashmap
                        PredictionEntry newPred = new PredictionEntry(newMZs, tmp.intensities,
                                tmp.fragNums, tmp.charges, tmp.fragmentIonTypes);
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

    public static void main(String[] args) throws FileNotFoundException {
        Constants.numThreads = 11;
        DiannSpeclibReader dslr = new DiannSpeclibReader(
                "Z:/yangkl/troubleshooting_issues/github1435/PublicSearch.266.HLAI.FDR3.pekin/MSBooster");
        dslr = new DiannSpeclibReader("C:/Users/yangkl/Downloads/proteomics/" +
                "diannspeclib/spectraRT.predicted.bin");
    }
}
