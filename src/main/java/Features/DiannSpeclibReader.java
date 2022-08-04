package Features;

import umich.ms.fileio.exceptions.FileParsingException;

import java.io.*;
import java.nio.ByteBuffer;
import java.nio.ByteOrder;
import java.sql.SQLException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.Map;
import java.util.concurrent.ExecutionException;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;

public class DiannSpeclibReader implements SpectralPredictionMapper{
    final ArrayList<String> filenames;
    HashMap<String, PredictionEntry> allPreds = new HashMap<>();

    //convert int flag to fragment ion type
    private static HashMap<Integer, String> makeFlagTOion() {
        HashMap<Integer, String> map = new HashMap<>();
        map.put(0, "b");
        map.put(1, "y");
        return map;
    }
    HashMap<Integer, String> flagTOion = makeFlagTOion();

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
                System.out.println("Error: no prediction file available at: " + bFile);
                System.exit(-1);
            }
            String textFile = bFile.substring(0, splitDot) + "tsv"; //enforces tsv naming convention
            if (! new File(textFile).exists()) {
                System.out.println("Error: no prediction file available at: " + textFile);
                System.exit(-1);
            }

            try{
                InputStream is = new FileInputStream(bFile);
                BufferedReader TSVReader = new BufferedReader(new FileReader(textFile));

                int len; //holds length of bytes
                byte[] buffer1 = new byte[12];
                String[] line = TSVReader.readLine().split("\t"); //header

                while ((len = is.read(buffer1)) != -1) {
                    line = TSVReader.readLine().split("\t");
                    MassCalculator mc = new MassCalculator(new PeptideFormatter(line[0], line[1], "diann").base, line[1]);

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
                        String ionType = flagTOion.get(flag);
                        float fragMZ = mc.calcMass(fragNum, ionType, charge);

                        //add to arrays
                        mzs[i] = fragMZ;
                        intensities[i] = intensity;
                        fragNums[i] = fragNum;
                        flags[i] = flag;
                        fragmentIonTypes[i] = ionType;
                        charges[i] = charge;
                    }

                    //add to hashmap
                    PredictionEntry newPred = new PredictionEntry();
                    newPred.setMzs(mzs);
                    newPred.setIntensities(intensities);
                    newPred.setFragNums(fragNums);
                    newPred.setFlags(flags);
                    newPred.setCharges(charges);
                    newPred.setRT(iRT);
                    newPred.setIM(IM);
                    newPred.setFragmentIonTypes(fragmentIonTypes);
                    allPreds.put(mc.fullPeptide, newPred);
                }
                is.close();

                if (TSVReader.readLine() != null) {
                    System.out.println("Prediction file is missing some entries. Please rerun MSBooster");
                    System.exit(-1);
                }
                TSVReader.close();

                //comment the next section out if doing if for TMT alignment
                //repeat this process with full peptides
                textFile = bFile.substring(0, splitDot - 1) + "_full.tsv";
                TSVReader = new BufferedReader(new FileReader(textFile));
                String l;

                while ((l = TSVReader.readLine()) != null) {
                    line = l.split("\t");
                    //check if diann to base results in same base peptide
                    PeptideFormatter pf = new PeptideFormatter(
                            new PeptideFormatter(line[0], line[1], "base").diann, line[1], "diann");

                    if (! pf.base.equals(line[0])) {
                        //get predictionEntry
                        PredictionEntry tmp = allPreds.get(pf.baseCharge);
                        MassCalculator mc = new MassCalculator(line[0], line[1]);
                        float[] newMZs = new float[tmp.mzs.length];
                        for (int i = 0; i < newMZs.length; i++) {
                            newMZs[i] = mc.calcMass(tmp.fragNums[i], flagTOion.get(tmp.flags[i]), tmp.charges[i]);
                        }

                        //add to hashmap
                        PredictionEntry newPred = new PredictionEntry();
                        newPred.setMzs(newMZs);
                        newPred.setIntensities(tmp.intensities);
                        newPred.setRT(tmp.RT);
                        newPred.setIM(tmp.IM);
                        newPred.setFragmentIonTypes(tmp.fragmentIonTypes);
                        newPred.setFragNums(tmp.fragNums);
                        newPred.setFlags(tmp.flags);
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

    public HashMap<String, PredictionEntry> getPreds() { return allPreds; }

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

    public static void main(String[] args)
            throws IOException, InterruptedException, ExecutionException, FileParsingException, SQLException {
        ExecutorService executorService = Executors.newFixedThreadPool(11);
        SpectralPredictionMapper spm = SpectralPredictionMapper.createSpectralPredictionMapper("C:/Users/yangkl/OneDriveUmich/proteomics/" +
                        "RTalignTMT/spectraRT.predicted.bin",
                executorService);
        executorService.shutdown();
    }
}
