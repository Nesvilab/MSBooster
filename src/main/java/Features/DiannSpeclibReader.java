package Features;

import java.io.*;
import java.nio.ByteBuffer;
import java.nio.ByteOrder;
import java.util.ArrayList;
import java.util.HashMap;

public class DiannSpeclibReader implements SpectralPredictionMapper{
    final ArrayList<String> filenames;
    HashMap<String, PredictionEntry> allPreds = new HashMap<>();

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
            String textFile = bFile.substring(0, splitDot) + "tsv"; //enforces tsv naming convention

            try (
                    InputStream is = new FileInputStream(bFile);
                    BufferedReader TSVReader = new BufferedReader(new FileReader(textFile));
            ) {
                int len; //holds length of bytes
                byte[] buffer1 = new byte[12];
                String[] line = TSVReader.readLine().split("\t"); //header

                while ((len = is.read(buffer1)) != -1) {
                    line = TSVReader.readLine().split("\t");
                    MassCalculator mc = new MassCalculator(line[0], line[1]);

                    //get data for precursor
                    int numFrags = ByteBuffer.wrap(buffer1, 0, 4).order(ByteOrder.LITTLE_ENDIAN).getInt();
                    float iRT = ByteBuffer.wrap(buffer1, 4, 4).order(ByteOrder.LITTLE_ENDIAN).getFloat();
                    float IM = ByteBuffer.wrap(buffer1, 8, 4).order(ByteOrder.LITTLE_ENDIAN).getFloat();

                    //arrays for hashmap
                    float[] mzs = new float[numFrags];
                    float[] intensities = new float[numFrags];

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
                        float fragMZ = mc.calcMass(fragNum, flag, charge);

                        //add to arrays
                        mzs[i] = fragMZ;
                        intensities[i] = intensity;
                    }

                    //add to hashmap
                    PredictionEntry newPred = new PredictionEntry();
                    newPred.setMzs(mzs);
                    newPred.setIntensities(intensities);
                    newPred.setRT(iRT);
                    newPred.setIM(IM);
                    allPreds.put(mc.fullPeptide, newPred);
                }

                //changing to deepccs values
//                BufferedReader newTSVReader = new BufferedReader(new FileReader("C:/Users/kevin/Downloads/proteomics/timsTOF/" +
//                        "DeepCCS_pred_diannName.csv"));
//                String line1 = "";
//                while ((line1 = newTSVReader.readLine()) != null) {
//                    String[] line1Split = line1.split(",");
//                    MassCalculator mc = new MassCalculator(line1Split[0], line1Split[1]);
//                    allPredIMs.put(mc.fullPeptide, Float.valueOf(line1Split[2]));
//                }
//                newTSVReader.close();

            } catch (IOException e) {
                e.printStackTrace();
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

    public static void main(String[] args) throws IOException {
        //DiannSpeclibReader d = new DiannSpeclibReader("C:/Users/kevin/OneDriveUmich/proteomics/preds/cptacDiann.predicted.bin");
        long startTime = System.nanoTime();
        SpectralPredictionMapper spm = SpectralPredictionMapper.createSpectralPredictionMapper(
                "C:/Users/kevin/Downloads/proteomics/wide/spectraRT.predicted.bin");
        long endTime = System.nanoTime();
        long duration = (endTime - startTime);
        System.out.println("Spectra/RT/IM prediction loading took " + duration / 1000000 +" milliseconds");
    }
}
