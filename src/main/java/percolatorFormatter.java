import umich.ms.fileio.exceptions.FileParsingException;

import java.io.*;
import java.util.Arrays;
import java.util.HashMap;
import java.util.Map;
import java.util.TreeMap;

public class percolatorFormatter {

    public static String percolatorPepFormat(String[] columns) {
        String pep = columns[31];
        pep = pep.substring(2, pep.length() - 2);
        StringBuilder mods = new StringBuilder();

        //n term acetylation
        if (pep.charAt(0) == 'n') {
            pep = pep.replace("n[42.0106]", "");
            mods.append("0,Acetyl[AnyN-term];");
        }

        pep = pep.replace("C[57.0215]", "c");
        pep = pep.replace("M[15.9949]", "m");
        TreeMap<Integer, String> modsMap = new TreeMap<>();

        //carbamidometylation
        while (pep.contains("c")) {
            int pos = pep.indexOf("c") + 1;
            pep = pep.replaceFirst("c", "C");
            modsMap.put(pos, ",Carbamidomethyl[C];");
        }

        //methionine oxidation
        while (pep.contains("m")) {
            int pos = pep.indexOf("m") + 1;
            pep = pep.replaceFirst("m", "M");
            modsMap.put(pos, ",Oxidation[M];");
        }

        for (Map.Entry<Integer, String> entry : modsMap.entrySet()) {
            mods.append(entry.getKey()).append(entry.getValue());
        }

        //charge
        String[] charges = Arrays.copyOfRange(columns, 24, 28);
        int charge = Arrays.asList(charges).indexOf("1") + 1;

        return pep + "|" + mods + "|" + charge;
    }

    public static void main(String[] args) throws IOException, FileParsingException {
        long startTime = System.nanoTime();

        //load in predicted spectra
        System.out.println("loading in all predicted spectra");
        mgfFileReader predictedSpectra = new mgfFileReader("C:/Users/kevin/OneDrive/IdeaProjects/" +
                "MSFraggerDIA_postprocess/preds/");
        HashMap<String, double[]> allPredMZs = predictedSpectra.allPredMZs;
        HashMap<String, double[]> allPredIntensities = predictedSpectra.allPredIntensities;

        //load experimental spectra from mzML
        //String mzmlDirectoryName = "C:/Users/yangkl/Documents/msfragger_dia_ms1/files/mzml/";
        String mzmlDirectoryName = "C:/Users/kevin/OneDriveUmich/proteomics/mzml/";
        File mzmlDirectory = new File(mzmlDirectoryName);
        String[] mzmlfiles = mzmlDirectory.list();
        Arrays.sort(mzmlfiles);

        //load pin files
        String pinDirectoryName = "C:/Users/kevin/OneDriveUmich/proteomics/new_results/pin/";
        File pinDirectory = new File(pinDirectoryName);
        String[] pinfiles = pinDirectory.list();
        Arrays.sort(pinfiles);

        boolean usedHeader = false; //so we only write new header once
        //write new file
        //BufferedWriter myWriter = new BufferedWriter(new FileWriter("windows_percolator.tsv"));
        BufferedWriter myWriter = new BufferedWriter(new FileWriter("windows_percolator.tsv"));

        for (int i = 0; i < mzmlfiles.length; i++) {
            System.out.println("loading in mzml file " + i);
            mzMLReader mzMLscans = new mzMLReader(mzmlDirectoryName + mzmlfiles[i]);
            //get mzFreq
            double[] mzFreqs = mzMLscans.getMzFreq();

            try { //https://www.javatpoint.com/how-to-read-file-line-by-line-in-java
                File file = new File(pinDirectoryName + pinfiles[i]);    //creates a new file instance
                FileReader fr = new FileReader(file);   //reads the file
                BufferedReader br = new BufferedReader(fr);  //creates a buffering character input stream
                String line;

                //read header
                String[] header = br.readLine().split("\t", 33);

                //save old scanNum (never -1, so initiated this way)
                int oldScanNum = -1;
                double[] expMZs = new double[0];
                double[] expIntensities = new double[0];

                System.out.println("getting predictions and similarity calculations for each row");
                while ((line = br.readLine()) != null) {
                    String[] columns = line.split("\t", 33);

                    //get experimental
                    int scanNum = Integer.parseInt(columns[2]);
                    String pep = percolatorPepFormat(columns);

                    ///////////////////////////////////////////////////////
                    //delete unless removing modified peptides
//                    String newMods = pep.split("\\|", 3)[1];
//                    if (newMods.contains("M") || newMods.contains("A")) {
//                        continue;
//                    }
                    ///////////////////////////////////////////////////////

                    if (scanNum != oldScanNum) {
                        expMZs = mzMLscans.getMZ(scanNum);
                        expIntensities = mzMLscans.getIntensity(scanNum);
                        oldScanNum = scanNum;
                    }

                    //get predicted
                    double[] predMZs = allPredMZs.get(pep);
                    double[] predIntensities = allPredIntensities.get(pep);

                    if (predMZs != null) {
                        //calculate similarity
                        spectrumComparison specAngle = new spectrumComparison(expMZs, expIntensities, predMZs, predIntensities);
                        TreeMap<String, Double> sims = new TreeMap<>(specAngle.getAllSimilarities(mzFreqs));

                        //write header to new file
                        if (!usedHeader) {

                            String[] pre = Arrays.copyOfRange(header, 0, 31);
                            String[] post = Arrays.copyOfRange(header, 31, 33);
                            StringBuilder newHeader = new StringBuilder("");
                            for (String s : pre) {
                                newHeader.append(s).append("\t");
                            }

                            for (String s : sims.keySet()) {
                                newHeader.append(s).append("\t");
                            }

                            for (String s : post) {
                                newHeader.append(s).append("\t");
                            }

                            myWriter.write(newHeader + "\n");
                            usedHeader = true;
                        }

                        //write row to new file
                        String[] pre = Arrays.copyOfRange(columns, 0, 31);
                        String[] post = Arrays.copyOfRange(columns, 31, 33);
                        StringBuilder newRow = new StringBuilder("");
                        for (String s : pre) {
                            newRow.append(s).append("\t");
                        }

                        for (Map.Entry<String, Double> e : sims.entrySet()) {
                            newRow.append(e.getValue()).append("\t");
                        }

                        for (String s : post) {
                            newRow.append(s).append("\t");
                        }

                        myWriter.write(newRow + "\n");

                    } else {
                        System.out.println(pep);
                    }
                }
                br.close();    //closes the stream and release the resources

            } catch (IOException e) {
                e.printStackTrace();
            }
        }

        myWriter.close();
        long stopTime = System.nanoTime();
        System.out.println((stopTime - startTime) / 1000000000.0);
    }
}
