import org.apache.commons.lang.ArrayUtils;
import umich.ms.fileio.exceptions.FileParsingException;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.lang.reflect.Method;
import java.util.Arrays;
import java.util.HashMap;

public class combinedFiles {

    public static void createSimilarityFile(String predFile, String pepXmlFiles, String mzmlFiles,
                                            String similarityMeasure, String outFile)
            throws IOException, FileParsingException {
        long startTime = System.nanoTime(); // timer: https://stackoverflow.com/questions/3382954/measure-execution-time-for-a-java-method

        ////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        //get predictions

        //load in predicted spectra
        mgfFileReader predictedSpectra = new mgfFileReader(predFile);

        //get predicted MZ's and intensities in form of HashMaps
        HashMap<String, double[]> predictedMZs = predictedSpectra.getMzDict();
        HashMap<String, double[]> predictedIntensities = predictedSpectra.getIntensityDict();

        ////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        //get paths for pepXML and mzML files
        String pepXMLDirectoryName = pepXmlFiles;
        File pepXMLDirectory = new File(pepXMLDirectoryName);
        String[] pepXMLfiles = pepXMLDirectory.list((d, name) -> name.endsWith(".pepXML"));
        Arrays.sort(pepXMLfiles);

        String mzmlDirectoryName = mzmlFiles;
        File mzmlDirectory = new File(mzmlDirectoryName);
        String[] mzmlfiles = mzmlDirectory.list((d, name) -> name.endsWith(".mzML"));
        Arrays.sort(mzmlfiles);

        ////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        String[] allSpectralAngles = new String[0];

        for (int j = 0; j < mzmlfiles.length; j++) {
            System.out.println("file" + j);

            //read in specific files
            String xmlPath = pepXMLfiles[j];
            String mzmlPath = mzmlfiles[j];

            pepXMLReader xmlReader = new pepXMLReader(pepXMLDirectoryName + xmlPath);

            //get all scan numbers from pepXML
            int[] scanNums = xmlReader.getScanNumbers();

            //get peptide names compatible with predicted HashMaps
            String[] xmlPeptides = xmlReader.getXMLpeptides();

            //get target or decoys
            int[] td = xmlReader.getTargetOrDecoy();

            //get expectation scores
            String[] eScores = xmlReader.getEScore();

            //load experimental spectra from mzML
            mzMLReader mzMLscans = new mzMLReader(mzmlDirectoryName + mzmlPath);

            //get mzFreq
            double[] mzFreqs = mzMLscans.getMzFreq();

            //iterate through all pepXMLhits
            int pepXMLhits = scanNums.length;

            //save spectral angles
            String[] sa = new String[pepXMLhits];

            for (int i = 0; i < pepXMLhits; i++) {
                //get scanNumber and peptide of current hit
                int num = scanNums[i];
                String pep = xmlPeptides[i];

                //get experimental spectra into
                double[] expMZs = mzMLscans.getMZ(num);
                double[] expIntensities = mzMLscans.getIntensity(num);

                //get predicted spectra info
                double[] predMZs = predictedMZs.get(pep);
                double[] predIntensities = predictedIntensities.get(pep);

                //calculate spectral angle
                try {
                    spectrumComparison specAngle = new spectrumComparison(expMZs, expIntensities,
                            predMZs, predIntensities);

                    //adapt for weighted or multiple sims, refer to pepXMLModifier
                    Method method = specAngle.getClass().getMethod(similarityMeasure);
                    double sim = (double) method.invoke(specAngle);

                    //get expectation score and target or decoy
                    String eScore = eScores[i];
                    int tdSingle = td[i];

                    sa[i] = pep + "\t" + sim + "\t" +
                            eScore + "\t" + tdSingle;

                } catch(Exception e) {
                    System.out.println(pep); //probably was not supported by pDeep2 prediction (ex. amino acid U)
                }

                //calculate deltaRT

            }
            //add to full list of spectral angles of peptides
            allSpectralAngles = (String[]) ArrayUtils.addAll(allSpectralAngles, sa);

        }
        ////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        //will need to remove null
        //send to file for analysis/data viz in python/R
        try {
            FileWriter myWriter = new FileWriter(outFile);
            //need to adapt for multiple similarities
            myWriter.write("peptide" + "\t" + similarityMeasure + "\t" +
                    "eScore" + "\t" + "target=1/decoy=0\n");

            for (String s : allSpectralAngles) {
                if (s != null) {
                    myWriter.write(s + "\n");
                } else {
                    System.out.println("found empty line");
                }
            }

            myWriter.close();
            System.out.println("Successfully wrote to the file.");
        } catch (IOException e) {
            System.out.println("An error occurred.");
            e.printStackTrace();
        }

        long stopTime = System.nanoTime();
        System.out.println((stopTime - startTime) / 1000000000.0);
    }

    public static void main(String[] args) throws IOException, FileParsingException {
//        createSimilarityFile("C:/Users/kevin/Downloads/proteomics/HeLa_preds_allmods.mgf",
//                "C:/Users/kevin/OneDriveUmich/proteomics/pepxml/rank1/",
//                "C:/Users/kevin/OneDriveUmich/proteomics/mzml/",
//                "brayCurtis",
//                "C:/Users/kevin/Downloads/proteomics/pDeep2_allMods_similarities.tsv");
        createSimilarityFile("C:/Users/kevin/Downloads/proteomics/HeLa_preds2_pretrained.mgf",
                "C:/Users/kevin/OneDriveUmich/proteomics/pepxml/rank1/",
                "C:/Users/kevin/OneDriveUmich/proteomics/mzml/",
                "brayCurtis",
                "C:/Users/kevin/Downloads/proteomics/pDeep2_pretrained_similarities.tsv");
    }
}
