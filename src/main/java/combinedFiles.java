import org.apache.commons.lang.ArrayUtils;
import umich.ms.fileio.exceptions.FileParsingException;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.Arrays;
import java.util.HashMap;

public class combinedFiles {

    public static void main(String[] args) throws FileParsingException, IOException {
        long startTime = System.nanoTime(); // timer: https://stackoverflow.com/questions/3382954/measure-execution-time-for-a-java-method

        ////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        //get predictions

        //load in predicted spectra
        mgfFileReader predictedSpectra = new mgfFileReader("rank4_preds.txt");

        //get predicted MZ's and intensities in form of HashMaps
        HashMap<String, double[]> predictedMZs = predictedSpectra.getMzDict();
        HashMap<String, double[]> predictedIntensities = predictedSpectra.getIntensityDict();

        ////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        //get paths for pepXML and mzML files
        String pepXMLDirectoryName = "C:/Users/yangkl/OneDriveUmich/proteomics/pepxml/rank4/";
        File pepXMLDirectory = new File(pepXMLDirectoryName);
        String[] pepXMLfiles = pepXMLDirectory.list();
        Arrays.sort(pepXMLfiles);

        String mzmlDirectoryName = "C:/Users/yangkl/Documents/msfragger_dia_ms1/files/mzml/";
        File mzmlDirectory = new File(mzmlDirectoryName);
        String[] mzmlfiles = mzmlDirectory.list();
        Arrays.sort(mzmlfiles);

        ////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        String[] allSpectralAngles = new String[0];

        //constants
        double binwidth = 1;

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
            double[] mzFreqs = mzMLscans.getMzFreq(binwidth, 1);

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
                            predMZs, predIntensities, 20);
                    double sim0 = specAngle.cosineSimilarity();
                    double sim1 = specAngle.spectralContrastAngle();
                    double sim2 = specAngle.euclideanDistance();
                    double sim3 = specAngle.brayCurtis();
                    double sim4 = specAngle.pearsonCorr();

                    //weighted
                    double[] weights = specAngle.getWeights(mzFreqs, binwidth);
                    double sim5 = specAngle.weightedCosineSimilarity(weights);
                    double sim6 = specAngle.weightedSpectralContrastAngle(weights);
                    double sim7 = specAngle.weightedEuclideanDistance(weights);
                    double sim8 = specAngle.weightedBrayCurtis(weights);
                    double sim9 = specAngle.weightedPearsonCorr(weights);

                    //dot product
                    double sim10 = specAngle.dotProduct();
                    double sim11 = specAngle.weightedDotProduct(weights);

                    //get expectation score and target or decoy
                    String eScore = eScores[i];
                    int tdSingle = td[i];

                    sa[i] = pep + "\t" +
                            sim0 + "\t" + sim1 + "\t" +
                            sim2 + "\t" + sim3 + "\t" +
                            sim4 + "\t" + sim5 + "\t" +
                            sim6 + "\t" + sim7 + "\t" +
                            sim8 + "\t" + sim9 + "\t" +
                            sim10 + "\t" + sim11 + "\t" +
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
            FileWriter myWriter = new FileWriter("rank4_analysis.tsv");
            myWriter.write("peptide" + "\t" + "cosine_similarity" + "\t" +
                    "contrast_angle" + "\t" + "euclidean" + "\t" +
                    "bray-curtis" + "\t" + "pearson" + "\t" + "weight_cosine" + "\t" +
                    "weight_spectral_contrast" + "\t" + "weight_euclidean" + "\t" +
                    "weight_bray-curtis" + "\t" + "weight_pearson" + "\t" +
                    "dot_product" + "\t" + "weight_dot" + "\t" +
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
}
