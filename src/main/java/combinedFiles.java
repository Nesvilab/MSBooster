import umich.ms.fileio.exceptions.FileParsingException;

import java.io.FileWriter;
import java.io.IOException;
import java.util.HashMap;

public class combinedFiles {
    public static void main(String[] args) throws FileParsingException, IOException {
        long startTime = System.nanoTime(); // timer: https://stackoverflow.com/questions/3382954/measure-execution-time-for-a-java-method

        //read in pepXML file
        pepXMLReader xmlReader = new pepXMLReader("23aug2017_hela_serum_timecourse_4mz_narrow_1_rank1.pepXML");

        //create file for pDeep2 prediction
        //xmlReader.createPDeepFile("narrow1rank1.txt"); //actually a tsv file

        //got predictions for all peptides in pepXML
        //python predict.py
        //{'nce': 0.27, 'instrument': 'QE', 'input': 'narrow1rank1.txt', 'output': 'narrow1rank1_pred.txt'}

        //load in predicted spectra
        mgfFileReader predictedSpectra = new mgfFileReader("narrow1rank1_pred.txt");

        //get predicted MZ's and intensities in form of HashMaps
        HashMap<String, double[]> predictedMZs = predictedSpectra.getMzDict();
        HashMap<String, double[]> predictedIntensities = predictedSpectra.getIntensityDict();

        //get all scan numbers from pepXML
        int[] scanNums = xmlReader.getScanNumbers();

        //get peptide names compatible with predicted HashMaps
        String[] xmlPeptides = xmlReader.getXMLpeptides();

        //get target or decoys
        int[] td = xmlReader.getTargetOrDecoy();

        //load experimental spectra from mzML
        mzMLReader mzMLscans = new mzMLReader("23aug2017_hela_serum_timecourse_4mz_narrow_1.mzML");

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

            //get target or decoy
            int tdSingle = td[i];

            //calculate spectral angle
            spectrumComparison specAngle = new spectrumComparison(expMZs, expIntensities,
                    predMZs, predIntensities, 20);
            //System.out.println(specAngle.cosineSimilarity());
            double sim = specAngle.cosineSimilarity();
            if (sim < 0 || sim > 1) {
                System.out.println(sim);
            }
            sa[i] = pep + "\t" + sim + "\t" + tdSingle;

            //calculate deltaRT

        }

        //send to file for analysis/data viz in python/R
        try {
            FileWriter myWriter = new FileWriter("narrow1_rank1_spectralAngles.tsv");
            myWriter.write("peptide" + "\t" + "spectral_angle" + "\t" + "target=1/decoy=0\n");

            for (String s : sa) {
                myWriter.write(s + "\n");
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
