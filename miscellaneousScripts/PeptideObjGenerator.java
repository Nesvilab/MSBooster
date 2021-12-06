package Features;

import umich.ms.fileio.exceptions.FileParsingException;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.lang.reflect.InvocationTargetException;
import java.util.ArrayList;
import java.util.Map;

public class PeptideObjGenerator {

    public static ArrayList<peptideObjLite> generatePeptideObj(String predFile, String pepXmlFiles, String mzmlFiles,
                                            String similarityMeasure)
            throws IOException, FileParsingException, NoSuchMethodException, IllegalAccessException, InvocationTargetException {
        long startTime = System.nanoTime(); // timer: https://stackoverflow.com/questions/3382954/measure-execution-time-for-a-java-method

        ////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        //set up arraylist to hold peptide objects
        ArrayList<peptideObjLite> pepObjContainer = new ArrayList<>();

        ////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        //get predictions by loading in predicted spectra
        SpectralPredictionMapper predictedSpectra = SpectralPredictionMapper.createSpectralPredictionMapper(predFile);

        ////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        //get paths for pepXML and mzML files
        File pepXMLDirectory = new File(pepXmlFiles);
        String[] pepXMLfiles = pepXMLDirectory.list((d, name) -> name.endsWith(".pepXML")); //not recursive

        //at some point, need to differentiate between single mzml and multiple
        String[] pepXMLroot = new String[pepXMLfiles.length];
        for (int i = 0; i < pepXMLfiles.length; i++) {
            pepXMLroot[i] = pepXMLfiles[i].split("_rank")[0];
        }

        File mzmlDirectory = new File(mzmlFiles);
        String[] mzmlfiles = mzmlDirectory.list((d, name) -> name.endsWith(".mzML"));
        for (String mFile : mzmlfiles) {
            mzMLReader mReader = new mzMLReader(mzmlFiles + "/" + mFile);
            mReader.createScanNumObjects(); //initiate
            System.out.println(mFile);

            //add all pepXML files from this mzml
            for (int pepIdx = 0; pepIdx < pepXMLfiles.length; pepIdx++) {
                if ((pepXMLroot[pepIdx] + ".mzML").equals(mFile)) {
                    pepXMLReader pReader = new pepXMLReader(pepXmlFiles + "/" + pepXMLfiles[pepIdx]);

                    //get rank of pepXML file
                    String[] ranks = pReader.pathname.split("\\.");
                    String rank = ranks[ranks.length - 2];
                    ranks = rank.split("rank");
                    int rankNum = Integer.parseInt(ranks[ranks.length - 1]);

                    //set peptide objects
                    mReader.setPepxmlEntries(pReader, rankNum, predictedSpectra); //make sure ranks are set in order

                    //iterate through the peptide objects again to add additional information
                    double[] massDiffs = pReader.getMassDiff();
                    int[] ntts = pReader.getNTT();
                    int[] nmcs = pReader.getNMCs();
                    int[] isoErrs = pReader.getIsotopeErrs();
                    int[] scanNums = pReader.getScanNumbers();

                    int iteration = 0;
                    for (int i : scanNums) {
                        peptideObj p = mReader.scanNumberObjects.get(i).getPeptideObject(rankNum);
//                        p.massDiff = massDiffs[iteration];
//                        p.ntt = ntts[iteration];
//                        p.nmc = nmcs[iteration];
//                        p.isoError = isoErrs[iteration];
                        p.setScore(similarityMeasure);
                        iteration++;
                    }

                    //concurrent modification exception
                    //pepXMLMap.remove(entry.getKey());
                }
            }

            //calculate betas and assigning delta RTs
            mReader.setBetas(predictedSpectra, Constants.RTregressionSize);
            mReader.normalizeRTs();

            //implement double counting here. Recalculate similarity score

            //retain peptide objects by converting them to lower memory peptideObjLite
            for (Map.Entry<Integer, mzmlScanNumber> entry : mReader.scanNumberObjects.entrySet()) {
                for (peptideObj entry2 : entry.getValue().peptideObjects) {
                    pepObjContainer.add(new peptideObjLite(entry2));
                }
            }
        }

        ////////////////////////////////////////////////////////////////////////////////////////////////////////////////

        System.out.println(pepObjContainer.size());
        long stopTime = System.nanoTime();
        System.out.println((stopTime - startTime) / 1000000000.0);
        return pepObjContainer;
    }

    public static void peptideObjWriter(String outfile, ArrayList<peptideObjLite> pObjs, String similarityMeasure)
            throws IOException {
        FileWriter myWriter = new FileWriter(outfile);
        myWriter.write("name\trank\tscanNum\tmzml\ttd\tescore\tdeltaRT\tmassdiff\tntt\tnmc\tisoError\tsim\n");
        for (peptideObjLite pObj : pObjs) {
            myWriter.write(pObj.name + "\t" + pObj.rank + "\t" + pObj.scanNum + "\t" + pObj.mzml + "\t" +
                    pObj.targetORdecoy + "\t" + pObj.escore + "\t" + pObj.deltaRT + "\t" + pObj.massDiff + "\t" +
                    pObj.ntt + "\t" + pObj.nmc + "\t" + pObj.isoError + "\t" + pObj.scores.get(similarityMeasure) + "\n");
        }
        myWriter.close();
    }

    public static void main(String[] args) throws IOException, FileParsingException, NoSuchMethodException,
            IllegalAccessException, InvocationTargetException {
//        createSimilarityFile("C:/Users/kevin/Downloads/proteomics/HeLa_preds_allmods.mgf",
//                "C:/Users/kevin/OneDriveUmich/proteomics/pepxml/rank1/",
//                "C:/Users/kevin/OneDriveUmich/proteomics/mzml/",
//                "brayCurtis",
//                "C:/Users/kevin/Downloads/proteomics/pDeep2_allMods_similarities.tsv");
        ArrayList<peptideObjLite> pObjs =
                generatePeptideObj("C:/Users/kevin/Downloads/proteomics/wideWindowPreds.mgf",
                "C:/Users/kevin/Downloads/proteomics/wideWindow",
                "C:/Users/kevin/Downloads/proteomics/wideWindow",
                "brayCurtis");
        peptideObjWriter(Constants.outputFolder + "/pepObjs.tsv", pObjs, "brayCurtis");
    }
}
