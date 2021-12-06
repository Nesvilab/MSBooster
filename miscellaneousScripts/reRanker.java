package Features;

import umich.ms.fileio.exceptions.FileParsingException;
import umich.ms.fileio.filetypes.pepxml.jaxb.standard.MsmsPipelineAnalysis;
import umich.ms.fileio.filetypes.pepxml.jaxb.standard.SpectrumQuery;

import javax.xml.bind.JAXBException;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.lang.reflect.InvocationTargetException;

import static org.apache.commons.io.FileUtils.listFiles;

public class reRanker {
    HashMap<String, double[]> predMZs;
    HashMap<String, double[]> predInts;
    mzMLReader mzmlScans;
    String filePrefix;

    //all the info we need about pepXML files, without saving pepXML files
    int numRanks = 0;
    int[] scanNumbers;
    ArrayList<List<SpectrumQuery>> spectrumQueriesList = new ArrayList<>();
    pepXMLReader analysisTemplate;

    public reRanker(String predictionFile, String mzmlFile, String pepxmlFilesDirectory)
            throws IOException, FileParsingException {
        SpectralPredictionMapper spm = SpectralPredictionMapper.createSpectralPredictionMapper(predictionFile);
        mzmlScans = new mzMLReader(mzmlFile);
        mzmlScans.createScanNumObjects(); //can later iterate through objects

        String prefix = formatPrefix(mzmlFile);
        filePrefix = prefix.split("/")[prefix.split("/").length - 1];

        Collection<File> pepxmlFiles = listFiles(new File(pepxmlFilesDirectory), new String[]{"pepXML"}, false);
        for (File pepxmlFile : pepxmlFiles) {
            String fileName = pepxmlFile.getCanonicalPath();

            if (fileName.contains(filePrefix)) {
                System.out.println("reading in " + fileName);
                numRanks += 1;

                int rank = Integer.parseInt(fileName.substring(fileName.indexOf("rank") + 4).split("\\.")[0]);
                pepXMLReader xmlReader = new pepXMLReader(fileName);

                //get template for pepXML files to be written
                if (analysisTemplate == null) {
                    analysisTemplate = xmlReader;
                }

                //set pepXML entries
                mzmlScans.setPepxmlEntries(xmlReader, rank, spm); //make sure ranks set in order

                //obtain spectrumQueries
                List<SpectrumQuery> sq = xmlReader.getSpectrumQueries();
                try {
                    spectrumQueriesList.add(rank - 1, sq); //helpful if rank is lower than previous ranks
                } catch (Exception e) { //indexOutOfBounds just means rank is higher than current available indices
                    spectrumQueriesList.add(sq); //just add to end of arraylist
                }
            }
        }

        //obtain scanNumbers. First spectrum query list has all scanNums
        scanNumbers = pepXMLReader.getScanNumbers(spectrumQueriesList.get(0));
    }

    private static String formatPrefix(String x) {
        String[] periodSplits = x.split("\\.");
        StringJoiner joiner = new StringJoiner(".");
        for (int i = 0; i < periodSplits.length - 1; i++) {
            joiner.add(periodSplits[i]);
        }
        return joiner.toString();
    }

    public void writeRerankedPepxmlFiles(String similarityMeasure, String directoryName)
            throws NoSuchMethodException, FileParsingException, IllegalAccessException, InvocationTargetException,
            JAXBException {

        //set up indices for querying spectrumQueries
        int[] sqIndices = new int[numRanks];

        //save an MsmsPipelineAnalysis for each new ranked file to be written
        MsmsPipelineAnalysis[] newPepxmlFiles = new MsmsPipelineAnalysis[numRanks];
        for (int i = 0; i < numRanks; i++) {
            newPepxmlFiles[i] = analysisTemplate.freshMsmsPipelineAnalysis();
        }

        //have newSpectrumQueryLists to store
        ArrayList<List<SpectrumQuery>> newSQlists = new ArrayList<>();
        for (int i = 0; i < numRanks; i++) {
            newSQlists.add(new ArrayList<>());
        }

        //work with one scan at a time
        for (int scanNum : scanNumbers) {

            //get new ranks. Can also use newRanks to replace expect score
            Object[] tmpOrders = mzmlScans.getScanNumObject(scanNum).reRankOrder(similarityMeasure);
            ArrayList<Integer> newRanks = (ArrayList<Integer>) tmpOrders[0];

            //need a check if only 1 rank
            for (int i = 0; i < newRanks.size(); i++) {
                int rank = newRanks.get(i);

                //within the rank's spectrumQueries, get the next spectrumQuery
                SpectrumQuery sq = spectrumQueriesList.get(rank - 1).get(sqIndices[rank - 1]);

                //add PSM to new ranked pepXMLs, based on ranking
                //if I want to modify hit rank in pepxml, modify sq
                newSQlists.get(i).add(sq);

                //now use next PSM when we use spectrumQuery from old pepXML again
                sqIndices[rank - 1] += 1;
            }
        }

        //set spectrumQueryLists from newSQlists into MsmsPipelineAnalyses
        //can optimize this if I figure out how to clone the analysis template
        for (int i = 0; i < newPepxmlFiles.length; i++) {
            newPepxmlFiles[i].getMsmsRunSummary().get(0).getSpectrumQuery().addAll(newSQlists.get(i));
        }

        //create folder to write to
        File writeToDirectory = new File(directoryName);
        if (! writeToDirectory.exists()) {
            writeToDirectory.mkdirs();
        }

        //write to files
        for (int i = 0; i < newPepxmlFiles.length; i++) {
            String rankDirectory = directoryName + "/rank" + String.valueOf(i + 1) + "/";
            File rankFile = new File(rankDirectory);
            if (! rankFile.exists()) {
                rankFile.mkdir();
            }
            String finalName = rankDirectory + filePrefix + ".pepXML";
            System.out.println("writing to " + finalName);
            pepXMLWriter.writepepXML(newPepxmlFiles[i], finalName);
        }
    }

    //method to see how often a decoy in the first rank gets replaced with target at first rank after rescoring
    //modify so only those that have target and decoy get scored
    public void getRankChanges(String similarityMeasure, String outfile) throws NoSuchMethodException, FileParsingException,
            IllegalAccessException, InvocationTargetException, IOException {

        FileWriter myWriter = new FileWriter(Constants.outputFolder + outfile);
        myWriter.write("change" + "\t" + "oldScore" + "\t" + "newScore\n");

        //initialize counter
        HashMap<String, Integer> rankCounts = new HashMap<>();
        for (String key : new String[]{"11", "00", "10", "01"}) {
            rankCounts.put(key, 0);
        }

        //add to counter
        for (Map.Entry<Integer, mzmlScanNumber> scanNumObj : mzmlScans.scanNumberObjects.entrySet()) {

            mzmlScanNumber scanNumVal = scanNumObj.getValue();
            if (scanNumVal.peptideObjects == null) {
                //no PSMs associated with it
                continue;
            }

            if (scanNumVal.peptideObjects.size() < 2) {
                //not enough PSMs associated with it
                continue;
            }

            int[] OGranks = scanNumVal.targetDecoyOrder();

            //make sure we can make appropriate comparison
            boolean contains0 = false;
            boolean contains1 = false;
            for (int r : OGranks) {
                if (r == 0) {
                    contains0 = true;
                } else { //must be 1
                    contains1 = true;
                }
            }
            if (contains0 && contains1) {
                int newRank = (int) scanNumVal.targetDecoyOrder(similarityMeasure).get(0);

                String key = String.valueOf(OGranks[0]) + String.valueOf(newRank);
                int oldCount = rankCounts.get(key);
                rankCounts.put(key, oldCount + 1);

                //collect score diff
                if (key.equals("01") || key.equals("10")) {
                    double oldScore = scanNumVal.getPeptideObject(1).getScore(similarityMeasure);
                    double newScore = -1.0;
                    for (peptideObj pObj : scanNumVal.peptideObjects) {
                        double ns = pObj.getScore(similarityMeasure);
                        if (ns > newScore) {
                            newScore = ns;
                        }
                    }
                    myWriter.write(key + "\t" + oldScore + "\t" + newScore + "\n");
                }
            }
        }

        myWriter.close();
        System.out.println(rankCounts.entrySet());
    }

    public static void main(String[] args) throws IOException, FileParsingException, InvocationTargetException,
            NoSuchMethodException, JAXBException, IllegalAccessException {
//        reRanker x = new reRanker("C:/Users/kevin/Downloads/proteomics/wideWindowPreds.mgf",
//                "C:/Users/kevin/Downloads/proteomics/wideWindow/" +
//                        "23aug2017_hela_serum_timecourse_pool_wide_003.mzML",
//                "C:/Users/kevin/Downloads/proteomics/wideWindow");
//        x.writeRerankedPepxmlFiles("brayCurtis",
//                "C:/Users/kevin/Downloads/proteomics/wideWindow/reranked");
        reRanker x = new reRanker("C:/Users/kevin/Downloads/proteomics/wideWindowPreds.mgf",
                "C:/Users/kevin/OneDriveUmich/proteomics/mzml/wideWindow/" +
                        "23aug2017_hela_serum_timecourse_pool_wide_001_170829031834.mzml",
                "C:/Users/kevin/Downloads/proteomics/wideWindow");
        x.getRankChanges("cosineSimilarity", "changes.tsv");
    }
}
