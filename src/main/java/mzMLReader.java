import umich.ms.datatypes.LCMSDataSubset;
import umich.ms.datatypes.scan.IScan;
import umich.ms.datatypes.scan.StorageStrategy;
import umich.ms.datatypes.scancollection.impl.ScanCollectionDefault;
import umich.ms.fileio.exceptions.FileParsingException;
import umich.ms.fileio.filetypes.mzml.MZMLFile;

import java.io.IOException;
import java.lang.reflect.InvocationTargetException;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.*;

public class mzMLReader{
    final Path path;
    final ScanCollectionDefault scans;
    double[] mzFreqs; //this can never be changed once set. It would be difficult if mzFreqs could change, as weighted
                      //similarity measures might be calculated using different weights. If you want to use different
                      //weights, just make new mzmlReader object

    HashMap<Integer, mzmlScanNumber> scanNumberObjects = new HashMap<>();
    HashMap<Double, ArrayList<Integer>> windowStartDict = new HashMap<>();
    HashMap<String, ArrayList<Integer>> peptideToScanNums = new HashMap<>();

    public mzMLReader(String filename) throws FileParsingException {

        // Creating data source
        path = Paths.get(filename);
        MZMLFile source = new MZMLFile(path.toString());

        scans = new ScanCollectionDefault();
        // Softly reference spectral data, make it reclaimable by GC
        scans.setDefaultStorageStrategy(StorageStrategy.SOFT);
        // Set it to automatically re-parse spectra from the file if spectra were not
        // yet parsed or were reclaimed to make auto-loading work you'll need to use
        // IScan#fetchSpectrum() method instead of IScan#getSpectrum()
        scans.isAutoloadSpectra(true); // this is actually the default

        // Set our mzXML file as the data source for this scan collection
        scans.setDataSource(source);
        // Set number of threads for multi-threaded parsing.
        // null means use as many cores as reported by Runtime.getRuntime().availableProcessors()
        source.setNumThreadsForParsing(null); // this is actually the default
        // load the meta-data about the whole run, with forced parsing of MS1 spectra
        // as we have enabled auto-loading, then if we ever invoke IScan#fetchSpectrum()
        // on an MS2 spectrum, for which the spectrum has not been parsed, it will be
        // obtained from disk automatically. And because of Soft referencing, the GC
        // will be able to reclaim it.
        scans.loadData(LCMSDataSubset.MS1_WITH_SPECTRA);
    }

    //get experimental spectra
    public double[] getIntensity(int scanNum) throws FileParsingException {
        IScan scan = scans.getScanByNum(scanNum);
        return scan.fetchSpectrum().getIntensities();
    }

    public double[] getMZ(int scanNum) throws FileParsingException {
        IScan scan = scans.getScanByNum(scanNum);
        return scan.fetchSpectrum().getMZs();
    }

    //getting fragment frequency distribution
    //final output might be edited to optionally add something to multiply at the end
    //this might be useful since the mz intensities are on a different scale
    public double[] getMzFreq() {
        if (mzFreqs != null) {
            return this.mzFreqs;
        } else {
            int scanLimit = scans.getScanCount();
            int scanNum = 0;
            HashMap<Integer, int[]> mzCounts = new HashMap<>();

            while (scanNum < scanLimit) {
                try {
                    IScan ms2Scan = scans.getNextScanAtMsLevel(scanNum, 2);

                    //increase count by 1
                    double[] mzs = ms2Scan.fetchSpectrum().getMZs();
                    for (double mz : mzs) {
                        int binIndex = (int) Math.floor(mz / constants.binwidth);

                        int[] value = mzCounts.get(binIndex);
                        if (value == null) {
                            mzCounts.put(binIndex, new int[]{1});
                        } else {
                            value[0]++;
                        }
                    }

                    //increase scan number
                    scanNum = ms2Scan.getNum();
                } catch (Exception e) {
                    break;
                }
            }

            //find max binIndex
            int maxKey = Collections.max(mzCounts.keySet());

            //create list
            int[] countsList = new int[maxKey];
            for (Map.Entry<Integer, int[]> entry : mzCounts.entrySet()) {
                int binIndex = entry.getKey();
                int counts = entry.getValue()[0];

                countsList[binIndex - 1] = counts;
            }

            //sliding window average
            double[] averagedCountsList = new double[maxKey];
            for (int i = 0; i < maxKey; i++) {
                double newLeft = Math.max(0, i - constants.mzFreqWindow);
                double newRight = Math.min(maxKey, i + constants.mzFreqWindow + 1);
                double sum = 0;

                for (int j = (int) newLeft; j < newRight; j++) {
                    sum += countsList[j];
                }
                double avg = sum / (newRight - newLeft);
                if (avg > 0) {
                    averagedCountsList[i] = 1 / avg;
                } else {
                    averagedCountsList[i] = 0; //fragment never detected, just ignore
                }
            }
            this.mzFreqs = averagedCountsList;
            return averagedCountsList;
        }
    }

    public void createScanNumObjects() throws FileParsingException {
        HashMap<Integer, mzmlScanNumber> scanMap = new HashMap<>();
        int scanNum = -1;
        IScan scan = scans.getNextScanAtMsLevel(scanNum, 2);
        while (scan != null) {
            scanNum = scan.getNum();
            double windowStart = scan.getPrecursor().getMzRangeStart();
            scanMap.put(scanNum, new mzmlScanNumber(this, scanNum, windowStart));
            scan = scans.getNextScanAtMsLevel(scanNum, 2);

            //for later use in getting spectra from same window
            if (windowStartDict.containsKey(windowStart)) {
                windowStartDict.get(windowStart).add(scanNum);
            } else {
                windowStartDict.put(windowStart, new ArrayList<Integer>(Arrays.asList(scanNum)));
            }
        }
        scanNumberObjects = scanMap;
    }

    public mzmlScanNumber getScanNumObject(int scanNum) {
        return scanNumberObjects.get(scanNum);
    }

    //can consider method for setting single pepxml entry
    public void setPepxmlEntries(pepXMLReader xmlReader, int rank,
                                 HashMap<String, double[]> allPredMZs, HashMap<String, double[]> allPredIntensities) {
        String[] peptides = xmlReader.getXMLpeptides();
        int[] tdArray = xmlReader.getTargetOrDecoy();
        int[] scanNums = xmlReader.getScanNumbers();
        int iterations = scanNums.length;

        for (int i = 0; i < iterations; i++) {
            try {
                int scanNum = scanNums[i];
                String pep = peptides[i];
                scanNumberObjects.get(scanNum).setPeptideObject(pep, rank, tdArray[i],
                        allPredMZs, allPredIntensities);

                if (peptideToScanNums.containsKey(pep)) {
                    peptideToScanNums.get(pep).add(scanNum);
                } else {
                    peptideToScanNums.put(pep, new ArrayList<Integer>(Arrays.asList(scanNum)));
                }

            } catch(Exception e) {
                System.out.println("failed for " + peptides[i]);
            }
        }
    }

    public static void getRankChanges(String[] args) throws FileParsingException, IOException, NoSuchMethodException,
            IllegalAccessException, InvocationTargetException {
        //get predictions
        System.out.println("getting predictions");
        mgfFileReader mgf = new mgfFileReader("preds/");

        //initialize counter
        HashMap<String, Integer> rankCounts = new HashMap<>();
        for (String key : new String[]{"11", "00", "10", "01"}) {
            rankCounts.put(key, 0);
        }

        //iterate over pepxml files
        for (int window = 1; window < 7; window++) {
            System.out.println("working on file " + window);

            mzMLReader mzml = new mzMLReader("C:/Users/kevin/OneDriveUmich/proteomics/mzml/" +
                    "23aug2017_hela_serum_timecourse_4mz_narrow_" + window +".mzml");
            mzml.createScanNumObjects();

            //given all ranked pepxml files associated with mzml, add all entries
            for (int i = 1; i < 5; i++) {
                String path = "C:/Users/kevin/OneDriveUmich/proteomics/pepxml/rank" +
                        i + "/23aug2017_hela_serum_timecourse_4mz_narrow_" + window + "_rank" + i + ".pepXML";
                System.out.println("setting pepXML entries for " + path);
                pepXMLReader xmlReader = new pepXMLReader(path);

                mzml.setPepxmlEntries(xmlReader, i, mgf.allPredMZs, mgf.allPredIntensities);
            }

            //read one pepXML
            String path = "C:/Users/kevin/OneDriveUmich/proteomics/pepxml/rank1/" +
                    "23aug2017_hela_serum_timecourse_4mz_narrow_" + window + "_rank1.pepXML";
            System.out.println("counting rank changes for " + path);
            pepXMLReader xmlReader = new pepXMLReader(path);
            int[] nums = xmlReader.getScanNumbers();

            //add to counter
            for (int num : nums) {
                int ogRank = mzml.getScanNumObject(num).targetDecoyOrder()[0];
                ArrayList x = mzml.getScanNumObject(num).targetDecoyOrder("dotProduct");
                int newRank = (int) x.get(0);

                String key = String.valueOf(ogRank) + String.valueOf(newRank);
                int oldCount = rankCounts.get(key);
                rankCounts.put(key, oldCount + 1);
            }
        }

        System.out.println(rankCounts.entrySet());
    }

    public static void main(String[] args) throws FileParsingException, IOException {
        mzMLReader mzml = new mzMLReader("C:/Users/kevin/OneDriveUmich/proteomics/mzml/" +
                "23aug2017_hela_serum_timecourse_4mz_narrow_1.mzml");
        mzml.createScanNumObjects();

        mgfFileReader mgf = new mgfFileReader("C:/Users/kevin/Downloads/proteomics/pDeep3preds.mgf");
        pepXMLReader xmlReader = new pepXMLReader("C:/Users/kevin/OneDriveUmich/proteomics/pepxml/rank1/" +
                "23aug2017_hela_serum_timecourse_4mz_narrow_1_rank1.pepXML");

        mzml.setPepxmlEntries(xmlReader, 1, mgf.allPredMZs, mgf.allPredIntensities);
        System.out.println(mzml.peptideToScanNums);
    }
}
