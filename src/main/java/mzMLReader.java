import umich.ms.datatypes.LCMSDataSubset;
import umich.ms.datatypes.scan.IScan;
import umich.ms.datatypes.scan.StorageStrategy;
import umich.ms.datatypes.scancollection.impl.ScanCollectionDefault;
import umich.ms.fileio.exceptions.FileParsingException;
import umich.ms.fileio.filetypes.mzml.MZMLFile;

import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashMap;
import java.util.Map;

public class mzMLReader{
    Path path;
    ScanCollectionDefault scans;

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
    public double[] getMzFreq(double bin, int window) {
        int scanLimit = scans.getScanCount();
        int scanNum = 0;
        HashMap<Integer, int[]> mzCounts = new HashMap<>();

        while (scanNum < scanLimit) {
            try {
                IScan ms2Scan = scans.getNextScanAtMsLevel(scanNum, 2);

                //increase count by 1
                double[] mzs = ms2Scan.fetchSpectrum().getMZs();
                for (double mz : mzs) {
                    int binIndex = (int) Math.floor(mz / bin);

                    int[] value = mzCounts.get(binIndex);
                    if (value == null) {
                        mzCounts.put(binIndex, new int[] {1});
                    } else {
                        value[0]++;
                    }
                }

                //increase scan number
                scanNum = ms2Scan.getNum();
            } catch(Exception e) {
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
            double newLeft = Math.max(0, i - window);
            double newRight = Math.min(maxKey, i + window + 1);
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
        return averagedCountsList;
    }

    public static void main(String[] args) throws FileParsingException {
        mzMLReader x = new mzMLReader("C://Users/kevin/Documents/proteomics/mzml/" +
                "23aug2017_hela_serum_timecourse_4mz_narrow_1.mzml");
        System.out.println(Arrays.toString(x.getMzFreq(1, 1)));

    }

}
