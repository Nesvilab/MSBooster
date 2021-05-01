package Features;

import kotlin.jvm.functions.Function1;
import smile.stat.distribution.KernelDensity;
import umich.ms.datatypes.LCMSDataSubset;
import umich.ms.datatypes.scan.IScan;
import umich.ms.datatypes.scan.StorageStrategy;
import umich.ms.datatypes.scancollection.impl.ScanCollectionDefault;
import umich.ms.fileio.exceptions.FileParsingException;
import umich.ms.fileio.filetypes.mzml.MZMLFile;

import java.io.IOException;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.Map;

import static Features.floatUtils.doubleToFloat;

public class mzMLReader {
    //final Path path;
    final String pathStr;
    transient ScanCollectionDefault scans; //need to implement serializable
    double[] mzFreqs; //this can never be changed once set. It would be difficult if mzFreqs could change, as weighted
                      //similarity measures might be calculated using different weights. If you want to use different
                      //weights, just make new mzmlReader object

    HashMap<Integer, mzmlScanNumber> scanNumberObjects = new HashMap<>();
    private float[] betas;
    public ArrayList<Float>[] RTbins;
    private Function1<Double, Double> LOESS;
    public double[][] expAndPredRTs;

    public mzMLReader(String filename) throws FileParsingException {

        // Creating data source
        //path = Paths.get(filename); //
        Path path = Paths.get(filename);
        pathStr = path.toString();
        MZMLFile source = new MZMLFile(pathStr);

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

        createScanNumObjects(); //seems like we always need this anyway
        //this.getMzFreq(); only if we end up using weights
    }

    //for deserialized mzmlReaders, need to reload scans
    public void setScans() throws FileParsingException {
        MZMLFile source = new MZMLFile(pathStr);
        scans = new ScanCollectionDefault();
        scans.setDefaultStorageStrategy(StorageStrategy.SOFT);
        scans.isAutoloadSpectra(true);
        scans.setDataSource(source);
        source.setNumThreadsForParsing(null);
        scans.loadData(LCMSDataSubset.MS1_WITH_SPECTRA);
    }

    //get experimental spectra
    public float[] getIntensity(int scanNum) throws FileParsingException {
        IScan scan = scans.getScanByNum(scanNum);
//        if (Constants.basePeakNormalization) { //for getting average matched fragment intensity for DIANN-MSFragger comparison
//            double[] intensities = scan.fetchSpectrum().getIntensities();
//
//            //get max intensity
//            double maxInt = 0.0;
//            for (double d : intensities) {
//                if (d > maxInt) {
//                    maxInt = d;
//                }
//            }
//
//            for (int i = 0; i < intensities.length; i++) {
//                intensities[i] = intensities[i] / maxInt * 100000.0; //arbitrary base peak intensity of 100000
//            }
//        }
        return doubleToFloat(scan.fetchSpectrum().getIntensities());
    }

    public float[] getMZ(int scanNum) throws FileParsingException {
        IScan scan = scans.getScanByNum(scanNum);
        return doubleToFloat(scan.fetchSpectrum().getMZs());
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
                        int binIndex = (int) Math.floor(mz / Constants.binwidth);

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
                double newLeft = Math.max(0, i - Constants.mzFreqWindow);
                double newRight = Math.min(maxKey, i + Constants.mzFreqWindow + 1);
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
        int scanNum = -1;
        IScan scan = scans.getNextScanAtMsLevel(scanNum, 2);
        while (scan != null) {
            scanNum = scan.getNum();
            //double windowStart = scan.getPrecursor().getMzRangeStart();
            scanNumberObjects.put(scanNum, new mzmlScanNumber(this, scanNum, scan.getRt().floatValue()));
            scan = scans.getNextScanAtMsLevel(scanNum, 2);
        }
    }

    public mzmlScanNumber getScanNumObject(int scanNum) {
        return scanNumberObjects.get(scanNum);
    }

    //can consider method for setting single pepxml entry
    public void setPepxmlEntries(pepXMLReader xmlReader, int rank, SpectralPredictionMapper spm) throws AssertionError, Exception {
        String[] peptides = xmlReader.getXMLpeptides();
        int[] tdArray = xmlReader.getTargetOrDecoy();
        int[] scanNums = xmlReader.getScanNumbers();
        String[] escore = xmlReader.getEScore();

        int iterations = scanNums.length;

        for (int i = 0; i < iterations; i++) {
            String pep = peptides[i];
            scanNumberObjects.get(scanNums[i]).setPeptideObject(pep, rank, tdArray[i], escore[i],
                    spm.getMzDict(), spm.getIntensityDict(), spm.getRtDict());
        }
    }

    public void setPinEntries(pinReader pin, SpectralPredictionMapper spm) throws AssertionError, Exception {
        while(pin.next()) {
            scanNumberObjects.get(pin.getScanNum()).setPeptideObject(pin.getPep(), pin.getRank(), pin.getTD(), pin.getEScore(),
                    spm.getMzDict(), spm.getIntensityDict(), spm.getRtDict());
        }
        pin.close();
        pin.reset();
    }

    public void setBetas(SpectralPredictionMapper preds, int RTregressionSize) throws IOException {
        betas = RTFunctions.getBetas(this, preds, RTregressionSize);
    }
    public void setBetas() {
        betas = RTFunctions.getBetas(expAndPredRTs);
    }

    public float[] getBetas() {
        return betas;
    }

    //get normalized RTs for regression
    public void normalizeRTs() {
        assert betas != null;
        for (Map.Entry<Integer, mzmlScanNumber> entry : scanNumberObjects.entrySet()) {
            //might be problematic if actual entry does not change
            mzmlScanNumber scanNum = entry.getValue();
            scanNum.normalizedRT = RTFunctions.normalizeRT(betas, scanNum.RT);

            //now calculate deltaRTs
            for (peptideObj pep : scanNum.peptideObjects) {
                pep.deltaRT = Math.abs(scanNum.normalizedRT - pep.RT);
            }
        }
    }

    public void setRTbins(SpectralPredictionMapper preds) throws IOException {
        RTbins = RTFunctions.RTbins(this, preds);
    }

    public ArrayList<Float>[] getRTbins() {
        return RTbins;
    }

    //apply stats to every scan number baseed on RT bin
    public void propagateRTBinStats() {
        assert RTbins != null;
        float[][] binStats = RTFunctions.characterizeRTbins(RTbins);

        for (Map.Entry<Integer, mzmlScanNumber> entry : scanNumberObjects.entrySet()) {
            mzmlScanNumber scanNum = entry.getValue();

            //get stats based on experimental RT bin
            int idx = Math.round(scanNum.RT);
            float binMean = binStats[idx][0];
            float binStd = binStats[idx][1];

            //also set bin size, for use with uniform prior prob
            scanNum.RTbinSize = RTbins[idx].size();

            //now calculate deltaRTs
            for (peptideObj pep : scanNum.peptideObjects) {
                pep.deltaRTbin = Math.abs(binMean - pep.RT);
                pep.RTzscore = Math.abs(StatMethods.zscore(pep.RT, binMean, binStd));
            }
        }
    }

    public void setKernelDensities() {
        KernelDensity[] kernelDensities = RTFunctions.generateEmpiricalDist(RTbins);
        for (Map.Entry<Integer, mzmlScanNumber> entry : scanNumberObjects.entrySet()) {
            mzmlScanNumber scanNum = entry.getValue();

            //now calculate RT probabilities
            for (peptideObj pep : scanNum.peptideObjects) {
                pep.RTprob = RTFunctions.RTprobability(scanNum.RT, pep.RT, kernelDensities);
            }
        }
    }

    public void setLOESS(SpectralPredictionMapper preds, int RTregressionSize, double bandwidth, int robustIters) {
        expAndPredRTs = RTFunctions.getRTarrays(this, preds, RTregressionSize);
        LOESS = RTFunctions.LOESS(expAndPredRTs, bandwidth, robustIters);
    }

    public double predictLOESS(double expRT) { return LOESS.invoke(expRT); }

    public HashMap<Integer, mzmlScanNumber> getScanNumberObjects(){
        return this.scanNumberObjects;
    }

    public static void main(String[] args) throws FileParsingException, IOException {
    }
}
