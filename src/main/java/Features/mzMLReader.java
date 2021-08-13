package Features;

import kotlin.jvm.functions.Function1;
import smile.stat.distribution.KernelDensity;
import umich.ms.datatypes.LCMSDataSubset;
import umich.ms.datatypes.scan.IScan;
import umich.ms.datatypes.scan.StorageStrategy;
import umich.ms.datatypes.scan.impl.ScanDefault;
import umich.ms.datatypes.scan.props.PrecursorInfo;
import umich.ms.datatypes.scancollection.impl.ScanCollectionDefault;
import umich.ms.datatypes.spectrum.ISpectrum;
import umich.ms.datatypes.spectrum.impl.SpectrumDefault;
import umich.ms.fileio.exceptions.FileParsingException;
import umich.ms.fileio.filetypes.mzml.MZMLFile;

import java.io.IOException;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.*;
import java.util.concurrent.ExecutionException;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Future;

import static Features.StatMethods.LOESS;
import static Features.StatMethods.movingAverage;
import static Features.floatUtils.doubleToFloat;

public class mzMLReader {
    //final Path path;
    final String pathStr;
    public ScanCollectionDefault scans; //need to implement serializable
    double[] mzFreqs; //this can never be changed once set. It would be difficult if mzFreqs could change, as weighted
                      //similarity measures might be calculated using different weights. If you want to use different
                      //weights, just make new mzmlReader object

    HashMap<Integer, mzmlScanNumber> scanNumberObjects = new HashMap<>();
    List<Integer> scanNums;
    private float[] betas;
    public ArrayList<Float>[] RTbins = null;
    public float[][] RTbinStats;
    private Function1<Double, Double> RTLOESS;
    public ArrayList<Float>[][] IMbins = null;
    public float[][][] IMbinStats = new float[IMFunctions.numCharges][2 * Constants.IMbinMultiplier + 1][3];
    private ArrayList<Function1<Double, Double>> IMLOESS = new ArrayList<>();
    public HashMap<String, ArrayList<Float>> peptideIMs = new HashMap<>();
    public double[][] expAndPredRTs;
    private List<Future> futureList = new ArrayList<>(Constants.numThreads);

    public mzMLReader(String filename) throws FileParsingException {

        // Creating data source
        //path = Paths.get(filename); //
        Path path = Paths.get(filename);
        pathStr = path.toString();
        MZMLFile source = new MZMLFile(pathStr);

        scans = new ScanCollectionDefault(true); //combined this with line 52
        // Softly reference spectral data, make it reclaimable by GC
        scans.setDefaultStorageStrategy(StorageStrategy.SOFT);
        // Set it to automatically re-parse spectra from the file if spectra were not
        // yet parsed or were reclaimed to make auto-loading work you'll need to use
        // IScan#fetchSpectrum() method instead of IScan#getSpectrum()
        //scans.isAutoloadSpectra(true); // this is actually the default

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

        Constants.useIM = false;
        createScanNumObjects(); //seems like we always need this anyway

        scanNums = new ArrayList<>(scanNumberObjects.size());
        scanNums.addAll(scanNumberObjects.keySet());
        //this.getMzFreq(); only if we end up using weights
    }

    public mzMLReader(mgfFileReader mgf) throws FileParsingException { //uncalibrated mgf from MSFragger .d search
        pathStr = null;

        scans = new ScanCollectionDefault(true);
        scans.setDefaultStorageStrategy(StorageStrategy.SOFT);
        for (String key : mgf.getRtDict().keySet()) {
            //get scan number
            String[] dotSplit = key.split("\\.");
            int scanNum = Integer.parseInt(dotSplit[dotSplit.length - 2]);

            IScan scan = new ScanDefault(scanNum, mgf.getRtDict().get(key), 2, true);
            scan.setIm(Double.valueOf(mgf.getIMDict().get(key)));
            PrecursorInfo pi = new PrecursorInfo();
            pi.setMzTargetMono(mgf.getMassesDict().get(key));
            scan.setPrecursor(pi);
            ISpectrum spectrum = new SpectrumDefault(floatUtils.floatToDouble(mgf.getMzDict().get(key)),
                    floatUtils.floatToDouble(mgf.getIntensityDict().get(key)), null);
            scan.setSpectrum(spectrum, false);
            scans.addScan(scan);
        }
        mgf.clear();

        Constants.useIM = true;
        createScanNumObjects(); //seems like we always need this anyway

        scanNums = new ArrayList<>(scanNumberObjects.size());
        scanNums.addAll(scanNumberObjects.keySet());
        //this.getMzFreq(); only if we end up using weights
    }

    //for deserialized mzmlReaders, need to reload scans
//    public void setScans() throws FileParsingException {
//        MZMLFile source = new MZMLFile(pathStr);
//        scans = new ScanCollectionDefault();
//        scans.setDefaultStorageStrategy(StorageStrategy.SOFT);
//        scans.isAutoloadSpectra(true);
//        scans.setDataSource(source);
//        source.setNumThreadsForParsing(null);
//        scans.loadData(LCMSDataSubset.MS1_WITH_SPECTRA);
//    }

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
        long startTime = System.nanoTime();
        futureList.clear();

        int scanNum = -1;
        IScan scan = scans.getNextScanAtMsLevel(scanNum, 2);
        if (Constants.useIM) {
            while (scan != null) {
                scanNum = scan.getNum();
                //double windowStart = scan.getPrecursor().getMzRangeStart();
                scanNumberObjects.put(scanNum, new mzmlScanNumber(this, scanNum, scan.getRt().floatValue(),
                        scan.getIm().floatValue()));
                scan = scans.getNextScanAtMsLevel(scanNum, 2);
            }
        } else {
            while (scan != null) {
                scanNum = scan.getNum();
                //double windowStart = scan.getPrecursor().getMzRangeStart();
                scanNumberObjects.put(scanNum, new mzmlScanNumber(this, scanNum, scan.getRt().floatValue(), null));
                scan = scans.getNextScanAtMsLevel(scanNum, 2);
            }
        }
        scans.reset(); //free up memory. Could consider setting to null as well, but idk how to check if it's garbage collected
        long endTime = System.nanoTime();
        long duration = (endTime - startTime);
        System.out.println("createScanNumObjects took " + duration / 1000000 +" milliseconds");
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
                    spm.getMzDict(), spm.getIntensityDict(), spm.getRtDict(), spm.getIMDict());
        }
    }

    public void setPinEntries(pinReader pin, SpectralPredictionMapper spm) throws AssertionError, Exception {
        while(pin.next()) {
            scanNumberObjects.get(pin.getScanNum()).setPeptideObject(pin.getPep(), pin.getRank(), pin.getTD(), pin.getEScore(),
                    spm.getMzDict(), spm.getIntensityDict(), spm.getRtDict(), spm.getIMDict());
        }
        pin.close();
        pin.reset();
    }

    public void setBetas(SpectralPredictionMapper preds, int RTregressionSize) throws IOException {
        betas = RTFunctions.getBetas(this, RTregressionSize);
    }
    public void setBetas() {
        betas = RTFunctions.getBetas(expAndPredRTs);
    }

    public float[] getBetas() {
        return betas;
    }

    //get normalized RTs for regression
    public void normalizeRTs(ExecutorService executorService) throws ExecutionException, InterruptedException {
        if (betas == null) {
            System.out.println("why are betas null?");
        };
        futureList.clear();

        //iterate over this list of scan numbers
        for (int i = 0; i < Constants.numThreads; i++) {
            int start = (int) (scanNumberObjects.size() * (long) i) / Constants.numThreads;
            int end = (int) (scanNumberObjects.size() * (long) (i + 1)) / Constants.numThreads;
            futureList.add(executorService.submit(() -> {
                for (int j = start; j < end; j++) {
                    mzmlScanNumber msn = getScanNumObject(scanNums.get(j));
                    msn.normalizedRT = RTFunctions.normalizeRT(betas, msn.RT);

                    //now calculate deltaRTs
                    for (peptideObj pep : msn.peptideObjects) {
                        pep.deltaRT = Math.abs(msn.normalizedRT - pep.RT);
                    }
                }
            }));
        }
        for (Future future : futureList) {
            future.get();
        }
    }

    public void setRTbins() throws IOException {
        RTbins = RTFunctions.RTbins(this);
    }

    public void setIMbins() throws IOException {
        IMbins = IMFunctions.IMbins(this);
    }

    public ArrayList<Float>[] getRTbins() {
        return RTbins;
    }

    public void calculateBinStats(String mode) {
        if (mode.equals("RT")) {
            if (RTbins == null) {
                System.out.println("why are RTbins null?");
            }
            RTbinStats = StatMethods.characterizebins(RTbins, Constants.RTIQR);

            //smoothing with window of 1
            //could make it at statmethods method
            for (int window = 0; window < RTbinStats.length; window++) {
                float[] b = RTbinStats[window];
                if (b.length < 2) {
                    continue;
                }

                RTbinStats[window] = movingAverage(b, 1);
            }
        } else if (mode.equals("IM")) {
            if (IMbins == null) {
                System.out.println("why are IMbins null?");
            }
            for (int charge = 0; charge < IMbins.length; charge ++) {
                IMbinStats[charge] = StatMethods.characterizebins(IMbins[charge], Constants.IMIQR);

                //smoothing with window of 1
                //could make it at statmethods method
                for (int window = 0; window < IMbinStats[charge].length; window++) {
                    float[] b = IMbinStats[charge][window];
                    if (b.length < 2) {
                        continue;
                    }

                    IMbinStats[charge][window] = movingAverage(b, 1);
                }
            }
        }
    }

    public void setRTBinSizes(ExecutorService executorService) throws ExecutionException, InterruptedException {
        futureList.clear();

        //iterate over this list of scan numbers
        for (int i = 0; i < Constants.numThreads; i++) {
            int start = (int) (scanNumberObjects.size() * (long) i) / Constants.numThreads;
            int end = (int) (scanNumberObjects.size() * (long) (i + 1)) / Constants.numThreads;
            futureList.add(executorService.submit(() -> {
                for (int j = start; j < end; j++) {
                    mzmlScanNumber msn = getScanNumObject(scanNums.get(j));

                    //also set bin size, for use with uniform prior prob
                    int idx = (int) (msn.RT * Constants.RTbinMultiplier);
                    msn.RTbinSize = RTbins[idx].size();
                }
            }));
        }
        for (Future future : futureList) {
            future.get();
        }
    }

    //apply stats to every scan number baseed on RT bin
    public void calculateDeltaRTbinAndRTzscore(ExecutorService executorService) throws ExecutionException, InterruptedException {
        futureList.clear();

        //iterate over this list of scan numbers
        for (int i = 0; i < Constants.numThreads; i++) {
            int start = (int) (scanNumberObjects.size() * (long) i) / Constants.numThreads;
            int end = (int) (scanNumberObjects.size() * (long) (i + 1)) / Constants.numThreads;
            futureList.add(executorService.submit(() -> {
                for (int j = start; j < end; j++) {
                    mzmlScanNumber msn = getScanNumObject(scanNums.get(j));

                    //get stats based on experimental RT bin
                    int idx = (int) (msn.RT * Constants.RTbinMultiplier);
                    float binMean = RTbinStats[idx][0];
                    float binStd = RTbinStats[idx][1];

                    //also set bin size, for use with uniform prior prob
                    //msn.RTbinSize = RTbins[idx].size();

                    //now calculate deltaRTs
                    for (peptideObj pep : msn.peptideObjects) {
                        pep.deltaRTbin = Math.abs(binMean - pep.RT);
                        pep.RTzscore = Math.abs(StatMethods.zscore(pep.RT, binMean, binStd));
                    }
                }
            }));
        }
        for (Future future : futureList) {
            future.get();
        }
    }

    public void calculateDeltaRTLOESSnormalized(ExecutorService executorService) throws ExecutionException, InterruptedException {
        futureList.clear();

        //iterate over this list of scan numbers
        for (int i = 0; i < Constants.numThreads; i++) {
            int start = (int) (scanNumberObjects.size() * (long) i) / Constants.numThreads;
            int end = (int) (scanNumberObjects.size() * (long) (i + 1)) / Constants.numThreads;
            futureList.add(executorService.submit(() -> {
                for (int j = start; j < end; j++) {
                    mzmlScanNumber msn = getScanNumObject(scanNums.get(j));

                    //get stats based on experimental RT bin
                    int idx = (int) (msn.RT * Constants.RTbinMultiplier);
                    float binStd = RTbinStats[idx][1];
                    float binIqr = RTbinStats[idx][2];

                    //now calculate deltaRTs
                    double LOESSRT = RTLOESS.invoke((double) msn.RT);
                    for (peptideObj pep : msn.peptideObjects) {
                        //pep.deltaRTLOESSnormalized = Math.abs(LOESSRT - pep.RT) / binStd;
                        pep.deltaRTLOESSnormalized = Math.abs(LOESSRT - pep.RT) / binIqr;
                    }
                }
            }));
        }
        for (Future future : futureList) {
            future.get();
        }
    }

    public void setIMBinSizes(ExecutorService executorService) throws ExecutionException, InterruptedException {
        futureList.clear();

        //iterate over this list of scan numbers
        for (int i = 0; i < Constants.numThreads; i++) {
            int start = (int) (scanNumberObjects.size() * (long) i) / Constants.numThreads;
            int end = (int) (scanNumberObjects.size() * (long) (i + 1)) / Constants.numThreads;
            futureList.add(executorService.submit(() -> {
                for (int j = start; j < end; j++) {
                    mzmlScanNumber msn = getScanNumObject(scanNums.get(j));

                    //also set bin size, for use with uniform prior prob
                    int idx = (int) (msn.IM * Constants.IMbinMultiplier);
                    if (msn.peptideObjects.size() > 0) {
                        msn.IMbinSize = IMbins[msn.getPeptideObject(1).charge - 1][idx].size();
                    }
                }
            }));
        }
        for (Future future : futureList) {
            future.get();
        }
    }

    public void calculateDeltaIMLOESSnormalized(ExecutorService executorService) throws ExecutionException, InterruptedException {
        futureList.clear();

        //iterate over this list of scan numbers
        for (int i = 0; i < Constants.numThreads; i++) {
            int start = (int) (scanNumberObjects.size() * (long) i) / Constants.numThreads;
            int end = (int) (scanNumberObjects.size() * (long) (i + 1)) / Constants.numThreads;
            futureList.add(executorService.submit(() -> {
                for (int j = start; j < end; j++) {
                    mzmlScanNumber msn = getScanNumObject(scanNums.get(j));

                    //get stats based on experimental RT bin
                    int idx = (int) (msn.IM * Constants.IMbinMultiplier);

                    //now calculate deltaRTs
                    for (peptideObj pep : msn.peptideObjects) {
                        float binIqr = IMbinStats[pep.charge - 1][idx][2];
                        double LOESSIM = IMLOESS.get(pep.charge - 1).invoke((double) msn.IM);
                        //pep.deltaIMLOESSnormalized = Math.abs(LOESSIM - pep.IM) / binStd;
                        pep.deltaIMLOESSnormalized = Math.abs(LOESSIM - pep.IM) / binIqr;
                    }
                }
            }));
        }
        for (Future future : futureList) {
            future.get();
        }
    }

    public void setKernelDensities(ExecutorService executorService, String mode) throws ExecutionException, InterruptedException {
        if (mode.equals("RT")) {
            KernelDensity[] kernelDensities = RTFunctions.generateEmpiricalDist(RTbins);

            long startTime = System.nanoTime();
            futureList.clear();

            //iterate over this list of scan numbers
            for (int i = 0; i < Constants.numThreads; i++) {
                int start = (int) (scanNumberObjects.size() * (long) i) / Constants.numThreads;
                int end = (int) (scanNumberObjects.size() * (long) (i + 1)) / Constants.numThreads;
                futureList.add(executorService.submit(() -> {
                    for (int j = start; j < end; j++) {
                        mzmlScanNumber msn = getScanNumObject(scanNums.get(j));

                        for (peptideObj pep : msn.peptideObjects) {
                            pep.RTprob = StatMethods.probability(msn.RT * Constants.RTbinMultiplier, pep.RT, kernelDensities);
                        }
                    }
                }));
            }
            for (Future future : futureList) {
                future.get();
            }
            long endTime = System.nanoTime();
            long duration = (endTime - startTime);
            System.out.println("Calculating RT probabilities took " + duration / 1000000000 + " seconds");
        } else if (mode.equals("IM")) {
            KernelDensity[][] kernelDensities = new KernelDensity[IMFunctions.numCharges][];
            for (int c = 0; c < IMFunctions.numCharges; c++) {
                kernelDensities[c] = RTFunctions.generateEmpiricalDist(IMbins[c]);
            }

//            long startTime = System.nanoTime();
            futureList.clear();

            //iterate over this list of scan numbers
            for (int i = 0; i < Constants.numThreads; i++) {
                int start = (int) (scanNumberObjects.size() * (long) i) / Constants.numThreads;
                int end = (int) (scanNumberObjects.size() * (long) (i + 1)) / Constants.numThreads;
                futureList.add(executorService.submit(() -> {
                    for (int j = start; j < end; j++) {
                        mzmlScanNumber msn = getScanNumObject(scanNums.get(j));

                        for (peptideObj pep : msn.peptideObjects) {
                            pep.IMprob = StatMethods.probability(msn.IM * Constants.IMbinMultiplier, pep.IM, kernelDensities[pep.charge - 1]);
                        }
                    }
                }));
            }
            for (Future future : futureList) {
                future.get();
            }
//            long endTime = System.nanoTime();
//            long duration = (endTime - startTime);
//            System.out.println("Calculating IM probabilities took " + duration / 1000000000 + " seconds");
        }
    }

    public void setLOESS(int regressionSize, double bandwidth, int robustIters, String mode) {
        if (mode.equals("RT")) {
            expAndPredRTs = RTFunctions.getRTarrays(this, regressionSize);
            RTLOESS = LOESS(expAndPredRTs, bandwidth, robustIters);
        } else if (mode.equals("IM")) {
            double[][][] expAndPredIMs = IMFunctions.getIMarrays(this, regressionSize);
            for (double[][] bins : expAndPredIMs) {
                if (bins[0] == null) {
                    continue;
                }
                IMLOESS.add(LOESS(bins, bandwidth, robustIters));
            }
        } else {
            throw new IllegalArgumentException("only choose RT or IM");
        }
    }

    public void predictRTLOESS(ExecutorService executorService) throws ExecutionException, InterruptedException {
        long startTime = System.nanoTime();
        futureList.clear();

        //iterate over this list of scan numbers
        for (int i = 0; i < Constants.numThreads; i++) {
            int start = (int) (scanNumberObjects.size() * (long) i) / Constants.numThreads;
            int end = (int) (scanNumberObjects.size() * (long) (i + 1)) / Constants.numThreads;
            futureList.add(executorService.submit(() -> {
                for (int j = start; j < end; j++) {
                    mzmlScanNumber msn = getScanNumObject(scanNums.get(j));
                    double LOESSRT = RTLOESS.invoke((double) msn.RT);

                    for (peptideObj pep : msn.peptideObjects) {
                        pep.deltaRTLOESS = Math.abs(LOESSRT - pep.RT);
                    }
                }
            }));
        }
        for (Future future : futureList) {
            future.get();
        }
        long endTime = System.nanoTime();
        long duration = (endTime - startTime);
        System.out.println("Calculating deltaRTLOESS took " + duration / 1000000000 +" seconds");
    }

    public void predictIMLOESS(ExecutorService executorService) throws ExecutionException, InterruptedException {
        long startTime = System.nanoTime();
        futureList.clear();

        //iterate over this list of scan numbers
        for (int i = 0; i < Constants.numThreads; i++) {
            int start = (int) (scanNumberObjects.size() * (long) i) / Constants.numThreads;
            int end = (int) (scanNumberObjects.size() * (long) (i + 1)) / Constants.numThreads;
            futureList.add(executorService.submit(() -> {
                for (int j = start; j < end; j++) {
                    mzmlScanNumber msn = getScanNumObject(scanNums.get(j));

                    for (peptideObj pep : msn.peptideObjects) {
                        double LOESSIM = IMLOESS.get(pep.charge - 1).invoke((double) msn.IM);
                        pep.deltaIMLOESS = Math.abs(LOESSIM - pep.IM);

//                        //min
//                        ArrayList<Float> potentialIMs = peptideIMs.get(pep.name);
//                        ArrayList<Double> deltas = new ArrayList<>(potentialIMs.size());
//                        for (int k = 0; k < potentialIMs.size(); k++) {
//                            deltas.add(Math.abs(IMLOESS.get(pep.charge - 1).invoke((double) potentialIMs.get(k)) - pep.IM));
//                            pep.deltaIMLOESS = Collections.min(deltas);
//                        }

//                        //median
//                        ArrayList<Float> potentialIMs = peptideIMs.get(pep.name);
//                        double med = StatMethods.median(potentialIMs);
//                        double LOESSIM = IMLOESS.get(pep.charge - 1).invoke(med);
//                        pep.deltaIMLOESS = Math.abs(LOESSIM - pep.IM);

//                        //(weighted) mean
//                        ArrayList<Float> potentialIMs = peptideIMs.get(pep.name);
//                        double mean = StatMethods.mean(potentialIMs);
//                        double LOESSIM = IMLOESS.get(pep.charge - 1).invoke(mean);
//                        pep.deltaIMLOESS = Math.abs(LOESSIM - pep.IM);
                    }
                }
            }));
        }
        for (Future future : futureList) {
            future.get();
        }
        long endTime = System.nanoTime();
        long duration = (endTime - startTime);
        System.out.println("Calculating deltaIMLOESS took " + duration / 1000000000 +" seconds");
    }

    public HashMap<Integer, mzmlScanNumber> getScanNumberObjects(){
        return this.scanNumberObjects;
    }

    public static void main(String[] args) throws Exception {
        //code below generates files for analysis/data exploration in jupyter notebook
//        mgfFileReader mgf = new mgfFileReader("C:/Users/kevin/Downloads/proteomics/timsTOF/" +
//                "20180819_TIMS2_12-2_AnBr_SA_200ng_HeLa_50cm_120min_100ms_11CT_3_A1_01_2769_uncalibrated.mgf");
//        System.out.println("done loading mgf");
//        mzMLReader mzml = new mzMLReader(mgf);
//        SpectralPredictionMapper spm = new DiannSpeclibReader("C:/Users/kevin/Downloads/proteomics/timsTOF/DIANN.predicted.bin");
//        //detectMap dm = new detectMap("C:/Users/kevin/Downloads/proteomics/timsTOF/detect_Predictions.txt");
//        System.out.println("Loading pin");
//        pinReader pin = new pinReader("C:/Users/kevin/Downloads/proteomics/timsTOF/" +
//                "20180819_TIMS2_12-2_AnBr_SA_200ng_HeLa_50cm_120min_100ms_11CT_3_A1_01_2769.pin");
//        System.out.println("Setting pin");
//        mzml.setPinEntries(pin, spm);
//        mzml.setBetas(spm, 5000);
//        Constants.numThreads = 11;
//        ExecutorService executorService = Executors.newFixedThreadPool(Constants.numThreads);
//        mzml.setLOESS(Constants.IMregressionSize, Constants.bandwidth, Constants.robustIters, "RT");
//        mzml.predictRTLOESS(executorService);
//        mzml.setLOESS(Constants.IMregressionSize, Constants.bandwidth, Constants.robustIters, "IM");
//        mzml.predictIMLOESS(executorService);
//        executorService.shutdown();
//        BufferedWriter writer = new BufferedWriter(new FileWriter("C:/Users/kevin/Downloads/proteomics/timsTOF/IMdccs.tsv"));
//        writer.write("expIM\tpredIM\tlabel\tevalue\tscanNum\tcharge\texpRT\tpredRT\tdeltaIMloess\tpep\tmass\n");
//        for (mzmlScanNumber msn : mzml.scanNumberObjects.values()) {
//            if (msn.peptideObjects.size() == 0) {
//                continue;
//            }
//            peptideObj pObj = msn.getPeptideObject(1);
//            Float expIM = msn.IM;
//            Float predIM = pObj.IM;
//            int td = pObj.targetORdecoy;
//            String escore = pObj.escore;
//            int scanNum = pObj.scanNum;
//            int charge = Integer.parseInt(pObj.name.split("\\|")[2]);
//            float expRT = msn.RT;
//            float predRT = pObj.RT;
//            double deltaIMloess = pObj.deltaIMLOESS;
//            String pep = pObj.name;
//            double mass = mzml.scans.getScanByNum(msn.scanNum).getPrecursor().getMzTargetMono(); //comment out scans.reset() in createScanNumObjects
//            writer.write(expIM + "\t" + predIM + "\t" + td + "\t" + escore + "\t" + scanNum + "\t" + charge + "\t"
//                    + expRT + "\t" + predRT + "\t" + deltaIMloess + "\t" + pep +  "\t" + mass + "\n");
//        }
//        writer.close();
//
//        //plot points of loess
//        for (int i = 1; i < 6; i++) {
//            BufferedWriter writer = new BufferedWriter(new FileWriter("C:/Users/kevin/Downloads/proteomics/timsTOF/IM" + i + ".tsv"));
//            Function1<Double, Double> loess = mzml.IMLOESS.get(i - 1);
//            double exp = 0.679;
//            while (exp < 1.488) {
//                writer.write(exp + "\t" + loess.invoke(exp) + "\n");
//                exp += 0.001;
//            }
//            writer.close();
//        }
//
//        BufferedWriter writer = new BufferedWriter(new FileWriter("C:/Users/kevin/Downloads/proteomics/timsTOF/RTLOESS.tsv"));
//        double exp = 50;
//        while (exp < 6800) {
//            writer.write(exp + "\t" + mzml.RTLOESS.invoke(exp) + "\n");
//            exp += 1;
//        }
//        writer.close();
//        mgfFileReader mgf = new mgfFileReader("C:/Users/kevin/Downloads/proteomics/timsTOF/" +
//                "20180819_TIMS2_12-2_AnBr_SA_200ng_HeLa_50cm_120min_100ms_11CT_3_A1_01_2769_uncalibrated.mgf");
//        System.out.println("done loading mgf");
//        mzMLReader mzml = new mzMLReader(mgf);
        long startTime = System.nanoTime();
        mzMLReader mzml = new mzMLReader("C:/Users/kevin/OneDriveUmich/proteomics/mzml/wideWindow/23aug2017_hela_serum_timecourse_pool_wide_003.mzML");
        long endTime = System.nanoTime();
        long duration = (endTime - startTime);
        System.out.println("total mzml processing took " + duration / 1000000 +" milliseconds");
    }
}
