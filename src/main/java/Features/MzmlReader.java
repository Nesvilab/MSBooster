/*
 * This file is part of MSBooster.
 *
 * MSBooster is free software: you can redistribute it and/or modify it under
 * the terms of the GNU Lesser General Public License as published by the Free
 * Software Foundation, either version 3 of the License, or (at your
 * option) any later version.
 *
 * MSBooster is distributed in the hope that it will be useful, but WITHOUT
 * ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
 * FITNESS FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public License
 * for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with MSBooster. If not, see <https://www.gnu.org/licenses/>.
 */

package Features;

import kotlin.jvm.functions.Function1;
import org.checkerframework.checker.units.qual.A;
import umich.ms.datatypes.LCMSDataSubset;
import umich.ms.datatypes.scan.IScan;
import umich.ms.datatypes.scan.StorageStrategy;
import umich.ms.datatypes.scancollection.impl.ScanCollectionDefault;
import umich.ms.fileio.exceptions.FileParsingException;
import umich.ms.fileio.filetypes.mzml.MZMLFile;
//import umontreal.ssj.gof.KernelDensity;
import umontreal.ssj.probdist.EmpiricalDist;

import java.io.IOException;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.*;
import java.util.concurrent.ConcurrentHashMap;
import java.util.concurrent.ExecutionException;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Future;

import static Features.FloatUtils.doubleToFloat;
import static Features.StatMethods.*;

public class MzmlReader {
    //final Path path;
    final String pathStr;
    public ScanCollectionDefault scans; //need to implement serializable
    double[] mzFreqs; //this can never be changed once set. It would be difficult if mzFreqs could change, as weighted
                      //similarity measures might be calculated using different weights. If you want to use different
                      //weights, just make new mzmlReader object

    public TreeMap<Integer, MzmlScanNumber> scanNumberObjects = new TreeMap<>();
    List<Integer> scanNums;
    private float[] betas;
    public ArrayList<Float>[] RTbins = null;
    public float[][] RTbinStats;
    public HashMap<String, Function1<Double, Double>> RTLOESS = new HashMap<>();
    public HashMap<String, Function1<Double, Double>> RTLOESS_realUnits = new HashMap<>();
    public int unifPriorSize;
    public float unifProb;
    public int[] unifPriorSizeIM;
    public float[] unifProbIM;
    public ArrayList<Float>[][] IMbins = null;
    public float[][][] IMbinStats = new float[IMFunctions.numCharges][2 * Constants.IMbinMultiplier + 1][3];
    private ArrayList<Function1<Double, Double>> IMLOESS = new ArrayList<>();
    public HashMap<String, ArrayList<Float>> peptideIMs = new HashMap<>();
    public HashMap<String, double[][]> expAndPredRTs;
    public List<Future> futureList = new ArrayList<>(Constants.numThreads);

    public MzmlReader(String filename) throws FileParsingException, ExecutionException, InterruptedException {
        // Creating data source
        //path = Paths.get(filename); //
        System.out.println("Processing " + filename);
        Path path = Paths.get(filename);
        pathStr = path.toString();
        MZMLFile source = new MZMLFile(pathStr);
        source.setExcludeEmptyScans(true);

        scans = new ScanCollectionDefault(true); //combined this with line 52
        // Softly reference spectral data, make it reclaimable by GC
        scans.setDefaultStorageStrategy(StorageStrategy.STRONG);
        // Set it to automatically re-parse spectra from the file if spectra were not
        // yet parsed or were reclaimed to make auto-loading work you'll need to use
        // IScan#fetchSpectrum() method instead of IScan#getSpectrum()
        //scans.isAutoloadSpectra(true); // this is actually the default

        // Set our mzXML file as the data source for this scan collection
        scans.setDataSource(source);
        // Set number of threads for multi-threaded parsing.
        // null means use as many cores as reported by Runtime.getRuntime().availableProcessors()
        source.setNumThreadsForParsing(Constants.numThreads); // this is actually the default
        // load the meta-data about the whole run, with forced parsing of MS1 spectra
        // as we have enabled auto-loading, then if we ever invoke IScan#fetchSpectrum()
        // on an MS2 spectrum, for which the spectrum has not been parsed, it will be
        // obtained from disk automatically. And because of Soft referencing, the GC
        // will be able to reclaim it.
        scans.loadData(LCMSDataSubset.MS2_WITH_SPECTRA);
        Constants.useIM = false;
        scanNumberObjects = new TreeMap<>();
        createScanNumObjects(); //seems like we always need this anyway

        scanNums = new ArrayList<>(scanNumberObjects.size());
        scanNums.addAll(scanNumberObjects.keySet());
        //this.getMzFreq(); only if we end up using weights
    }

    public MzmlReader(MgfFileReader mgf) throws FileParsingException, ExecutionException, InterruptedException { //uncalibrated mgf from MSFragger .d search
        pathStr = mgf.filenames.get(0);
        System.out.println("Processing " + pathStr);

        //Constants.useIM = true;
        //scanNumberObjects = mgf.scanNumberObjects;
        scanNumberObjects.putAll(mgf.scanNumberObjects);
        mgf.clear();

        scanNums = new ArrayList<>(scanNumberObjects.size());
        scanNums.addAll(scanNumberObjects.keySet());
        //this.getMzFreq(); only if we end up using weights
    }

    //get experimental spectra

    //getting fragment frequency distribution
    //final output might be edited to optionally add something to multiply at the end
    //this might be useful since the mz intensities are on a different scale
//    public double[] getMzFreq() {
//        if (mzFreqs != null) {
//            return this.mzFreqs;
//        } else {
//            int scanLimit = scans.getScanCount();
//            int scanNum = 0;
//            HashMap<Integer, int[]> mzCounts = new HashMap<>();
//
//            while (scanNum < scanLimit) {
//                try {
//                    IScan ms2Scan = scans.getNextScanAtMsLevel(scanNum, 2);
//
//                    //increase count by 1
//                    double[] mzs = ms2Scan.fetchSpectrum().getMZs();
//                    for (double mz : mzs) {
//                        int binIndex = (int) Math.floor(mz / Constants.binwidth);
//
//                        int[] value = mzCounts.get(binIndex);
//                        if (value == null) {
//                            mzCounts.put(binIndex, new int[]{1});
//                        } else {
//                            value[0]++;
//                        }
//                    }
//
//                    //increase scan number
//                    scanNum = ms2Scan.getNum();
//                } catch (Exception e) {
//                    break;
//                }
//            }
//
//            //find max binIndex
//            int maxKey = Collections.max(mzCounts.keySet());
//
//            //create list
//            int[] countsList = new int[maxKey];
//            for (Map.Entry<Integer, int[]> entry : mzCounts.entrySet()) {
//                int binIndex = entry.getKey();
//                int counts = entry.getValue()[0];
//
//                countsList[binIndex - 1] = counts;
//            }
//
//            //sliding window average
//            double[] averagedCountsList = new double[maxKey];
//            for (int i = 0; i < maxKey; i++) {
//                double newLeft = Math.max(0, i - Constants.mzFreqWindow);
//                double newRight = Math.min(maxKey, i + Constants.mzFreqWindow + 1);
//                double sum = 0;
//
//                for (int j = (int) newLeft; j < newRight; j++) {
//                    sum += countsList[j];
//                }
//                double avg = sum / (newRight - newLeft);
//                if (avg > 0) {
//                    averagedCountsList[i] = 1 / avg;
//                } else {
//                    averagedCountsList[i] = 0; //fragment never detected, just ignore
//                }
//            }
//            this.mzFreqs = averagedCountsList;
//            return averagedCountsList;
//        }
//    }

    public void createScanNumObjects() throws FileParsingException {
        //for checking resolution
        boolean hasFTMS = false;
        boolean hasITMS = false;

        //get all scan nums
        for (IScan scan : scans.getMapNum2scan().values()) {
            if (scan.getMsLevel() != 1) {
                if (scan.getFilterString() != null) {
                    if (!hasFTMS) {
                        if (scan.getFilterString().contains("FTMS")) {
                            hasFTMS = true;
                        }
                    }
                    if (!hasITMS) {
                        if (scan.getFilterString().contains("ITMS")) {
                            hasITMS = true;
                        }
                    }
                }

                MzmlScanNumber msn = new MzmlScanNumber(scan);
                scanNumberObjects.put(scan.getNum(), msn);
            }
        }

        //what happens with resolution
        //needs to be updated for each mzml
        if (hasFTMS && ! hasITMS) { //only high resoltuion
            Constants.ppmTolerance = Constants.highResppmTolerance;
        } else if (hasITMS) { //low resolution, or both high and low are present so default to low
            Constants.ppmTolerance = Constants.lowResppmTolerance;
        }

        scans.reset(); //free up memory. Could consider setting to null as well, but idk how to check if it's garbage collected
        //long endTime = System.nanoTime();
        //long duration = (endTime - startTime);
        //System.out.println("createScanNumObjects took " + duration / 1000000 +" milliseconds");
    }

    public MzmlScanNumber getScanNumObject(int scanNum) {
        return scanNumberObjects.get(scanNum);
    }

    //can consider method for setting single pepxml entry
//    public void setPepxmlEntries(pepXMLReader xmlReader, int rank, SpectralPredictionMapper spm) throws AssertionError, Exception {
//        String[] peptides = xmlReader.getXMLpeptides();
//        int[] tdArray = xmlReader.getTargetOrDecoy();
//        int[] scanNums = xmlReader.getScanNumbers();
//        String[] escore = xmlReader.getEScore();
//
//        int iterations = scanNums.length;
//
//        for (int i = 0; i < iterations; i++) {
//            String pep = peptides[i];
//            scanNumberObjects.get(scanNums[i]).setPeptideObject(pep, rank, tdArray[i], escore[i],
//                    spm.getPreds());
//        }
//    }

    class setScanNumPepObj implements Runnable {
        private ArrayList<String> lines = new ArrayList<>();
        private int scanNum;
        private ArrayList<PeptideFormatter> peps = new ArrayList<>();
        private ArrayList<Integer> ranks = new ArrayList<>();
        private ArrayList<Integer> tds = new ArrayList<>();
        private ArrayList<String> escores = new ArrayList<>();
        private final ConcurrentHashMap<String, PredictionEntry> allPreds;
        private int specIdx;
        private int pepIdx;
        private int rankIdx;
        private int labelIdx;
        private int eScoreIdx;
        private boolean calcEvalue;
        public setScanNumPepObj(String scanNum, ConcurrentHashMap<String, PredictionEntry> allPreds,
                                int specIdx, int pepIdx, int rankIdx, int labelIdx, int eScoreIdx, boolean calcEvalue) {
            this.scanNum = Integer.parseInt(scanNum);
            this.allPreds = allPreds;
            this.specIdx = specIdx;
            this.pepIdx = pepIdx;
            this.rankIdx = rankIdx;
            this.labelIdx = labelIdx;
            this.eScoreIdx = eScoreIdx;
            this.calcEvalue = calcEvalue;
        }

        //TODO move this to parallel part
        public void add(String line) {
            lines.add(line);
        }

        @Override
        public void run() {
            for (String line : lines) {
                String[] row = line.split("\t");
                String[] periodSplit = row[specIdx].split("\\.");
                peps.add(new PeptideFormatter(row[pepIdx],
                        periodSplit[periodSplit.length - 1].split("_")[0], "pin"));
                try {
                    ranks.add(Integer.parseInt(row[rankIdx]));
                } catch (Exception e) {
                    String[] specIdxSplit = row[specIdx].split("_");
                    ranks.add(Integer.parseInt(specIdxSplit[specIdxSplit.length - 1]));
                }
                tds.add(Math.max(0, Integer.parseInt(row[labelIdx])));
                if (calcEvalue) {
                    escores.add(String.valueOf(Math.exp(15.0 - Double.parseDouble(row[eScoreIdx]))));
                } else {
                    escores.add(String.valueOf(Math.pow(10, Double.parseDouble(row[eScoreIdx]))));
                }
            }

            for (int i = 0; i < peps.size(); i++) {
                scanNumberObjects.get(scanNum).setPeptideObject(peps.get(i), ranks.get(i),
                        tds.get(i), escores.get(i), allPreds, true);
            }
        }
    }
    public void setPinEntries(PinReader pin, SpectralPredictionMapper spm, ExecutorService executorService)
            throws AssertionError, Exception {
        //TODO: multithread?
        ConcurrentHashMap<String, PredictionEntry> allPreds = spm.getPreds();
        ProgressReporter pr = new ProgressReporter(pin.length);
        futureList.clear();

        String currentScanNum = "0";
        setScanNumPepObj task = null;
        int limit = pin.scanNumIdx + 2;
        while (pin.next(false)) {
            try {
                pr.progress();
                //get scanNum as string
                String scanNum = pin.line.split("\t", limit)[pin.scanNumIdx];
                if (!Objects.equals(scanNum, currentScanNum)) {
                    //send previous setScanNumPepObj
                    if (task != null) {
                        executorService.execute(task);
                    }

                    //make new setScanNumPepObj
                    currentScanNum = scanNum;
                    task = new setScanNumPepObj(currentScanNum, allPreds,
                            pin.specIdx, pin.pepIdx, pin.rankIdx, pin.labelIdx, pin.eScoreIdx, pin.calcEvalue);
                }
                //add to it
                task.add(pin.line);
            } catch (Exception e) {
                e.printStackTrace();
                System.exit(1);
            }
        }
        try {
            executorService.execute(task);
        } catch (Exception e) {
            e.printStackTrace();
            System.exit(1);
        }

        //set RT filter
        float maxRT = pin.getRT();
        pin.close();
        pin.reset();

        if (Constants.realMinuteFilter == 10000f) {
            if (Constants.percentRTgradientFilter != 100f) {
                Constants.realMinuteFilter = Constants.percentRTgradientFilter / 100 *
                        maxRT;
                System.out.println("Setting minute filter to " + Constants.realMinuteFilter);
            }
        }
    }

    public void setBetas(SpectralPredictionMapper preds, int RTregressionSize) throws IOException {
        betas = RTFunctions.getBetas(this, RTregressionSize);
    }
    public void setBetas() {
        betas = RTFunctions.getBetas(expAndPredRTs.get(""));
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
                    MzmlScanNumber msn = getScanNumObject(scanNums.get(j));
                    msn.normalizedRT = RTFunctions.normalizeRT(betas, msn.RT);

                    //now calculate deltaRTs
                    for (PeptideObj pep : msn.peptideObjects) {
                        if (pep == null) {
                            break;
                        }
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
            RTbinStats = characterizebins(RTbins, Constants.RTIQR);

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
                IMbinStats[charge] = characterizebins(IMbins[charge], Constants.IMIQR);

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
                    MzmlScanNumber msn = getScanNumObject(scanNums.get(j));

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
                    MzmlScanNumber msn = getScanNumObject(scanNums.get(j));

                    //get stats based on experimental RT bin
                    int idx = (int) (msn.RT * Constants.RTbinMultiplier);
                    float binMean = RTbinStats[idx][0];
                    float binStd = RTbinStats[idx][1];

                    //also set bin size, for use with uniform prior prob
                    //msn.RTbinSize = RTbins[idx].size();

                    //now calculate deltaRTs
                    for (PeptideObj pep : msn.peptideObjects) {
                        if (pep == null) {
                            break;
                        }
                        pep.deltaRTbin = Math.abs(binMean - pep.RT);
                        pep.RTzscore = Math.abs(zscore(pep.RT, binMean, binStd));
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
                HashMap<String, Double> LOESSRT = new HashMap<>();
                for (int j = start; j < end; j++) {
                    MzmlScanNumber msn = getScanNumObject(scanNums.get(j));

                    //get stats based on experimental RT bin
                    int idx = (int) (msn.RT * Constants.RTbinMultiplier);
//                    float binStd = RTbinStats[idx][1];
                    float binIqr = RTbinStats[idx][2];

                    //now calculate deltaRTs
                    for (String mass : RTLOESS.keySet()) {
                        LOESSRT.put(mass, RTLOESS.get(mass).invoke((double) msn.RT));
                    }
                    for (PeptideObj pep : msn.peptideObjects) {
                        if (pep == null) {
                            break;
                        }
                        double finalDelta = Double.MAX_VALUE;
                        for (String mass : LOESSRT.keySet()) {
                            if (pep.name.contains(mass)) {
                                double delta = Math.abs(LOESSRT.get(mass) - pep.RT);
                                if (delta < finalDelta) {
                                    finalDelta = delta;
                                }
                            }
                        }
                        pep.deltaRTLOESSnormalized = finalDelta / binIqr;
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
                    MzmlScanNumber msn = getScanNumObject(scanNums.get(j));

                    //also set bin size, for use with uniform prior prob
                    int idx = (int) (msn.IM * Constants.IMbinMultiplier);
                    if (msn.peptideObjects[0] != null) {
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
                    MzmlScanNumber msn = getScanNumObject(scanNums.get(j));

                    //get stats based on experimental RT bin
                    int idx = (int) (msn.IM * Constants.IMbinMultiplier);

                    //now calculate deltaRTs
                    for (PeptideObj pep : msn.peptideObjects) {
                        if (pep == null) {
                            break;
                        }
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
            EmpiricalDist[] kernelDensities = RTFunctions.generateEmpiricalDist(RTbins);

            //long startTime = System.nanoTime();
            futureList.clear();

            //iterate over this list of scan numbers
            for (int i = 0; i < Constants.numThreads; i++) {
                int start = (int) (scanNumberObjects.size() * (long) i) / Constants.numThreads;
                int end = (int) (scanNumberObjects.size() * (long) (i + 1)) / Constants.numThreads;
                futureList.add(executorService.submit(() -> {
                    for (int j = start; j < end; j++) {
                        MzmlScanNumber msn = getScanNumObject(scanNums.get(j));

                        for (PeptideObj pep : msn.peptideObjects) {
                            if (pep == null) {
                                break;
                            }
                            pep.RTprob = probability(msn.RT * Constants.RTbinMultiplier, pep.RT, kernelDensities);
                        }
                    }
                }));
            }
            for (Future future : futureList) {
                future.get();
            }
            //long endTime = System.nanoTime();
            //long duration = (endTime - startTime);
            //System.out.println("Calculating RT probabilities took " + duration / 1000000 + " milliseconds");
        } else if (mode.equals("IM")) {
            EmpiricalDist[][] kernelDensities = new EmpiricalDist[IMFunctions.numCharges][];
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
                        MzmlScanNumber msn = getScanNumObject(scanNums.get(j));

                        for (PeptideObj pep : msn.peptideObjects) {
                            if (pep == null) {
                                break;
                            }
                            pep.IMprob = probability(msn.IM * Constants.IMbinMultiplier, pep.IM, kernelDensities[pep.charge - 1]);
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

    public void setLOESS(int regressionSize, String bandwidth, int robustIters, String mode) {
        if (mode.equals("RT")) {
            expAndPredRTs = RTFunctions.getRTarrays(this, regressionSize);
            ArrayList<String> masses = new ArrayList<>();
            if (Constants.RTmassesForCalibration.isEmpty()) {
                masses.add("");
            } else {
                masses.addAll(Arrays.asList(Constants.RTmassesForCalibration.split(",")));
                masses.add("others");
            }
            for (String mass : masses) {
                String[] bandwidths = bandwidth.split(",");
                float[] bestBandwidths = new float[Constants.regressionSplits];
                double[][] rts = expAndPredRTs.get(mass);
                if (rts[0].length < 50) {
                    RTLOESS.put(mass, null);
                    RTLOESS_realUnits.put(mass, null);
                    continue;
                }

                //divide into train and test sets
                ArrayList<double[][][]> splits = trainTestSplit(rts);

                for (int Nsplit = 0; Nsplit < splits.size(); Nsplit++) {
                    System.out.println("Iteration " + (Nsplit + 1));
                    double bestMSE = Double.MAX_VALUE;
                    float bestBandwidth = 1f;

                    double[][][] split = splits.get(Nsplit);
                    double[][] train = split[0];
                    double[][] test = split[1];

                    for (String b : bandwidths) { //for bandwidth in grid search
                        float floatb = Float.valueOf(b);

                        //get the loess model
                        try {
                            Function1<Double, Double> loess = LOESS(train, floatb, robustIters);

                            //calculate MSE by comparing calibrated expRT to predRT
                            double[] calibratedRTs = new double[test[0].length];
                            for (int i = 0; i < calibratedRTs.length; i++) {
                                double rt = test[0][i];
                                calibratedRTs[i] = loess.invoke(rt);
                            }
                            double mse = meanSquaredError(calibratedRTs, test[1]);

                            //choose best model
                            if (mse < bestMSE) {
                                bestMSE = mse;
                                bestBandwidth = floatb;
                            }
                        } catch (Exception e) {
                            System.out.println("Bandwidth " + floatb + " failed. Moving on");
                        } //bandwidth too small?
                    }
                    bestBandwidths[Nsplit] = bestBandwidth;
                }
                float finalBandwidth = Float.parseFloat(String.format("%.4f", mean(bestBandwidths)));
                System.out.println("Best average bandwidth for mass " + mass + " from grid search of " +
                        Constants.rtBandwidth + " after " + Constants.regressionSplits + " iterations is " + finalBandwidth);

                //final model trained on all data
                while (true) {
                    try {
                        RTLOESS.put(mass, LOESS(rts, finalBandwidth, robustIters));
                        double[][] reverseRts = new double[2][];
                        reverseRts[0] = rts[1];
                        reverseRts[1] = rts[0];
                        RTLOESS_realUnits.put(mass, LOESS(reverseRts, finalBandwidth, robustIters));
                        break;
                    } catch (Exception e) {
                        if (finalBandwidth == 1) {
                            System.out.println("Regression still not possible with bandwidth 1. Setting RT score to 0");
                            RTLOESS.put(mass, null);
                            RTLOESS_realUnits.put(mass, null);
                            break;
                        }
                        finalBandwidth = Math.min(finalBandwidth * 2, 1);
                        System.out.println("Regression failed, retrying with double the bandwidth: " + finalBandwidth);
                    }
                }
            }
        } else if (mode.equals("IM")) {
            double[][][] expAndPredIMs = IMFunctions.getIMarrays(this, regressionSize);
            for (double[][] bins : expAndPredIMs) {
                if (bins[0] == null) {
                    continue;
                }
                IMLOESS.add(LOESS(bins, Float.parseFloat(bandwidth), robustIters)); //needs more smoothing?
            }
        } else {
            throw new IllegalArgumentException("only choose RT or IM");
        }
    }

    //this assumes min delta RT is the best method, but could also be average across mass calibratioin curves
    public void predictRTLOESS(ExecutorService executorService) throws ExecutionException, InterruptedException {
        //long startTime = System.nanoTime();
        futureList.clear();

        //iterate over this list of scan numbers
        for (int i = 0; i < Constants.numThreads; i++) {
            int start = (int) (scanNumberObjects.size() * (long) i) / Constants.numThreads;
            int end = (int) (scanNumberObjects.size() * (long) (i + 1)) / Constants.numThreads;
            futureList.add(executorService.submit(() -> {
                HashMap<String, Double> LOESSRT = new HashMap<>();
                for (int j = start; j < end; j++) {
                    MzmlScanNumber msn = getScanNumObject(scanNums.get(j));
                    for (String mass : RTLOESS.keySet()) {
                        //if null
                        if (RTLOESS.get(mass) == null) {
                            continue;
                        }
                        LOESSRT.put(mass, RTLOESS.get(mass).invoke((double) msn.RT));
                    }
                    for (PeptideObj pep : msn.peptideObjects) {
                        if (pep == null) {
                            break;
                        }
                        double finalDelta = Double.MAX_VALUE;
                        boolean isNone = true;
                        for (String mass : LOESSRT.keySet()) {
                            String[] masses = mass.split("/");
                            for (String minimass : masses) {
                                if (pep.name.contains(minimass)) {
                                    isNone = false;
                                    double rt = LOESSRT.get(minimass);
                                    double delta = Math.abs(rt - pep.RT);
                                    if (delta < finalDelta) {
                                        finalDelta = delta;
                                        pep.calibratedRT = rt;
                                        pep.predRTrealUnits = RTLOESS_realUnits.get(minimass).invoke((double) pep.RT);
                                        pep.deltaRTLOESS_real = Math.abs(msn.RT - pep.predRTrealUnits);
                                    }
                                }
                            }
                        }
                        if (isNone) {
                            if (LOESSRT.isEmpty()) {
                                finalDelta = 0;
                                pep.deltaRTLOESS_real = 0;
                            } else {
                                pep.calibratedRT = LOESSRT.get("others");
                                finalDelta = Math.abs(pep.calibratedRT - pep.RT);
                                pep.predRTrealUnits = RTLOESS_realUnits.get("others").invoke((double) pep.RT);
                                pep.deltaRTLOESS_real = Math.abs(msn.RT - pep.predRTrealUnits);
                            }
                        }
                        pep.deltaRTLOESS = finalDelta;
                    }
                }
            }));
        }
        for (Future future : futureList) {
            future.get();
        }
    }

    public void predictIMLOESS(ExecutorService executorService) throws ExecutionException, InterruptedException {
        //long startTime = System.nanoTime();
        futureList.clear();

        //iterate over this list of scan numbers
        for (int i = 0; i < Constants.numThreads; i++) {
            int start = (int) (scanNumberObjects.size() * (long) i) / Constants.numThreads;
            int end = (int) (scanNumberObjects.size() * (long) (i + 1)) / Constants.numThreads;
            futureList.add(executorService.submit(() -> {
                for (int j = start; j < end; j++) {
                    MzmlScanNumber msn = getScanNumObject(scanNums.get(j));

                    for (PeptideObj pep : msn.peptideObjects) {
                        if (pep == null) {
                            break;
                        }
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
        //long endTime = System.nanoTime();
        //long duration = (endTime - startTime);
        //System.out.println("Calculating deltaIMLOESS took " + duration / 1000000 +" milliseconds");
    }

    public void clear() {
        scanNumberObjects.clear();
        scanNums.clear();
        IMLOESS.clear();
        peptideIMs.clear();
    }
}
