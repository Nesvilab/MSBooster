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
import umich.ms.datatypes.LCMSDataSubset;
import umich.ms.datatypes.scan.IScan;
import umich.ms.datatypes.scan.StorageStrategy;
import umich.ms.datatypes.scancollection.impl.ScanCollectionDefault;
import umich.ms.fileio.exceptions.FileParsingException;
import umich.ms.fileio.filetypes.mzml.MZMLFile;
//import umontreal.ssj.gof.KernelDensity;
import umontreal.ssj.probdist.EmpiricalDist;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
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
    public Function1<Double, Double> RTLOESS;
    public ArrayList<Float>[][] IMbins = null;
    public float[][][] IMbinStats = new float[IMFunctions.numCharges][2 * Constants.IMbinMultiplier + 1][3];
    private ArrayList<Function1<Double, Double>> IMLOESS = new ArrayList<>();
    public HashMap<String, ArrayList<Float>> peptideIMs = new HashMap<>();
    public double[][] expAndPredRTs;
    private List<Future> futureList = new ArrayList<>(Constants.numThreads);

    public mzMLReader(String filename) throws FileParsingException, ExecutionException, InterruptedException {
        // Creating data source
        //path = Paths.get(filename); //
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
        scanNumberObjects = new HashMap<>(scans.getScanCount());
        createScanNumObjects(); //seems like we always need this anyway

        scanNums = new ArrayList<>(scanNumberObjects.size());
        scanNums.addAll(scanNumberObjects.keySet());
        //this.getMzFreq(); only if we end up using weights
    }

    public mzMLReader(mgfFileReader mgf) throws FileParsingException, ExecutionException, InterruptedException { //uncalibrated mgf from MSFragger .d search
        pathStr = mgf.filenames.get(0);

        //Constants.useIM = true;
        //scanNumberObjects = mgf.scanNumberObjects;
        for (Map.Entry<Integer, mzmlScanNumber> entry : mgf.scanNumberObjects.entrySet()) {
            scanNumberObjects.put(entry.getKey(), entry.getValue());
        }
        mgf.clear();

        scanNums = new ArrayList<>(scanNumberObjects.size());
        scanNums.addAll(scanNumberObjects.keySet());
        //this.getMzFreq(); only if we end up using weights
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

    public void createScanNumObjects() throws FileParsingException, ExecutionException, InterruptedException {
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

                mzmlScanNumber msn = new mzmlScanNumber(scan);
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
    }

    public mzmlScanNumber getScanNumObject(int scanNum) {
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

    public void setPinEntries(pinReader pin, SpectralPredictionMapper spm) throws AssertionError, Exception {
        //TODO: multithread?
        HashMap<String, PredictionEntry> allPreds = spm.getPreds();
        while (pin.next()) {
            try {
                scanNumberObjects.get(pin.getScanNum()).setPeptideObject(pin.getPep(), pin.getRank(), pin.getTD(), pin.getEScore(),
                        allPreds);
            } catch (Exception e) {
                e.printStackTrace();
                System.exit(1);
            }
        }
        pin.close();
        pin.reset();
    }

    //could consider not limiting ourselves to peaks predicted by model,
    //but we want something quick to compute as a sanity check
    public HashMap<String, ArrayList<spectrumComparison>> findMS2replicability(
            HashMap<String, ArrayList<spectrumComparison>> peptidoforms) throws IOException {
        //create hashmap with key:peptidoform and value:arraylist of spectrumComparisons
        //filtered by e value
        for (mzmlScanNumber msn : scanNumberObjects.values()) {
            for (peptideObj pobj : msn.peptideObjects) {
                if (Double.parseDouble(pobj.escore) < Constants.MS2escore) {
                    ArrayList<spectrumComparison> arrayList = new ArrayList<>();
                    if (peptidoforms.containsKey(pobj.name)) {
                        arrayList = peptidoforms.get(pobj.name);
                    }
                    if (pobj.spectralSimObj.matchedIntensities != null) {
                        arrayList.add(pobj.spectralSimObj);
                    } else {
                        arrayList.add(pobj.spectralSimObj.spectrumComparisons.get(1));
                        //TODO expand to multiple similarity comparisons
                    }
                    peptidoforms.put(pobj.name, arrayList);
                }
            }
        }
        return peptidoforms;
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
            EmpiricalDist[] kernelDensities = RTFunctions.generateEmpiricalDist(RTbins);

            //long startTime = System.nanoTime();
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
                IMLOESS.add(LOESS(bins, bandwidth * 2, robustIters)); //needs more smoothing?
            }
        } else {
            throw new IllegalArgumentException("only choose RT or IM");
        }
    }

    public void predictRTLOESS(ExecutorService executorService) throws ExecutionException, InterruptedException {
        //long startTime = System.nanoTime();
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
                        pep.calibratedRT = LOESSRT;
                    }
                }
            }));
        }
        for (Future future : futureList) {
            future.get();
        }
        //long endTime = System.nanoTime();
        //long duration = (endTime - startTime);
        //System.out.println("Calculating deltaRTLOESS took " + duration / 1000000 +" milliseconds");
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
