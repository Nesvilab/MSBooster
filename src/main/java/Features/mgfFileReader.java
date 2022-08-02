package Features;

import umich.ms.fileio.exceptions.FileParsingException;

import java.io.*;
import java.nio.charset.StandardCharsets;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.List;
import java.util.concurrent.*;

public class mgfFileReader implements SpectralPredictionMapper{
    //mgfFileReader can handle both single files and entire directories

    ArrayList<String> filenames = new ArrayList<>();
    private ConcurrentHashMap<String, PredictionEntry> allPreds = new ConcurrentHashMap<>();
    HashMap<String, PredictionEntry> allPredsHashMap = new HashMap<>();
    public ConcurrentHashMap<Integer, mzmlScanNumber> scanNumberObjects = new ConcurrentHashMap<>();
    private List<Future> futureList = new ArrayList<>(Constants.numThreads);

    private String returnString(char endChar, byte[] myData, int startInt) {
        int addInt = 0;
        while (!(myData[startInt + addInt] == endChar)) {
            addInt++;
        }
        byte[] byteArray = Arrays.copyOfRange(myData, startInt, startInt + addInt);
        return new String(byteArray);
    }

    private String returnString(char[] endChar, byte[] myData, int startInt) {
        int addInt = 0;
        while (!(myData[startInt + addInt] == endChar[0]) && !(myData[startInt + addInt] == endChar[1])) {
            addInt++;
        }
        byte[] byteArray = Arrays.copyOfRange(myData, startInt, startInt + addInt);
        return new String(byteArray);
    }

    private int returnAdd(char endChar, byte[] myData, int startInt) {
        int addInt = 0;
        while (!(myData[startInt + addInt] == endChar)) {
            addInt++;
        }
        return addInt;
    }

    //this version for uncalibrated mgf acting as mzml
    public mgfFileReader(String file, boolean createScanNumObjects, ExecutorService executorService)
            throws IOException, FileParsingException, ExecutionException, InterruptedException {
        try {
            //add name
            filenames.add(file);

            //load data
            File myFile = new File(file);
            ArrayList<String> allPaths = new ArrayList<>();

            //potentially split data files if mgf is too big
            boolean delFiles = false;
            if (myFile.length() > Integer.MAX_VALUE) {
                int numSplits = (int) (myFile.length() / Integer.MAX_VALUE) + 2; //+2? to be safe

                delFiles = true;

                BufferedReader br = new BufferedReader(new FileReader(myFile));
                ArrayList<Integer> startLines = new ArrayList<>();
                String line;
                int lineNum = 0;
                while ((line = br.readLine()) != null) {
                    if (line.contains("BEGIN IONS")) {
                        startLines.add(lineNum);
                    }
                    lineNum += 1;
                }
                int[] splitPoints = new int[numSplits + 1];
                for (int i = 0; i < numSplits; i++) {
                    int indexer = i * startLines.size() / numSplits;
                    splitPoints[i] = startLines.get(indexer);
                }
                splitPoints[numSplits] = lineNum;
                br.close();

                //generate smaller sub files
                br = new BufferedReader(new FileReader(myFile));
                lineNum = 0;
                line = "";
                for (int split : Arrays.copyOfRange(splitPoints, 1, splitPoints.length)) {
                    String name = file.substring(0, file.length() - 4) + "tmp" + split + ".mgf";
                    BufferedWriter bw = new BufferedWriter(new FileWriter(name));
                    bw.write(line + "\n");
                    while ((line = br.readLine()) != null && lineNum < split) {
                        bw.write(line + "\n");
                        lineNum += 1;
                    }
                    lineNum += 1;
                    allPaths.add(name);
                    bw.close();
                }
            } else {
                allPaths.add(file);
            }

            for (String filePath : allPaths) {
                myFile = new File(filePath);

                byte[] data = new byte[(int) myFile.length()];
                DataInputStream in = new DataInputStream(new BufferedInputStream(new FileInputStream(myFile), 1 << 24));
                in.read(data);
                in.close();
                if (delFiles) {
                    myFile.delete();
                }

                //removing comments
                //can we assume that comments means it's from timsTOF?
                for (byte b : data) {
                    if (b == 35) {
                        data = (new String(data, StandardCharsets.UTF_8)).replaceAll("#[^\r\n]*[\r\n]", "").trim()
                                .getBytes(StandardCharsets.UTF_8);
                        break;
                    }
                }
                //remove windows carriage return
                for (int i = 0; i < data.length; i++) {
                    if (data[i] == 13) {
                        if (data[i + 1] == 10) {
                            data = (new String(data, StandardCharsets.UTF_8)).replaceAll("\r\n", "\n").trim()
                                    .getBytes(StandardCharsets.UTF_8);
                            break;
                        }
                    }
                }

                //find where specific lines are
                int[] chunks = null;
                int ogNumThreads = Constants.numThreads;
                try {
                    chunks = new int[Constants.numThreads + 1];
                    chunks[0] = 0;
                    chunks[Constants.numThreads] = data.length;
                    int jump = data.length / Constants.numThreads;
                    for (int l = 1; l < Constants.numThreads; l++) {
                        int start = jump * l;
                        while (true) {
                            if (data[start] == '\n') {
                                if (new String(Arrays.copyOfRange(data, start + 1, start + 11)).equals("BEGIN IONS")) {
                                    chunks[l] = start + 12; //skip first BEGIN IONS
                                    break;
                                }
                            }
                            start++;
                        }
                    }
                } catch (Exception e) { //short mgf
                    Constants.numThreads = 1;
                    chunks = new int[Constants.numThreads + 1];
                    chunks[0] = 0;
                    chunks[1] = data.length;
                }

                //parallelize
                for (int i = 0; i < Constants.numThreads; i++) {
                    int finalI = i;
                    byte[] finalData = data;
                    int[] finalChunks = chunks;
                    futureList.add(executorService.submit(() -> {
                        int scanNum = 0;
                        StringBuilder sb = new StringBuilder();
                        float RT = 0;
                        float IM = 0;
                        String charge;
                        ArrayList<Float> intensities = new ArrayList<>(12000);
                        ArrayList<Float> mzs = new ArrayList<>(12000);
                        int start = finalChunks[finalI];
                        int end = finalChunks[finalI + 1];
                        String line = "";

                        while (start < end - 11) {
                            switch (finalData[start]) {
                                case 'T': //TITLE
                                    start += 6;
                                    line = returnString('\n', finalData, start);
                                    sb.append(line).append("|");

                                    //timstof
                                    if (line.contains("Cmpd")) {
                                        String[] dotSplit = line.split(",");
                                        scanNum = Integer.parseInt(dotSplit[0].split(" ")[1]);

                                        //read in 1/K0
                                        for (String s : dotSplit) {
                                            if (s.startsWith("1/K0")) {
                                                IM = Float.parseFloat(s.split("=")[1]);
                                            }
                                        }
                                    } else {
                                        String[] dotSplit = line.split("\\.");
                                        try { //for create scan num obj
                                            scanNum = Integer.parseInt(dotSplit[dotSplit.length - 2]);
                                        } catch (Exception ignored) { }
                                    }
                                    start += line.length() + 1;
                                    break;
                                case 'C': //CHARGE
                                    if (finalData[start + 1] == 'H') {
                                        start += 7;
                                    }
                                    charge = returnString('\n', finalData, start);
                                    charge = charge.replace("+", "");
                                    charge = charge.replace("-", "");
                                    sb.append(charge);
                                    start += returnAdd('\n', finalData, start) + 1;
                                    break;
                                case 'P': //PEPMASS or PEPTIDE
                                    start += 8;
                                    start += returnAdd('\n', finalData, start) + 1;
                                    break;
                                case 'R': //RTINSECONDS
                                    if (finalData[start + 1] == 'T') {
                                        start += 12;
                                        line = returnString('\n', finalData, start);
                                        RT = Float.parseFloat(line);
                                    } else {
                                        line = returnString('\n', finalData, start);
                                    }
                                    start += line.length() + 1;
                                    break;
                                case 'E': // END IONS
                                    start += 9;
                                    //do create scanNumObj
                                    float[] mzArray = new float[mzs.size()];
                                    for (int h = 0; h < mzs.size(); h++) {
                                        mzArray[h] = mzs.get(h);
                                    }
                                    float[] intArray = new float[intensities.size()];
                                    for (int h = 0; h < intensities.size(); h++) {
                                        intArray[h] = intensities.get(h);
                                    }
                                    if (createScanNumObjects) { //act as mzml
                                        scanNumberObjects.put(scanNum, new mzmlScanNumber(scanNum, mzArray, intArray, RT, IM));
                                    } else { //act as predictions
                                        PredictionEntry newPred = new PredictionEntry();
                                        newPred.setMzs(mzArray);
                                        newPred.setIntensities(intArray);
                                        newPred.setRT(RT);
                                        newPred.setIM(IM);
                                        //convert title to base format
                                        String basePep = sb.toString();
                                        if (Constants.spectraRTPredModel.contains("PredFull")) {
                                            String[] basePepSplit = basePep.split("\\|");
                                            basePep = new PeptideFormatter(basePepSplit[0],
                                                    basePepSplit[1], "predfull").baseCharge;
                                        } else if (Constants.spectraRTPredModel.contains("pDeep")) {
                                            basePep = new PeptideFormatter(basePep,
                                                    basePep.split("\\|")[2], "pdeep3").baseCharge;
                                        }
                                        allPreds.put(basePep, newPred);
                                        sb.setLength(0);
                                    }

                                    //reset for next peptide/PSM
                                    mzs.clear();
                                    intensities.clear();
                                    break;
                                case 'B': // BEGIN IONS
                                    start += 11;
                                    break;
                                case '1': // 1/K0
                                    if (finalData[start + 1] == '/') {
                                        start += 5;
                                        line = returnString('\n', finalData, start);
                                        IM = Float.parseFloat(line);
                                    } else {
                                        line = returnString(new char[]{' ', '\t'}, finalData, start);
                                        mzs.add(Float.parseFloat(line));
                                        start += line.length() + 1;

                                        line = returnString('\n', finalData, start);
                                        line = line.split(" ")[0]; //for pdeep2 predictions, included space and fragment ion type
                                        intensities.add(Float.parseFloat(line.split("\t")[0]));
                                    }
                                    start += line.length() + 1;
                                    break;
                                case '2': //if number, skip switch
                                case '3':
                                case '4':
                                case '5':
                                case '6':
                                case '7':
                                case '8':
                                case '9':
                                case '0':
                                    line = returnString(new char[]{' ', '\t'}, finalData, start);
                                    mzs.add(Float.parseFloat(line));
                                    start += line.length() + 1;

                                    line = returnString('\n', finalData, start);
                                    line = line.split(" ")[0]; //for pdeep2 predictions, included space and fragment ion type
                                    intensities.add(Float.parseFloat(line.split("\t")[0]));
                                    start += line.length() + 1;
                                    break;
                                default: //USERNAME
                                    line = returnString('\n', finalData, start);
                                    start += line.length() + 1;
                                    break;
                            }
                        }

                        //last try to see if any remaining entries to add
                        float[] mzArray = new float[mzs.size()];
                        for (int h = 0; h < mzs.size(); h++) {
                            mzArray[h] = mzs.get(h);
                        }
                        float[] intArray = new float[intensities.size()];
                        for (int h = 0; h < intensities.size(); h++) {
                            intArray[h] = intensities.get(h);
                        }
                        if (createScanNumObjects) { //act as mzml
                            if (! scanNumberObjects.containsKey(scanNum)) {
                                scanNumberObjects.put(scanNum, new mzmlScanNumber(scanNum, mzArray, intArray, RT, IM));
                            }
                        } else { //act as predictions
                            PredictionEntry newPred = new PredictionEntry();
                            newPred.setMzs(mzArray);
                            newPred.setIntensities(intArray);
                            newPred.setRT(RT);
                            newPred.setIM(IM);
                            if (mzArray.length != 0) {
                                //convert title to base format
                                String basePep = sb.toString();
                                if (Constants.spectraRTPredModel.contains("PredFull")) {
                                    String[] basePepSplit = basePep.split("\\|");
                                    basePep = new PeptideFormatter(basePepSplit[0],
                                            basePepSplit[1], "predfull").baseCharge;
                                } else if (Constants.spectraRTPredModel.contains("pDeep")) {
                                    basePep = new PeptideFormatter(basePep,
                                            basePep.split("\\|")[2], "pdeep3").baseCharge;
                                }
                                allPreds.put(basePep, newPred);
                            }
                        }

                        //reset for next peptide/PSM
                        mzs.clear();
                        intensities.clear();
                    }));
                }
                for (Future future : futureList) {
                    future.get();
                }
                Constants.numThreads = ogNumThreads;
            }
        } catch (Exception e) {
            e.printStackTrace();
            System.exit(1);
        }
    }

    public HashMap<String, PredictionEntry> getPreds() {
        if (allPredsHashMap.size() == 0) {
            allPredsHashMap = new HashMap<>(allPreds);
            allPreds = null; //no longer need concurrency
        }
        return allPredsHashMap;
    }

    public float getMaxPredRT() {
        float maxRT = 0f;
        for (PredictionEntry entry : allPreds.values()) {
            if (entry.RT > maxRT) {
                maxRT = entry.RT;
            }
        }
        return maxRT;
    }

    public void clear() {
        allPreds.clear();
        scanNumberObjects.clear();
        futureList.clear();
    }

    public static void main(String[] args) throws IOException, FileParsingException, ExecutionException, InterruptedException {
        Constants.numThreads = 11;
        ExecutorService executorService = Executors.newFixedThreadPool(Constants.numThreads);
        long startTime = System.nanoTime();
        mgfFileReader mgf = new mgfFileReader("C:/Users/kevin/Downloads/proteomics/File_2.mgf", true,
                executorService);
        long endTime = System.nanoTime();
        long duration = (endTime - startTime);
        System.out.println("loading took " + duration / 1000000 +" milliseconds");
        executorService.shutdown();
    }
}
