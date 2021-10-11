package Features;

import umich.ms.fileio.exceptions.FileParsingException;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.lang.reflect.InvocationTargetException;
import java.util.concurrent.ExecutionException;

public class RandomThings {

    public RandomThings(){

    }

    public static void main(String[] args) throws IOException, FileParsingException, NoSuchMethodException, IllegalAccessException, InvocationTargetException, InterruptedException, ExecutionException {
        BufferedReader br = new BufferedReader(new FileReader(new File("C:/Users/yangkl/Downloads/proteomics/" +
                "timstof/20180819_TIMS2_12-2_AnBr_SA_200ng_HeLa_50cm_120min_100ms_11CT_3_A1_01_2769_quant.csv")));
        System.out.println(br.readLine());
        System.out.println(br.readLine());
        System.out.println(br.readLine().split(",")[11].equals(""));
        br.close();
        //        ExecutorService executorService = Executors.newFixedThreadPool(Constants.numThreads);
//        List<Future> futureList = new ArrayList<>(Constants.numThreads);
//
//        for (int i = 0; i < Constants.numThreads; i++) {
//            int start = (int) (1000000 * (long) i) / Constants.numThreads;
//            int end = (int) (1000000 * (long) (i + 1)) / Constants.numThreads;
//            futureList.add(executorService.submit(() -> {
//                for (int j = start; j < end; j++) {
//                    int[] arr = new int[1];
//                    arr[0] = 1;
//                }
//            }));
//        }
//
//        long startTime = System.nanoTime();
//        for (Future future : futureList) {
//            future.get();
//        }
//
//        long endTime = System.nanoTime();
//        long duration = (endTime - startTime);
//        System.out.println("process took " + duration / 1000000 +" milliseconds");
//        executorService.shutdown();

    }

    //this is for creating tsv for creating histograms for scores between target and decoy in python
//    public static void main(String[] args) throws FileParsingException, IOException {
//        String[] allSpectralAngles = new String[0];
//
//        //read in specific files
//        mgfFileReader predictedSpectra = new mgfFileReader("C:/Users/kevin/Downloads/proteomics/wideWindowPreds.mgf");
//        HashMap<String, double[]> predictedMZs = predictedSpectra.getMzDict();
//        HashMap<String, double[]> predictedIntensities = predictedSpectra.getMzDict();
//
//        String mzmlPath = "C:/Users/kevin/Downloads/proteomics/wideWindow/23aug2017_hela_serum_timecourse_pool_wide_002.mzML";
//        //load experimental spectra from mzML
//        mzMLReader mzMLscans = new mzMLReader(mzmlPath);
//
//        //get mzFreq
//        double[] mzFreqs = mzMLscans.getMzFreq();
//
//        File pepXMLDirectory = new File("C:/Users/kevin/Downloads/proteomics/wideWindow/");
//        String[] pepXMLfiles = pepXMLDirectory.list((d, name) -> name.endsWith(".pepXML")
//                && name.startsWith("23aug2017_hela_serum_timecourse_pool_wide_002")); //not recursive
//
//        for (String xmlPath : pepXMLfiles) {
//            System.out.println(xmlPath);
//            pepXMLReader xmlReader = new pepXMLReader("C:/Users/kevin/Downloads/proteomics/wideWindow/" + xmlPath);
//
//            //get all scan numbers from pepXML
//            int[] scanNums = xmlReader.getScanNumbers();
//
//            //get peptide names compatible with predicted HashMaps
//            String[] xmlPeptides = xmlReader.getXMLpeptides();
//
//            //get target or decoys
//            int[] td = xmlReader.getTargetOrDecoy();
//
//            //get expectation scores
//            String[] eScores = xmlReader.getEScore();
//
//            //iterate through all pepXMLhits
//            int pepXMLhits = scanNums.length;
//
//            //save spectral angles
//            String[] sa = new String[pepXMLhits];
//
//            for (int i = 0; i < pepXMLhits; i++) {
//                //get scanNumber and peptide of current hit
//                int num = scanNums[i];
//                String pep = xmlPeptides[i];
//
//                //get experimental spectra into
//                double[] expMZs = mzMLscans.getMZ(num);
//                double[] expIntensities = mzMLscans.getIntensity(num);
//
//                //get predicted spectra info
//                double[] predMZs = predictedMZs.get(pep);
//                double[] predIntensities = predictedIntensities.get(pep);
//
//                //calculate spectral angle
//                try {
//                    spectrumComparison specAngle = new spectrumComparison(expMZs, expIntensities,
//                            predMZs, predIntensities);
//
//                    //adapt for weighted or multiple sims, refer to pepXMLModifier
//                    HashMap<String, Double> sims = specAngle.getAllSimilarities(mzFreqs);
//
//                    //get expectation score and target or decoy
//                    String eScore = eScores[i];
//                    int tdSingle = td[i];
//
//                    sa[i] = pep + "\t" + sims.get("cosine") + "\t" + sims.get("weightCosine") + "\t" +
//                            sims.get("contrast") + "\t" + sims.get("weightContrast") + "\t" +
//                            sims.get("euclidean") + "\t" + sims.get("weightEuclidean") + "\t" +
//                            sims.get("bray-curtis") + "\t" + sims.get("weightBray-curtis") + "\t" +
//                            sims.get("pearson") + "\t" + sims.get("weightPearson") + "\t" +
//                            sims.get("dot") + "\t" + sims.get("weightDot") + "\t" +
//                            eScore + "\t" + tdSingle;
//
//                } catch (Exception e) {
//                    System.out.println(pep); //probably was not supported by pDeep2 prediction (ex. amino acid U)
//                }
//
//            }
//            //add to full list of spectral angles of peptides
//            allSpectralAngles = (String[]) ArrayUtils.addAll(allSpectralAngles, sa);
//        }
//
//
//        ////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//        //will need to remove null
//        //send to file for analysis/data viz in python/R
//        try {
//            FileWriter myWriter = new FileWriter("C:/Users/kevin/Downloads/proteomics/wideWindow/sims.tsv");
//            //need to adapt for multiple similarities
//            myWriter.write("peptide" + "\t" + "cosineSimilarity" + "\t" + "weightedCosineSimilarity" + "\t" +
//                    "spectralContrastAngle" + "\t" + "weightedSpectralContrastAngle" + "\t" +
//                    "euclideanDistance" + "\t" + "weightedEuclideanDistance" + "\t" +
//                    "brayCurtis" + "\t" + "weightedBrayCurtis" + "\t" +
//                    "pearsonCorr" + "\t" + "weightedPearsonCorr" + "\t" +
//                    "dotProduct" + "\t" + "weightedDotProduct" + "\t" +
//                    "eScore" + "\t" + "target=1/decoy=0\n");
//
//            for (String s : allSpectralAngles) {
//                if (s != null) {
//                    myWriter.write(s + "\n");
//                } else {
//                    System.out.println("found empty line");
//                }
//            }
//
//            myWriter.close();
//            System.out.println("Successfully wrote to the file.");
//        } catch (IOException e) {
//            System.out.println("An error occurred.");
//            e.printStackTrace();
//        }
//    }

//    public static void RTbins(mzMLReader mzml, mgfFileReader preds) throws IOException {
//        //get max index
//        int maxKey = Collections.max(mzml.scanNumberObjects.keySet());
//        int numBins = (int) mzml.scanNumberObjects.get(maxKey).RT + 1;
//
//        //set up arrays for predRT, assuming minute long bins
//        ArrayList<Float>[][] predRTfloor = new ArrayList[2][numBins];
//        for (int row = 0; row < 2; row++) {
//            for (int col = 0; col < numBins; col++) {
//                predRTfloor[row][col] = new ArrayList<Float>();
//            }
//        }
//
//        ArrayList<Float>[][] predRTround = new ArrayList[2][numBins + 1];
//        for (int row = 0; row < 2; row++) {
//            for (int col = 0; col < numBins + 1; col++) {
//                predRTround[row][col] = new ArrayList<Float>();
//            }
//        }
//
//        //write experimental RT for python plotting 2d kde
//        FileWriter[] files = new FileWriter[2];
//        files[0] = new FileWriter("C:/Users/kevin/Downloads/proteomics/scanNumberRTdecoy.txt");
//        files[1] = new FileWriter("C:/Users/kevin/Downloads/proteomics/scanNumberRTtarget.txt");
//
//        //iterate through scanNumbers
//        for (int scanNum : mzml.scanNumberObjects.keySet()) {
//            mzmlScanNumber scanNumObj = mzml.getScanNumObject(scanNum);
//            float rt = scanNumObj.RT; //experimental RT for this scan
//            int floor = (int) rt;
//            int round = (int) Math.round(rt);
//
//            //iterate through PSMs
//            for (int i = 1; i < scanNumObj.peptideObjects.size() + 1; i++) {
//                peptideObj pep = scanNumObj.getPeptideObject(i);
//                //identify if target or decoy
//                int td = pep.targetORdecoy;
//
//                //testing out adding more info for decoys
//                if (td == 0) {
//                    predRTfloor[0][floor].add(pep.RT);
//                    predRTround[0][round].add(pep.RT);
//                    files[0].write(rt + "\n");
//                    continue;
//                }
//                //end this section
//
//                //decide instances to add based on eScore (subject to change)
//                int instances = -1 * (int) Math.ceil(Math.log10(Double.parseDouble(pep.escore)));
//                if (instances > 0) { //add directly to distribution now
//                    for (int j = 0; j < instances; j++) {
//                        predRTfloor[td][floor].add(pep.RT);
//                        predRTround[td][round].add(pep.RT);
//                        files[td].write(rt + "\n");
//                    }
//                }
//            }
//        }
//        files[0].close();
//        files[1].close();
//
//        //now we have our distributions
//        FileWriter w = new FileWriter("C:/Users/kevin/Downloads/proteomics/empDistDecoy.txt");
//        for (ArrayList a : predRTfloor[0]) {
//            for (Object d : a) {
//                w.write(d + "\t");
//            }
//            w.write("\n");
//        }
//        w.close();
//
//
//    }

//for creating diannFrags and pDeep3Frags for comparing number fragments nad similarity predictions (predictionCompare.ipynb)
//SpectralPredictionMapper spm = SpectralPredictionMapper.createSpectralPredictionMapper
//        ("C:/Users/kevin/OneDriveUmich/proteomics/preds/cptacDiann.predicted.bin");
//    //        SpectralPredictionMapper spm = new mgfFileReader("C:/Users/kevin/OneDriveUmich/proteomics/preds/cptacPreds.mgf");
//    FileWriter myWriter = new FileWriter("C:/Users/kevin/Downloads/proteomics/cptac/diannFrags.tsv");
//    ////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
//    mzMLReader mzml = new mzMLReader("C:/Users/kevin/OneDriveUmich/proteomics/mzml/cptac/" +
//            "CPTAC_CCRCC_W_JHU_LUMOS_C3L-01665_T.mzML");
//
//        for (int i = 1; i < 17; i++) {
//        pepXMLReader pepxml = new pepXMLReader("C:/Users/kevin/Downloads/proteomics/cptac/2021-2-21/" +
//                "CPTAC_CCRCC_W_JHU_LUMOS_C3L-01665_T_rank" + i + ".pepXML");
//        mzml.setPepxmlEntries(pepxml, i, spm);
//    }
//
//        for (int i : mzml.scanNumberObjects.keySet()) {
//        mzmlScanNumber msn = mzml.getScanNumObject(i);
//
//        for (peptideObj pObj : msn.peptideObjects) {
//            int numPred = pObj.spectralSimObj.predMZs.length;
//            int numMatch = 0;
//            for (float f : pObj.spectralSimObj.matchedIntensities) {
//                if (f != 0f) {
//                    numMatch += 1;
//                }
//            }
//            myWriter.write(numPred + "\t" + numMatch + "\t" + pObj.targetORdecoy + "\t" +
//                    pObj.spectralSimObj.brayCurtis() + "\n");
//        }
//    }
//        myWriter.close();
//    public static void peptideRTForPython() throws FileParsingException, IOException {
//        String prefix = "23aug2017_hela_serum_timecourse_4mz_narrow_1";
//        String outfile = "C:/Users/kevin/Downloads/proteomics/narrowWindow/RT2-21.csv";
//
//        mzMLReader mzml = new mzMLReader("C:/Users/kevin/OneDriveUmich/proteomics/mzml/" +
//                "narrowWindow/" + prefix + ".mzML");
//
//        //mgfFileReader mgf = new mgfFileReader("C:/Users/kevin/OneDriveUmich/proteomics/preds/narrowPDeepPreds.mgf");
//        SpectralPredictionMapper spm = SpectralPredictionMapper.createSpectralPredictionMapper("C:/Users/kevin/OneDriveUmich/proteomics/preds/narrowPDeepPreds.mgf");
//        //end modify
//
//        for (int rank = 1; rank < 5; rank++) {
//            System.out.println(rank);
//            pepXMLReader xmlReader = new pepXMLReader("C:/Users/kevin/Downloads/proteomics/narrowWindow/" +
//                    "rank" + rank + "/" + prefix + ".pepXML");
//
//            mzml.setPepxmlEntries(xmlReader, rank, spm);
//        }
//
//        //write csv file, x and y columns for real RT and predicted
//        FileWriter myWriter = new FileWriter(outfile);
//        for (mzmlScanNumber s : mzml.scanNumberObjects.values()) {
//            if (s.peptideObjects.size() > 0) {
//                double RT = s.RT;
//                for (peptideObj p : s.peptideObjects) {
//    //                    if (Double.parseDouble(p.escore) > 0.000001) {
//    //                        break;
//    //                    }
//                    myWriter.write(RT + "," + p.RT + "," + p.escore + "\n");
//                }
//            }
//        }
//        myWriter.close();
//    }
//get list of peptides we did and did not identify
//for each row in DAINN.tsv, check which list it is in
//    CsvParserSettings settings = new CsvParserSettings();
//        RowListProcessor rowProcessor = new RowListProcessor();
//            settings.getFormat().setLineSeparator("\n");
//            settings.setHeaderExtractionEnabled(true);
//            settings.setProcessor(rowProcessor);
//        CsvParser parser = new CsvParser(settings);
//            parser.parse(new File("C:/Users/kevin/Downloads/proteomics/cptac/2021-2-21/DIANNGroupsTop128.csv"));
//        String[] headers = rowProcessor.getHeaders();
//        List<String[]> allRows = rowProcessor.getRows();
//        int idIdx = ArrayUtils.indexOf(headers, "Precursor.Id");
//        int idScanNum = ArrayUtils.indexOf(headers, "MS2.Scan");
//        int group = ArrayUtils.indexOf(headers, "list");
//
//        //load mzmlreader
//        mzMLReader mzml = new mzMLReader("C:/Users/kevin/OneDriveUmich/proteomics/mzml/cptac/" +
//                "CPTAC_CCRCC_W_JHU_LUMOS_C3L-01665_T.mzML");
//
//        //writer
//        FileWriter myWriter = new FileWriter("C:/Users/kevin/Downloads/proteomics/cptac/2021-2-21/averageIntensities128.csv");
//
//            try {
//            for (String[] row : allRows) {
//                //format peptide: for precursor.id, replace () with [], and get charge from last character
//                String peptide = row[idIdx].replaceAll("\\(", "[").replaceAll("\\)", "]");
//
//                //get fragments MZs for the peptide of the row
//                MassCalculator mc = new MassCalculator(peptide.substring(0, peptide.length() - 1),
//                        peptide.substring(peptide.length() - 1));
//                float[] predMZs = mc.calcAllMasses();
//                float[] predInts = new float[predMZs.length];
//                for (int i = 0; i < predMZs.length; i++) {
//                    predInts[i] = 100f;
//                }
//
//                //get experimental mz list from mzmlreader
//                int scanNum = Integer.parseInt(row[idScanNum]);
//    //                float[] expMZs = mzml.getScanNumObject(scanNum).getExpMZs();
//    //                float[] expInts = mzml.getScanNumObject(scanNum).getExpIntensities();
//                float[] expMZs = mzml.getMZ(scanNum);
//                float[] expInts = mzml.getIntensity(scanNum);
//
//                //generate spectrumComparison object, with predicted fragments with intensities of 100 each
//                spectrumComparison sc = new spectrumComparison(expMZs, expInts, predMZs, predInts);
//
//                //see what fraction of fragments were matched and average intensity of matched fragments
//                int matched = 0;
//                float matchedInts = 0f;
//                ArrayList<Float> intArray = new ArrayList<Float>();
//                for (float f : sc.matchedIntensities) {
//                    if (f != 0f) {
//                        matched += 1;
//                        intArray.add(f);
//                    }
//                    matchedInts += f;
//                }
//                float avgInt = matchedInts / (float) matched;
//                float fraction = (float) matched / (float) predMZs.length;
//                float maxInt = 0f;
//                float median = 0f;
//                if (Float.isNaN(avgInt)) {
//                    avgInt = 0f;
//                } else {
//                    Collections.sort(intArray);
//                    int arraySize = intArray.size();
//                    maxInt = intArray.get(arraySize - 1);
//                    if (intArray.size() % 2 == 0)
//                        median = (intArray.get((arraySize / 2) - 1) + intArray.get((arraySize / 2))) / 2;
//                    else
//                        median = intArray.get((arraySize - 1) / 2);
//                }
//                //write to file
//                myWriter.write(fraction + "," + avgInt + "," + row[group] + "," + mc.charge + "," +
//                        median + "," + maxInt + "\n");
//            }
//            myWriter.close();
//            System.out.println("done");
//        } catch (Exception e) {
//            System.out.println(e);
//            myWriter.close();
//        }
    //    private static void collisions() throws IOException, FileParsingException{
//        //constants
//        //int window = 3;
//        int maxRank = 10;
//        float fragTol = 10f; //ppm
//
//        //read in peptideIDs.tsv and save to hashmap
//        //also define counter here
//        int[] counter = new int[3];
//        HashSet[] pepSet = {new HashSet(), new HashSet(), new HashSet()};
//        HashMap<String, Integer> pepIDs = new HashMap<String, Integer>();
//        try (BufferedReader br = new BufferedReader(new FileReader("peptideIDs.tsv"))) {
//            String line;
//            while ((line = br.readLine()) != null) {
//                String[] info = line.split("\t", -1);
//                int id = Integer.parseInt(info[1]);
//                pepIDs.put(info[0], id);
//                counter[id] += 1; //add to counter
//            }
//        }
//        System.out.println(Arrays.toString(counter));
//
//        //read in predictions
//        System.out.println("loading predictions");
//        mgfFileReader mgf = new mgfFileReader("C:/Users/kevin/Downloads/proteomics/wideWindowPredsAllRanks.mgf");
//        HashMap<String, float[]> predMZs = mgf.getMzDict();
//
//        //keep track of how many times there is fragment collision
//        double collisions = 0.0;
//        double total = 0.0;
//        double lostPeps = 0.0; //peptides lost from interact.pepxml
//
//        //write to file
//        //FileWriter myWriter = new FileWriter("collisions.txt");
//
//        //get mzmlReaders ready
//        for (int window = 1; window < 4; window++) {
//            System.out.println("loading mzml");
//            mzMLReader m = new mzMLReader("C:/Users/kevin/OneDriveUmich/proteomics/mzml/wideWindow/" +
//                    "23aug2017_hela_serum_timecourse_pool_wide_00" + window + ".mzML");
//
//            //fill with pepxml files
//            System.out.println("loading pepxml");
//            for (int rank = 1; rank < maxRank + 1; rank++) {
//                pepXMLReader p = new pepXMLReader("C:/Users/kevin/Downloads/proteomics/wideWindow/encyclopedia_params/rank"
//                        + rank + "/23aug2017_hela_serum_timecourse_pool_wide_00" + window + ".pepXML");
//                m.setPepxmlEntries(p, rank, mgf);
//            }
//
//
//            System.out.println("starting collision counting");
//
//            //for each scanNum, get predicted MZs for each peptide and check for overlapping peaks
//            for (mzmlScanNumber msn : m.scanNumberObjects.values()) {
//                String[] names = msn.getNames(true);
//
//                if (names.length < 2) { //ms2 scans with no or only 1 peptide IDs
//                    continue;
//                }
//                total += 1.0;
//                //arraylist of mz values
//                ArrayList<Float> mzs = new ArrayList<>();
//                for (String s : names) {
//                    try {
//                        float[] allMZs = mgf.allPredMZs.get(s);
//
//                        //remove zero intensity fragments
//                        float[] ints = msn.getPeptideObject(s).spectralSimObj.matchedIntensities;
//                        for (int k = ints.length - 1; k > -1; k--){
//                            if (ints[k] == 0f){
//                                allMZs = ArrayUtils.remove(allMZs, k);
//                            }
//                        }
//
//                        Arrays.sort(allMZs);
//                        //remove fragments that are too close to each other
//                        //add largest fragment and others that qualify
//                        mzs.add(allMZs[allMZs.length - 1]);
//                        for (int i = allMZs.length - 2; i > -1; i--) {
//                            double lowerBound = allMZs[i + 1] * (1000000 - fragTol) / 1000000;
//                            if (allMZs[i] < lowerBound) {
//                                mzs.add(allMZs[i]);
//                            }
//                        }
//                    } catch (Exception ignored) {
//                    }
//                }
//
//                //continue if no matched fragment mzs
//                if (mzs.size() == 0){
//                    continue;
//                }
//
//                //sort
//                Collections.sort(mzs);
//                Collections.reverse(mzs);
//
//                //calculate lower bound of each fragment's error
//                double[] lowerBounds = new double[mzs.size() - 1];
//                int i = 0;
//                boolean notCollided = true; //change when no need to add to collide anymore
//                ArrayList<Float> collideFrags = new ArrayList<>(); //keep track of colliding fragments
//                for (float d : mzs.subList(0, lowerBounds.length)) {
//                    float lowerBound = d * (1000000f - fragTol) / 1000000f;
//                    i++;
//
//                    //compare to see if adjacent fragment is greater than lower bound
//                    //note: assumes that this doesn't occur within a peptide (peptide has fragments close to each other)
//                    float otherFrag = mzs.get(i);
//                    if (lowerBound <= otherFrag) {
//                        if (notCollided) {
//                            collisions += 1.0;
//                            notCollided = false;
//                        }
//
//                        //add fragments to study later
//                        collideFrags.add(otherFrag);
//                        //collideFrags.add(mzs.get(i - 1));
//                        collideFrags.add(d);
//
//                        //break;
//                    }
//                }
//
//                //now need to find peptide rank entries that were involved in collision
//                //names will be helpful here to focus on just targets
//                ArrayList<Integer> results = new ArrayList<>(); //holds in pairs of (rank, id)
//                for (String pep : names) {
//                    float[] thisMZs = predMZs.get(pep);
//                    if (thisMZs == null){ //probably has U
//                        continue;
//                    }
//                    //check if anything in intersection
//                    Set<Float> collideFragsSet = new HashSet<>(collideFrags);
//                    Set<Float> thisMZsSet = new HashSet<>(Arrays.asList(ArrayUtils.toObject(thisMZs)));
//                    collideFragsSet.retainAll(thisMZsSet);
//                    if (collideFragsSet.size() != 0) {
//                        //if so, collect rank and id
//                        int collideRank = msn.getPeptideObject(pep).rank;
//                        int pID = pepIDs.get(pep.split("\\|", 2)[0]);
//
//                        //add to results (this is for writing file
//                        results.add(collideRank);
//                        results.add(pID);
//
//                        //percent unique appearances, to see what percent of each pID type is involved in at least 1 collision
//                        pepSet[pID].add(pep);
//                        if (pID == 1) { //lost peptide
//                            //System.out.println(pep);
//                            lostPeps += 1;
//                        }
//
//                        //what happens when we change intensity of shared fragment?
//                        //we have the shared fragments
//                        //get predicted intensities for peptide
//                        //we have exp intensities we can manipulate and make new spectrumComparison obj
//                        //msn.getPeptideObject(pep).spectralSimObj.matchedIntensities;
//                        //with bray curtis, expect experimental to be bigger than predicted (unit normalized)
//                        //can also try with cosine sim
//                    }
//                }
//
//                //write in columns
////                if (results.size() > 0) {
////                    for (int number : results) {
////                        myWriter.write(number + "\t");
////                    }
////                    myWriter.write("\n");
////                }
//            }
//        }
//        System.out.println(collisions / total);
//        System.out.println(lostPeps / total); //different normalization? Only count when lost peptide is present?
//        for (int p = 0; p < 3; p++){
//            System.out.println((double) pepSet[p].size() / (double) counter[p]);
//        }
//        //myWriter.close();
//    }
//
//    private static void chiSquare() throws IOException, FileParsingException {
//        //constants
//        //int window = 3;
//        int maxRank = 10;
//        float fragTol = 10f; //ppm
//        int increased = 0;
//        int decreased = 0;
//
//        //read in peptideIDs.tsv and save to hashmap
//        //also define counter here
//        //int[] counter = new int[3];
//        //HashSet[] pepSet = {new HashSet(), new HashSet(), new HashSet()};
//        HashMap<String, Integer> pepIDs = new HashMap<String, Integer>();
//        try (BufferedReader br = new BufferedReader(new FileReader("peptideIDs.tsv"))) {
//            String line;
//            while ((line = br.readLine()) != null) {
//                String[] info = line.split("\t", -1);
//                int id = Integer.parseInt(info[1]);
//                pepIDs.put(info[0], id);
//            }
//        }
//
//        //read in predictions
//        System.out.println("loading predictions");
//        mgfFileReader mgf = new mgfFileReader("C:/Users/kevin/Downloads/proteomics/wideWindowPredsAllRanks.mgf");
//        HashMap<String, float[]> predMZs = mgf.getMzDict();
//        HashMap<String, float[]> predInts = mgf.getIntensityDict();
//
//        //get mzmlReaders ready
//        for (int window = 1; window < 4; window++) {
//            System.out.println("loading mzml");
//            mzMLReader m = new mzMLReader("C:/Users/kevin/OneDriveUmich/proteomics/mzml/wideWindow/" +
//                    "23aug2017_hela_serum_timecourse_pool_wide_00" + window + ".mzML");
//
//            //fill with pepxml files
//            System.out.println("loading pepxml");
//            for (int rank = 1; rank < maxRank + 1; rank++) {
//                pepXMLReader p = new pepXMLReader("C:/Users/kevin/Downloads/proteomics/wideWindow/encyclopedia_params/rank"
//                        + rank + "/23aug2017_hela_serum_timecourse_pool_wide_00" + window + ".pepXML");
//                m.setPepxmlEntries(p, rank, mgf);
//            }
//
//
//            System.out.println("starting collision counting");
//
//            //for each scanNum, get predicted MZs for each peptide and check for overlapping peaks
//            for (mzmlScanNumber msn : m.scanNumberObjects.values()) {
//
//                //names that are interact peptides
//                String[] names = msn.getNames(true);
//                for (int i = names.length - 1; i > -1; i--){
//                    if (pepIDs.get(names[i].split("\\|")[0]) == 1){
//                        names = (String[]) ArrayUtils.remove(names, i);
//                    }
//                }
//
//                //names that are decoy
////                String[] names = msn.getNames(false);
//
//
//                if (names.length < 2) { //ms2 scans with no or only 1 peptide IDs
//                    continue;
//                }
//
//                //arraylist of mz values
//                ArrayList<Float> mzs = new ArrayList<>();
//                for (String s : names) {
//                    try {
//                        float[] allMZs = mgf.allPredMZs.get(s);
//                        Arrays.sort(allMZs);
//                        //remove fragments that are too close to each other
//                        //add largest fragment and others that qualify
//                        mzs.add(allMZs[allMZs.length - 1]);
//                        for (int i = allMZs.length - 2; i > -1; i--) {
//                            double lowerBound = allMZs[i + 1] * (1000000 - fragTol) / 1000000;
//                            if (allMZs[i] < lowerBound) {
//                                mzs.add(allMZs[i]);
//                            }
//                        }
//                    } catch (Exception ignored) {
//                    }
//                }
//
//                //sort
//                Collections.sort(mzs);
//                Collections.reverse(mzs);
//
//                //calculate lower bound of each fragment's error
//                double[] lowerBounds = new double[mzs.size() - 1];
//                int i = 0;
//                ArrayList<Float> collideFrags = new ArrayList<>(); //keep track of colliding fragments
//                for (float d : mzs.subList(0, lowerBounds.length)) {
//                    float lowerBound = d * (1000000f - fragTol) / 1000000f;
//                    i++;
//
//                    //compare to see if adjacent fragment is greater than lower bound
//                    //note: assumes that this doesn't occur within a peptide (peptide has fragments close to each other)
//                    float otherFrag = mzs.get(i);
//                    if (lowerBound <= otherFrag) {
//                        //add fragments to study later
//                        collideFrags.add(otherFrag);
//                        //collideFrags.add(mzs.get(i - 1));
//                        collideFrags.add(d);
//
//                        //break;
//                        System.out.println(Arrays.toString(names));
//                    }
//                }
//
//                //now need to find peptide rank entries that were involved in collision
//                //names will be helpful here to focus on just targets
//                for (String pep : names) {
//
//                    //for decoys vs all
////                    if (msn.getPeptideObject(pep).targetORdecoy == 1){
////                        continue;
////                    }
//
//                    float[] thisMZs = predMZs.get(pep);
//                    if (thisMZs == null){ //probably has U
//                        continue;
//                    }
//                    //check if anything in intersection
//                    Set<Float> collideFragsSet = new HashSet<>(collideFrags);
//                    Set<Float> thisMZsSet = new HashSet<>(Arrays.asList(ArrayUtils.toObject(thisMZs)));
//                    collideFragsSet.retainAll(thisMZsSet);
//                    if (collideFragsSet.size() != 0) {
//                        //what happens when we change intensity of shared fragment?
//                        //we have the shared fragments
//                        //get predicted intensities for peptide
//                        //we have exp intensities we can manipulate and make new spectrumComparison obj
//                        //msn.getPeptideObject(pep).spectralSimObj.matchedIntensities;
//                        //with bray curtis, expect experimental to be bigger than predicted (unit normalized)
//                        //can also try with cosine sim
//
//                        spectrumComparison og = new spectrumComparison(msn.getExpMZs(), msn.getExpIntensities(),
//                                predMZs.get(pep), predInts.get(pep));
//                        spectrumComparison notOg = new spectrumComparison(msn.getExpMZs(), msn.getExpIntensities(),
//                                predMZs.get(pep), predInts.get(pep));
//
//                        for (int j = 0; j < thisMZs.length; j++) { //arbitrarily decrease int for all colliding frags
//                            if (collideFragsSet.contains(thisMZs[j])) {
//                                if (notOg.matchedIntensities[j] > 0) {
//                                    notOg.matchedIntensities[j] -= 1;
//                                }
//                            }
//                        }
//
//                        if (og.brayCurtis() > notOg.brayCurtis()){
//                            decreased += 1;
//                        } else if (og.brayCurtis() < notOg.brayCurtis()) {
//                            increased += 1;
//                        }
//
//                    }
//                }
//            }
//        }
//        System.out.println(increased);
//        System.out.println(decreased);
//    }
}

