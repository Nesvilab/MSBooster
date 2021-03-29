package Features;

import org.expasy.mzjava.core.io.ms.spectrum.MgfReader;
import org.expasy.mzjava.core.ms.peaklist.PeakList;
import org.expasy.mzjava.core.ms.spectrum.MsnSpectrum;
import umich.ms.fileio.exceptions.FileParsingException;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;

import static Features.floatUtils.doubleToFloat;

public class mgfFileReader implements SpectralPredictionMapper{
    //mgfFileReader can handle both single files and entire directories

    final ArrayList<String> filenames;
    private HashMap<String, float[]> allPredMZs = new HashMap<>();
    private HashMap<String, float[]> allPredIntensities = new HashMap<>();
    private HashMap<String, Float> allPredRTs = new HashMap<>();

    public mgfFileReader(String files) throws IOException {
        File predsDirectory = new File(files);
        String[] predsFiles = predsDirectory.list();
        filenames = new ArrayList<String>();

        if (predsFiles == null) { //if user provided a file, not a directory
            filenames.add(files);
        } else { //user provides directory
            for (String predsFile : predsFiles) {
                if (predsFile.contains(".mgf")) {
                    filenames.add(files + File.separator + predsFile);
                }
            }
        }

        for (String fname : filenames) {
            MgfReader reader = new MgfReader(new File(fname), PeakList.Precision.FLOAT);
            reader.acceptUnsortedSpectra();

            // hasNext() returns true if there is more spectrum to read
            while (reader.hasNext()) {

                // next() returns the next spectrum or throws an IOException is something went wrong
                MsnSpectrum spectrum = reader.next();

                // do some stuff with your spectrum
                String peptide = spectrum.getComment();
                double[] intensities = new double[spectrum.size()];
                double[] mzs = new double[spectrum.size()];
                spectrum.getIntensities(intensities);
                spectrum.getMzs(mzs);

                //some fragments have zero intensity
                ArrayList<Integer> zeroFrags = new ArrayList<>();
                for (int i = 0; i < intensities.length; i++) {
                    if (intensities[i] == 0.0) {
                        zeroFrags.add(i);
                    }
                }

                if (zeroFrags.size() == 0) {
                    allPredMZs.put(peptide, doubleToFloat(mzs));
                    allPredIntensities.put(peptide, doubleToFloat(intensities));
                } else { //some empty frags
                    double[] newIntensities = new double[spectrum.size() - zeroFrags.size()];
                    double[] newMzs = new double[spectrum.size() - zeroFrags.size()];

                    int j = 0;
                    int k = 0;
                    int exclude = zeroFrags.get(j);
                    for (int i = 0; i < intensities.length; i++) {
                        if (i == exclude) {
                            j += 1;
                            try {
                                exclude = zeroFrags.get(j);
                            } catch (Exception e) {
                                exclude = -1; //no more empty frags
                            }
                        } else {
                            newIntensities[k] = intensities[i];
                            newMzs[k] = mzs[i];
                            k += 1;
                        }
                    }
                    allPredMZs.put(peptide, doubleToFloat(newMzs));
                    allPredIntensities.put(peptide, doubleToFloat(newIntensities));
                }

                //retention time, without determining number of sig figs
                float rt = (float) spectrum.getRetentionTimes().getFirst().getTime();
                allPredRTs.put(peptide, rt);
            }
            reader.close();
        }
    }

    public HashMap<String, float[]> getMzDict() { return allPredMZs; }

    public HashMap<String, float[]> getIntensityDict() { return allPredIntensities; }

    public HashMap<String, Float> getRtDict() { return allPredRTs; }

    public HashMap<String, Float> getIMDict() {
        System.out.println("no IM predictions");
        return null;
    }

    public float getMaxPredRT() { return Collections.max(allPredRTs.values()); }

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

    public static void main(String[] args) throws IOException, FileParsingException {

    }
}
