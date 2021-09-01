package Features;

import org.apache.commons.lang3.ArrayUtils;
import umich.ms.fileio.exceptions.FileParsingException;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;

public class mgfFileReader implements SpectralPredictionMapper{
    //mgfFileReader can handle both single files and entire directories

    final ArrayList<String> filenames;
    HashMap<String, PredictionEntry> allPreds = new HashMap<>();
    public HashMap<Integer, mzmlScanNumber> scanNumberObjects = new HashMap<>();

    //this version if loading pDeep3 predictions
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

        String line;
        String title = null;
        //int charge;
        double pepmass;
        //double precursorIntensity;
        float RT = 0;
        float IM = 0;
        ArrayList<Float> intensities = new ArrayList<>();
        ArrayList<Float> mzs = new ArrayList<>();
        String[] lineSplit;
        float[] finalMzs;
        float[] finalIntensities;

        for (String fname : filenames) {
//            MgfReader reader = new MgfReader(new File(fname), PeakList.Precision.FLOAT);
//            reader.acceptUnsortedSpectra();
            BufferedReader reader = new BufferedReader(new FileReader(fname));
//
//            // hasNext() returns true if there is more spectrum to read
//            while (reader.hasNext()) {
            while ((line = reader.readLine()) != null) {
//
//                // next() returns the next spectrum or throws an IOException is something went wrong
//                MsnSpectrum spectrum = reader.next();
//
//                // do some stuff with your spectrum

                if (line.contains("=")) {
                    lineSplit = line.split("=");
                    switch (lineSplit[0]) {
                        case "TITLE":
                            title = lineSplit[1];
                            break;
                        case "CHARGE":
//                            charge = Integer.parseInt(lineSplit[1].replace("+", ""));
                            break;
                        case "PEPMASS":
//                            lineSplit = lineSplit[1].split(" ");
//                            pepmass = Double.parseDouble(lineSplit[0]);
//                            allMasses.put(title, pepmass);
                            break;
                        case "RTINSECONDS":
                            RT = Float.parseFloat(lineSplit[1]);
                            break;
                        case "1/K0":
                            IM = Float.parseFloat(lineSplit[1]);
                            break;
                    }
                } else {
                    lineSplit = line.split(" ");
                    try { //fragment ions
                        mzs.add(Float.parseFloat(lineSplit[0]));
                        intensities.add(Float.parseFloat(lineSplit[1]));
                    } catch (Exception e) {
                        if (lineSplit[0].equals("END")) {
                            //unload
                            //some fragments have zero intensity
                            ArrayList<Integer> zeroFrags = new ArrayList<>();
                            for (int i = 0; i < intensities.size(); i++) {
                                if (intensities.get(i) == 0.0) {
                                    zeroFrags.add(i);
                                }
                            }

                            if (zeroFrags.size() == 0) {
                                finalMzs = ArrayUtils.toPrimitive(mzs.toArray(new Float[0]), 0.0F);
                                finalIntensities = ArrayUtils.toPrimitive(intensities.toArray(new Float[0]), 0.0F);
                            } else { //some empty frags
                                float[] newIntensities = new float[intensities.size() - zeroFrags.size()];
                                float[] newMzs = new float[intensities.size() - zeroFrags.size()];

                                int j = 0;
                                int k = 0;
                                int exclude = zeroFrags.get(j);
                                for (int i = 0; i < intensities.size(); i++) {
                                    if (i == exclude) {
                                        j += 1;
                                        try {
                                            exclude = zeroFrags.get(j);
                                        } catch (Exception e1) {
                                            exclude = -1; //no more empty frags
                                        }
                                    } else {
                                        newIntensities[k] = intensities.get(i);
                                        newMzs[k] = mzs.get(i);
                                        k += 1;
                                    }
                                }
                                finalMzs = newMzs;
                                finalIntensities = newIntensities;
                            }

                            PredictionEntry newPred = new PredictionEntry();
                            newPred.setMzs(finalMzs);
                            newPred.setIntensities(finalIntensities);
                            newPred.setRT(RT);
                            newPred.setIM(IM);
                            allPreds.put(title, newPred);

                            //reset for next peptide/PSM
                            mzs.clear();
                            intensities.clear();
                        }
                    }
                }
            }
            reader.close();
        }
    }

    //this version for uncalibrated mgf acting as mzml
    public mgfFileReader(String files, boolean createScanNumObjects) throws IOException, FileParsingException {
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

        String line;
        String title = null;
        int scanNum = 0;
        float RT = 0;
        float IM = 0;
        ArrayList<Float> intensities = new ArrayList<>();
        ArrayList<Float> mzs = new ArrayList<>();
        String[] lineSplit;
        for (String fname : filenames) {
            BufferedReader reader = new BufferedReader(new FileReader(fname));

            // hasNext() returns true if there is more spectrum to read
            while ((line = reader.readLine()) != null) {
                if (line.contains("=")) {
                    lineSplit = line.split("=");
                    switch (lineSplit[0]) {
                        case "TITLE":
                            title = lineSplit[1];
                            String[] dotSplit = title.split("\\.");
                            scanNum = Integer.parseInt(dotSplit[dotSplit.length - 2]);
                            break;
                        case "CHARGE":
                            break;
                        case "PEPMASS":
                            break;
                        case "RTINSECONDS":
                            RT = Float.parseFloat(lineSplit[1]);
                            break;
                        case "1/K0":
                            IM = Float.parseFloat(lineSplit[1]);
                            break;
                    }
                } else {
                    lineSplit = line.split(" ");
                    try { //fragment ions
                        mzs.add(Float.parseFloat(lineSplit[0]));
                        intensities.add(Float.parseFloat(lineSplit[1]));
                    } catch (Exception e) {
                        if (lineSplit[0].equals("END")) {
                            //unload
                            //some fragments have zero intensity
                            ArrayList<Integer> zeroFrags = new ArrayList<>();
                            for (int i = 0; i < intensities.size(); i++) {
                                if (intensities.get(i) == 0.0) {
                                    zeroFrags.add(i);
                                }
                            }

                            if (zeroFrags.size() == 0) {
                                scanNumberObjects.put(scanNum, new mzmlScanNumber(scanNum, ArrayUtils.toPrimitive(mzs.toArray(new Float[0]), 0.0F),
                                        ArrayUtils.toPrimitive(intensities.toArray(new Float[0]), 0.0F), RT, IM));
                            } else { //some empty frags
                                float[] newIntensities = new float[intensities.size() - zeroFrags.size()];
                                float[] newMzs = new float[intensities.size() - zeroFrags.size()];

                                int j = 0;
                                int k = 0;
                                int exclude = zeroFrags.get(j);
                                for (int i = 0; i < intensities.size(); i++) {
                                    if (i == exclude) {
                                        j += 1;
                                        try {
                                            exclude = zeroFrags.get(j);
                                        } catch (Exception e1) {
                                            exclude = -1; //no more empty frags
                                        }
                                    } else {
                                        newIntensities[k] = intensities.get(i);
                                        newMzs[k] = mzs.get(i);
                                        k += 1;
                                    }
                                }
                                scanNumberObjects.put(scanNum, new mzmlScanNumber(scanNum, newMzs,
                                        newIntensities, RT, IM));
                            }

                            //reset for next peptide/PSM
                            mzs.clear();
                            intensities.clear();
                        }
                    }
                }
            }
            reader.close();
        }
    }

    public HashMap<String, PredictionEntry> getPreds() { return allPreds; }

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
    }

    public static void main(String[] args) throws IOException, FileParsingException {
        //mgfFileReader mgf = new mgfFileReader("C:/Users/kevin/OneDriveUmich/proteomics/preds/cptacPreds.mgf");
        mgfFileReader mgf = new mgfFileReader("C:/Users/kevin/Downloads/proteomics/timsTOF/" +
                "20180819_TIMS2_12-2_AnBr_SA_200ng_HeLa_50cm_120min_100ms_11CT_3_A1_01_2769_uncalibrated.mgf");
    }
}
