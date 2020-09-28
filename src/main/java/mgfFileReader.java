import org.expasy.mzjava.core.io.ms.spectrum.MgfReader;
import org.expasy.mzjava.core.ms.peaklist.PeakList;
import org.expasy.mzjava.core.ms.spectrum.MsnSpectrum;

import java.io.File;
import java.io.IOException;
import java.util.Arrays;
import java.util.HashMap;
import java.util.Map;

public class mgfFileReader {
    //mgfFileReader can handle both single files and entire directories

    final String[] filenames;
    HashMap<String, double[]> allPredMZs = new HashMap<>();
    HashMap<String, double[]> allPredIntensities = new HashMap<>();

    public mgfFileReader(String files) throws IOException {
        File predsDirectory = new File(files);
        String[] predsFiles = predsDirectory.list();

        if (predsFiles == null) { //if user provided a file, not a directory
            filenames = new String[] {files};
        } else { //user provides directory
            String[] newPredsFiles = new String[predsFiles.length];
            for (int i = 0; i < predsFiles.length; i++) {
                newPredsFiles[i] = files + '/' + predsFiles[i];
            }
            filenames = newPredsFiles;
        }

        this.createDicts();
    }

    private void createDicts() throws IOException {
        for (String fname : filenames) {
            MgfReader reader = new MgfReader(new File(fname), PeakList.Precision.DOUBLE);
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

                //if we decide to filter out low intensity fragments, do it here

                allPredMZs.put(peptide, mzs);
                allPredIntensities.put(peptide, intensities);
            }
            reader.close();
        }
    }

    public HashMap<String, double[]> getMzDict() throws IOException {
        return allPredMZs;
    }

    public HashMap<String, double[]> getIntensityDict() throws IOException {
        return allPredIntensities;
    }

    public static void main(String[] args) throws IOException {
        mgfFileReader x = new mgfFileReader("C:/Users/kevin/Downloads/proteomics/test2.mgf");
        //mgfFileReader x = new mgfFileReader("preds/");
        for (Map.Entry<String, double[]> e : x.allPredMZs.entrySet()) {
            System.out.println(Arrays.toString(e.getValue()));
        }
    }
}
