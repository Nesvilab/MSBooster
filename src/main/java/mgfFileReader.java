import org.expasy.mzjava.core.io.ms.spectrum.MgfReader;
import org.expasy.mzjava.core.ms.peaklist.PeakList;
import org.expasy.mzjava.core.ms.spectrum.MsnSpectrum;

import java.io.File;
import java.io.IOException;
import java.util.HashMap;


public class mgfFileReader {
    String fname;
    HashMap<String, double[]> intensityDict;
    HashMap<String, double[]> mzDict;
    boolean alreadyCreated;

    public mgfFileReader(String filename) {
        fname = filename;
        intensityDict = new HashMap<>();
        mzDict = new HashMap<>();
        alreadyCreated = false;
    }

    public void createDicts() throws IOException {
        String mgfFilename = fname;
        MgfReader reader = new MgfReader(new File(mgfFilename), PeakList.Precision.DOUBLE);

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
            intensityDict.put(peptide, intensities);
            mzDict.put(peptide, mzs);
        }

        reader.close();
        alreadyCreated = true;
    }

    public HashMap<String, double[]> getIntensityDict() throws IOException {
        if (! alreadyCreated) {
            this.createDicts();
        }
        return intensityDict;
    }

    public HashMap<String, double[]> getMzDict() throws IOException {
        if (! alreadyCreated) {
            this.createDicts();
        }
        return mzDict;
    }
}
