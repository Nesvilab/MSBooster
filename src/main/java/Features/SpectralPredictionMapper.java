package Features;

import java.io.IOException;
import java.util.HashMap;

public interface SpectralPredictionMapper {

    public HashMap<String, float[]> getMzDict() throws IOException;

    public HashMap<String, float[]> getIntensityDict() throws IOException;;

    public HashMap<String, Float> getRtDict() throws IOException;

    public HashMap<String, Float> getIMDict() throws IOException;

    public float getMaxPredRT();

    public static SpectralPredictionMapper createSpectralPredictionMapper(String file) throws IOException {
        if (Constants.predFileFormat.equals("bin") || Constants.predFileFormat.equals("DIA-NN")) {
            return new DiannSpeclibReader(file);
        } else if (Constants.predFileFormat.equals("mgf") || Constants.predFileFormat.equals("pDeep3")) {
            return new mgfFileReader(file);
        } else {
            System.out.println("not valid prediction file format");
            return null;
        }
    }
}
