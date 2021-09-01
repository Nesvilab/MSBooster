package Features;

import java.io.IOException;
import java.util.HashMap;

public interface SpectralPredictionMapper {

    public HashMap<String, PredictionEntry> getPreds() throws IOException;
    public void clear();

    public float getMaxPredRT();

    public static SpectralPredictionMapper createSpectralPredictionMapper(String file) throws IOException {
        //detecting file extension
        String[] extensionSplit = file.split("\\.");
        String extension = extensionSplit[extensionSplit.length - 1];

        if (extension.equals("bin")) {
            return new DiannSpeclibReader(file);
        } else if (extension.equals("mgf")) {
            return new mgfFileReader(file);
        } else if (extension.equals("msp")) {
            return new MspReader(file);
        } else {
            System.out.println("not valid prediction file format");
            return null;
        }
    }
}
