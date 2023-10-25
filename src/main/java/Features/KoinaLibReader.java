package Features;

import java.io.IOException;
import java.util.concurrent.ConcurrentHashMap;

public class KoinaLibReader implements SpectralPredictionMapper {
    ConcurrentHashMap<String, PredictionEntry> allPreds = new ConcurrentHashMap<>();

    public KoinaLibReader() {}

    public ConcurrentHashMap<String, PredictionEntry> getPreds() {return allPreds;}

    public void clear() {
        allPreds.clear();
    }

    public float getMaxPredRT() {
        return 0;
    }
}
