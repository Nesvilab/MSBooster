package Features;

public class PredictionEntry {
    float[] mzs;
    float[] intensities;
    float RT;
    float IM;

    public PredictionEntry(float[] mzs, float[] intensities, float RT, float IM) {
        this.mzs = mzs;
        this.intensities = intensities;
        this.RT = RT;
        this.IM = IM;
    }
}
