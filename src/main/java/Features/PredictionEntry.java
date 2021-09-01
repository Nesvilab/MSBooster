package Features;

public class PredictionEntry {
    float[] mzs;
    float[] intensities;
    float RT;
    float IM;
    float detectability;
    int counter;

    //TODO: make all set methods
    public PredictionEntry() {}

    public void setMzs(float[] mzs) {this.mzs = mzs;}

    public void setIntensities(float[] intensities) {this.intensities = intensities;}

    public void setRT(float RT) {this.RT = RT;}

    public void setIM(float IM) {this.IM = IM;}

    public void setDetectability(float detectability) {this.detectability = detectability;}

    public void setCounter(int counts) {this.counter = counts;}
}
