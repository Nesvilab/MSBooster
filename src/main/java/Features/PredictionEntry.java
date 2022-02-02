package Features;

public class PredictionEntry {
    float[] mzs;
    float[] intensities;
    int[] fragNums;
    int[] flags;
    int[] charges;
    String[] fragmentIonTypes;
    MassCalculator massCalculator;
    float RT;
    float IM;
    float detectability;
    int counter;

    public PredictionEntry() {}

    public void setMzs(float[] mzs) {this.mzs = mzs;}

    public void setIntensities(float[] intensities) {this.intensities = intensities;}

    public void setFragNums(int[] fragNums) {this.fragNums = fragNums;}

    public void setFlags(int[] flags) {this.flags = flags;}

    public void setCharges(int[] charges) {this.charges = charges;}

    public void setRT(float RT) {this.RT = RT;}

    public void setIM(float IM) {this.IM = IM;}

    public void setDetectability(float detectability) {this.detectability = detectability;}

    public void setCounter(int counts) {this.counter = counts;}

    public void setFragmentIonTypes(String[] ions) {this.fragmentIonTypes = ions;}
}
