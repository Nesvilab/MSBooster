package Features;

import java.util.stream.IntStream;

public class PredictionEntry {
    float[] mzs;
    float[] intensities;
    int[] fragNums;
    int[] flags;
    int[] charges;
    String[] fragmentIonTypes;
    float RT;
    float IM;
    float detectability;
    int counter;

    public PredictionEntry() {}

    public PredictionEntry(float[] mzs, float[] intensities, boolean needToSort) {
        if (needToSort) {
            this.mzs = new float[mzs.length];
            this.intensities = new float[intensities.length];

            int[] sortedIndices = IntStream.range(0, mzs.length)
                    .boxed().sorted((k, j) -> Float.compare(mzs[k], mzs[j]))
                    .mapToInt(ele -> ele).toArray();

            for (int i = 0; i < sortedIndices.length; i++) {
                this.mzs[i] = mzs[sortedIndices[i]];
                this.intensities[i] = intensities[sortedIndices[i]];
            }
        } else {
            setMzs(mzs);
            setIntensities(intensities);
        }
    }

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
