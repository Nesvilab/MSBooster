package Features;

import java.io.IOException;
import java.util.ArrayList;
import java.util.Comparator;
import java.util.TreeSet;
import java.util.stream.IntStream;

public class IMFunctions {
    public static double[][][] getIMarrays(mzMLReader mzml, int IMregressionSize) {
        final int numCharges = 7;
        double[][][] finalIMs = new double[numCharges][2][];

        ArrayList<ArrayList<Float>> expIMs = new ArrayList<ArrayList<Float>>(numCharges);
        ArrayList<ArrayList<Float>> predIMs = new ArrayList<ArrayList<Float>>(numCharges);
        ArrayList<ArrayList<Float>> eScores = new ArrayList<ArrayList<Float>>(numCharges); //for sorting
        for (int i = 0; i < numCharges; i++) {
            expIMs.add(new ArrayList<>());
            predIMs.add(new ArrayList<>());
            eScores.add(new ArrayList<>());
        }

        //collect RTs and escores

        for (int scanNum : new TreeSet<Integer>(mzml.scanNumberObjects.keySet())) {
            mzmlScanNumber scanNumObj = mzml.getScanNumObject(scanNum);
            float im = scanNumObj.IM; //experimental RT for this scan

            for (int i = 1; i < scanNumObj.peptideObjects.size() + 1; i++) {
                peptideObj pep = scanNumObj.getPeptideObject(i);
                int charge = pep.charge - 1;
                float e = Float.parseFloat(pep.escore);
                if (e > Constants.IMescoreCutoff) {
                    break;
                }
                expIMs.get(charge).add(im);
                predIMs.get(charge).add(pep.IM);
                eScores.get(charge).add(e);
            }
        }

        //get top peptides based on eScore
        //https://stackoverflow.com/questions/4859261/get-the-indices-of-an-array-after-sorting
        //also consider taking them all, if want more samples for regression/z scoring
        //then divide into bins with constant size (higher precursor density in middle of RT)

        //if negative, use all
        //can consider e score cutoff in constants
        for (int c = 0; c < numCharges; c++) {
            int sizeLimit = expIMs.get(c).size();
            if (sizeLimit == 0) { //no PSMs with this charge
                continue;
            }

            if (IMregressionSize > 0 && IMregressionSize <= sizeLimit) {
                //get top e score PSMs
                int[] sortedIndices = IntStream.range(0, eScores.get(c).size())
                        .boxed().sorted(Comparator.comparing(eScores.get(c)::get))
                        .mapToInt(ele -> ele).toArray();

                double[][] IMs = new double[2][IMregressionSize];
                for (int i = 0; i < IMregressionSize; i++) {
                    int idx = sortedIndices[i];
                    IMs[0][i] = expIMs.get(c).get(idx);
                    IMs[1][i] = predIMs.get(c).get(idx);
                }

                //sort by experimental IM
                int[] sortedIndices2 = IntStream.range(0, IMregressionSize)
                        .boxed().sorted(Comparator.comparingDouble(k -> IMs[0][k]))
                        .mapToInt(ele -> ele).toArray();

                double[][] IMs2 = new double[2][IMregressionSize];
                for (int i = 0; i < IMregressionSize; i++) {
                    int idx = sortedIndices2[i];
                    IMs2[0][i] = IMs[0][idx];
                    IMs2[1][i] = IMs[1][idx];
                }

                 finalIMs[c] = IMs2;

            } else {
                //sort by experimental IM
                int[] sortedIndices = IntStream.range(0, expIMs.get(c).size())
                        .boxed().sorted(Comparator.comparing(expIMs.get(c)::get))
                        .mapToInt(ele -> ele).toArray();

                double[][] IMs = new double[2][sizeLimit];
                for (int i = 0; i < sizeLimit; i++) {
                    int idx = sortedIndices[i];
                    IMs[0][i] = expIMs.get(c).get(idx);
                    IMs[1][i] = predIMs.get(c).get(idx);
                }

                finalIMs[c] = IMs;
            }
        }

        return finalIMs;
    }

    public static ArrayList<Float>[] IMbins(mzMLReader mzml) throws IOException {
        //hard coded as 2, but if there are higher IM values, this can change
        int numBins = 2 * Constants.IMbinMultiplier;

        ArrayList<Float>[] predIMround = new ArrayList[numBins + 1];
        for (int col = 0; col < numBins + 1; col++) {
            predIMround[col] = new ArrayList<Float>();
        }

        //iterate through scanNumbers
        for (int scanNum : mzml.scanNumberObjects.keySet()) {
            mzmlScanNumber scanNumObj = mzml.getScanNumObject(scanNum);
            int round = (int) (scanNumObj.IM *  Constants.IMbinMultiplier); //experimental RT for this scan, assume in minutes

            //iterate through PSMs
            for (int i = 1; i < scanNumObj.peptideObjects.size() + 1; i++) {
                peptideObj pep = scanNumObj.getPeptideObject(i);

                int instances = Math.max(1, -1 * (int) Math.ceil(Math.log10(Double.parseDouble(pep.escore)))); //this version avoids empty bins
                for (int j = 0; j < instances; j++) {
                    predIMround[round].add(pep.IM);
                }
            }
        }

        return predIMround;
    }
}
