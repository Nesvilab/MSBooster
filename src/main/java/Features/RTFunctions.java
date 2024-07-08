/*
 * This file is part of MSBooster.
 *
 * MSBooster is free software: you can redistribute it and/or modify it under
 * the terms of the GNU Lesser General Public License as published by the Free
 * Software Foundation, either version 3 of the License, or (at your
 * option) any later version.
 *
 * MSBooster is distributed in the hope that it will be useful, but WITHOUT
 * ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
 * FITNESS FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public License
 * for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with MSBooster. If not, see <https://www.gnu.org/licenses/>.
 */

package Features;

import static utils.Print.printError;
import static utils.Print.printInfo;

import umich.ms.fileio.exceptions.FileParsingException;
import umontreal.ssj.probdist.EmpiricalDist;

import java.io.IOException;
import java.util.*;
import java.util.stream.IntStream;

public class RTFunctions {
    //given experimental RT, what should it be on predicted scale?
    public static float normalizeRT(float[] betas, float expRT) {
        return betas[0] + betas[1] * expRT;
    }

//    public static float[] normalizeRT(float[] betas, float[] expRTs) {
//        float[] newRTs = new float[expRTs.length];
//        for (int i = 0; i < expRTs.length; i++) {
//            newRTs[i] = normalizeRT(betas, expRTs[i]);
//        }
//        return newRTs;
//    }
//
//    public static float[] calculateDeltaRT(float[] newRTs, float[] predRTs) {
//        //newRTs can also be array of mean predRTs from bins
//        float[] deltaRTs = new float[predRTs.length];
//        for (int i = 0; i < predRTs.length; i++) {
//            deltaRTs[i] = Math.abs(newRTs[i] - predRTs[i]);
//        }
//        return deltaRTs;
//    }

    //given mzmlreader, get all peptide objects and their RTs
    //get predicted RTs
    //assumes peptide objects already set
    //TODO: only supporting getting PSMs for supported PTMs
    public static float[] getBetas(MzmlReader mzml, int RTregressionSize) throws FileParsingException {
        double[][] RTs = LoessUtilities.getArrays(mzml, RTregressionSize, "RT", 0).get("");
        return StatMethods.linearRegression(RTs[0], RTs[1]);
    }

    public static float[] getBetas(double[][] RTs) {
        return StatMethods.linearRegression(RTs[0], RTs[1]);
    }

    //function that returns double[] of bin boundaries, with mean and var od each
    //exp on x axis, pred on y
    //to do: test for normal dist nature of each bin (Shapiro Wilk test)
    public static ArrayList<Float>[] RTbins(MzmlReader mzml) throws IOException, FileParsingException {
        //get max index
        float maxRT = -1;
        for (int num : mzml.getScanNums()) {
            MzmlScanNumber msn = mzml.getScanNumObject(num);
            if (msn.RT > maxRT) {
                maxRT = msn.RT;
            }
        }
        int numBins = (int) (maxRT * Constants.RTbinMultiplier) + 1;

        ArrayList<Float>[] predRTround = new ArrayList[numBins + 1];
        for (int col = 0; col < numBins + 1; col++) {
            predRTround[col] = new ArrayList<Float>();
        }

        //iterate through scanNumbers
        for (int scanNum : mzml.getScanNums()) {
            MzmlScanNumber scanNumObj = mzml.getScanNumObject(scanNum);
            int round = (int) (scanNumObj.RT *  Constants.RTbinMultiplier); //experimental RT for this scan, assume in minutes

            //iterate through PSMs
            for (PeptideObj pep : scanNumObj.peptideObjects) {
                if (pep == null) {
                    continue;
                }

                int instances = Math.max(1, -1 * (int) Math.ceil(Math.log10(Double.parseDouble(pep.escore)))); //this version avoids empty bins
                for (int j = 0; j < instances; j++) {
                    predRTround[round].add(pep.RT);
                }
            }
        }

        return predRTround;
    }

    private static EmpiricalDist generateEmpiricalDist(ArrayList<Float> bin) { //could also use umontreal version
        int binSize = bin.size();

        if (binSize > 1) {
            //check that this is sorted
            float maxFloat = bin.get(0);
            boolean needToSort = false;
            for (float f : bin) {
                if (f >= maxFloat) {
                    maxFloat = f;
                } else {
                    needToSort = true;
                    break;
                }
            }
            if (needToSort) {
                Collections.sort(bin);
            }
            return new EmpiricalDist(FloatUtils.floatToDouble(bin));
        } else if (binSize == 1) {
            return new EmpiricalDist(new double[]{(double) bin.get(0), (double) bin.get(0)});
        } else {
            return null; //empty bin never accessed
        }
        //get probability at point with KernelDensity.p();
    }

    public static EmpiricalDist[] generateEmpiricalDist(ArrayList<Float>[] bins) {
        EmpiricalDist[] empDists = new EmpiricalDist[bins.length];

        int i = 0;
        for (ArrayList<Float> bin : bins) {
            empDists[i] = generateEmpiricalDist(bin);
            i++;
        }

        return empDists;
    }

//    public static double RTprobability(float expRT, float predRT, KernelDensity[] bins) {
//        //get right bin to search
//        KernelDensity kd = bins[Math.round(expRT)];
//
//        //check probability at point
//        try {
//            return kd.p(predRT);
//        } catch (Exception e) { //nothing in bin
//            return 0;
//        }
//    }
//
//    public static float RTprobabilityWithUniformPrior(int unifPriorSize, float unifProb,
//                                                       int binSize, float empiricalProb) {
//        float w1 = (float) unifPriorSize / (float) (unifPriorSize + binSize);
//        float w2 = (float) binSize / (float) (unifPriorSize + binSize);
//
//        return w1 * unifProb + w2 * empiricalProb;
//    }
}
