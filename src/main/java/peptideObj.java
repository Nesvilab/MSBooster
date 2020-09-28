import umich.ms.fileio.exceptions.FileParsingException;

import java.io.IOException;
import java.lang.reflect.InvocationTargetException;
import java.lang.reflect.Method;
import java.util.HashMap;

public class peptideObj {
    final String name;
    final int rank;
    final mzmlScanNumber scanNumObj;
    final int targetORdecoy;
    spectrumComparison spectralSimObj;
    HashMap<String, Double> scores = new HashMap<>();

    public peptideObj(mzmlScanNumber scanNumObj, String name, int rank, int targetORdecoy,
                      double[] predMZs, double[] predIntensities) {
        this.name = name;
        this.rank = rank;
        this.scanNumObj = scanNumObj;
        this.targetORdecoy = targetORdecoy;
        this.spectralSimObj = new spectrumComparison(scanNumObj.expMZs, scanNumObj.expIntensities,
                predMZs, predIntensities);
    }

    public void setScore(String similarityMeasure) throws FileParsingException, NoSuchMethodException,
            InvocationTargetException, IllegalAccessException {

        //get similarity
        //check that similarityMeasure is valid

        double sim = 0;
        //only if need weights
        if (similarityMeasure.substring(0, 8).equals("weighted")) {
            Method method = spectralSimObj.getClass().getMethod(similarityMeasure, double[].class);
            double[] weights = spectralSimObj.getWeights(scanNumObj.mzmlScans.getMzFreq());
            sim = (double) method.invoke(spectralSimObj, weights);
        } else {
            Method method = spectralSimObj.getClass().getMethod(similarityMeasure);
            sim = (double) method.invoke(spectralSimObj);
        }
        scores.put(similarityMeasure, sim);
    }

    public double getScore(String similarityMeasure) throws NoSuchMethodException, FileParsingException,
            IllegalAccessException, InvocationTargetException {
        if (! scores.containsKey(similarityMeasure)) {
            setScore(similarityMeasure);
        }
        return scores.get(similarityMeasure);
    }

    public static void main(String[] args) throws FileParsingException, IOException, NoSuchMethodException,
            IllegalAccessException, InvocationTargetException {
//        //can check if getMZFreq needs to be run once or multiple times
//        mzMLReader x = new mzMLReader("C:/Users/kevin/OneDriveUmich/proteomics/mzml/" +
//                "23aug2017_hela_serum_timecourse_4mz_narrow_1.mzml");
//        x.createScanNumObjects();
//
//        //load in predicted spectra
//        mgfFileReader predictedSpectra = new mgfFileReader("preds/rank1_preds.txt");
//
//        //get predicted MZ's and intensities in form of HashMaps
//        HashMap<String, double[]> predictedMZs = predictedSpectra.getMzDict();
//        HashMap<String, double[]> predictedIntensities = predictedSpectra.getIntensityDict();
//
//        x.scanNumberObjects.get(5047).setPeptideObject("NGFGKNGFLQSR||3", 1,
//                predictedMZs, predictedIntensities);
//        x.scanNumberObjects.get(5047).getPeptideObject(1).setScore("weightedBrayCurtis");
//        System.out.println(Arrays.toString(x.mzFreqs));
//        x.scanNumberObjects.get(5047).getPeptideObject(1).setScore("weightedEuclideanDistance");

    }
}
