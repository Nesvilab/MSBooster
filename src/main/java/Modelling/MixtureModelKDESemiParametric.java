package Modelling;

import Features.peptideObjLite;
import umich.ms.fileio.exceptions.FileParsingException;

import java.io.IOException;
import java.lang.reflect.InvocationTargetException;
import java.util.ArrayList;

import static Features.PeptideObjGenerator.generatePeptideObj;

public class MixtureModelKDESemiParametric {
    ArrayList<Double> tExpect;
    ArrayList<Double> dExpect;
    ArrayList<Double> tDeltaRT;
    ArrayList<Double> dDeltaRT;
    ArrayList<Double> tMassDiff;
    ArrayList<Double> dMassDiff;
    ArrayList<Double> tNTT;
    ArrayList<Double> dNTT;
    ArrayList<Double> tNMC;
    ArrayList<Double> dNMC;
    ArrayList<Double> tIsoError;
    ArrayList<Double> dIsoError;
    ArrayList<Double> tSim;
    ArrayList<Double> dSim;

    public MixtureModelKDESemiParametric(ArrayList<peptideObjLite> pObjs) {
        System.out.println("Starting");
        int pObjSize = pObjs.size();
        dExpect = new ArrayList<>(pObjSize);
        tExpect = new ArrayList<>(pObjSize);

        for (peptideObjLite pObj : pObjs) {
            if (pObj.targetORdecoy == 0) {
                dExpect.add(Double.valueOf(pObj.escore));
            } else {
                tExpect.add(Double.valueOf(pObj.escore));
            }
        }
    }

    //generate instance from file

    public static void main(String[] args) throws InvocationTargetException, NoSuchMethodException,
            FileParsingException, IllegalAccessException, IOException {
        MixtureModelKDESemiParametric model = new MixtureModelKDESemiParametric(
                generatePeptideObj("C:/Users/kevin/Downloads/proteomics/wideWindowPreds.mgf",
                "C:/Users/kevin/Downloads/proteomics/wideWindow",
                "C:/Users/kevin/Downloads/proteomics/wideWindow",
                "brayCurtis"));
    }
}
