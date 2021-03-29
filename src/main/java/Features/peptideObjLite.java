package Features;

import java.util.HashMap;

public class peptideObjLite {
    final String name;
    final int rank;
    final int scanNum;
    final String mzml;
    final public int targetORdecoy;
    final public String escore;
    //final double RT;
    public double deltaRT;
    public double massDiff;
    public double ntt;
    public double nmc;
    public double isoError;
    public HashMap<String, Double> scores = new HashMap<>();

    public peptideObjLite(peptideObj pObj) {
        this.name = pObj.name;
        this.rank = pObj.rank;
        this.scanNum = pObj.scanNum;
        this.mzml = pObj.mzml;
        this.targetORdecoy = pObj.targetORdecoy;
        this.escore = pObj.escore;
        //this.RT = pObj.RT;
        this.deltaRT = pObj.deltaRT;
        this.massDiff = pObj.massDiff;
        this.ntt = pObj.ntt;
        this.nmc = pObj.nmc;
        this.isoError = pObj.isoError;
        this.scores = pObj.scores;

    }
}
