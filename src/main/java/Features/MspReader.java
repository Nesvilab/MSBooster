package Features;

import java.io.*;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashMap;

public class MspReader implements SpectralPredictionMapper {
    final ArrayList<String> filenames;
    private HashMap<String, float[]> allPredMZs = new HashMap<>();
    private HashMap<String, float[]> allPredIntensities = new HashMap<>();
    private HashMap<String, Float> allPredRTs = new HashMap<>();

    private String prositPepFormat(String pep) {
        String[] pepSplit = pep.split("/");
        StringBuilder sb = new StringBuilder();
        sb.append(pepSplit[0]).append("|"); //base pep
        if (!pepSplit[2].equals("")) {
            String[] mods = pepSplit[2].split("; ");
            for (String mod : mods) {
                String[] modSplit = mod.split("@");
                String aa = modSplit[1].substring(0, 1);
                String position = modSplit[1].substring(1, modSplit.length);

                sb.append(position).append(",").append(modSplit[0]).append("[").append(aa).append("];");
            }
        }
        sb.append("|").append(pepSplit[3]); //charge

        return sb.toString();
    }

    public MspReader(String files) throws FileNotFoundException {
        File predsDirectory = new File(files);
        String[] predsFiles = predsDirectory.list();
        filenames = new ArrayList<String>();

        if (predsFiles == null) { //if user provided a file, not a directory
            filenames.add(files);
        } else { //user provides directory
            for (String predsFile : predsFiles) {
                if (predsFile.contains(".msp")) {
                    filenames.add(files + File.separator + predsFile);
                }
            }
        }

        String line = null;
        for (String f : filenames) {
            try (
                    BufferedReader msp = new BufferedReader(new FileReader(f));
            ) {

                while ((line = msp.readLine()) != null) { //Name
                    msp.readLine(); //MW
                    line = msp.readLine(); //comment
                    String[] lineSplit = line.split(" ");

                    //peptide
                    StringBuilder sb = new StringBuilder();
                    for (String l : Arrays.copyOfRange(lineSplit, 4, lineSplit.length - 1)) {
                        sb.append(l);
                    }
                    String pep = prositPepFormat(sb.toString().split("=")[1]);

                    //RT
                    float RT = Float.parseFloat(lineSplit[lineSplit.length - 1].split("=")[1]);

                    //fragments
                    line = msp.readLine(); //Num peaks
                    int numPeaks = Integer.parseInt(line.split(" ")[2]);
                    float[] mzs = new float[numPeaks];
                    float[] ints = new float[numPeaks];

                    for (int i = 0; i < numPeaks; i++) {
                        line = msp.readLine();
                        lineSplit = line.split("\t");
                        mzs[i] = Float.parseFloat(lineSplit[0]);
                        ints[i] = Float.parseFloat(lineSplit[1]);
                    }

                    allPredMZs.put(pep, mzs);
                    allPredIntensities.put(pep, ints);
                    allPredRTs.put(pep, RT);
                }
            } catch (IOException e) {
                e.printStackTrace();
            }
        }
    }

    public HashMap<String, float[]> getMzDict() { return allPredMZs; }

    public HashMap<String, float[]> getIntensityDict() { return allPredIntensities; }

    public HashMap<String, Float> getRtDict() { return allPredRTs; }

    public HashMap<String, Float> getIMDict() {
        System.out.println("no IM predictions");
        return null;
    }

    public float getMaxPredRT() { return Collections.max(allPredRTs.values()); }

    public static void main(String[] args) throws IOException {
        //MspReader m = new MspReader("C:/Users/kevin/OneDriveUmich/proteomics/preds/cptacProsit.msp");
        SpectralPredictionMapper m = SpectralPredictionMapper.createSpectralPredictionMapper(
                "C:/Users/kevin/OneDriveUmich/proteomics/preds/cptacProsit.msp");
    }
}
