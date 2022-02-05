package Features;

import umich.ms.fileio.exceptions.FileParsingException;

import java.io.*;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.concurrent.ExecutionException;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;

public class MspReader implements SpectralPredictionMapper {
    final ArrayList<String> filenames;
    HashMap<String, PredictionEntry> allPreds = new HashMap<>();

    //convert int flag to fragment ion type
    private static HashMap<Integer, String> makeFlagTOion() {
        HashMap<Integer, String> map = new HashMap<>();
        map.put(0, "b");
        map.put(1, "y");
        return map;
    }
    HashMap<Integer, String> flagTOion = makeFlagTOion();

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

        String line;
        for (String f : filenames) {
            try (BufferedReader msp = new BufferedReader(new FileReader(f))) {

                while ((line = msp.readLine()) != null) { //Name
                    StringBuilder sb = new StringBuilder();
                    String[] lineSplit = line.split(": ")[1].split("/");
                    sb.append(lineSplit[0]).append("|");
                    String charge = lineSplit[1];

                    msp.readLine(); //MW
                    line = msp.readLine(); //comment
                    lineSplit = line.split(" ");

                    String pep;
                    float RT = 0f;

                    for (String s : lineSplit) {
                        String[] commentContents = s.split("=");
                        switch (commentContents[0]) {
                            case "Mods":
                                sb.append(commentContents[1]);
                                break;
                            case "iRT":
                                RT = Float.parseFloat(commentContents[1]);
                                break;
                        }
                    }
                    pep = new PeptideFormatter(sb.toString(), charge, "msp").baseCharge;

                    //fragments
                    //if it keeps reporting 174 peaks, set that as initial arraylist size
                    line = msp.readLine(); //Num peaks
                    int numPeaks = Integer.parseInt(line.split(" ")[2]);
                    ArrayList<Float> mzsList = new ArrayList<>(numPeaks);
                    ArrayList<Float> intsList = new ArrayList<>(numPeaks);
                    ArrayList<Integer> fragNumsList = new ArrayList<>(numPeaks);
                    ArrayList<Integer> flagsList = new ArrayList<>(numPeaks);
                    ArrayList<Integer> chargesList = new ArrayList<>(numPeaks);

                    for (int i = 0; i < numPeaks; i++) {
                        line = msp.readLine();
                        lineSplit = line.split("\t");
                        float tmpInt = Float.parseFloat(lineSplit[1]);
                        if (! (tmpInt == 0f)) {
                            mzsList.add(Float.parseFloat(lineSplit[0]));
                            intsList.add(tmpInt);
                        }
                        String fragment = lineSplit[2].substring(1, lineSplit[2].length() - 1).split("/")[0];
                        String[] chargeSplit = fragment.split("\\^");
                        if (chargeSplit.length == 1) {
                            chargesList.add(1);
                        } else {
                            chargesList.add(Integer.parseInt(chargeSplit[1]));
                        }
                        if (fragment.charAt(0) == 'y') {
                            flagsList.add(1);
                        } else {
                            flagsList.add(0);
                        }
                        fragNumsList.add(Integer.parseInt(chargeSplit[0].substring(1)));
                    }

                    float[] mzs = new float[mzsList.size()];
                    float[] ints = new float[intsList.size()];
                    int[] fragNums = new int[intsList.size()];
                    int[] flags = new int[intsList.size()];
                    int[] charges = new int[intsList.size()];
                    for (int i = 0; i < mzs.length; i++) {
                        mzs[i] = mzsList.get(i);
                        ints[i] = intsList.get(i);
                        fragNums[i] = fragNumsList.get(i);
                        flags[i] = flagsList.get(i);
                        charges[i] = chargesList.get(i);
                    }

                    PredictionEntry newPred = new PredictionEntry();
                    newPred.setMzs(mzs);
                    newPred.setIntensities(ints);
                    newPred.setRT(RT);
                    newPred.setIM(0f);
                    newPred.setFragNums(fragNums);
                    newPred.setFlags(flags);
                    newPred.setCharges(charges);
                    allPreds.put(pep, newPred);
                }

                //read in extra files
                BufferedReader TSVReader = new BufferedReader(new FileReader(
                        f.substring(0, f.length() - 4) + "_full.tsv"));
                while ((line = TSVReader.readLine()) != null) {
                    String[] lineSplit2 = line.split("\t");
                    //check if diann to base results in same base peptide
                    PeptideFormatter pf = new PeptideFormatter(
                            new PeptideFormatter(lineSplit2[0], lineSplit2[1], "base").prosit,
                            lineSplit2[1], "prosit");
                    if (! pf.base.equals(lineSplit2[0])) {
                        //get predictionEntry
                        PredictionEntry tmp = allPreds.get(pf.baseCharge);

                        if (tmp == null) { //valid reasons to be empty
                            if ((pf.stripped.length() > 20)) {
                                continue;
                            }
                            if (pf.stripped.contains("O") || pf.stripped.contains("U") ||
                                    pf.stripped.contains("Z") || pf.stripped.contains("B") ||
                                    pf.stripped.contains("X")) {
                                continue;
                            }
                            if (Integer.parseInt(pf.charge) > 6) {
                                continue;
                            }
                        }

                        MassCalculator mc = new MassCalculator(lineSplit2[0], lineSplit2[1]);
                        float[] newMZs = new float[tmp.mzs.length];
                        for (int i = 0; i < newMZs.length; i++) {
                            newMZs[i] = mc.calcMass(tmp.fragNums[i], flagTOion.get(tmp.flags[i]), tmp.charges[i]);
                        }

                        //add to hashmap
                        PredictionEntry newPred = new PredictionEntry();
                        newPred.setMzs(newMZs);
                        newPred.setIntensities(tmp.intensities);
                        newPred.setRT(tmp.RT);
                        newPred.setIM(0f);
                        allPreds.put(mc.fullPeptide, newPred);
                    }
                }
            } catch (IOException e) {
                e.printStackTrace();
                System.exit(1);
            }
        }
    }

    public HashMap<String, PredictionEntry> getPreds() { return allPreds; }

    public float getMaxPredRT() {
        float maxRT = 0f;
        for (PredictionEntry entry : allPreds.values()) {
            if (entry.RT > maxRT) {
                maxRT = entry.RT;
            }
        }
        return maxRT;
    }

    public void clear() {
        allPreds.clear();
    }

    public static void main(String[] args) throws IOException, InterruptedException, ExecutionException, FileParsingException {
        //MspReader m = new MspReader("C:/Users/kevin/Downloads/proteomics/newHLA/msfragger3.5rc9/myPrositLib31.msp");
        ExecutorService executorService = Executors.newFixedThreadPool(11);
        SpectralPredictionMapper spm = SpectralPredictionMapper.createSpectralPredictionMapper(
                "C:/Users/kevin/Downloads/proteomics/newHLA/msfragger3.5rc9/myPrositLib31.msp", executorService);
        executorService.shutdown();
    }
}
