package External;

import Features.*;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.HashMap;
import java.util.HashSet;
import java.util.LinkedList;
import java.util.concurrent.ConcurrentHashMap;
import java.util.concurrent.ScheduledThreadPoolExecutor;

public class KoinaMethods {
    //these fields are shared regardless of which model is called
    public HashSet<String> peptideSet = new HashSet<>();
    public HashMap<String, LinkedList<Integer>> scanNums = new HashMap<>();
    public HashMap<String, LinkedList<String>> peptides = new HashMap<>();
    public KoinaMethods() {}

    public void getTopPeptides(PinMzmlMatcher pmMatcher) throws IOException {
        //need to collect top 1000 peptides for calibration
        //approximate by doing a subset per pin
        int numTopPSMs = (int) Math.ceil((float) Constants.numPSMsToCalibrate /
                (float) pmMatcher.pinFiles.length);

        for (int j = 0; j < pmMatcher.pinFiles.length; j++) {
            File pinFile = pmMatcher.pinFiles[j];
            PinReader pinReader = new PinReader(pinFile.getAbsolutePath());
            LinkedList[] topPSMs = pinReader.getTopPSMs(numTopPSMs);
            peptideSet.addAll(topPSMs[0]);
            scanNums.put(pmMatcher.mzmlFiles[j].getName(), topPSMs[1]);
            peptides.put(pmMatcher.mzmlFiles[j].getName(), topPSMs[0]);
            if (Constants.instrument.isEmpty()) {
                pinReader.attachMzml(pmMatcher.mzmlReaders[j]);
                Constants.instrument = pinReader.getInstrument();
            }
        }
    }

    public HashSet<String> writeFullPeptideFile(String filePath, String currentModel) throws IOException {
        FileWriter myWriter = new FileWriter(filePath);
        HashSet<String> allHits = new HashSet<>();
        for (String peptide : peptideSet) {
            String[] peptFormats = peptide.split(",");

            String stripped = peptFormats[2];

            if ((currentModel.contains("Prosit") || currentModel.contains("ms2pip") || currentModel.contains("Deeplc"))
                    && stripped.contains("U")) { // no peptides with U
                continue;
            }
            if (currentModel.contains("ms2pip") && stripped.length() > 30) { //peptide has length limit
                continue;
            }
            if (currentModel.contains("Prosit") && currentModel.contains("TMT") && stripped.length() > 30) {
                continue;
            }

            String pep = peptFormats[1].replace("UniMod", "UNIMOD"); //need diann here
            if (pep.contains("[TMT]")) {
                pep = pep.replace("[TMT]", "[UNIMOD:737]");
            }

            if (pep.startsWith("[")) { //this is the case for all n term mods //TODO deal with c term mods
                int splitpoint = pep.indexOf("]");
                if (currentModel.contains("Prosit")) {
                    if (currentModel.contains("TMT") && pep.startsWith("[UNIMOD:737]")) {
                        pep = pep.substring(0, splitpoint + 1) + "-" + pep.substring(splitpoint + 1);
                    } else {
                        pep = pep.substring(splitpoint + 1);
                    }
                } else {
                    pep = pep.substring(0, splitpoint + 1) + "-" + pep.substring(splitpoint + 1);
                }
            }
            if (currentModel.contains("Prosit") && currentModel.contains("TMT")) {
                pep = pep.replace("S[UNIMOD:737]", "S");
                if (!pep.startsWith("[")) {
                    pep = "[UNIMOD:737]-" + pep;
                }
            }
            String[] baseCharge = peptFormats[0].split("\\|");
            allHits.add(pep + "," + baseCharge[1]);

            //need to generate full.tsv
            myWriter.write(baseCharge[0] + "\t" + baseCharge[1] + "\n");
        }
        myWriter.close();

        return allHits;
    }

    public ConcurrentHashMap<String, PredictionEntry> getKoinaPredictions(
            HashSet<String> allHits, String model, int NCE, String folder, String fulltsv) {
        ScheduledThreadPoolExecutor executorService = new ScheduledThreadPoolExecutor(1);

        if (Constants.FragmentationType.isEmpty()) {
            System.out.println("Setting fragmentation type to HCD. " +
                    "You can specify this with '--FragmentationType' via the command line " +
                    "or 'FragmentationType=' in the param file.");
            Constants.FragmentationType = "HCD";
        }

        HashSet<String> hits = new HashSet<>();
        for (String s : allHits) {
            hits.add(s + "," + NCE + "," + Constants.instrument + "," + Constants.FragmentationType);
        }

        JSONWriter jw = new JSONWriter(model, hits);

        String jsonFolder = "";
        try {
            jsonFolder = jw.write(true, folder, executorService);
        } catch (Exception e) {
            e.printStackTrace();
            System.exit(1);
        }

        //send predictions to Koina
        KoinaLibReader klr = new KoinaLibReader();
        KoinaModelCaller kmc = new KoinaModelCaller();
        kmc.callModel(model, klr, jsonFolder, executorService, false, false);
        executorService.shutdown();
        try {
            kmc.assignMissingPeptidePredictions(klr, fulltsv);
        } catch (IOException e) {
            throw new RuntimeException(e);
        }
        return klr.allPreds;
    }
}
