package Features;

import com.univocity.parsers.tsv.TsvWriter;
import com.univocity.parsers.tsv.TsvWriterSettings;
import org.apache.commons.lang3.ArrayUtils;
import umich.ms.fileio.exceptions.FileParsingException;

import java.io.File;
import java.io.IOException;
import java.util.*;
import java.util.stream.IntStream;

import static org.apache.commons.io.FileUtils.listFiles;

public class percolatorFormatter {

    public static String percolatorPepFormat(String[] columns, int pepIdx, int specIDidx) {
        String pep = columns[pepIdx];
        pep = pep.substring(2, pep.length() - 2);
        StringBuilder mods = new StringBuilder();

        //n term acetylation
        if (pep.charAt(0) == 'n') {
            pep = pep.replace("n[42.0106]", "");
            mods.append("0,Acetyl[AnyN-term];");
        }

        pep = pep.replace("C[57.0215]", "c");
        pep = pep.replace("M[15.9949]", "m");
        TreeMap<Integer, String> modsMap = new TreeMap<>();

        //carbamidometylation
        while (pep.contains("c")) {
            int pos = pep.indexOf("c") + 1;
            pep = pep.replaceFirst("c", "C");
            modsMap.put(pos, ",Carbamidomethyl[C];");
        }

        //methionine oxidation
        while (pep.contains("m")) {
            int pos = pep.indexOf("m") + 1;
            pep = pep.replaceFirst("m", "M");
            modsMap.put(pos, ",Oxidation[M];");
        }

        for (Map.Entry<Integer, String> entry : modsMap.entrySet()) {
            mods.append(entry.getKey()).append(entry.getValue());
        }

        //charge
        String[] charges = columns[specIDidx].split("_");
        String chargeStr = charges[charges.length - 2];
        int charge = Integer.parseInt(chargeStr.substring(chargeStr.length() - 1));

        return pep + "|" + mods + "|" + charge;
    }

    //set mgf or detectFile as null if not applicable
    //baseNames is the part before mzml or pin extensions
    public static void editPin(String pinDirectory, String mzmlDirectory, String mgf, String detectFile,
                               String[] features, String outfile)
            throws IOException {
        List<String> featuresList = Arrays.asList(features);
        //remove features from array for multiple protein formatting
        if (featuresList.contains("detectability")) {
            int idx = ArrayUtils.indexOf(features, "detectability");
            features = ArrayUtils.remove(features, idx);
        }

        //booleans for future determination of what to do
        boolean needsMGF = false;

        //get names of mzml files
        //check if file or directory
        String[] mzmlFiles;
        if (mzmlDirectory.substring(mzmlDirectory.length() - 4).toLowerCase().equals("mzml")) {
            File f = new File(mzmlDirectory);
            mzmlFiles = new String[]{f.getName()};
            mzmlDirectory = f.getAbsoluteFile().getParent();
        } else {
            Collection<File> mzmlFilesCollection = listFiles(new File(mzmlDirectory), new String[]{"mzML"}, false);
            mzmlFiles = new String[mzmlFilesCollection.size()];
            int FileIdx = 0;
            for (File f : mzmlFilesCollection) {
                mzmlFiles[FileIdx] = f.getName();
                FileIdx++;
            }
            Arrays.sort(mzmlFiles);
        }

        //check that corresponding pin files exist
        String[] pinFiles;
        if (pinDirectory.substring(pinDirectory.length() - 3).toLowerCase().equals("pin")) {
            File f = new File(pinDirectory);
            pinFiles = new String[]{f.getName()};
            pinDirectory = f.getAbsoluteFile().getParent();
        } else {
            Collection<File> pinFilesCollection = listFiles(new File(pinDirectory), new String[]{"pin"}, false);
            HashSet<String> pinFilesSet = new HashSet<>();
            for (File f : pinFilesCollection) {
                pinFilesSet.add(f.getName());
            }

            pinFiles = new String[mzmlFiles.length];
            for (int i = 0; i < mzmlFiles.length; i++) {
                pinFiles[i] = mzmlFiles[i].substring(0, mzmlFiles[i].length() - 4) + "pin";
                if (!pinFilesSet.contains(pinFiles[i])) {
                    throw new AssertionError("mzML file must have corresponding pin file. " +
                            pinFiles[i] + " does not exist");
                }
            }
        }

        //load predicted spectra
        SpectralPredictionMapper predictedSpectra = null;

        //Special preparations dependent on features we require
        //only time this isn't needed is detect
        //could consider an mgf constant
        if (featuresList.contains("deltaRTlinear") || featuresList.contains("deltaRTbins") ||
                featuresList.contains("RTzscore") || featuresList.contains("RTprobability") ||
                featuresList.contains("RTprobabilityUnifPrior") || featuresList.contains("brayCurtis") ||
                featuresList.contains("cosineSimilarity") || featuresList.contains("spectralContrastAngle") ||
                featuresList.contains("euclideanDistance") || featuresList.contains("pearsonCorr") ||
                featuresList.contains("dotProduct") || featuresList.contains("deltaRTLOESS")) {
            System.out.println("Loading predicted spectra");
            predictedSpectra = SpectralPredictionMapper.createSpectralPredictionMapper(mgf);
            needsMGF = true;
        }

        //detectability
        detectMap dm = null;
        if (detectFile != null) {
            dm = new detectMap(detectFile);
        }

        try {
            //////////////////////////////iterate through pin and mzml files//////////////////////////////////////////
            for (int i = 0; i < pinFiles.length; i++) {
                long startTime = System.nanoTime();
                String newOutfile = outfile + pinFiles[i];
                TsvWriter writer = new TsvWriter(new File(newOutfile), new TsvWriterSettings());
                //load mzml file
                System.out.println("Loading " + mzmlFiles[i]);
                mzMLReader mzml = new mzMLReader(mzmlDirectory + File.separator + mzmlFiles[i]);

                //load pin file, which already includes all ranks
                System.out.println("Loading " + pinFiles[i]);
                pinReader pin = new pinReader(pinDirectory + File.separator + pinFiles[i]);

                //add header to written tsv
                ArrayList<String> newHeader = new ArrayList<>();
                newHeader.addAll(Arrays.asList(pin.header));
                newHeader.addAll(pin.pepIdx, featuresList); //add features before Peptide
                newHeader.remove("detectability");
                writer.writeHeaders(newHeader);

                //Special preparations dependent on features we require
                if (needsMGF) {
                    System.out.println("Loading PSMs onto mzml object");
                    //TODO: can we detect before this how many ranks there are?
                    mzml.setPinEntries(pin, predictedSpectra);
                    System.out.println("Done loading PSMs onto mzml object");
                }
                if (featuresList.contains("deltaRTLOESS")) {
                    System.out.println("Generating LOESS regression");
                    mzml.setLOESS(predictedSpectra, Constants.RTregressionSize, Constants.bandwidth, Constants.robustIters);
                }
                if (featuresList.contains("deltaRTlinear")) {
                    System.out.println("Calculating delta RT linear");
                    if (mzml.expAndPredRTs != null) { mzml.setBetas();
                    } else { mzml.setBetas(predictedSpectra, Constants.RTregressionSize);
                    }
                    //System.out.println(Arrays.toString(mzml.getBetas())); //print beta 0 and 1
                    mzml.normalizeRTs();
                }
                if (featuresList.contains("deltaRTbins") || featuresList.contains("RTzscore") ||
                        featuresList.contains("RTprobability") || featuresList.contains("RTprobabilityUnifPrior")) {
                    System.out.println("Generating RT bins");
                    mzml.setRTbins(predictedSpectra);
                    mzml.propagateRTBinStats(); //may need to make more specific
                }
                int unifPriorSize = 0; //only used if using uniform prior
                float unifProb = 0;
                if (featuresList.contains("RTprobabilityUnifPrior")) {
                    //decide how many uniform points to add
                    int[] binSizes = new int[mzml.RTbins.length];
                    for (int bin = 0; bin < mzml.RTbins.length; bin++) {
                        binSizes[bin] = mzml.RTbins[bin].size();
                    }
                    Arrays.sort(binSizes);
                    int cutoff = (int) Math.floor(((double) mzml.RTbins.length / 100.0) * Constants.uniformPriorPercentile);
                    unifPriorSize = binSizes[cutoff];

                    //also need uniform probability with right bound of max predicted RT
                    unifProb = 1.0f / predictedSpectra.getMaxPredRT();
                }
                if (featuresList.contains("RTprobability") || featuresList.contains("RTprobabilityUnifPrior")) {
                    System.out.println("Generating empirical densities");
                    mzml.setKernelDensities();
                }

                //for storing detects and whether peptides are present
                HashMap<String, float[]> detects = new HashMap<String, float[]>();
                HashMap<String, float[]> presence = new HashMap<String, float[]>();
                if (featuresList.contains("detectFractionGreater")) {
                    //get all peptides present in pin
                    HashSet<String> allPeps = pin.getAllPep();

                    //load fasta
                    System.out.println("Loading fasta");
                    FastaReader fasta = new FastaReader(Constants.fasta);

                    System.out.println("Loading detectabilities for unique peptides from each protein");
                    for (Map.Entry<String, ArrayList<String>> e : fasta.protToPep.entrySet()) {
                        float[] protDetects = new float[e.getValue().size()]; //for storing initial detect order

                        //store detect
                        for (int pep = 0; pep < e.getValue().size(); pep++) {
                            protDetects[pep] = dm.getDetectability(e.getValue().get(pep));
                        }

                        //save sorted detect in detects hashmap
                        int[] sortedIndices = IntStream.range(0, protDetects.length)
                                .boxed().sorted((k, j) -> Float.compare(protDetects[k], protDetects[j]))
                                .mapToInt(ele -> ele).toArray();
                        float[] sortedDetect = new float[protDetects.length];
                        for (int j = 0; j < sortedDetect.length; j++) {
                            sortedDetect[j] = protDetects[sortedIndices[j]];
                        }
                        detects.put(e.getKey(), sortedDetect);

                        //check which peptides present
                        float[] protPresence = new float[sortedDetect.length];
                        for (int j = 0; j < sortedIndices.length; j++) {
                            if (allPeps.contains(e.getValue().get(j))) {
                                protPresence[j] = 1f;
                            } else {
                                protPresence[j] = 0f;
                            }
                        }
                        presence.put(e.getKey(), protPresence);
                    }
                }

                System.out.println("Getting predictions for each row");
                while (pin.next()) {
                    //peptide name
                    String pep = pin.getPep();

                    //trying filtering out low detectability
                    if (featuresList.contains("detectability")) {
                        if (dm.getDetectability(pep) < Constants.detectThreshold) {
                            continue;
                        }
                    }

                    //get entry
                    peptideObj pepObj = null;
                    if (needsMGF) {
                        pepObj = mzml.scanNumberObjects.get(pin.getScanNum()).getPeptideObject(pep);
                    }

                    //write everything we already have, not including extra protein columns
                    String[] row = pin.getRow();
                    for (int j = 0; j < pin.header.length; j++) {
                        writer.addValue(pin.header[j], row[j]);
                    }
                    //add extra protein columns
                    if (pin.getRow().length > pin.header.length) {
                        for (int j = pin.header.length; j < pin.getRow().length; j++) {
                            writer.addValue(j + features.length, row[j]);
                        }
                    }

                    //switch case
                    for (String feature : features) {
                        switch (feature) {
                            case "detectFractionGreater":
                                float d = dm.getDetectability(pep);
                                //for each protein, get the position of pep's detect and see how many peptides with greater detect are present
                                //take max (proxy for protein that actually generated peptide)
                                String[] r = pin.getRow();
                                String[] prots = Arrays.copyOfRange(r, pin.pepIdx + 1, r.length);
                                float maxFraction = 0f;
                                for (String prot : prots) { //if more than one, this peptide is shared among proteins
                                    String protAbr = prot.split("\\|")[1];
                                    float[] arr = detects.get(protAbr);
                                    if (arr == null) { //no peptides qualify from this protein
                                        continue;
                                    }
                                    int idx = Arrays.binarySearch(arr, d);
                                    if (idx < 0) { //not found
                                        idx = (-1 * idx) - 1;
                                    } else {
                                        idx += 1; //don't want to include itself in calculation
                                    }
                                    float[] presenceArr = Arrays.copyOfRange(presence.get(protAbr), idx, arr.length);
                                    float total = 0f;
                                    for (float j : presenceArr) {
                                        total += j;
                                    }
                                    float fraction = (total + Constants.detectFractionGreaterNumerator) /
                                            (presenceArr.length + Constants.detectFractionGreaterDenominator); //customizable prior
                                    if (fraction > maxFraction) {
                                        maxFraction = fraction;
                                    }
                                }
                                writer.addValue("detectFractionGreater", maxFraction);
                                break;
                            case "deltaRTlinear":
                                writer.addValue("deltaRTlinear", pepObj.deltaRT);
                                break;
                            case "deltaRTbins":
                                writer.addValue("deltaRTbins", pepObj.deltaRTbin);
                                break;
                            case "deltaRTLOESS":
                                writer.addValue("deltaRTLOESS",
                                        Math.abs(mzml.predictLOESS(pepObj.scanNumObj.RT) - pepObj.RT));
                                break;
                            case "RTzscore":
                                writer.addValue("RTzscore", pepObj.RTzscore);
                                break;
                            case "RTprobability":
                                writer.addValue("RTprobability", pepObj.RTprob);
                                break;
                            case "RTprobabilityUnifPrior":
                                float prob = RTFunctions.RTprobabilityWithUniformPrior(unifPriorSize, unifProb,
                                        pepObj.scanNumObj.RTbinSize, (float) pepObj.RTprob);
                                writer.addValue("RTprobabilityUnifPrior", prob);
                                break;
                            case "brayCurtis":
                                writer.addValue("brayCurtis", pepObj.spectralSimObj.brayCurtis());
                                break;
                            case "cosineSimilarity":
                                writer.addValue("cosineSimilarity", pepObj.spectralSimObj.cosineSimilarity());
                                break;
                            case "spectralContrastAngle":
                                writer.addValue("spectralContrastAngle", pepObj.spectralSimObj.spectralContrastAngle());
                                break;
                            case "euclideanDistance":
                                writer.addValue("euclideanDistance", pepObj.spectralSimObj.euclideanDistance());
                                break;
                            case "pearsonCorr":
                                writer.addValue("pearsonCorr", pepObj.spectralSimObj.pearsonCorr());
                                break;
                            case "dotProduct":
                                writer.addValue("dotProduct", pepObj.spectralSimObj.dotProduct());
                                break;
                        }
                    }
                    //flush values to output
                    writer.writeValuesToRow();
                }
                pin.close();
                long endTime = System.nanoTime();

                long duration = (endTime - startTime);
                System.out.println("Pin editing took " + duration / 1000000000 +" seconds");
                writer.close();
                System.out.println("Edited pin file at " + newOutfile);
            }
        } catch (IOException | FileParsingException e) {
            e.printStackTrace();
        } catch (Exception e) {
            e.printStackTrace();
        }
    }

    public static void main(String[] args) throws IOException {
        //CHANGE PPM TO 10 if wide, narrow

        //CHANGE PPM TO 20 if cptac
        editPin("C:/Users/kevin/Downloads/proteomics/wide/",
                "C:/Users/kevin/OneDriveUmich/proteomics/mzml/wideWindow/",
                "C:/Users/kevin/Downloads/proteomics/wide/spectraRT.predicted.bin",
                "C:/Users/kevin/OneDriveUmich/proteomics/preds/detectwideAll_Predictions.txt",
                ("detectability,brayCurtis").split(","),
                "C:/Users/kevin/Downloads/proteomics/wide/edited_");
    }
}
