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

package writers;

import allconstants.Constants;
import jakarta.xml.bind.JAXBException;
import mainsteps.MainClass;
import mainsteps.PinMzmlMatcher;
import peptideptmformatting.PeptideFormatter;
import peptideptmformatting.PeptideSkipper;
import readers.datareaders.MzmlReader;
import readers.datareaders.PinReader;
import umich.ms.fileio.exceptions.FileParsingException;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;
import java.util.concurrent.ConcurrentHashMap;
import java.util.concurrent.ExecutionException;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Future;

import static allconstants.ModelCollections.KoinaModels;
import static utils.Print.printError;
import static utils.Print.printInfo;

public class PeptideFileCreator {

    //TODO: remove psmFormat
    public static void createPeptideFile(PinMzmlMatcher pmm, String outfile, String modelFormat)
            throws IOException, InterruptedException, ExecutionException, FileParsingException { //pepXML or pin
        //diff versions based on submitting File[] or pinReader
        long startTime = System.nanoTime();
        File[] pinFiles = pmm.pinFiles;
        File[] mzmlFiles = pmm.mzmlFiles;

        //read in mzml files if needed
        for (int i = 0; i < mzmlFiles.length; i++) {
            if (pmm.mzmlReaders[i] == null) {
                pmm.mzmlReaders[i] = new MzmlReader(mzmlFiles[i].getCanonicalPath());
            }
        }

        //read in pin files
        //filter out redundant peptides
        //this step can reduce number of predictions needed to 1/3, decreasing prediction time
        ConcurrentHashMap<String, Integer> allHits = new ConcurrentHashMap<>();

        printInfo("Creating input file for " + modelFormat);
        List<Future> futureList = new ArrayList<>();
        ExecutorService executorService = MainClass.executorService;
        for (int i = 0; i < pinFiles.length; i++) {
            int finalI = i;
            futureList.add(executorService.submit(() -> {
                try {
                    File f = pinFiles[finalI];
                    File mzmlf = mzmlFiles[finalI];

                    String fileName = f.getCanonicalPath();
                    PinReader pin = new PinReader(fileName);
                    String[] hitsToAdd = new String[0];
                    switch (modelFormat) {
                        case "pDeep2":
                            hitsToAdd = pin.createPDeep2List();
                            break;
                        case "pDeep3":
                            hitsToAdd = pin.createPDeep3List();
                            break;
                        case "DeepMSPeptide": //ignores charge and mods
                            hitsToAdd = pin.createDeepMSPeptideList();
                            break;
                        case "DeepMSPeptideAll": //ignores charge and mods
                            hitsToAdd = pin.createDeepMSPeptideList();
                            break;
                        case "Diann":
                            hitsToAdd = pin.createDiannList();
                            break;
                        case "Prosit":
                            hitsToAdd = pin.createPrositList();
                            break;
                        case "PrositTMT":
                            hitsToAdd = pin.createPrositTMTList();
                            break;
                        case "createFull":
                            hitsToAdd = pin.createFull();
                            break;
                        case "alphapeptdeep":
                            hitsToAdd = pin.createAlphapeptdeepList();
                            break;
                    }
                    if (KoinaModels.contains(modelFormat)) {
                        hitsToAdd = pin.createJSON(modelFormat);
                    }
                    pin.close();

                    for (String hit : hitsToAdd) {
                        allHits.put(hit, 0);
                    }
                } catch (Exception e) {
                    e.printStackTrace();
                    System.exit(1);
                }
            }));
        }
        for (Future future : futureList) {
            future.get();
        }

        HashSet<String> hSetHits = new HashSet<>(allHits.keySet());
        printInfo(hSetHits.size() + " PSMs for prediction");

        if ((KoinaModels.contains(modelFormat)) &&
                (modelFormat.contains("AlphaPept"))) {
            //want to see what unimods were assigned
            HashSet<String> unimodCodes = new HashSet<>();
            for (String hit : hSetHits) {
                String[] unimodSplit = hit.split(",")[0].split("UNIMOD:");
                if (unimodSplit.length > 1) {
                    for (int i = 1; i < unimodSplit.length; i++) {
                        String split = unimodSplit[i];
                        unimodCodes.add(split.split("]")[0]);
                    }
                }
            }
            printInfo("AlphaPept using UniMod codes " + unimodCodes);
        }

        //write to file
        try {
            String filename = "";
            if (KoinaModels.contains(modelFormat)) {
                JSONWriter jw = new JSONWriter(modelFormat, hSetHits, true);
                if (!Constants.usedKoina) {
                    filename = jw.write(true,
                            Constants.outputDirectory + File.separator + "jsonFiles",
                            MainClass.executorService);
                    Constants.usedKoina = true;
                } else {
                    filename = jw.write(false,
                            Constants.outputDirectory + File.separator + "jsonFiles",
                            MainClass.executorService);
                }
            } else {
                //TODO: outfile name with prediction model in name
                FileWriter myWriter = new FileWriter(outfile);
                switch (modelFormat) {
                    case "pDeep2":
                        printInfo("Writing pDeep2 input file");
                        myWriter.write("peptide" + "\t" + "modification" + "\t" + "charge\n");
                        break;
                    case "pDeep3":
                        printInfo("Writing pDeep3 input file");
                        myWriter.write("raw_name" + "\t" + "scan" + "\t" + "peptide" + "\t" +
                                "modinfo" + "\t" + "charge\n");
                        break;
                    case "fineTune":
                        myWriter.write("raw_name" + "\t" + "scan" + "\t" + "peptide" + "\t" +
                                "modinfo" + "\t" + "charge" + "\t" + "RTInSeconds\n");
                        break;
                    case "DeepMSPeptide":
                        printInfo("Writing DeepMSPeptide input file");
                        //hSetHits.removeIf(hSetHit -> hSetHit.contains("B") || hSetHit.contains("X") || hSetHit.contains("Z"));
                        break; //no header
                    case "DeepMSPeptideAll":
                        printInfo("Writing DeepMSPeptideAll input file");
                        //hSetHits.removeIf(hSetHit -> hSetHit.contains("B") || hSetHit.contains("X") || hSetHit.contains("Z"));
                        break; //no header
                    case "Diann":
                        printInfo("Writing DIA-NN input file");
                        myWriter.write("peptide" + "\t" + "charge\n");
                        break;
                    case "PredFull":
                        printInfo("Writing PredFull input file");
                        myWriter.write("Peptide" + "\t" + "Charge" + "\t" + "Type" + "\t" + "NCE\n");
                        break;
                    case "Prosit":
                        printInfo("Writing Prosit input file");
                        myWriter.write("modified_sequence,collision_energy,precursor_charge\n");
                        break;
                    case "PrositTMT":
                        printInfo("Writing Prosit TMT input file");
                        myWriter.write("modified_sequence,collision_energy,precursor_charge,fragmentation\n");
                        break;
                    case "alphapeptdeep":
                        printInfo("Writing alphapeptdeep input file");
                        myWriter.write("sequence,mods,mod_sites,charge,nce,instrument,base\n");
                        break;
                }
                for (String hSetHit : hSetHits) {
                    myWriter.write(hSetHit + "\n");
                }
                myWriter.close();

                if (Constants.modelSplit && modelFormat.equals("alphapeptdeep")) {
                    PepXMLDivider pxd = new PepXMLDivider(Constants.modelSplitNum);
                    pxd.dividePinPepxml(pmm.pinFiles, Constants.spectraRTPrefix + ".csv");
                }
            }

            long endTime = System.nanoTime();
            long duration = (endTime - startTime);
            printInfo(modelFormat + " input file generation took " + duration / 1000000 +" milliseconds");
            if (KoinaModels.contains(modelFormat)) {
                printInfo("Input files in " + filename);
                Constants.JsonDirectory = filename;
            } else {
                printInfo("Input file at " + outfile);
            }
            //return fasta; //save fasta for later
        } catch (IOException | JAXBException e) {
            printError("An error occurred");
            e.printStackTrace();
            System.exit(1);
            //return null;
        }
    }

    //give it a list of peptides to write instead of all peptides in pin file
    //writes spectraRT_full.tsv that has base peptide format
    public static void createPartialFile(String filePath, String currentModel,
                                  ArrayList<PeptideFormatter> peptideFormatters) throws IOException {
        FileWriter myWriter = new FileWriter(filePath);

        for (PeptideFormatter pf : peptideFormatters) {
            if (PeptideSkipper.skipPeptide(pf.getStripped(), pf.getCharge(), currentModel)) {
                continue;
            }

            myWriter.write(pf.getBase() + "\t" + pf.getCharge() + "\n");
        }

        myWriter.close();
    }
}
