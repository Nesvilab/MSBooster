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

import jakarta.xml.bind.JAXBException;
import org.apache.commons.lang3.ArrayUtils;
import umich.ms.fileio.exceptions.FileParsingException;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.*;
import java.util.concurrent.ConcurrentHashMap;
import java.util.concurrent.ExecutionException;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Future;

import static org.apache.commons.io.FileUtils.listFiles;

public class PeptideFileCreator {

    //TODO: remove psmFormat
    public static void createPeptideFile(PinMzmlMatcher pmm, String outfile, String modelFormat)
            throws IOException, InterruptedException, ExecutionException, FileParsingException { //pepXML or pin
        //diff versions based on submitting File[] or pinReader
        long startTime = System.nanoTime();
        File[] pinFiles = pmm.pinFiles;
        File[] mzmlFiles = pmm.mzmlFiles;

        //read in mzml files if needed
        if ((Constants.useKoina && !modelFormat.equals("Diann")) || modelFormat.contains("Prosit") ||
                modelFormat.equals("PredFull") || modelFormat.equals("alphapeptdeep")) {
            for (int i = 0; i < mzmlFiles.length; i++) {
                if (Constants.NCE.equals("") && pmm.mzmlReaders[i] == null) {
                    pmm.mzmlReaders[i] = new MzmlReader(mzmlFiles[i].getCanonicalPath());
                }
            }
        }

        //read in pin files
        //filter out redundant peptides
        //this step can reduce number of predictions needed to 1/3, decreasing prediction time
        ConcurrentHashMap<String, Integer> allHits = new ConcurrentHashMap<>();

        System.out.println("Creating input file for " + modelFormat);
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
                            if (Constants.NCE.equals("")) {
                                System.out.println("If mzml file is available, will read in NCE from there");
                                Constants.FragmentationType = "HCD";
                            }
                            hitsToAdd = pin.createPrositList(mzmlf, pmm);
                            break;
                        case "PrositTMT":
                            if (Constants.FragmentationType.equals("")) {
                                System.out.println("Missing information for Prosit file generation. " +
                                        "Please provide FragmentationType (HCD or CID) in " +
                                        "parameter file via --paramsList or as arguments on command line." +
                                        "For now, setting as HCD.");
                                Constants.FragmentationType = "HCD";
                            }
                            if (Constants.NCE.equals("")) {
                                System.out.println("If mzml file is available, will read in NCE from there");
                            }
                            hitsToAdd = pin.createPrositTMTList(mzmlf, pmm);
                            break;
                        case "createFull":
                            hitsToAdd = pin.createFull();
                            break;
                        case "PredFull":
                            if (Constants.FragmentationType.equals("")) {
                                System.out.println("Missing fragmentation type for PredFull file generation. " +
                                        "You can provide FragmentationType (HCD or ETD) in the " +
                                        "parameter file via --paramsList or as arguments on command line. " +
                                        "For now, setting as HCD");
                                Constants.FragmentationType = "HCD";
                            }
                            if (Constants.NCE.equals("")) {
                                System.out.println("If mzml file is available, will read in NCE from there");
                            }
                            hitsToAdd = pin.createPredFullList(mzmlf, pmm);
                            break;
                        case "alphapeptdeep":
                            if (Constants.instrument.equals("")) {
                                System.out.println("Missing instrument for alphapeptdeep file generation. " +
                                        "You can provide an instrument in the " +
                                        "parameter file via --paramsList or as arguments on command line. " +
                                        "For now, automatically setting it based on mzml metadata. " +
                                        "The following instruments are allowed, along with the model mode that " +
                                        "alphapeptdeep converts them to (user input: alphapeptdeep mode):\n" +
                                        "    Lumos: Lumos\n" +
                                        "    QE: QE\n" +
                                        "    timsTOF: timsTOF\n" +
                                        "    SciexTOF: SciexTOF\n" +
                                        "    Fusion: Lumos\n" +
                                        "    Eclipse: Lumos\n" +
                                        "    Velos: Lumos\n" +
                                        "    Elite: Lumos\n" +
                                        "    OrbitrapTribrid: Lumos\n" +
                                        "    ThermoTribrid: Lumos\n" +
                                        "    QE+: QE\n" +
                                        "    QEHF: QE\n" +
                                        "    QEHFX: QE\n" +
                                        "    Exploris: QE\n" +
                                        "    Exploris480: QE" +
                                        "    ThermoTOF: ThermoTOF" +
                                        "    Astral: ThermoTOF");
                                Constants.instrument = "Lumos";
                            }
                            if (Constants.NCE.equals("")) {
                                System.out.println("If mzml file is available, will read in NCE from there");
                                Constants.FragmentationType = "HCD";
                            }
                            hitsToAdd = pin.createAlphapeptdeepList(mzmlf, pmm);
                            break;
                    }
                    if (Constants.useKoina && !modelFormat.equals("Diann")) {
                        hitsToAdd = pin.createJSON(mzmlf, pmm, modelFormat);
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
        System.out.println(hSetHits.size() + " PSMs for prediction");

        //write to file
        try {
            String filename = "";
            if (Constants.useKoina && !modelFormat.equals("Diann")) {
                JSONWriter jw = new JSONWriter(modelFormat, hSetHits);
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
                        System.out.println("Writing pDeep2 input file");
                        myWriter.write("peptide" + "\t" + "modification" + "\t" + "charge\n");
                        break;
                    case "pDeep3":
                        System.out.println("Writing pDeep3 input file");
                        myWriter.write("raw_name" + "\t" + "scan" + "\t" + "peptide" + "\t" +
                                "modinfo" + "\t" + "charge\n");
                        break;
                    case "fineTune":
                        myWriter.write("raw_name" + "\t" + "scan" + "\t" + "peptide" + "\t" +
                                "modinfo" + "\t" + "charge" + "\t" + "RTInSeconds\n");
                        break;
                    case "DeepMSPeptide":
                        System.out.println("Writing DeepMSPeptide input file");
                        //hSetHits.removeIf(hSetHit -> hSetHit.contains("B") || hSetHit.contains("X") || hSetHit.contains("Z"));
                        break; //no header
                    case "DeepMSPeptideAll":
                        System.out.println("Writing DeepMSPeptideAll input file");
                        //hSetHits.removeIf(hSetHit -> hSetHit.contains("B") || hSetHit.contains("X") || hSetHit.contains("Z"));
                        break; //no header
                    case "Diann":
                        System.out.println("Writing DIA-NN input file");
                        myWriter.write("peptide" + "\t" + "charge\n");
                        break;
                    case "PredFull":
                        System.out.println("Writing PredFull input file");
                        myWriter.write("Peptide" + "\t" + "Charge" + "\t" + "Type" + "\t" + "NCE\n");
                        break;
                    case "Prosit":
                        System.out.println("Writing Prosit input file");
                        myWriter.write("modified_sequence,collision_energy,precursor_charge\n");
                        break;
                    case "PrositTMT":
                        System.out.println("Writing Prosit TMT input file");
                        myWriter.write("modified_sequence,collision_energy,precursor_charge,fragmentation\n");
                        break;
                    case "alphapeptdeep":
                        System.out.println("Writing alphapeptdeep input file");
                        myWriter.write("sequence,mods,mod_sites,charge,nce,instrument,base\n");
                        break;
                }
                for (String hSetHit : hSetHits) {
                    myWriter.write(hSetHit + "\n");
                }
                myWriter.close();

                if (Constants.modelSplit && modelFormat.equals("alphapeptdeep")) {
                    PepXMLDivider pxd = new PepXMLDivider(Constants.modelSplitNum);
                    pxd.dividePinPepxml(pmm.pinFiles, Constants.spectraRTPredInput);
                }
            }

            long endTime = System.nanoTime();
            long duration = (endTime - startTime);
            System.out.println(modelFormat + " input file generation took " + duration / 1000000 +" milliseconds");
            if (Constants.useKoina && !modelFormat.equals("Diann")) {
                System.out.println("Input files in " + filename);
                Constants.JsonDirectory = filename;
            } else {
                System.out.println("Input file at " + outfile);
            }
            System.out.println();
            //return fasta; //save fasta for later
        } catch (IOException | JAXBException e) {
            System.out.println("An error occurred");
            e.printStackTrace();
            System.exit(1);
            //return null;
        }
    }
}
