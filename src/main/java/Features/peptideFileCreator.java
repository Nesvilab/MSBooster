package Features;

import org.apache.commons.lang3.ArrayUtils;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.*;

import static org.apache.commons.io.FileUtils.listFiles;

public class peptideFileCreator {
    //private static String[] acceptableFormats = new String[] {"pDeep2", "pDeep3", "prosit"};

    public static HashSet<String> getUniqueHits(String[] allHits) {
        //remove duplicates from allHits
        //can reduce number of hits to a third
        HashSet<String> hSetHits = new HashSet<>();
        Collections.addAll(hSetHits, allHits);
        System.out.println(hSetHits.size() + " unique peptides from " + allHits.length + " PSMs");
        return hSetHits;
    }

    //infile is pepXML file locations
    public static void createPeptideFile(String infile, String outfile, String modelFormat, String psmFormat)
            throws IOException { //pepXML or pin
        //long startTime = System.nanoTime();

        //file or directory
        Collection<File> x = new ArrayList<File>();
        File newFile = new File(infile);
        if (newFile.isFile()) {
            x.add(newFile);
        } else { //directory
            x = listFiles(new File(infile), new String[]{psmFormat}, false);
        }

        File[] infileArray = new File[x.size()];
        int i = 0;
        for (File f : x) {
            infileArray[i] = f;
            i++;
        }

        createPeptideFile(infileArray, outfile, modelFormat);
    }

    //TODO: remove psmFormat
    public static void createPeptideFile(File[] x, String outfile, String modelFormat)
            throws IOException { //pepXML or pin
        //diff versions based on submitting File[] or pinReader
        long startTime = System.nanoTime();
        //read in pin files
        String[] allHits = new String[0];

        for (File f : x) {
            String fileName = f.getCanonicalPath();
            pinReader pin = new pinReader(fileName);
            String[] hitsToAdd = new String[0];
            switch (modelFormat) {
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
                    if (Constants.FragmentationType == null || Constants.NCE == null) {
                        System.out.println("Missing information for Prosit file generation. " +
                                "Please provide FragmentationType (HCD or CID) and NCE (normalized collision energy) in " +
                                "parameter file via --paramsList or as arguments on command line.");
                        System.exit(-1);
                    }
                    hitsToAdd = pin.createPrositList();
                    break;
                case "createFull":
                    hitsToAdd = pin.createFull();
                    break;
                case "PredFull":
                    if (Constants.FragmentationType == null || Constants.NCE == null) {
                        System.out.println("Missing information for PredFull file generation. " +
                                "Please provide FragmentationType (HCD or ETD) and NCE (normalized collision energy) in " +
                                "parameter file via --paramsList or as arguments on command line.");
                        System.exit(-1);
                    }
                    hitsToAdd = pin.createPredFullList();
                    break;
            }
            pin.close();

            allHits = ArrayUtils.addAll(allHits, hitsToAdd);
        }

        //filter out redundant peptides
        //this step can reduce number of predictions needed to 1/3, decreasing prediction time
        HashSet<String> hSetHits = getUniqueHits(allHits);

        //write to file
        try {
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
                    myWriter.write("Peptide" + "\t" + "Charge" + "\t" + "Type" + "\t" + "NCE\n");
                    break;
                case "Prosit":
                    myWriter.write("modified_sequence,collision_energy,precursor_charge\n");
                    break;
            }
            for (String hSetHit : hSetHits) {
                myWriter.write(hSetHit + "\n");
            }

            long endTime = System.nanoTime();
            long duration = (endTime - startTime);
            System.out.println(modelFormat + " input file generation took " + duration / 1000000 +" milliseconds");
            myWriter.close();
            System.out.println("Input file at  " + outfile);
            //return fasta; //save fasta for later
        } catch (IOException e) {
            System.out.println("An error occurred");
            e.printStackTrace();
            System.exit(1);
            //return null;
        }
    }

    public static void main(String[] args) throws IOException {
        createPeptideFile("C:/Users/kevin/Downloads/20190627_QX0_AnBr_SA_BPP_DDA_M01_01.pin",
                "C:/Users/kevin/Downloads/20190627_QX0_AnBr_SA_BPP_DDA_M01_01_length63.tsv",
                "Diann", "pin");
    }
}
