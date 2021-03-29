package Features;

import org.apache.commons.lang.ArrayUtils;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.Collection;
import java.util.Collections;
import java.util.HashSet;

import static org.apache.commons.io.FileUtils.listFiles;

public class peptideFileCreator {
    //private static String[] acceptableFormats = new String[] {"pDeep2", "pDeep3", "prosit"};

    public static HashSet<String> getUniqueHits(String[] allHits) {
        //remove duplicates from allHits
        //can reduce number of hits to a third
        HashSet<String> hSetHits = new HashSet<>();
        Collections.addAll(hSetHits, allHits);
        System.out.println("before filtering: " + allHits.length +
                ", after filtering: " + hSetHits.size());
        return hSetHits;
    }

    public static void createPeptideFile(String infile, String outfile, String modelFormat) throws IOException {
        Collection<File> x = listFiles(new File(infile), new String[]{"pepXML"}, true);

        //read in pepXML files
        String[] allHits = new String[0];
        for (File f : x) {
            String fileName = f.getCanonicalPath();
            System.out.println(fileName);
            pepXMLReader xmlReader = new pepXMLReader(fileName);
            String[] hitsToAdd = new String[0];
            switch (modelFormat) {
                case "pDeep3":
                    hitsToAdd = xmlReader.createPDeep3List();
                    break;
                case "DeepMSPeptide": //ignores charge and mods
                    hitsToAdd = xmlReader.createDeepMSPeptideList();
                    break;
                case "Diann":
                    hitsToAdd = xmlReader.createDiannList();
                    break;
            }

            allHits = (String[]) ArrayUtils.addAll(allHits, hitsToAdd);
        }

        //filter out redundant peptides
        //this step can reduce number of predictions needed to 1/3, decreasing prediction time
        HashSet<String> hSetHits = getUniqueHits(allHits);

        //write to file
        try {
            FileWriter myWriter = new FileWriter(outfile);
            switch (modelFormat) {
                case "prosit":
                    System.out.println("writing prosit");
                    myWriter.write("modified_sequence" + "," + "collision_energy" + "," + "precursor_charge\n");
                    //instances of null being added because no n-acetyl
                    hSetHits.remove(null);
                    break;
                case "pDeep2":
                    System.out.println("writing pDeep2");
                    myWriter.write("peptide" + "\t" + "modification" + "\t" + "charge\n");
                    break;
                case "pDeep3":
                    System.out.println("writing pDeep3");
                    myWriter.write("raw_name" + "\t" + "scan" + "\t" + "peptide" + "\t" +
                            "modinfo" + "\t" + "charge\n");
                    break;
                case "fineTune":
                    myWriter.write("raw_name" + "\t" + "scan" + "\t" + "peptide" + "\t" +
                            "modinfo" + "\t" + "charge" + "\t" + "RTInSeconds\n");
                    break;
                case "DeepMSPeptide":
                    break; //no header
                case "Diann":
                    System.out.println("writing Diann");
                    myWriter.write("peptide" + "\t" + "charge\n");
                    break;
            }

            for (String hSetHit : hSetHits) {
                myWriter.write(hSetHit + "\n");
            }

            myWriter.close();
            System.out.println("Successfully wrote to the file.");
        } catch (IOException e) {
            System.out.println("An error occurred.");
            e.printStackTrace();
        }
    }

    public static void main(String[] args) throws IOException {
        createPeptideFile("C:/Users/kevin/Downloads/proteomics/narrow/",
                "C:/Users/kevin/OneDriveUmich/proteomics/preds/narrowPDeep3.tsv",
                "pDeep3");
    }
}
