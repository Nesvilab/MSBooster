import org.apache.commons.lang.ArrayUtils;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.Collection;
import java.util.Collections;
import java.util.HashSet;

import static org.apache.commons.io.FileUtils.listFiles;

public class FileCreator {

    public static void createPDeepFile(String[] allHits, String outfile) {
        //remove duplicates from allHits
        //can reduce number of hits to a third
        HashSet<String> hSetHits = new HashSet<>();
        Collections.addAll(hSetHits, allHits);
        System.out.println(allHits.length);
        System.out.println(hSetHits.size());

        //write to file
        try {
            FileWriter myWriter = new FileWriter(outfile);
            myWriter.write("peptide" + "\t" + "modification" + "\t" + "charge\n");

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
        //get all files to be analyzed
        //FileCreator x = new FileCreator("C:/Users/kevin/OneDriveUmich/proteomics/pepxml/rank4/");
        Collection<File> x = listFiles(new File("C:/Users/kevin/OneDriveUmich/proteomics/pepxml/"),
                new String[]{"pepXML"}, true);

        //read in pepXML files
        String[] allHits = new String[0];
        for (File f : x) {
            String fileName = f.getCanonicalPath();
            System.out.println(fileName);
            pepXMLReader xmlReader = new pepXMLReader(fileName);
            //xmlReader.createPDeepList();
            xmlReader.createPDeepListNoMods();
            allHits = (String[]) ArrayUtils.addAll(allHits, xmlReader.allHitsPDeep);
        }

        //create file for pDeep2 prediction
        createPDeepFile(allHits, "C:/Users/kevin/Downloads/proteomics/peptides_for_pDeep_noMods.tsv"); //actually a tsv file

        //got predictions for all peptides in pepXML
        //python predict.py
        //{'nce': 0.27, 'instrument': 'QE', 'input': 'narrow1rank1.txt', 'output': 'narrow1rank1_pred.txt'}
    }
}
