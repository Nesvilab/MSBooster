package Features;

import java.io.File;
import java.io.FilenameFilter;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;

//deprecated, as we no longer have to rearrange them
public class organization {
    public static void sortRankPepxml(String directory, int numRanks) throws IOException {

        File f = new File(directory);

        //create directories
        for (int rank = 1; rank < numRanks + 1; rank++) {
            File newDirectory = new File(f.getAbsolutePath() + File.separator + "rank" + rank);
            if (!newDirectory.exists()) {
                newDirectory.mkdirs();
            }
        }

        // This filter will only include files ending with .py
        FilenameFilter filter = new FilenameFilter() {
            @Override
            public boolean accept(File f, String name) {
                return name.endsWith(".pepXML");
            }
        };

        // This is how to apply the filter
        String[] pathnames = f.list(filter);
        for (String name : pathnames) {
            String[] parts = name.split("_");
            String direc = parts[parts.length - 1].split("\\.")[0];
            //String newName = String.join("_", Arrays.copyOfRange(parts, 0, parts.length - 1)) + ".pepXML";
            Files.move(Paths.get(f.getAbsolutePath() + File.separator + name),
                    Paths.get(f.getAbsolutePath() + File.separator + direc + File.separator + name));
        }

    }

    public static void main(String[] args) throws IOException {
        //sortRankPepxml("C:/Users/kevin/Downloads/proteomics/cptac/2021-2-21/", 16);
    }
}
