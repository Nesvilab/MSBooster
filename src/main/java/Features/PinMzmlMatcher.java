package Features;

import java.io.File;
import java.io.IOException;
import java.util.Arrays;
import java.util.Collection;
import java.util.HashSet;

import static org.apache.commons.io.FileUtils.listFiles;

public class PinMzmlMatcher {
    public File[] mzmlFiles;
    public File[] pinFiles;
    public File mzmlDirectory;
    public File pinDirectory;
    public PinMzmlMatcher(String mzmlDirectory, String pinDirectory) throws IOException {

        //get names of mzml files
        //check if file or directory
        if (mzmlDirectory.substring(mzmlDirectory.length() - 4).toLowerCase().equals("mzml")) {
            mzmlFiles = new File[]{new File(mzmlDirectory)};
            this.mzmlDirectory = new File(mzmlFiles[0].getAbsoluteFile().getParent());
        } else {
            this.mzmlDirectory = new File(mzmlDirectory);
            Collection<File> mzmlFilesCollection = listFiles(this.mzmlDirectory, new String[]{"mzML"}, false);
            mzmlFiles = new File[mzmlFilesCollection.size()];
            int FileIdx = 0;
            for (File f : mzmlFilesCollection) {
                mzmlFiles[FileIdx] = f;
                FileIdx++;
            }
            Arrays.sort(mzmlFiles);
        }

        //check that corresponding pin files exist
        if (pinDirectory.substring(pinDirectory.length() - 3).toLowerCase().equals("pin")) {
            //check that only one mzml file provided
            if (mzmlFiles.length > 1) {
                throw new IllegalArgumentException("If only one pin file provided, only one mzML file can be provided");
            }

            pinFiles = new File[]{new File(pinDirectory)};
            this.pinDirectory = new File(pinFiles[0].getAbsoluteFile().getParent());
        } else {
            this.pinDirectory = new File(pinDirectory);
            Collection<File> pinFilesCollection = listFiles(this.pinDirectory, new String[]{"pin"}, false);
            HashSet<String> pinFilesSet = new HashSet<String>();
            for (File f : pinFilesCollection) {
                pinFilesSet.add(f.getName());
            }

            pinFiles = new File[mzmlFiles.length];
            for (int i = 0; i < mzmlFiles.length; i++) {
                String name = mzmlFiles[i].getName();
                name = name.substring(0, name.length() - 4) + "pin";

                if (!pinFilesSet.contains(name)) {
                    throw new AssertionError("mzML file must have corresponding pin file. " +
                            pinFiles[i] + " does not exist");
                }
                pinFiles[i] = new File(pinDirectory + File.separator + name);
            }
        }
    }
}
