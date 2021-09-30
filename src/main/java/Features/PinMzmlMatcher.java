package Features;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collection;
import java.util.HashMap;

import static org.apache.commons.io.FileUtils.listFiles;

public class PinMzmlMatcher {
    public File[] mzmlFiles;
    public File[] pinFiles;

    public PinMzmlMatcher(String mzmlDirectory, String pinDirectory) throws IOException {
        //get pin files
        String[] allPinDirectories = pinDirectory.split(" ");
        ArrayList<String> pinFileList = new ArrayList<>();
        for (String directory : allPinDirectories) {
            if (directory.substring(directory.length() - 3).toLowerCase().equals("pin")) { //single file
                pinFileList.add(directory);
            } else { //directory, but not recursive
                Collection<File> pinFilesCollection = listFiles(new File(pinDirectory), new String[]{"pin"}, false);
                for (File file : pinFilesCollection) {
                    pinFileList.add(file.getCanonicalPath());
                }
            }
        }

        //look for corresponding mzml files
        String[] allMzmlDirectories = mzmlDirectory.split(" ");
        HashMap<String, File> mzmlFileMap = new HashMap<>();
        for (String directory : allMzmlDirectories) {
            File f = new File(directory);
            if (directory.substring(directory.length() - 4).toLowerCase().equals("mzml")) {
                mzmlFileMap.put(f.getName().substring(0, f.getName().length() - 5), f);
            } else if (directory.substring(directory.length() - 3).toLowerCase().equals("mgf")) {
                mzmlFileMap.put(f.getName().substring(0, f.getName().length() - 4), f);
            } else { //directory
                Collection<File> mzmlFilesCollection = listFiles(f, new String[]{"mzML"}, false);
                for (File file : mzmlFilesCollection) {
                    mzmlFileMap.put(file.getName().substring(0, file.getName().length() - 5), file);
                }

                mzmlFilesCollection = listFiles(f, new String[]{"mgf"}, false);
                for (File file : mzmlFilesCollection) {
                    //editing if calibrated or uncalibrated in name
                    String mgfName = file.getName();
                    if (mgfName.contains("_uncalibrated.mgf")) {
                        mgfName = mgfName.substring(0, mgfName.length() - 17);
                    } else if (mgfName.contains("_calibrated.mgf")) {
                        mgfName = mgfName.substring(0, mgfName.length() - 15);
                    }
                    mzmlFileMap.put(mgfName, file);
                }
            }
        }

        //add files to array
        pinFileList.sort(String::compareToIgnoreCase);
        pinFiles = new File[pinFileList.size()];
        mzmlFiles = new File[mzmlFileMap.size()];
        for (int i = 0; i < pinFiles.length; i++) {
            pinFiles[i] = new File(pinFileList.get(i));
            String baseName = pinFiles[i].getName().substring(0, pinFiles[i].getName().length() - 4);
            if (mzmlFileMap.containsKey(baseName)) {
                mzmlFiles[i] = mzmlFileMap.get(baseName);
            } else { //no matching mzml file
                throw new IllegalArgumentException("No matching mzml/mgf file for " + baseName + ".pin, " +
                        "please check that provided mzml directories contain proper mzml files.");
            }
        }
    }
}
