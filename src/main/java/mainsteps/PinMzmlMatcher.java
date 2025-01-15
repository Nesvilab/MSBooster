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

package mainsteps;

import allconstants.Constants;
import readers.datareaders.MzmlReader;
import readers.datareaders.PinReader;
import umich.ms.datatypes.LCMSDataSubset;
import umich.ms.datatypes.scan.IScan;
import umich.ms.fileio.exceptions.FileParsingException;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.concurrent.ExecutionException;
import java.util.stream.Collectors;

public class PinMzmlMatcher {
    public File[] mzmlFiles;
    public File[] pinFiles;
    public MzmlReader[] mzmlReaders;

    public PinMzmlMatcher(String mzmlDirectory, String pinDirectory)
            throws IOException, FileParsingException, ExecutionException, InterruptedException {
        //get pin files
        String[] allPinDirectories = pinDirectory.split(" ");
        ArrayList<String> pinFileList = new ArrayList<>();
        for (String directory : allPinDirectories) {
            if (directory.substring(directory.length() - 3).equalsIgnoreCase("pin")) { //single file
                if (! directory.contains("_" + Constants.editedPin)) {
                    pinFileList.add(directory);
                }
            } else { //directory, but not recursive
                List<File> pinFilesCollection = Files.list(Paths.get(directory))
                    .filter(Files::isRegularFile)
                    .filter(p -> p.getFileName().toString().toLowerCase().endsWith(".pin"))
                    .map(Path::toFile).collect(Collectors.toList());
                for (File file : pinFilesCollection) {
                    if (! file.getCanonicalPath().contains("_" + Constants.editedPin)) {
                        pinFileList.add(file.getCanonicalPath());
                    }
                }
            }
        }

        //look for corresponding mzml files
        //can maybe allow uncalibrated.pin to look for uncalibrated.mzml
        String[] allMzmlDirectories = mzmlDirectory.split(" ");
        HashMap<String, File> mzmlFileMap = new HashMap<>();
        for (String directory : allMzmlDirectories) {
            File f = new File(directory);
            if (f.isFile()) {
                if (directory.substring(directory.length() - 17).equalsIgnoreCase("uncalibrated.mzml")) {
                    mzmlFileMap.put(f.getName().substring(0, f.getName().length() - 18), f);
                    mzmlFileMap.put(f.getName().substring(0, f.getName().length() - 5), f);
                } else if (directory.substring(directory.length() - 4).equalsIgnoreCase("mzml")) {
                    mzmlFileMap.put(f.getName().substring(0, f.getName().length() - 5), f);
                } else if (directory.substring(directory.length() - 16).equalsIgnoreCase("uncalibrated.mgf")) {
                    mzmlFileMap.put(f.getName().substring(0, f.getName().length() - 17), f);
                } else if (directory.substring(directory.length() - 3).equalsIgnoreCase("mgf")) {
                    mzmlFileMap.put(f.getName().substring(0, f.getName().length() - 4), f);
                }
            } else { //directory
                List<File> mzmlFilesCollection = Files.list(f.toPath())
                    .filter(Files::isRegularFile)
                    .filter(p -> p.getFileName().toString().toLowerCase().endsWith(".mzml"))
                    .map(Path::toFile).collect(Collectors.toList());
                for (File file : mzmlFilesCollection) {
                    if (file.getName().toLowerCase().contains("_uncalibrated.mzml")) {
                        mzmlFileMap.put(file.getName().substring(0, file.getName().length() - 18), file);
                    }
                    mzmlFileMap.put(file.getName().substring(0, file.getName().length() - 5), file);
                }

                mzmlFilesCollection = Files.list(f.toPath())
                    .filter(Files::isRegularFile)
                    .filter(p -> p.getFileName().toString().toLowerCase().endsWith(".mgf"))
                    .map(Path::toFile).collect(Collectors.toList());
                for (File file : mzmlFilesCollection) {
                    //editing if calibrated or uncalibrated in name
                    String mgfName = file.getName();
                    if (mgfName.toLowerCase().contains("_uncalibrated.mgf")) {
                        mgfName = mgfName.substring(0, mgfName.length() - 17);
                    }
                    if (! mzmlFileMap.containsKey(mgfName.substring(0, mgfName.length() - 4))) { //only want mzml
                        mzmlFileMap.put(mgfName.substring(0, mgfName.length() - 4), file);
                    }
                }
            }
        }

        //add files to array
        pinFileList.sort(String::compareToIgnoreCase);
        pinFiles = new File[pinFileList.size()];
        mzmlFiles = new File[pinFileList.size()];
        for (int i = 0; i < pinFiles.length; i++) {
            pinFiles[i] = new File(pinFileList.get(i));
            String baseName = pinFiles[i].getName().substring(0, pinFiles[i].getName().length() - 4);
            if (mzmlFileMap.containsKey(baseName)) {
                mzmlFiles[i] = mzmlFileMap.get(baseName);
            } else { //no matching mzml file
                throw new IllegalArgumentException("No matching mzml/mgf file for " + baseName + ".pin, " +
                        "please check that provided mzml/mgf directories contain proper mzml/mgf files.");
            }
        }

        mzmlReaders = new MzmlReader[mzmlFiles.length];
        loadMzmlReaders();
        mzmlReaders[0].detectMassSpecType();
        mzmlReaders[0].setFragmentationType();
        mzmlReaders[0].setNCE();
    }

    private void loadMzmlReaders() throws IOException, FileParsingException, ExecutionException, InterruptedException {
        for (int j = 0; j < mzmlReaders.length; j++) {
            if (mzmlReaders[j] == null) {
                MzmlReader mzml = new MzmlReader(mzmlFiles[j].getCanonicalPath());
                if (j == 0) { //initialize first ms2 scan in first mzml so we can read in metadata
                    PinReader pin = new PinReader(pinFiles[j].getCanonicalPath());
                    if (pin.next(true)) {
                        mzml.getScanNumObject(pin.getScanNum());
                    } else { //empty pin
                        mzml.scans.loadData(LCMSDataSubset.STRUCTURE_ONLY);
                        IScan iscan = mzml.scans.getNextScanAtMsLevel(0, 2);
                        mzml.getScanNumObject(iscan.getNum());
                    }
                    pin.close();
                }
                mzmlReaders[j] = mzml;
            }
        }
    }
}
