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

import static utils.Print.printInfo;

import allconstants.Constants;
import readers.datareaders.MzmlReader;
import readers.datareaders.PinReader;
import umich.ms.fileio.exceptions.FileParsingException;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Set;
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
        setFragmentationType();
        setNCE();
    }

    private void loadMzmlReaders() throws IOException, FileParsingException, ExecutionException, InterruptedException {
        boolean needToRemove = false;
        for (int j = 0; j < mzmlReaders.length; j++) {
            if (mzmlReaders[j] == null) {
                MzmlReader mzml = new MzmlReader(mzmlFiles[j].getCanonicalPath());
                //set hasITMS and ppmTolerance
                PinReader pin = new PinReader(pinFiles[j].getCanonicalPath());
                if (pin.next(true)) {
                    mzml.getScanNumObject(pin.getScanNum());
                    mzmlReaders[j] = mzml;
                } else { //if any empty pins, remove
                    printInfo(pin.name + " is empty. Will not process this file.");
                    needToRemove = true;
                }
                pin.close();
            }
        }

        //if any empty pins, remove
        if (needToRemove) {
            ArrayList<File> mzmlFilesList = new ArrayList<>();
            ArrayList<File> pinFilesList = new ArrayList<>();
            ArrayList<MzmlReader> mzmlReadersList = new ArrayList<>();

            for (int j = 0; j < mzmlReaders.length; j++) {
                if (mzmlReaders[j] != null) {
                    mzmlFilesList.add(mzmlFiles[j]);
                    mzmlFiles[j] = null;
                    pinFilesList.add(pinFiles[j]);
                    pinFiles[j] = null;
                    mzmlReadersList.add(mzmlReaders[j]);
                    mzmlReaders[j] = null;
                }
            }

            mzmlReaders = new MzmlReader[mzmlReadersList.size()];
            pinFiles = new File[mzmlReadersList.size()];
            mzmlFiles = new File[mzmlReadersList.size()];
            for (int i = 0; i < mzmlReadersList.size(); i++) {
                mzmlReaders[i] = mzmlReadersList.get(i);
                pinFiles[i] = pinFilesList.get(i);
                mzmlFiles[i] = mzmlFilesList.get(i);
            }
        }
    }

    private void setFragmentationType() {
        if (Constants.FragmentationType.isEmpty()) {
            try {
                Set<String> fragTypes = mzmlReaders[0].
                        getScanNumObject(mzmlReaders[0].getScanNums().first()).NCEs.keySet();
                if (fragTypes.contains("HCD")) {
                    Constants.FragmentationType = "HCD";
                    printInfo("Fragmentation type detected: " + Constants.FragmentationType);
                } else if (fragTypes.contains("CID")) {
                    Constants.FragmentationType = "CID";
                    printInfo("Fragmentation type detected: " + Constants.FragmentationType);
                } else {
                    printInfo("No fragmentation type detected. Setting fragmentation type to HCD. " +
                            "You can specify this with '--FragmentationType' via the command line " +
                            "or 'FragmentationType=' in the param file.");
                    Constants.FragmentationType = "HCD";
                }
            } catch (Exception e) {
                printInfo("No fragmentation type detected. Setting fragmentation type to HCD. " +
                        "You can specify this with '--FragmentationType' via the command line " +
                        "or 'FragmentationType=' in the param file.");
                Constants.FragmentationType = "HCD";
            }
        }
    }

    private void setNCE() {
        if (Constants.NCE.isEmpty()) {
            try {
                Float NCE = mzmlReaders[0].getScanNumObject(mzmlReaders[0].getScanNums().first()).
                        NCEs.get(Constants.FragmentationType);
                if (NCE != null) {
                    Constants.NCE = String.valueOf(NCE);
                    printInfo("NCE detected: " + Constants.NCE);
                } else {
                    printInfo("No NCE detected. Setting NCE to 25. " +
                            "You can specify this with '--NCE' via the command line " +
                            "or 'NCE=' in the param file.");
                    Constants.NCE = "25";
                }
            } catch (Exception e) {
                printInfo("No NCE detected. Setting NCE to 25. " +
                        "You can specify this with '--NCE' via the command line " +
                        "or 'NCE=' in the param file.");
                Constants.NCE = "25";
            }
        }
    }
}
