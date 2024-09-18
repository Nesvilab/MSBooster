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

package readers.predictionreaders;

import static utils.Print.printError;
import static utils.Print.printInfo;

import java.io.File;
import java.io.IOException;
import java.sql.SQLException;
import java.util.Map;
import java.util.concurrent.ExecutionException;
import java.util.concurrent.ExecutorService;

import readers.MgfFileReader;
import predictions.PredictionEntry;
import predictions.PredictionEntryHashMap;
import umich.ms.fileio.exceptions.FileParsingException;

public interface LibraryPredictionMapper {
    //TODO: how many more methods can be made default?
    PredictionEntryHashMap getPreds() throws IOException;
    void setPreds(PredictionEntryHashMap preds);
    void clear();

    //TODO: multithread all of these
    static LibraryPredictionMapper createLibraryPredictionMapper(String file, String model,
                                                                 ExecutorService executorService)
            throws IOException, InterruptedException, ExecutionException, FileParsingException, SQLException {
        //detecting file extension
        String[] extensionSplit = file.split("\\.");
        String extension = extensionSplit[extensionSplit.length - 1];

        switch (extension) {
            case "bin":
                return new DiannSpeclibReader(file);
            case "mgf":
                //for pdeep and alphapeptdeep
                return new MgfFileReader(file, false, executorService, model);
            case "msp":
                return new MspReader(file);
            case "dlib":
                return new DlibReader(file);
            default:
                printError(extension + " is not a valid prediction file format");
                return null;
        }
    }

    static LibraryPredictionMapper createLibraryPredictionMapper(String file, File[] pinFiles,
                                                                 ExecutorService executorService)
            throws IOException, InterruptedException, ExecutionException, FileParsingException {
        return new PredFullSpeclibReader(file, false, pinFiles, executorService);
    }

    //merge the allPreds library into another one
    //could move this to PredictionEntryHashMap
    //need three methods
    //1. merge libraries as is (merge libraries)
    //2. for the final library with full peptide entries, get each property from model-specific peptide entry and transfer it to full peptide entry (assign missing peptide predictions)
    default void mergeLibraries(PredictionEntryHashMap library1, String mode) throws IOException {
        //switch case for spectra, RT, IM, aux
        PredictionEntryHashMap library2 = getPreds();

        for (Map.Entry<String, PredictionEntry> entry : library2.entrySet()) {
            String peptide = entry.getKey();
            PredictionEntry pe2 = entry.getValue();
            if (library1.containsKey(peptide)) {
                PredictionEntry pe1 = library1.get(peptide);
                switch (mode) {
                    //for RT and IM, just replace that
                    case "RT":
                        pe1.setRT(pe2.getRT());
                        break;
                    case "IM":
                        pe1.setIM(pe2.getIM());
                        break;
                    //for spectra and aux, probably need two separate vectors to be saved
                    case "spectra":
                        pe1.mzs = pe2.mzs;
                        pe1.intensities = pe2.intensities;
                        pe1.fragNums = pe2.fragNums;
                        pe1.charges = pe2.charges;
                        pe1.fragmentIonTypes = pe2.fragmentIonTypes;
                        pe1.flags = pe2.flags;
                        pe1.filtered = pe2.filtered;
                        break;
                    case "auxSpectra": //TODO: look below for extra code to implement
                        break;
                }

                library1.put(peptide, pe1); //entry in library1, but getting property for it from library2
            } else {
                library1.put(peptide, pe2); //if entry missing in library1, add it from library2
            }
        }

        //delete old one at end and save library
        clear();

        //Special preparations dependent on features we require
        //get all possible keys from both preds1 and preds2
//        Set<String> totalKeyset = new HashSet<String>();
//        allPreds = predictedSpectra.getPreds();
//        totalKeyset.addAll(allPreds.keySet());
//        totalKeyset.addAll(predictedSpectra2.getPreds().keySet());
//
//        //check what fragment ion types have been predicted by model 1
//        HashSet<String> model1FragmentIonTypes = new HashSet<>();
//        for (PredictionEntry pe : allPreds.values()) {
//            model1FragmentIonTypes.addAll(Arrays.asList(pe.fragmentIonTypes));
//        }

        //        for (String key : totalKeyset) {
//            PredictionEntry pe = allPreds.get(key);
//            if (pe == null) { //missing in prosit/diann
//                allPreds.put(key, predictedSpectra2.getPreds().get(key));
//            } else { //add non-y/b ions
//                ArrayList<Float> mzs = new ArrayList<>();
//                ArrayList<Float> intensities = new ArrayList<>();
//                ArrayList<String> fragTypes = new ArrayList<>();
//
//                if (Constants.addNonYb) {
//                    float maxIntensity = Constants.modelMaxIntensity.get(modelSplit[0]);
//                    float maxIntensityMZ = Float.NaN;
//
//                    //add original peaks
//                    for(int i = 0; i < pe.mzs.length; i++) {
//                        float mz = pe.mzs[i];
//                        float intensity = pe.intensities[i];
//                        String fragType = pe.fragmentIonTypes[i];
//
//                        if (intensity == maxIntensity) {
//                            maxIntensityMZ = mz;
//                        }
//
//                        mzs.add(mz);
//                        intensities.add(intensity);
//                        fragTypes.add(fragType);
//                    }
//
//                    float minMZ = maxIntensityMZ - Constants.DaTolerance;
//                    float maxMZ = maxIntensityMZ + Constants.DaTolerance;
//
//                    //add new peaks
//                    //Scale so that max intensity fragment of diann has same intensity as matched fragment in predfull
//                    //TODO: multiply pe2 intensity by (diann max intensity / predfull intensity of matching fragment)
//                    PredictionEntry pe2 = predictedSpectra2.getPreds().get(key);
//                    //if null, convert to base format
//
//                    if ((!Objects.isNull(pe2)) && (!Objects.isNull(pe2.fragmentIonTypes))) {
//                        float matchedFragInt = Constants.modelMaxIntensity.get(modelSplit[1]);
//                        for (int i = 0; i < pe2.mzs.length; i++) {
//                            float potentialMZ = pe2.mzs[i];
//                            float potentialInt = pe2.intensities[i];
//                            if ((potentialMZ >= minMZ) & (potentialMZ <= maxMZ) & (potentialInt > matchedFragInt)) {
//                                matchedFragInt = potentialInt;
//                            }
//                        }
//
//                        for (int i = 0; i < pe2.fragmentIonTypes.length; i++) {
//                            if (!model1FragmentIonTypes.contains(pe2.fragmentIonTypes[i])) {
//                                mzs.add(pe2.mzs[i]);
//                                intensities.add(pe2.intensities[i] * maxIntensity / matchedFragInt); //putting intensities on same scale
//                                fragTypes.add(pe2.fragmentIonTypes[i]);
//                            }
//                        }
//                    }
//
//                    //convert back to array
//                    float[] mzArray = new float[mzs.size()];
//                    float[] intArray = new float[intensities.size()];
//                    String[] typeArray = new String[fragTypes.size()];
//                    for (int i = 0; i < mzArray.length; i++) {
//                        mzArray[i] = mzs.get(i);
//                        intArray[i] = intensities.get(i);
//                        typeArray[i] = fragTypes.get(i);
//                    }
//
//                    PredictionEntry newPe = new PredictionEntry(mzArray, intArray,
//                            pe.getFragNums(), pe.getCharges(), typeArray, new int[0]);
//                    allPreds.put(key, newPe);
//                } else { //retain predfull intensities, just add RT from other model
//                    //but if predfull has missing entry, use other model instead
//                    PredictionEntry pe2 = predictedSpectra2.getPreds().get(key);
//                    if (!Objects.isNull(pe2)) {
//                        pe2.setRT(pe.getRT());
//                        allPreds.put(key, pe2);
//                    }
//                }
//            }
//        }
    }
}
