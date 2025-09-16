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

import predictions.PredictionEntryHashMap;
import readers.MgfFileReader;

import java.io.File;
import java.io.IOException;
import java.util.concurrent.ExecutorService;

import static utils.Print.printError;

public interface LibraryPredictionMapper {
    //TODO: how many more methods can be made default?
    PredictionEntryHashMap getPreds() throws IOException;
    void setPreds(PredictionEntryHashMap preds);
    void clear();

    //TODO: multithread all of these
    static LibraryPredictionMapper createLibraryPredictionMapper(String file, String model,
                                                                 ExecutorService executorService)
            throws Exception {
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
            case "tsv":
                return new LibraryTsvReader(file, "unimod.obo");
            default:
                printError(extension + " is not a valid prediction file format");
                return null;
        }
    }

    static LibraryPredictionMapper createLibraryPredictionMapper(String file, File[] pinFiles,
                                                                 ExecutorService executorService)
            throws Exception {
        return new PredFullSpeclibReader(file, false, pinFiles, executorService);
    }

//old notes for merging predfull library
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
