/*
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or
 * implied. See the License for the specific language governing
 * permissions and limitations under the License.
 */

package readers.predictionreaders;

import predictions.PredictionEntryHashMap;
import readers.MgfFileReader;
import readers.datareaders.PinReader;
import utils.Print;
import writers.LibraryTsvReducer;

import java.io.File;
import java.io.IOException;
import java.util.HashSet;
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
        return createLibraryPredictionMapper(file, model, executorService, new File[0]);
    }

    //assume diann predicted.bin
    static LibraryPredictionMapper createLibraryPredictionMapper(String file) throws Exception {
        return new DiannSpeclibReader(file);
    }

    static LibraryPredictionMapper createLibraryPredictionMapper(String file, String model,
                                                                 ExecutorService executorService,
                                                                 File[] pinFiles)
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
                //TODO: if we ever encounter a situation where other prediction formats have OOM, extend this idea
                //write a reduced library tsv first
                LibraryTsvReducer ltr = new LibraryTsvReducer(file, pinFiles);
                ltr.writeReducedLib();

                if (ltr.needReduced) {
                    //get set of precursors in pin files, so only have to read in relevant entries
                    Print.printInfo("Collecting precursors from pin file to only load relevant precursors");
                    //TODO: another hashset to store previously checked peptides, not peptide formatters
                    HashSet<String> allPrecursors = new HashSet<>();
                    HashSet<String> allPrecursorsUnprocessed = new HashSet<>();
                    for (File pinFile : pinFiles) {
                        PinReader pinReader = new PinReader(pinFile.getAbsolutePath());
                        while (pinReader.next(true)) {
                            String unprocessedPeptide = pinReader.getColumn("Peptide");
                            if (!allPrecursors.contains(unprocessedPeptide)) {
                                allPrecursors.add(pinReader.getPep().getBaseCharge());
                                allPrecursorsUnprocessed.add(unprocessedPeptide);
                            }
                        }
                    }

                    return new LibraryTsvReader(ltr.newLib, executorService, "unimod.obo", allPrecursors);
                } else {
                    return new LibraryTsvReader(file, executorService, "unimod.obo");
                }
            case "parquet":
            case "pq": {
                //FragCast (or any DIA-NN/Spectronaut) library written as Parquet. Collect the pin
                //precursors so only relevant entries are kept (empty when no pins -> keep everything).
                HashSet<String> parquetPrecursors = new HashSet<>();
                for (File pinFile : pinFiles) {
                    PinReader pinReader = new PinReader(pinFile.getAbsolutePath());
                    while (pinReader.next(true)) {
                        parquetPrecursors.add(pinReader.getPep().getBaseCharge());
                    }
                }
                return new ParquetSpeclibReader(file, executorService, parquetPrecursors);
            }
            default:
                printError(file + ": ." + extension + " is not a valid prediction file format. " +
                        "Use one of .bin, .mgf, .msp, .dlib, .tsv, .parquet or .pq.");
                return null;
        }
    }
}
