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

import allconstants.Constants;
import allconstants.FragmentIonConstants;
import com.google.common.collect.Multimap;
import peptideptmformatting.PeptideFormatter;
import predictions.PredictionEntry;
import predictions.PredictionEntryHashMap;
import umich.ms.fileio.filetypes.library.LibraryTsv;
import umich.ms.fileio.filetypes.library.Transition;
import utils.Multithreader;
import utils.ProgressReporter;

import java.nio.file.Paths;
import java.util.*;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Future;

import static peptideptmformatting.PTMhandler.setUnimodObo;
import static utils.Print.printInfo;

public class LibraryTsvReader implements LibraryPredictionMapper {

    public List<String> filenames = new ArrayList<>();
    private PredictionEntryHashMap allPreds = new PredictionEntryHashMap();
    protected PredictionEntryHashMap allPredsHashMap = new PredictionEntryHashMap();

    public LibraryTsvReader(String file, ExecutorService executorService, String format) throws Exception {
        this(file, executorService, format, new HashSet<>());
    }

    public LibraryTsvReader(String file, ExecutorService executorService, String format, HashSet<String> allowedPrecursors) throws Exception {
        filenames.add(file);
        printInfo("Reading " + file);
        setUnimodObo();
        LibraryTsv libraryTsv = new LibraryTsv(Constants.unimodObo);
        Multimap<String, Transition> transitions = libraryTsv.read(Paths.get(file));

        try {
            Set<String> ignoredFragmentIonTypesSet = FragmentIonConstants.makeIgnoredFragmentIonTypes();

            List<Future> futureList = new ArrayList<>(Constants.numThreads);
            ProgressReporter pr = new ProgressReporter(transitions.size());
            Multithreader mt = new Multithreader(transitions.size(), Constants.numThreads);
            String[] transitionKeys = transitions.keySet().toArray(new String[0]);

            for (int j = 0; j < Constants.numThreads; j++) {
                int finalJ = j;
                futureList.add(executorService.submit(() -> {
                    for (int k = mt.indices[finalJ]; k < mt.indices[finalJ + 1]; k++) {
                        Collection<Transition> tt = transitions.get(transitionKeys[k]);
                        for (Transition t : tt) {
                            String peptide = t.peptide.getUnimodPeptide(false, libraryTsv.massSiteUnimodTable,
                                            null, null, null, '[', ']')
                                    .replaceFirst("^n", "");
                            String charge = t.peptideCharge + "";
                            String basePep = new PeptideFormatter(peptide, charge, format).getBaseCharge(); //format is unimod.obo

                            //only read in precursors in the pin files
                            if (!allowedPrecursors.isEmpty()) {
                                if (!allowedPrecursors.contains(basePep)) {
                                    continue;
                                }
                            }

                            float[] mzArray = new float[t.fragments.length];
                            float[] intArray = new float[t.fragments.length];
                            String[] fragmentIonTypes = new String[t.fragments.length];
                            int[] fragNums = new int[t.fragments.length];
                            int[] fragCharges = new int[t.fragments.length];
                            for (int i = 0; i < t.fragments.length; i++) { //might need to decrease array length if removing some fragments
                                if (!ignoredFragmentIonTypesSet.contains(t.fragments[i].type + "")) {
                                    mzArray[i] = t.fragments[i].mz;
                                    intArray[i] = t.fragments[i].intensity;
                                    fragmentIonTypes[i] = t.fragments[i].type + "";
                                    fragNums[i] = t.fragments[i].ordinal;
                                    fragCharges[i] = t.fragments[i].charge;
                                }
                            }

                            PredictionEntry newPred = new PredictionEntry(mzArray, intArray, fragNums, fragCharges, fragmentIonTypes);
                            newPred.setRT(t.normalizedRetentionTime);
                            newPred.setIM(t.precursorIonMobility);
                            allPreds.put(basePep, newPred);
                        }
                        pr.progress();
                    }
                }));
            }

            for (Future future : futureList) {
                future.get();
            }
        } catch (Exception e) {
            e.printStackTrace();
            System.exit(1);
        }
    }

    public PredictionEntryHashMap getPreds() {
        if (allPredsHashMap.isEmpty()) {
            allPredsHashMap.putAll(allPreds);
            allPreds.clear(); //no longer need concurrency
        }
        return allPredsHashMap;
    }
    public void setPreds(PredictionEntryHashMap preds) {
        allPreds = preds;
    }

    public void clear() {
        allPreds.clear();
        allPredsHashMap.clear();
    }
}
