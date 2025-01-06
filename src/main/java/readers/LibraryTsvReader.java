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

package readers;

import allconstants.Constants;
import allconstants.FragmentIonConstants;
import com.google.common.collect.Multimap;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.Collection;
import java.util.List;
import java.util.Set;
import peptideptmformatting.PeptideFormatter;
import predictions.PredictionEntry;
import predictions.PredictionEntryHashMap;
import readers.predictionreaders.LibraryPredictionMapper;
import umich.ms.fileio.filetypes.library.LibraryTsv;
import umich.ms.fileio.filetypes.library.Transition;

public class LibraryTsvReader implements LibraryPredictionMapper {

    public List<String> filenames = new ArrayList<>();
    private PredictionEntryHashMap allPreds = new PredictionEntryHashMap();
    protected PredictionEntryHashMap allPredsHashMap = new PredictionEntryHashMap();

    public LibraryTsvReader(String file, String format) throws Exception {
        filenames.add(file);
        LibraryTsv libraryTsv = new LibraryTsv(Constants.numThreads, Constants.unimodObo);
        Multimap<String, Transition> transitions = libraryTsv.read(Paths.get(file));

        try {
            Set<String> ignoredFragmentIonTypesSet = FragmentIonConstants.makeIgnoredFragmentIonTypes();

            for (String k : transitions.keySet()) {
                Collection<Transition> tt = transitions.get(k);
                for (Transition t : tt) {
                    String peptide = t.peptide.getUnimodPeptide(false, libraryTsv.massSiteUnimodTable, null, null, null, '[', ']').replaceFirst("^n", "");
                    String charge = t.peptideCharge + "";
                    String basePep = new PeptideFormatter(peptide, charge, format).getBaseCharge();

                    float[] mzArray = new float[t.fragments.length];
                    float[] intArray = new float[t.fragments.length];
                    String[] fragmentIonTypes = new String[t.fragments.length];
                    for (int i = 0; i < t.fragments.length; i++) { //might need to decrease array length if removing some fragments
                        if (!ignoredFragmentIonTypesSet.contains(t.fragments[i].type + "")) {
                            mzArray[i] = t.fragments[i].mz;
                            intArray[i] = t.fragments[i].intensity;
                            fragmentIonTypes[i] = t.fragments[i].type + "";
                        }
                    }

                    PredictionEntry newPred = new PredictionEntry(mzArray, intArray, new int[0], new int[0], fragmentIonTypes, new int[0]);
                    newPred.setRT(t.normalizedRetentionTime);
                    newPred.setIM(t.precursorIonMobility);
                    allPreds.put(basePep, newPred);
                }
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
