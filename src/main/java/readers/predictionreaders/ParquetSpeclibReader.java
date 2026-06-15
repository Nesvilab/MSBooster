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
import blue.strategic.parquet.Hydrator;
import blue.strategic.parquet.HydratorSupplier;
import blue.strategic.parquet.ParquetReader;
import peptideptmformatting.PeptideFormatter;
import predictions.PredictionEntry;
import predictions.PredictionEntryHashMap;
import utils.Multithreader;
import utils.ProgressReporter;

import java.io.File;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashSet;
import java.util.Iterator;
import java.util.List;
import java.util.Set;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Future;
import java.util.regex.Pattern;
import java.util.stream.Stream;

import static peptideptmformatting.PTMhandler.setUnimodObo;
import static utils.Print.printInfo;

/**
 * Reads a FragCast (or any DIA-NN/Spectronaut) spectral library written as a Parquet file
 * ({@code build-library --format parquet}). It produces the exact same {@link PredictionEntryHashMap}
 * as {@link LibraryTsvReader} does for the equivalent TSV: precursors are grouped on consecutive
 * {@code (ModifiedPeptideSequence, PrecursorCharge)} rows, the key is the canonical base+charge
 * obtained by rewriting {@code (UniMod:N)} to {@code [UniMod:N]} and running it through
 * {@link PeptideFormatter} with the {@code unimod.obo} format, and per-fragment fields map 1:1.
 */
public class ParquetSpeclibReader implements LibraryPredictionMapper {
    private PredictionEntryHashMap allPreds = new PredictionEntryHashMap();
    private PredictionEntryHashMap allPredsHashMap = new PredictionEntryHashMap();

    //FragCast writes mods as (UniMod:N) parentheses; rewrite to [UniMod:N] so the unimod.obo
    //PeptideFormatter (which expects brackets) yields the same canonical key as the TSV reader.
    private static final Pattern UNIMOD_PARENS = Pattern.compile("\\(UniMod:(\\d+)\\)");

    //the subset of the 19-column library schema we actually need (column projection -> faster reads)
    private static final String MODSEQ = "ModifiedPeptideSequence", CHARGE = "PrecursorCharge",
            PRODMZ = "ProductMz", INTENSITY = "LibraryIntensity", RT = "NormalizedRetentionTime",
            IM = "PrecursorIonMobility", FTYPE = "FragmentType", FCHARGE = "FragmentCharge",
            FSERIES = "FragmentSeriesNumber";
    private static final List<String> PROJECTION = Arrays.asList(
            MODSEQ, CHARGE, PRODMZ, INTENSITY, RT, IM, FTYPE, FCHARGE, FSERIES);

    //one parquet row (a single transition / fragment)
    private static final class Row {
        String modSeq = "";
        int charge;
        double prodMz, intensity, rt, im;
        String fType = "";
        int fCharge, fSeries;
    }

    //all fragment rows for one precursor (consecutive in the file)
    private static final class Group {
        String modSeq;
        int charge;
        float rt, im;
        final ArrayList<Row> rows = new ArrayList<>();
    }

    public ParquetSpeclibReader(String file, ExecutorService executorService,
                                HashSet<String> allowedPrecursors) throws Exception {
        printInfo("Reading " + file);
        setUnimodObo();

        //1) stream parquet rows (in file order) and group consecutive rows by (modSeq, charge)
        ArrayList<Group> groups = readGroups(file);

        //2) build predictions in parallel, mirroring LibraryTsvReader's per-fragment logic
        if (groups.isEmpty()) {
            return;
        }
        Set<String> ignoredFragmentIonTypesSet = FragmentIonConstants.makeIgnoredFragmentIonTypes();
        ProgressReporter pr = new ProgressReporter(groups.size());
        Multithreader mt = new Multithreader(groups.size(), Constants.numThreads);
        List<Future> futureList = new ArrayList<>(Constants.numThreads);
        for (int j = 0; j < Constants.numThreads; j++) {
            int finalJ = j;
            futureList.add(executorService.submit(() -> {
                for (int k = mt.indices[finalJ]; k < mt.indices[finalJ + 1]; k++) {
                    Group g = groups.get(k);
                    String bracketPep = UNIMOD_PARENS.matcher(g.modSeq).replaceAll("[UniMod:$1]");
                    String basePep = new PeptideFormatter(bracketPep, g.charge + "", "unimod.obo").getBaseCharge();

                    //only keep precursors present in the pin files (empty set = keep everything)
                    if (allowedPrecursors.isEmpty() || allowedPrecursors.contains(basePep)) {
                        int nf = g.rows.size();
                        float[] mzArray = new float[nf];
                        float[] intArray = new float[nf];
                        String[] fragmentIonTypes = new String[nf];
                        int[] fragNums = new int[nf];
                        int[] fragCharges = new int[nf];
                        for (int i = 0; i < nf; i++) {
                            Row row = g.rows.get(i);
                            String type = row.fType.isEmpty() ? "" : row.fType.substring(0, 1);
                            if (!ignoredFragmentIonTypesSet.contains(type)) {
                                mzArray[i] = (float) row.prodMz;
                                intArray[i] = (float) row.intensity;
                                fragmentIonTypes[i] = type;
                                fragNums[i] = row.fSeries;
                                fragCharges[i] = row.fCharge;
                            }
                        }
                        PredictionEntry pred = new PredictionEntry(mzArray, intArray, fragNums, fragCharges, fragmentIonTypes);
                        pred.setRT(g.rt);
                        pred.setIM(g.im);
                        allPreds.put(basePep, pred);
                    }
                    pr.progress();
                }
            }));
        }
        for (Future future : futureList) {
            future.get();
        }
    }

    private ArrayList<Group> readGroups(String file) throws Exception {
        ArrayList<Group> groups = new ArrayList<>();
        Hydrator<Row, Row> hydrator = new Hydrator<Row, Row>() {
            @Override public Row start() { return new Row(); }
            @Override public Row add(Row r, String name, Object v) {
                switch (name) {
                    case MODSEQ: r.modSeq = str(v); break;
                    case CHARGE: r.charge = intv(v); break;
                    case PRODMZ: r.prodMz = dbl(v); break;
                    case INTENSITY: r.intensity = dbl(v); break;
                    case RT: r.rt = dbl(v); break;
                    case IM: r.im = dbl(v); break;
                    case FTYPE: r.fType = str(v); break;
                    case FCHARGE: r.fCharge = intv(v); break;
                    case FSERIES: r.fSeries = intv(v); break;
                    default: break;
                }
                return r;
            }
            @Override public Row finish(Row r) { return r; }
        };

        Group cur = null;
        try (Stream<Row> stream = ParquetReader.streamContent(new File(file),
                HydratorSupplier.constantly(hydrator), PROJECTION)) {
            Iterator<Row> it = stream.iterator();
            while (it.hasNext()) {
                Row r = it.next();
                if (cur == null || cur.charge != r.charge || !cur.modSeq.equals(r.modSeq)) {
                    cur = new Group();
                    cur.modSeq = r.modSeq;
                    cur.charge = r.charge;
                    cur.rt = (float) r.rt;
                    cur.im = (float) r.im;
                    groups.add(cur);
                }
                cur.rows.add(r);
            }
        }
        return groups;
    }

    @Override
    public PredictionEntryHashMap getPreds() {
        if (allPredsHashMap.isEmpty()) {
            allPredsHashMap.putAll(allPreds);
            allPreds.clear(); //no longer need concurrency
        }
        return allPredsHashMap;
    }

    @Override
    public void setPreds(PredictionEntryHashMap preds) {
        allPreds = preds;
    }

    @Override
    public void clear() {
        allPreds.clear();
        allPredsHashMap.clear();
    }

    private static double dbl(Object o) {
        if (o instanceof Number) {
            return ((Number) o).doubleValue();
        }
        String s = str(o);
        return s.isEmpty() ? 0d : Double.parseDouble(s);
    }

    private static int intv(Object o) {
        if (o instanceof Number) {
            return ((Number) o).intValue();
        }
        String s = str(o);
        return s.isEmpty() ? 0 : Integer.parseInt(s);
    }

    private static String str(Object o) {
        return o == null ? "" : o.toString();
    }
}
