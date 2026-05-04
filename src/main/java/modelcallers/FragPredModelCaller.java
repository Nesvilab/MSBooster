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

package modelcallers;

import allconstants.Constants;
import com.fragpred.cli.ModelDir;
import com.fragpred.properties.api.PropertyType;
import com.fragpred.properties.data.ParsedPeptide;
import com.fragpred.properties.featurize.ChargeOneHot;
import com.fragpred.properties.featurize.PeptideEncoder;
import com.fragpred.properties.library.Calibration;
import com.fragpred.properties.models.OnnxFragRtImPredictor;
import com.fragpred.properties.models.OnnxFragSpecPredictor;
import features.spectra.MassCalculator;
import peptideptmformatting.PeptideFormatter;
import predictions.PredictionEntry;
import predictions.PredictionEntryHashMap;

import java.io.IOException;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collection;
import java.util.Comparator;
import java.util.List;
import java.util.concurrent.ExecutionException;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.Future;
import java.util.concurrent.TimeUnit;

import static utils.Print.printError;
import static utils.Print.printInfo;

/**
 * In-process FragPred prediction (RT, CCS/IM, and MS2 fragment intensities via
 * ONNX Runtime). Builds a {@link PredictionEntryHashMap} directly from in-memory
 * peptide records — no intermediate files are read or written.
 *
 * <p>Input is a collection of {@code peptide<tab>charge} strings with the
 * peptides in UniMod bracket syntax (e.g. {@code AAA[UniMod:35]K}). The optional
 * {@code fullRecords} argument is the corresponding base-format peptide list and
 * is used to alias predictions for peptides whose base format does not
 * round-trip through the UniMod conversion (mods that get stripped during
 * formatting).
 */
public class FragPredModelCaller {

    // Mirror FragPred LibraryBuilder.Options defaults — apply the same fragment
    // filters here so MSBooster's top-fragment selection sees the same candidate
    // set the file-based library would.
    private static final int FRAG_TOP_N = 20;
    private static final double FRAG_MIN_MZ = 200.0;
    private static final double FRAG_MIN_REL_INTENSITY = 0.01;
    private static final int FRAG_MIN_SIZE = 2;

    // FragPred model limits: PeptideEncoder uses MAX_LEN-residue arrays and the
    // standard 20 amino acids, ChargeOneHot covers charges 1..6, and LibraryBuilder
    // rejects plain length <3 or >MAX_LEN. Records outside these bounds get
    // truncated/zero-filled tensors and produce nonsense RT/IM/MS2, so we drop
    // them before encoding.
    private static final int FRAG_MIN_PLAIN_LEN = 3;
    private static final int FRAG_MAX_PLAIN_LEN = PeptideEncoder.MAX_LEN;
    private static final int FRAG_MIN_CHARGE = 1;
    private static final int FRAG_MAX_CHARGE = ChargeOneHot.N_CHARGES;
    private static final String FRAG_AA_LIST = PeptideEncoder.AA_LIST;

    /**
     * Run FragPred on a set of peptide records and return predictions in-memory.
     *
     * @param records      {@code peptide<tab>charge} strings (UniMod bracket syntax).
     * @param fullRecords  optional collection of {@code basePeptide<tab>charge} strings
     *                     used to alias predictions when {@code base→fragpred→base}
     *                     does not round-trip. Pass {@code null} to skip the alias pass.
     * @param verbose      print progress messages.
     * @return {@link PredictionEntryHashMap} keyed by {@code base|charge}.
     */
    public static PredictionEntryHashMap predict(Collection<String> records,
                                                  Collection<String> fullRecords,
                                                  boolean verbose) {
        long startTime = System.nanoTime();
        if (verbose) {
            printInfo("Generating FragPred predictions");
        }
        PredictionEntryHashMap allPreds;
        try {
            // FragPred ships its own model resolver — looks under
            // data/pretrained_models/ relative to the cwd or to the FragPred jar,
            // and falls back to the bundled ONNX resources inside the FragPred jar.
            Path modelDir = ModelDir.resolveDir(null);
            allPreds = runPrediction(records, modelDir, verbose);

            // When the user-supplied base peptide round-trips through
            // base→fragpred→base to a different base, copy the prediction
            // (with re-computed m/z) under the original base.
            if (fullRecords != null && !fullRecords.isEmpty()) {
                applyAliases(fullRecords, allPreds, verbose);
            }

            if (verbose) {
                printInfo("Done generating FragPred predictions");
            }
        } catch (Exception e) {
            e.printStackTrace();
            printError("FragPred prediction failed");
            System.exit(1);
            return null; // unreachable
        }

        if (verbose) {
            long endTime = System.nanoTime();
            printInfo("Model running took " + (endTime - startTime) / 1_000_000 + " milliseconds");
        }
        return allPreds;
    }

    private static PredictionEntryHashMap runPrediction(java.util.Collection<String> records,
                                                        Path modelDir, boolean verbose)
            throws IOException, ExecutionException, InterruptedException {
        // Each record is "peptide<tab>charge" (UniMod bracket syntax).
        List<String> peptidesIn = new ArrayList<>(records.size());
        List<String> charges = new ArrayList<>(records.size());
        for (String rec : records) {
            int tab = rec.indexOf('\t');
            if (tab < 0) continue;
            peptidesIn.add(rec.substring(0, tab));
            charges.add(rec.substring(tab + 1).trim());
        }

        int n = peptidesIn.size();
        PredictionEntryHashMap allPreds = new PredictionEntryHashMap(n);
        if (n == 0) {
            if (verbose) {
                printInfo("FragPred: predicting 0 peptides");
            }
            return allPreds;
        }

        // FragPred's PeptideEncoder.parse auto-detects three syntaxes — UniMod
        // parentheses (UniMod:N), delta-mass [mass], and inline short tags (ac)/(ox).
        // We hand it the parentheses-UniMod form because that is FragPred's preferred
        // input and uses its UnimodTable lookup directly (avoiding the narrower
        // DeltaMassTable). PeptideFormatter already has a converter to that exact
        // syntax — getLibrarytsv() — which produces e.g. AAAM(UniMod:35)K from the
        // base [mass] form. Keep the [mass] base form alongside for MassCalculator.
        // Filter out records outside FragPred's supported bounds — those would be
        // silently truncated/zero-filled by the encoder/charge one-hot.
        List<ParsedPeptide> parsedKept = new ArrayList<>(n);
        List<Integer> chargeKept = new ArrayList<>(n);
        List<Integer> plainLenKept = new ArrayList<>(n);
        List<MassCalculator> mcKept = new ArrayList<>(n);
        int dropped = 0;
        for (int i = 0; i < n; i++) {
            String peptide = peptidesIn.get(i);
            String chargeStr = charges.get(i);
            int charge = parseIntSafe(chargeStr, -1);
            if (charge < FRAG_MIN_CHARGE || charge > FRAG_MAX_CHARGE) {
                dropped++;
                continue;
            }
            PeptideFormatter pf = new PeptideFormatter(peptide, chargeStr, "fragpred");
            String base = pf.getBase();
            int plainLen = countPlainResidues(base);
            if (plainLen < 0 || plainLen < FRAG_MIN_PLAIN_LEN || plainLen > FRAG_MAX_PLAIN_LEN) {
                dropped++;
                continue;
            }
            ParsedPeptide pp = PeptideEncoder.parse(pf.getLibrarytsv());
            parsedKept.add(pp);
            chargeKept.add(charge);
            plainLenKept.add(plainLen);
            mcKept.add(new MassCalculator(base, chargeStr));
        }

        int kept = parsedKept.size();
        if (verbose) {
            if (dropped > 0) {
                printInfo("FragPred: predicting " + kept + " peptides (skipped " + dropped
                        + " outside model limits: length " + FRAG_MIN_PLAIN_LEN + ".." + FRAG_MAX_PLAIN_LEN
                        + ", charge " + FRAG_MIN_CHARGE + ".." + FRAG_MAX_CHARGE
                        + ", residues " + FRAG_AA_LIST + ")");
            } else {
                printInfo("FragPred: predicting " + kept + " peptides");
            }
        }
        if (kept == 0) return allPreds;

        ParsedPeptide[] parsed = parsedKept.toArray(new ParsedPeptide[0]);
        int[] chargeArr = new int[kept];
        int[] plainLens = new int[kept];
        MassCalculator[] mcs = mcKept.toArray(new MassCalculator[0]);
        for (int i = 0; i < kept; i++) {
            chargeArr[i] = chargeKept.get(i);
            plainLens[i] = plainLenKept.get(i);
        }

        int threads = Math.max(1, Constants.numThreads == null ? 1 : Constants.numThreads);
        int batchSize = 4;
        int specBatchSize = 4;

        float[] rtZ = new float[kept];
        float[] imZ = new float[kept];
        float[][][] specOut = new float[kept][][];

        // ONNX Runtime sessions are heavyweight; create once per task and run batches in parallel.
        ExecutorService pool = Executors.newFixedThreadPool(threads);
        try {
            try (OnnxFragRtImPredictor rt = new OnnxFragRtImPredictor(
                    PropertyType.RT,
                    ModelDir.resolveModel(ModelDir.DEFAULT_RT_NAME, modelDir), 1)) {
                runScalar(rt, parsed, chargeArr, batchSize, pool, rtZ);
            }
            try (OnnxFragRtImPredictor im = new OnnxFragRtImPredictor(
                    PropertyType.IM,
                    ModelDir.resolveModel(ModelDir.DEFAULT_IM_NAME, modelDir), 1)) {
                runScalar(im, parsed, chargeArr, batchSize, pool, imZ);
            }
            try (OnnxFragSpecPredictor spec = new OnnxFragSpecPredictor(
                    ModelDir.resolveModel(ModelDir.DEFAULT_SPEC_NAME, modelDir), 1)) {
                runSpec(spec, parsed, chargeArr, specBatchSize, pool, specOut);
            }
        } finally {
            pool.shutdown();
            try { pool.awaitTermination(5, TimeUnit.MINUTES); } catch (InterruptedException ignored) {}
        }

        // FragPred's RT/IM ONNX heads emit z-scores in normalized space; LibraryBuilder
        // calibrates them per-charge before writing a library. We mirror that here with
        // the default ballpark fit (RT = 40z + 50, IM = 0.05z + 0.75) so that downstream
        // consumers see physical RT minutes and 1/K0 values rather than raw z-scores.
        Calibration.PerChargeCalibration cal = Calibration.ballpark();

        // Build PredictionEntry objects and populate the hashmap.
        for (int i = 0; i < kept; i++) {
            float rt = (float) cal.applyRt(chargeArr[i], rtZ[i]);
            float im = (float) cal.applyIm(chargeArr[i], imZ[i]);
            PredictionEntry pe = buildEntry(specOut[i], mcs[i], chargeArr[i], plainLens[i],
                    rt, im);
            allPreds.put(mcs[i].fullPeptide, pe);
        }
        return allPreds;
    }

    /**
     * Count plain (mod-stripped) residues in a base-format peptide string. Letters
     * inside [...] mass brackets are ignored. Returns {@code -1} if any residue
     * outside {@link #FRAG_AA_LIST} appears outside a bracket — those peptides are
     * unsupported by FragPred's encoder (token 0 fallback) and must be skipped.
     */
    private static int countPlainResidues(String base) {
        int plainLen = 0;
        int depth = 0;
        for (int k = 0; k < base.length(); k++) {
            char c = base.charAt(k);
            if (c == '[') { depth++; continue; }
            if (c == ']') { if (depth > 0) depth--; continue; }
            if (depth > 0) continue;
            if (Character.isLetter(c)) {
                if (FRAG_AA_LIST.indexOf(Character.toUpperCase(c)) < 0) {
                    return -1;
                }
                plainLen++;
            }
        }
        return plainLen;
    }

    private static PredictionEntry buildEntry(float[][] spec, MassCalculator mc,
                                              int precursorCharge, int plainLen,
                                              float rt, float im) {
        if (spec == null || plainLen < 2) {
            return emptyEntry(rt, im);
        }
        int maxFragZ = Math.max(1, Math.min(3, precursorCharge - 1));
        int maxBonds = Math.min(spec.length, plainLen - 1);

        // Upper-bound m/z = protonated parent mass; LibraryBuilder filters fragments
        // at fragMz < precursorMz * charge, which simplifies to (mass + charge*proton).
        double precursorMz = mc.calcMassPrecursor(0, 0f, precursorCharge);
        double maxFragMz = precursorMz * precursorCharge;

        // Collect all candidate fragments that pass the per-fragment filters.
        // Mirrors LibraryBuilder.buildLinesForPrecursor.
        List<Cand> cands = new ArrayList<>();
        for (int b = 0; b < maxBonds; b++) {
            int fragNum = b + 1; // 1-based fragment series number
            if (fragNum < FRAG_MIN_SIZE) continue; // skip y1/b1
            float[] row = spec[b];
            for (int j = 0; j < OnnxFragSpecPredictor.N_ION_SLOTS; j++) {
                int fz = (j < 3) ? (j + 1) : (j - 2);
                if (fz > maxFragZ) continue;
                float intRel = row[j];
                if (intRel <= 0f) continue;
                String ionType = (j < 3) ? "y" : "b";
                float fragMZ = mc.calcMass(fragNum, ionType, fz, 0);
                if (!Float.isFinite(fragMZ)) continue;
                if (fragMZ < FRAG_MIN_MZ) continue;
                if (fragMZ >= maxFragMz) continue;
                cands.add(new Cand(ionType, fragNum, fz, fragMZ, intRel));
            }
        }
        if (cands.isEmpty()) {
            return emptyEntry(rt, im);
        }

        // Sort by sigmoid intensity descending and keep the topN.
        cands.sort(Comparator.comparingDouble((Cand c) -> -c.intensity));
        int kept = Math.min(FRAG_TOP_N, cands.size());
        float maxInt = cands.get(0).intensity;
        if (maxInt <= 0f) {
            return emptyEntry(rt, im);
        }

        // Drop anything below the relative-intensity floor within the topN window;
        // because cands are sorted, this trims a contiguous tail.
        float minAbsInt = (float) (FRAG_MIN_REL_INTENSITY * maxInt);
        while (kept > 0 && cands.get(kept - 1).intensity < minAbsInt) {
            kept--;
        }
        if (kept == 0) {
            return emptyEntry(rt, im);
        }

        float[] mzArr = new float[kept];
        float[] intArr = new float[kept];
        int[] fragNumArr = new int[kept];
        int[] fragChargeArr = new int[kept];
        String[] ionTypeArr = new String[kept];
        for (int i = 0; i < kept; i++) {
            Cand c = cands.get(i);
            mzArr[i] = c.mz;
            // Keep the historical 0–60000 scale used by the DIA-NN-style intensity
            // column so downstream code that treats this as relative intensity is
            // unaffected by the new filtering.
            intArr[i] = c.intensity / maxInt * 60000f;
            fragNumArr[i] = c.fragNum;
            fragChargeArr[i] = c.fragCharge;
            ionTypeArr[i] = c.ionType;
        }
        PredictionEntry pe = new PredictionEntry(mzArr, intArr, fragNumArr, fragChargeArr, ionTypeArr);
        pe.setRT(rt);
        pe.setIM(im);
        return pe;
    }

    private static PredictionEntry emptyEntry(float rt, float im) {
        PredictionEntry empty = new PredictionEntry(new float[0], new float[0],
                new int[0], new int[0], new String[0]);
        empty.setRT(rt);
        empty.setIM(im);
        return empty;
    }

    private static final class Cand {
        final String ionType;
        final int fragNum;
        final int fragCharge;
        final float mz;
        final float intensity;
        Cand(String ionType, int fragNum, int fragCharge, float mz, float intensity) {
            this.ionType = ionType;
            this.fragNum = fragNum;
            this.fragCharge = fragCharge;
            this.mz = mz;
            this.intensity = intensity;
        }
    }

    private static void applyAliases(Collection<String> fullRecords,
                                     PredictionEntryHashMap allPreds, boolean verbose) {
        int aliased = 0;
        for (String rec : fullRecords) {
            int tab = rec.indexOf('\t');
            if (tab < 0) continue;
            String peptideSeq = rec.substring(0, tab);
            String chargeStr = rec.substring(tab + 1).trim();
            PeptideFormatter pf = new PeptideFormatter(
                    new PeptideFormatter(peptideSeq, chargeStr, "base").getFragpred(),
                    chargeStr, "fragpred");
            if (pf.getBase().equals(peptideSeq)) continue;
            PredictionEntry tmp = allPreds.get(pf.getBaseCharge());
            if (tmp == null) continue;
            MassCalculator mc = new MassCalculator(peptideSeq, chargeStr);
            int nfrag = tmp.numFragments();
            float[] newMZs = new float[nfrag];
            for (int i = 0; i < nfrag; i++) {
                newMZs[i] = mc.calcMass(tmp.getFragNum(i), tmp.getIonTypeString(i),
                        tmp.getCharge(i),
                        tmp.isotopes.length > 0 ? tmp.isotopes[i] : 0);
            }
            PredictionEntry copy = new PredictionEntry(newMZs, tmp.getIntensities(),
                    tmp.getFragNums(), tmp.getCharges(), tmp.getFragmentIonTypes());
            copy.setRT(tmp.RT);
            copy.setIM(tmp.IM);
            allPreds.put(mc.fullPeptide, copy);
            aliased++;
        }
        if (verbose && aliased > 0) {
            printInfo("FragPred: aliased " + aliased + " entries from full peptide list");
        }
    }

    private static void runScalar(OnnxFragRtImPredictor predictor, ParsedPeptide[] parsed,
                                  int[] charges, int batch, ExecutorService pool,
                                  float[] outZ)
            throws ExecutionException, InterruptedException {
        int n = parsed.length;
        List<Future<?>> futures = new ArrayList<>();
        for (int s = 0; s < n; s += batch) {
            final int start = s;
            final int end = Math.min(s + batch, n);
            futures.add(pool.submit(() -> {
                List<ParsedPeptide> sub = Arrays.asList(parsed).subList(start, end);
                int[] cz = Arrays.copyOfRange(charges, start, end);
                float[] part = predictor.predictBatch(sub, cz);
                System.arraycopy(part, 0, outZ, start, end - start);
            }));
        }
        for (Future<?> f : futures) f.get();
    }

    private static void runSpec(OnnxFragSpecPredictor predictor, ParsedPeptide[] parsed,
                                int[] charges, int batch, ExecutorService pool,
                                float[][][] outSpec)
            throws ExecutionException, InterruptedException {
        int n = parsed.length;
        List<Future<?>> futures = new ArrayList<>();
        for (int s = 0; s < n; s += batch) {
            final int start = s;
            final int end = Math.min(s + batch, n);
            futures.add(pool.submit(() -> {
                List<ParsedPeptide> sub = Arrays.asList(parsed).subList(start, end);
                int[] cz = Arrays.copyOfRange(charges, start, end);
                float[] cce = new float[end - start];
                Arrays.fill(cce, 30f); // default normalized collision energy
                float[][][] part = predictor.predictBatch(sub, cz, cce);
                for (int i = 0; i < part.length; i++) {
                    outSpec[start + i] = part[i];
                }
            }));
        }
        for (Future<?> f : futures) f.get();
    }

    private static int parseIntSafe(String s, int def) {
        try { return Integer.parseInt(s.trim()); } catch (NumberFormatException e) { return def; }
    }
}
