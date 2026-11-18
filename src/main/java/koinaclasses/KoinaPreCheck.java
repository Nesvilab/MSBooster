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

package koinaclasses;

import peptideptmformatting.PTMhandler;
import peptideptmformatting.PeptideFormatter;
import peptideptmformatting.PeptideSkipper;

import java.util.ArrayList;
import java.util.Comparator;
import java.util.LinkedHashMap;
import java.util.LinkedHashSet;
import java.util.List;
import java.util.Map;
import java.util.concurrent.ConcurrentHashMap;
import java.util.concurrent.atomic.LongAdder;

import static utils.Print.printInfo;

/**
 * What a Koina model is sent of a run's PSMs, by its spec, before any request is written. A PSM the model cannot
 * take ({@link KoinaModelSpec#limits()}: residues, length, charge, free cysteines) is not sent, and transferKoinaPreds
 * scores it as unpredictable. A mod the model cannot encode ({@link KoinaModelSpec#modRule()},
 * {@link KoinaModelSpec#encodableMods()}, or a terminal mod on a terminal residue with a mod of its own for a format
 * that writes terminal mods on residues) is left out of the peptide sent, and its mass is added back to the fragments
 * that hold it when the predictions are moved onto the PSM's own peptide. These are the rules PeptideSkipper and
 * PeptideFormatter.getModel apply; this class counts what they did, so a run says it once per model
 * ({@link #report}), and ends when nothing is left to send ({@link #requirePredictable}).
 * <p>
 * The pin files may be read on several threads at once, so the counts are thread-safe.
 * <p>
 * It also holds the check of the search's isobaric label against a model ({@link #checkIsobaricLabel}).
 */
public final class KoinaPreCheck {
    private final KoinaModelSpec spec;
    private final LongAdder psms = new LongAdder();
    //PSMs not sent, by why
    private final ConcurrentHashMap<String, LongAdder> skipped = new ConcurrentHashMap<>();
    //PSMs sent with a mod left out, by that mod as its residue and mass, e.g. S[79.9663]
    private final ConcurrentHashMap<String, LongAdder> droppedMods = new ConcurrentHashMap<>();

    public KoinaPreCheck(KoinaModelSpec spec) {
        this.spec = spec;
    }

    public KoinaModelSpec spec() {
        return spec;
    }

    /**
     * The peptide of this PSM as the model is sent it, or null if the model cannot take the PSM. Counts the PSM: why
     * it is not sent, or the mods left out of the peptide sent, each once per PSM.
     * <p>
     * It asks PeptideSkipper and then PeptideFormatter.getModel, in that order, as the request writers always have,
     * so a model that remembers the UniMod codes it finds (PeptideFormatter.foundUnimods) finds the same ones.
     */
    public String modelPeptide(PeptideFormatter pf) {
        psms.increment();
        String reason = PeptideSkipper.skipReason(pf, spec);
        if (reason != null) {
            count(skipped, reason);
            return null;
        }
        List<String> dropped = new ArrayList<>(0);
        String peptide = pf.getModel(spec, dropped);
        for (String mod : new LinkedHashSet<>(dropped)) {
            count(droppedMods, mod);
        }
        return peptide;
    }

    private static void count(ConcurrentHashMap<String, LongAdder> counts, String key) {
        counts.computeIfAbsent(key, k -> new LongAdder()).increment();
    }

    /** The PSMs looked at. */
    public long psms() {
        return psms.sum();
    }

    /** The PSMs the model cannot take. */
    public long skippedPsms() {
        long sum = 0;
        for (LongAdder count : skipped.values()) {
            sum += count.sum();
        }
        return sum;
    }

    /** The PSMs not sent, by why, most first. */
    public Map<String, Long> skippedByReason() {
        return sorted(skipped);
    }

    /** The PSMs sent with a mod left out, by that mod (residue, or n or c for a terminus, and mass), most first. */
    public Map<String, Long> droppedModPsms() {
        return sorted(droppedMods);
    }

    private static Map<String, Long> sorted(ConcurrentHashMap<String, LongAdder> counts) {
        List<Map.Entry<String, Long>> entries = new ArrayList<>();
        for (Map.Entry<String, LongAdder> e : counts.entrySet()) {
            entries.add(Map.entry(e.getKey(), e.getValue().sum()));
        }
        entries.sort(Map.Entry.<String, Long>comparingByValue(Comparator.reverseOrder())
                .thenComparing(Map.Entry.comparingByKey()));
        Map<String, Long> sorted = new LinkedHashMap<>();
        for (Map.Entry<String, Long> e : entries) {
            sorted.put(e.getKey(), e.getValue());
        }
        return sorted;
    }

    /**
     * Prints, in one block, how many PSMs the model is sent, those it is not sent by why, and the mods left out of
     * the peptides sent with the PSMs that have them.
     *
     * @param precursors the distinct precursors sent
     */
    public void report(int precursors) {
        long all = psms();
        long skippedPsms = skippedPsms();
        printInfo(spec.name() + " is sent " + (all - skippedPsms) + " of " + all + " PSMs, as " + precursors
                + (precursors == 1 ? " precursor" : " precursors"));
        if (skippedPsms > 0) {
            printInfo("  " + skippedPsms + " PSMs not sent, as the model cannot take them: "
                    + counts(skippedByReason()));
        }
        if (!droppedMods.isEmpty()) {
            printInfo("  mods the model cannot encode, left out of the peptides sent (PSMs): "
                    + counts(droppedModPsms())
                    + (spec.kind() == KoinaModelSpec.Kind.MS2 ? "; the fragments that hold one are shifted by its mass"
                    : "; predicted as if unmodified there"));
        }
    }

    /**
     * Ends the run when the model is sent nothing: none of the PSMs is one it can take (D28).
     *
     * @param precursors the distinct precursors sent
     * @throws KoinaException with the reasons the PSMs are not sent, if there is no precursor to send
     */
    public void requirePredictable(int precursors) {
        if (precursors > 0) {
            return;
        }
        throw new KoinaException(psms() == 0 ? noneMessage()
                : noneMessage() + ", so it has nothing to predict. Choose another model");
    }

    /** Why the model is sent no precursor: there are no PSMs, or it cannot take any, by why. */
    public String noneMessage() {
        long all = psms();
        if (all == 0) {
            return "Koina model " + spec.name() + " has nothing to predict: there are no PSMs";
        }
        return "Koina model " + spec.name() + " cannot take any of the " + all + " PSMs: "
                + counts(skippedByReason());
    }

    //"reason (n), reason (m)"
    private static String counts(Map<String, Long> counts) {
        StringBuilder sb = new StringBuilder();
        for (Map.Entry<String, Long> e : counts.entrySet()) {
            if (sb.length() > 0) {
                sb.append(", ");
            }
            sb.append(e.getKey()).append(" (").append(e.getValue()).append(')');
        }
        return sb.toString();
    }

    /**
     * Ends the run if it would send peptides with the search's isobaric label (TMT or iTRAQ, from fragger.params) to
     * a model that cannot encode one ({@link KoinaModelSpec#takesIsobaricLabels()}), such as IM2Deep (D25). Any other
     * model, and a run without such a label, passes.
     *
     * @throws KoinaException naming models of the same kind that take the label
     */
    public static void checkIsobaricLabel(String model) {
        if (!takesSearchLabel(model)) {
            throw new KoinaException(isobaricLabelMessage(KoinaModelRegistry.require(model)));
        }
    }

    /**
     * Whether the run can send this model its peptides: true unless it is a Koina model that cannot encode isobaric
     * labels and the search has one.
     */
    public static boolean takesSearchLabel(String model) {
        KoinaModelSpec spec = KoinaModelRegistry.find(model);
        return spec == null || spec.takesIsobaricLabels() || !PTMhandler.hasIsobaricLabel();
    }

    /** Why a model that cannot encode isobaric labels cannot be used, and what to use instead. */
    public static String isobaricLabelMessage(KoinaModelSpec spec) {
        List<String> others = new ArrayList<>();
        for (KoinaModelSpec other : KoinaModelRegistry.all()) {
            if (other.kind() == spec.kind() && other.takesIsobaricLabels()) {
                others.add(other.name());
            }
        }
        return spec.name() + " cannot encode isobaric labels (TMT/iTRAQ), and the search in fragger.params has one ("
                + String.format("%.4f", PTMhandler.getTmtMass()) + ")."
                + (others.isEmpty() ? " Choose another model" : " Use " + String.join(" or ", others) + " instead");
    }
}
