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

package modelcallers;

import java.io.File;
import java.util.List;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

/**
 * What one FragCast fine-tuning run reported, and whether it left a weights file to predict with.
 *
 * <p>FragCast exits 0 on every one of its degenerate paths - an empty training set, a training set
 * truncated to nothing, a run that never left its learning-rate warmup - and prints a confident
 * "Wrote fine-tuned model" line each time. The exit code alone therefore says nothing, so this
 * class reads the numbers out of the log; {@link #unusableReason} then asks the one question that
 * decides anything, which is whether there is an ONNX file to point {@code --rt-onnx},
 * {@code --im-onnx} or {@code --spec-onnx} at. A fine-tuned model is never refused for being
 * disappointing: how much it helped is FragCast's business and the user's, and every number parsed
 * here is reported rather than judged.
 *
 * <p>One of those checks is less obvious than the rest:
 *
 * <ul>
 *   <li><b>Mean and standard deviation must be real when FragCast reports them.</b> Retention time
 *       and ion mobility are trained on standardized targets and the scaling is folded back into the
 *       exported model. Standardizing an empty set yields a mean of NaN and a standard deviation of
 *       zero, which FragCast writes into the weights without complaint, and every prediction that
 *       file makes is then NaN or the same constant. MS2 standardizes nothing, and an export that
 *       named no scale stated none to check, so neither is refused over a scale nobody reported.</li>
 * </ul>
 */
public class FragCastFineTuneOutcome {
    //every FragCast line is prefixed "[yyyy-MM-dd HH:mm:ss LEVEL] "; strip it before matching
    private static final Pattern LOG_PREFIX = Pattern.compile("^\\[[^\\]]*\\]\\s*");

    private static final Pattern RT_COUNTS =
            Pattern.compile("^Loaded train n=(\\d+).*eval n=(\\d+)");
    private static final Pattern SPEC_COUNTS =
            Pattern.compile("^Spec: train precursors=(\\d+) eval precursors=(\\d+)");
    //the error is printed to exactly three decimals with its unit appended and no separator,
    //so ion mobility reads "MAE=0.1491/K0" - a looser number pattern swallows the unit's
    //leading digit and reports an error that is not the one FragCast measured. The unit is captured
    //with it: "34.688 -> 12.400" says nothing on its own about whether those are minutes or seconds.
    private static final String MAE = "MAE=(-?\\d+\\.\\d{3})(\\S*)";
    private static final Pattern RT_BASELINE =
            Pattern.compile("^Baseline: Pearson=(\\S+) " + MAE);
    private static final Pattern RT_FINAL =
            Pattern.compile("^Fine-tuned: Pearson=(\\S+).*?" + MAE);
    private static final Pattern SPEC_BASELINE =
            Pattern.compile("^Spec baseline: cosine=(\\S+)");
    private static final Pattern SPEC_FINAL =
            Pattern.compile("^Spec fine-tuned: cosine=(\\S+)");
    private static final Pattern EXPORT_DENORM =
            Pattern.compile("^Wrote fine-tuned model \\(de-normalized to physical units, "
                    + "mu=(\\S+) sd=(\\S+)\\) -> (.*)$");
    private static final Pattern EXPORT_PLAIN =
            Pattern.compile("^Wrote fine-tuned model -> (.*)$");

    /**
     * FragCast's exit code for "this library carries nothing this task can be trained on" - a
     * library with no ion mobility in it, most often. Not a failure, and deliberately a number
     * rather than a line of log text, so the two programs agree on it by contract.
     */
    public static final int EXIT_NO_TRAINING_DATA = 3;

    public final String task;
    public final int exitCode;
    public final File savedOnnx;
    public final int trainCount;
    public final int evalCount;
    /** Correlation for retention time and ion mobility, cosine similarity for MS2. */
    public final double baselineMetric;
    public final double finalMetric;
    /**
     * Mean absolute error, in the target's own units; NaN for MS2, which reports no such thing.
     *
     * <p>Worth reading beside the correlation rather than instead of it: correlation is
     * scale-invariant, so it can improve while the predictions drift away from the measurements.
     * Neither number is a bar the model has to clear - both are printed by {@link #summary()} and
     * judged by nobody here.
     */
    public final double baselineMae;
    public final double finalMae;
    /** The unit FragCast appended to the error ({@code s}, {@code 1/K0}), empty if it named none. */
    public final String maeUnit;
    /**
     * The standardization FragCast folded back into the exported weights; NaN for MS2, and NaN for
     * an export that did not name it - which is what {@link #denormalizedExport} distinguishes.
     */
    public final double labelMean;
    public final double labelSd;
    /** Did the export line state the standardization? Only then are the two above real numbers. */
    public final boolean denormalizedExport;
    public final boolean exported;

    FragCastFineTuneOutcome(String task, int exitCode, File savedOnnx, int trainCount, int evalCount,
                            double baselineMetric, double finalMetric, double baselineMae,
                            double finalMae, String maeUnit, double labelMean, double labelSd,
                            boolean denormalizedExport, boolean exported) {
        this.task = task;
        this.exitCode = exitCode;
        this.savedOnnx = savedOnnx;
        this.trainCount = trainCount;
        this.evalCount = evalCount;
        this.baselineMetric = baselineMetric;
        this.finalMetric = finalMetric;
        this.baselineMae = baselineMae;
        this.finalMae = finalMae;
        this.maeUnit = maeUnit == null ? "" : maeUnit;
        this.labelMean = labelMean;
        this.labelSd = labelSd;
        this.denormalizedExport = denormalizedExport;
        this.exported = exported;
    }

    /**
     * Did FragCast find nothing to train this task on? Then there is no model to be had and nothing
     * went wrong - the caller skips the task rather than reporting a failure.
     */
    public boolean noTrainingData() {
        return exitCode == EXIT_NO_TRAINING_DATA;
    }

    /** The name of the metric this task reports, for a log line a user can act on. */
    public String metricName() {
        return "spec".equals(task) ? "cosine" : "Pearson";
    }

    /** Parse one run's output. {@code savedOnnx} is where the run was told to write its weights. */
    public static FragCastFineTuneOutcome parse(String task, int exitCode, File savedOnnx,
                                                List<String> output) {
        int trainCount = 0;
        int evalCount = 0;
        double baselineMetric = Double.NaN;
        double finalMetric = Double.NaN;
        double baselineMae = Double.NaN;
        double finalMae = Double.NaN;
        String maeUnit = "";
        double labelMean = Double.NaN;
        double labelSd = Double.NaN;
        boolean denormalizedExport = false;
        boolean exported = false;

        for (String raw : output) {
            String line = LOG_PREFIX.matcher(raw).replaceFirst("").trim();

            Matcher m = RT_COUNTS.matcher(line);
            if (m.find()) {
                trainCount = Integer.parseInt(m.group(1));
                evalCount = Integer.parseInt(m.group(2));
                continue;
            }
            m = SPEC_COUNTS.matcher(line);
            if (m.find()) {
                trainCount = Integer.parseInt(m.group(1));
                evalCount = Integer.parseInt(m.group(2));
                continue;
            }
            m = RT_BASELINE.matcher(line);
            if (m.find()) {
                baselineMetric = parseDouble(m.group(1));
                baselineMae = parseDouble(m.group(2));
                maeUnit = m.group(3);
                continue;
            }
            m = RT_FINAL.matcher(line);
            if (m.find()) {
                finalMetric = parseDouble(m.group(1));
                finalMae = parseDouble(m.group(2));
                maeUnit = m.group(3);
                continue;
            }
            m = SPEC_BASELINE.matcher(line);
            if (m.find()) {
                baselineMetric = parseDouble(m.group(1));
                continue;
            }
            m = SPEC_FINAL.matcher(line);
            if (m.find()) {
                finalMetric = parseDouble(m.group(1));
                continue;
            }
            m = EXPORT_DENORM.matcher(line);
            if (m.find()) {
                labelMean = parseDouble(m.group(1));
                labelSd = parseDouble(m.group(2));
                denormalizedExport = true;
                exported = true;
                continue;
            }
            if (EXPORT_PLAIN.matcher(line).find()) {
                exported = true;
            }
        }
        return new FragCastFineTuneOutcome(task, exitCode, savedOnnx, trainCount, evalCount,
                baselineMetric, finalMetric, baselineMae, finalMae, maeUnit, labelMean, labelSd,
                denormalizedExport, exported);
    }

    /**
     * Is there a fine-tuned weights file to predict with?
     *
     * <p>Every check below asks the same question - did this run leave an ONNX file behind that
     * {@code FragCastRtOnnx}, {@code FragCastImOnnx} or {@code FragCastSpecOnnx} can name - and none
     * of them looks at how well the model scored. A fine-tuned model is always used; the metrics
     * parsed alongside it are the run log, not a bar it has to clear. FragCast selects the best
     * epoch against its own held-out slice, and how good is good enough is a judgement it and the
     * user make, not one MSBooster makes on their behalf.
     *
     * @return {@code null} when the exported model can be used, otherwise why there is none to use
     */
    public String unusableReason() {
        if (exitCode != 0) {
            return "FragCast exited with code " + exitCode;
        }
        if (!exported) {
            return "FragCast reported no exported model";
        }
        if (savedOnnx == null || !savedOnnx.isFile() || savedOnnx.length() == 0) {
            return "no weights were written to " + savedOnnx;
        }
        if (!"spec".equals(task) && denormalizedExport) {
            //standardizing an empty target set gives a NaN mean and a zero deviation, which the
            //exporter folds into the weights without complaint; every prediction that file makes is
            //then NaN or the same constant. Only the de-normalized export line states the scale, so
            //a run that exported without one is left alone rather than assumed to be degenerate.
            if (!isFinite(labelMean) || !(labelSd > 0)) {
                return "FragCast standardized the targets to mean=" + labelMean + " sd=" + labelSd +
                        ", which corrupts the exported model";
            }
        }
        return null;
    }

    /** A one-line summary of what the run measured, for the run log. */
    public String summary() {
        String line = String.format("%s: train=%d eval=%d %s %.4f -> %.4f", task, trainCount,
                evalCount, metricName(), baselineMetric, finalMetric);
        //MS2 measures no error at all, so there is nothing to append for it
        if (isFinite(baselineMae) || isFinite(finalMae)) {
            line += String.format(" MAE %.3f -> %.3f%s", baselineMae, finalMae,
                    maeUnit.isEmpty() ? "" : " " + maeUnit);
        }
        return line;
    }

    private static double parseDouble(String text) {
        try {
            return Double.parseDouble(text);
        } catch (NumberFormatException e) {
            return Double.NaN; //FragCast prints a literal NaN for a degenerate run
        }
    }

    //Double.isFinite arrived in Java 8 but reads less clearly at the call sites above
    private static boolean isFinite(double value) {
        return !Double.isNaN(value) && !Double.isInfinite(value);
    }
}
