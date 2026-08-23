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

import allconstants.Constants;

import java.io.BufferedReader;
import java.io.File;
import java.io.IOException;
import java.io.InputStreamReader;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.List;
import java.util.Map;

import static utils.Print.printError;
import static utils.Print.printInfo;

/**
 * The one place a FragCast child process is started, shared by prediction
 * ({@link FragCastModelCaller}) and fine-tuning ({@link FragCastFineTuneCaller}).
 *
 * <p>It exists because the child's <i>environment</i> is part of FragCast's contract, and getting it
 * wrong fails quietly rather than loudly:
 *
 * <ul>
 *   <li>{@code FRAGCAST_MODEL_DIR} tells FragCast where its pretrained weights live. Without it a
 *       relocated binary falls back to the directory baked in when it was compiled.</li>
 *   <li>{@code OMP_NUM_THREADS}/{@code OPENBLAS_NUM_THREADS} are set to 1. FragCast parallelizes
 *       across peptides itself, so a multi-threaded BLAS underneath it oversubscribes every core.
 *       (FragCast also pins BLAS internally; setting them is belt and braces for other builds.)</li>
 *   <li>{@code RUST_LOG} is forced to {@code info}. FragCast announces a completed export with a
 *       "Wrote fine-tuned model" line printed at info level, and that line is the only evidence
 *       {@link FragCastFineTuneOutcome#unusableReason()} has that a model was written at all, so an
 *       inherited {@code RUST_LOG=warn} would hide it and a perfectly good fine-tuned model would be
 *       discarded for an export nobody saw.</li>
 *   <li>FragCast's diagnostic variables are removed. Each of {@code RT_DUMP}, {@code SPEC_DUMP},
 *       {@code SPEC_GRADCHECK} and friends makes it dump some state and return <em>successfully</em>
 *       without training or predicting, so one inherited from the user's shell would produce an
 *       exit code of 0 and no model.</li>
 * </ul>
 *
 * <p>All FragCast output — including the metrics — goes to stderr, so the two streams are merged and
 * every line is captured for the caller to parse. Nothing here calls {@code System.exit}: it returns
 * the exit code and lets the caller decide whether a failure is fatal.
 */
public class FragCastProcess {
    /**
     * FragCast's own log prefix: {@code [2026-08-23 12:19:31 INFO] }.
     *
     * <p>Echoing a line whole prints two timestamps and two level tags, the second of them wrong -
     * every line goes out under MSBooster's own {@code [INFO]}, so FragCast's warnings arrive
     * labelled as information. The one that matters says a library is too large to finish, and it
     * has to look like a warning. DIA-NN needs none of this because it stamps nothing itself.
     */
    private static final java.util.regex.Pattern LOG_PREFIX = java.util.regex.Pattern.compile(
            "^\\[\\d{4}-\\d{2}-\\d{2} \\d{2}:\\d{2}:\\d{2} (TRACE|DEBUG|INFO|WARN|ERROR)]\\s?(.*)$");

    /** The line as it should be printed, with FragCast's own prefix taken off when it has one. */
    static String message(String line) {
        java.util.regex.Matcher m = LOG_PREFIX.matcher(line);
        return m.matches() ? m.group(2) : line;
    }

    /**
     * Whether this line should print as a problem rather than as information. Only FragCast's own
     * levels decide: a line it did not stamp - a panic, an allocation failure - is passed through as
     * it came, and the exit code the caller reports is what says the run failed.
     */
    static boolean isProblem(String line) {
        java.util.regex.Matcher m = LOG_PREFIX.matcher(line);
        return m.matches() && ("WARN".equals(m.group(1)) || "ERROR".equals(m.group(1)));
    }

    /** Put one line of FragCast's output on the console, at its own level. */
    private static void echo(String line) {
        if (isProblem(line)) {
            printError(message(line));
        } else {
            printInfo(message(line));
        }
    }

    /** FragCast variables that make it return early, successfully, having done no real work. */
    private static final List<String> DIAGNOSTIC_VARS = Arrays.asList(
            "RT_DUMP", "SPEC_DUMP", "SPEC_GRADCHECK", "SPEC_HEAD_ONLY",
            "SPEC_DUMP_SA", "SPEC_DUMP_BASELINE", "SPEC_DUMP_FINAL",
            "FAST_SPEC_ONNX", "FAST_SPEC_FIXTURE");

    /** The outcome of one FragCast invocation: its exit code and everything it printed. */
    public static class Result {
        public final int exitCode;
        public final List<String> output;

        Result(int exitCode, List<String> output) {
            this.exitCode = exitCode;
            this.output = Collections.unmodifiableList(output);
        }

        public boolean succeeded() {
            return exitCode == 0;
        }

        /**
         * What this exit code means in words, or empty when nothing is known about it.
         *
         * <p>The codes a native process dies by say nothing to the person reading the log:
         * {@code -1073740791} is how a FASTA-scale library reports that it ran out of memory, and a
         * run that spent four hours getting there deserves better than the number. The DIA-NN caller
         * translates its own three the same way, and for the same reason.
         *
         * <p>The allocation failure is read off what FragCast printed rather than off the code alone,
         * because {@code 0xC0000409} is Rust's abort for any panic that cannot unwind - the message
         * is what separates "out of memory" from the rest of them.
         */
        public String diagnosis() {
            if (succeeded()) {
                return "";
            }
            for (String line : output) {
                if (line.contains("memory allocation of") && line.contains("failed")) {
                    return "FragCast ran out of memory. A library is held whole until it is written, " +
                            "so the peptide list is the thing to cut: narrow the charge range, or " +
                            "lower fragCastTopN. FragCast reports the size it projects when it starts.";
                }
            }
            switch (exitCode) {
                case 137:
                    return "FragCast was killed for using too much memory.";
                case -1073741515: //0xC0000135, STATUS_DLL_NOT_FOUND
                    return "FragCast could not start: a runtime DLL is missing. The Microsoft Visual " +
                            "C++ Redistributable is the usual one, at " +
                            "https://learn.microsoft.com/en-us/cpp/windows/latest-supported-vc-redist";
                case -1073741819: //0xC0000005, STATUS_ACCESS_VIOLATION
                    return "FragCast crashed with an access violation.";
                case -1073740791: //0xC0000409, what Rust aborts with
                    return "FragCast aborted. Its last lines above say why.";
                default:
                    return "";
            }
        }

        /** The last {@code n} lines, for an error message that points at the real failure. */
        public String tail(int n) {
            return String.join(System.lineSeparator(),
                    output.subList(Math.max(0, output.size() - n), output.size()));
        }
    }

    /**
     * Run FragCast to completion.
     *
     * @param command the full command line, {@code command.get(0)} being the executable
     * @param verbose whether to echo FragCast's output as it arrives
     */
    public static Result run(List<String> command, boolean verbose) throws IOException, InterruptedException {
        ProcessBuilder builder = new ProcessBuilder(command);
        Map<String, String> env = builder.environment();

        String modelDir = resolveModelDir();
        if (modelDir != null && !modelDir.isEmpty()) {
            env.put("FRAGCAST_MODEL_DIR", modelDir);
            if (verbose) {
                printInfo("Using FragCast model directory " + modelDir);
            }
        }
        env.put("OMP_NUM_THREADS", "1");
        env.put("OPENBLAS_NUM_THREADS", "1");
        env.put("RUST_LOG", "info");
        for (String var : DIAGNOSTIC_VARS) {
            env.remove(var);
        }

        if (verbose) {
            printInfo(String.join(" ", command));
        }
        builder.redirectErrorStream(true);
        Process process = builder.start();

        ArrayList<String> output = new ArrayList<>();
        try (BufferedReader reader = new BufferedReader(new InputStreamReader(process.getInputStream()))) {
            String line;
            while ((line = reader.readLine()) != null) {
                output.add(line);
                if (verbose) {
                    echo(line);
                }
            }
        } catch (IOException | RuntimeException e) {
            //nothing will read this process again, so it would sit there indefinitely
            process.destroy();
            throw e;
        }
        return new Result(process.waitFor(), output);
    }

    /**
     * Resolve the directory holding FragCast's pretrained ONNX weights, to be exported as
     * {@code FRAGCAST_MODEL_DIR}. An explicit {@link Constants#FragCastModelDir} param takes
     * precedence; otherwise the directory is derived from the FragCast executable's location, trying
     * both the layout MSBooster ships ({@code pretrained_models} beside or above the exe) and the
     * one FragCast's own repository uses ({@code data/pretrained_models}), so a developer build and
     * an installed build both resolve. Each candidate must actually contain
     * {@code FragCast-RT.onnx}. Returns {@code null} when nothing resolves, leaving FragCast to fall
     * back to its own bundled-model lookup.
     */
    public static String resolveModelDir() {
        if (Constants.FragCastModelDir != null && !Constants.FragCastModelDir.isEmpty()) {
            return Constants.FragCastModelDir;
        }
        if (Constants.FragCast == null || Constants.FragCast.isEmpty()) {
            return null;
        }
        File dir = new File(Constants.FragCast).getAbsoluteFile().getParentFile();
        //walk up from the executable: tools/FragCast/windows/fragcast.exe finds
        //tools/FragCast/pretrained_models, and target/debug/fragcast.exe finds data/pretrained_models
        for (int i = 0; i < 4 && dir != null; i++) {
            for (String relative : new String[] {"pretrained_models", "data/pretrained_models"}) {
                File candidate = new File(dir, relative);
                if (new File(candidate, "FragCast-RT.onnx").isFile()) {
                    return candidate.getPath();
                }
            }
            dir = dir.getParentFile();
        }
        return null;
    }
}
