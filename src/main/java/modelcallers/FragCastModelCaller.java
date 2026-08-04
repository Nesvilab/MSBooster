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

import java.io.*;
import java.nio.file.Files;
import java.util.ArrayList;
import java.util.List;

import static utils.Print.printError;
import static utils.Print.printInfo;

/**
 * Runs the FragCast Rust predictor (selectable as the FragCast model for spectra/RT/IM).
 *
 * <p>The input is a {@code peptide<TAB>charge} list (header {@code peptide\tcharge}) written by
 * {@link writers.PeptideFileCreator}, where each line is one precursor and each peptide carries its
 * modifications as delta masses (e.g. {@code C[57.0215]}); FragCast resolves those against its full
 * UniMod table, so every modification MSBooster knows is preserved. FragCast's {@code build-library}
 * reads the peptide and charge columns (the charge is mandatory per row). A single call predicts RT,
 * IM and MS2 at once and writes a 19-column DIA-NN/Spectronaut spectral library as Parquet, read back
 * by {@link readers.predictionreaders.ParquetSpeclibReader}.
 *
 * <p>{@code fast} selects FragCast's small/fast Spec model (its {@code --fast} flag,
 * {@code FragCast-Spec-Fast.onnx}) instead of the default Conformer; see
 * {@link allconstants.FragCastModels}. Only the MS2 intensities change — RT and IM are predicted by
 * the same weights either way, which is why a job naming both variants (fast spectra with Conformer
 * RT/IM) calls this once, under the spectra model's weights. The two variants still write to
 * different files: the name records which Spec weights produced a library, so findBest scoring them
 * against each other, and a rerun that switches variants while keeping its predictions, never read
 * or overwrite the other's.
 */
public class FragCastModelCaller {
    /**
     * The output library path for an input peptide file: {@code spectraRT.tsv} ->
     * {@code spectraRT.predicted.parquet}, or {@code spectraRT.fast.predicted.parquet} under
     * {@code fast}. The two Spec variants must not share a path — the name is what says which
     * weights wrote a library.
     */
    static String predictionFile(String inputFile, boolean fast) {
        return inputFile.substring(0, inputFile.length() - 4) //strip the trailing ".tsv"
                + (fast ? ".fast" : "") + ".predicted.parquet";
    }

    /**
     * The FragCast command line. Every flag here exists in FragCast's {@code build-library} task;
     * {@code --fast} swaps the Conformer Spec weights for the small/fast ones and is FragCast's only
     * way to choose between them (there is no per-property task to call).
     */
    static List<String> buildCommand(String inputFile, String predFileString, boolean fast) {
        List<String> command = new ArrayList<>();
        command.add(Constants.FragCast);
        command.add("--task");
        command.add("build-library");
        if (fast) {
            command.add("--fast");
        }
        command.add("--in");
        command.add(inputFile);
        command.add("--out");
        command.add(predFileString);
        command.add("--format");
        command.add("parquet");
        command.add("--threads");
        command.add(String.valueOf(Constants.numThreads));
        command.add("--top-n");
        command.add(String.valueOf(Constants.fragCastTopN));
        command.add("--min-frag-mz");
        command.add(String.valueOf(Constants.fragCastMinFragMz));
        command.add("--min-rel-intensity");
        command.add(String.valueOf(Constants.fragCastMinRelIntensity));
        command.add("--min-frag-size");
        command.add(String.valueOf(Constants.fragCastMinFragSize));
        return command;
    }

    public static String callModel(String inputFile, boolean verbose, boolean fast) {
        long startTime = System.nanoTime();
        String predFileString = predictionFile(inputFile, fast);
        try {
            if (Constants.FragCast == null) {
                printError("path to FragCast executable must be provided (set the FragCast parameter)");
                System.exit(1);
            }
            if (verbose) {
                printInfo("Generating FragCast predictions" + (fast ? " (small/fast Spec model)" : ""));
            }

            ProcessBuilder builder = new ProcessBuilder(buildCommand(inputFile, predFileString, fast));
            //Tell FragCast where its pretrained ONNX weights live via FRAGCAST_MODEL_DIR. An explicit
            //FragCastModelDir param wins; otherwise we derive a "pretrained_models" directory from the
            //executable's location (e.g. tools/FragCast/pretrained_models when the exe is in
            //tools/FragCast/windows). If neither resolves, FragCast falls back to its own lookup.
            String modelDir = resolveModelDir();
            if (modelDir != null && !modelDir.isEmpty()) {
                builder.environment().put("FRAGCAST_MODEL_DIR", modelDir);
                if (verbose) {
                    printInfo("Using FragCast model directory " + modelDir);
                }
            }
            //FragCast links the OpenMP build of OpenBLAS and parallelizes across peptides with its own
            //thread pool, so BLAS itself must stay single-threaded. The OpenMP pool is governed by
            //these env vars (read before the runtime initializes); without them every worker spawns a
            //full BLAS pool and oversubscribes the cores, which is orders of magnitude slower.
            builder.environment().put("OMP_NUM_THREADS", "1");
            builder.environment().put("OPENBLAS_NUM_THREADS", "1");
            if (verbose) {
                printInfo(String.join(" ", builder.command()));
            }
            builder.redirectErrorStream(true);
            Process process = builder.start();

            //print FragCast output while running
            try (BufferedReader reader = new BufferedReader(new InputStreamReader(process.getInputStream()))) {
                String line;
                while ((line = reader.readLine()) != null) {
                    if (verbose) {
                        printInfo(line);
                    }
                }
            }

            int termination = process.waitFor();
            if (termination != 0) {
                printError("Abnormal FragCast termination: " + termination + ", please run the " +
                        "following command from the command line for more information\n" +
                        String.join(" ", builder.command()));
                System.exit(1);
            }

            File predFile = new File(predFileString);
            if (Files.isReadable(predFile.toPath())) {
                if (verbose) {
                    printInfo("Done generating FragCast predictions");
                }
            } else {
                printError("Cannot find FragCast's output. Please rerun MSBooster");
                System.exit(1);
            }
        } catch (IOException | InterruptedException e) {
            e.printStackTrace();
            System.exit(1);
        }
        if (verbose) {
            long endTime = System.nanoTime();
            long duration = (endTime - startTime);
            printInfo("Model running took " + duration / 1000000 + " milliseconds");
        }

        return predFileString;
    }

    /**
     * Resolve the directory holding FragCast's pretrained ONNX weights, to be exported as
     * {@code FRAGCAST_MODEL_DIR}. An explicit {@link Constants#FragCastModelDir} param takes
     * precedence; otherwise a {@code pretrained_models} folder is derived from the FragCast
     * executable's location: first beside the exe, then one level up (so an exe at
     * {@code tools/FragCast/windows/fragcast.exe} finds {@code tools/FragCast/pretrained_models}).
     * Each candidate must actually contain {@code FragCast-RT.onnx}. Returns {@code null} when
     * nothing resolves, leaving FragCast to fall back to its own bundled-model lookup.
     */
    private static String resolveModelDir() {
        if (Constants.FragCastModelDir != null && !Constants.FragCastModelDir.isEmpty()) {
            return Constants.FragCastModelDir;
        }
        if (Constants.FragCast == null || Constants.FragCast.isEmpty()) {
            return null;
        }
        File exeDir = new File(Constants.FragCast).getAbsoluteFile().getParentFile();
        if (exeDir == null) {
            return null;
        }
        List<File> candidates = new ArrayList<>();
        candidates.add(new File(exeDir, "pretrained_models"));
        File parent = exeDir.getParentFile();
        if (parent != null) {
            candidates.add(new File(parent, "pretrained_models"));
        }
        for (File c : candidates) {
            if (new File(c, "FragCast-RT.onnx").isFile()) {
                return c.getPath();
            }
        }
        return null;
    }
}
