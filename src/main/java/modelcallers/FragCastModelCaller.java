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
 */
public class FragCastModelCaller {
    public static String callModel(String inputFile, boolean verbose) {
        long startTime = System.nanoTime();
        // spectraRT.tsv -> spectraRT.predicted.parquet (strip the trailing ".tsv")
        String predFileString = inputFile.substring(0, inputFile.length() - 4) + ".predicted.parquet";
        try {
            if (Constants.FragCast == null) {
                printError("path to FragCast executable must be provided (set the FragCast parameter)");
                System.exit(1);
            }
            if (verbose) {
                printInfo("Generating FragCast predictions");
            }

            List<String> command = new ArrayList<>();
            command.add(Constants.FragCast);
            command.add("--task");
            command.add("build-library");
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

            ProcessBuilder builder = new ProcessBuilder(command);
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
