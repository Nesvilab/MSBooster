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
import allconstants.FragCastWeights;

import java.io.File;
import java.io.IOException;
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
 *
 * <p>{@link FragCastWeights} extends that naming rule to custom per-property weights, which is how
 * a locally fine-tuned model is used for prediction (see {@link mainsteps.FragCastTransferLearner}).
 * The pretrained weights carry an empty tag, so a job that fine-tunes nothing writes exactly the
 * files it always did.
 */
public class FragCastModelCaller {
    /**
     * The output library path for an input peptide file: {@code spectraRT.tsv} ->
     * {@code spectraRT.predicted.parquet}, or {@code spectraRT.fast.predicted.parquet} under
     * {@code fast}, plus {@link FragCastWeights#fileTag()} for custom weights. Nothing predicted
     * under different weights may share a path — the name is what says which weights wrote a
     * library, so neither the two Spec variants nor a fine-tuned pass and the pretrained one can
     * overwrite, or be mistaken for, each other.
     */
    static String predictionFile(String inputFile, boolean fast, FragCastWeights weights) {
        return inputFile.substring(0, inputFile.length() - 4) //strip the trailing ".tsv"
                + (fast ? ".fast" : "") + weights.fileTag() + ".predicted.parquet";
    }

    /**
     * The FragCast command line. Every flag here exists in FragCast's {@code build-library} task;
     * {@code --fast} swaps the Conformer Spec weights for the small/fast ones and is FragCast's only
     * way to choose between them (there is no per-property task to call). A
     * {@code --rt-onnx}/{@code --im-onnx}/{@code --spec-onnx} flag follows for each property with
     * custom weights, appended last so the flags FragCast has always been given keep their positions.
     */
    static List<String> buildCommand(String inputFile, String predFileString, boolean fast,
                                     FragCastWeights weights) {
        List<String> command = new ArrayList<>();
        command.add(Constants.FragCast);
        command.add("--task");
        command.add("build-library");
        //--fast is shorthand for one bundled weights file, and FragCast refuses it alongside
        //--spec-onnx because the two would both be naming the Spec weights. Named weights win: the
        //ONNX itself declares which architecture it is, so a fine-tuned fast model still runs fast
        if (fast && !weights.supersedesFastFlag(true)) {
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
        //The decoy prefix is on the database token of a protein label (rev_sp|P1|A_HUMAN), so reading
        //the accession out of it drops the prefix and a decoy is labelled exactly like the target it
        //was reversed from. Naming it here is what lets FragCast keep it, as the AlphaPeptDeep server
        //keeps its own decoy_tag - the two backends must not disagree about which entries are decoys.
        command.add("--decoy-tag");
        command.add(Constants.decoyPrefix);
        //no --min-charge/--max-charge: every row of the list names its own charge, so FragCast
        //predicts exactly those precursors and the range would only be reported as overridden
        weights.appendFlags(command);
        return command;
    }

    /**
     * Predict one library and hand back what FragCast reported, without judging it.
     *
     * <p>{@link #callModel} is the rescoring path: a failure there is fatal, so it exits. This
     * returns the outcome instead, so a caller can treat a failed prediction as "do not trust this
     * model" rather than "abandon the run".
     *
     * <p>FragCast's own output is echoed as it arrives. This builds the run's product and is the
     * longest thing in it - hours, on a peptide list digested from a FASTA - and FragCast is the only
     * one of the two that knows how far along it is. Captured and shown only on failure, as it used
     * to be, the size of a run was reported solely in its post-mortem: a list of 15.8 million
     * precursors looked exactly like one of 61 thousand until the four hours were already spent.
     */
    public static FragCastProcess.Result predictLibrary(String inputFile, String outputFile,
                                                        boolean fast, FragCastWeights weights)
            throws IOException, InterruptedException {
        return FragCastProcess.run(buildCommand(inputFile, outputFile, fast, weights), true);
    }

    public static String callModel(String inputFile, boolean verbose, boolean fast, FragCastWeights weights) {
        long startTime = System.nanoTime();
        String predFileString = predictionFile(inputFile, fast, weights);
        try {
            if (Constants.FragCast == null) {
                printError("path to FragCast executable must be provided (set the FragCast parameter)");
                System.exit(1);
            }
            if (verbose) {
                printInfo("Generating FragCast predictions" + (fast ? " (small/fast Spec model)" : "") +
                        (weights.isBase() ? "" : " using " + weights));
            }

            List<String> command = buildCommand(inputFile, predFileString, fast, weights);
            FragCastProcess.Result result = FragCastProcess.run(command, verbose);
            if (!result.succeeded()) {
                printError("Abnormal FragCast termination: " + result.exitCode + ", please run the " +
                        "following command from the command line for more information\n" +
                        String.join(" ", command));
                if (!result.diagnosis().isEmpty()) {
                    printError(result.diagnosis());
                }
                //a quiet call captured its output instead of streaming it, so this is the only place
                //it can be seen; a verbose one already printed every line as it arrived
                if (!verbose) {
                    printError(result.tail(10));
                }
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
}
