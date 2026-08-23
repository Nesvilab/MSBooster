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

package transferlearn;

import allconstants.Constants;
import allconstants.FragCastModels;
import allconstants.FragCastWeights;
import mainsteps.FragCastTransferLearner;
import mainsteps.ParameterUtils;
import peptideptmformatting.PTMhandler;
import utils.Print;

import java.io.File;
import java.util.Arrays;
import java.util.HashMap;
import java.util.List;
import java.util.Locale;

import static allconstants.Constants.versionNumber;
import static utils.Print.printInfo;

/**
 * Transfer learning without a server: the local counterpart to {@link SecondPassWorkflow}.
 *
 * <p>{@link SecondPassWorkflow} sends a spectral library to a training server, waits for a job to
 * come back, and then sends peptides to a prediction server. This entry point does the same job with
 * neither: FragCast fine-tunes on the CPU of the machine it is running on, and predicts on it too.
 * No URL, no API key, no queue, nothing leaves the machine.
 *
 * <p>What it learns from is always a spectral library you supply with {@code --library}, as TSV, CSV
 * or Parquet, carrying the column names {@code FragCastPredictor}
 * documents. One library is all it needs: FragCast holds out its own slice of it to select the best
 * epoch against.
 *
 * <p>It fine-tunes and stops. The product is one zip holding every model it wrote, named in the
 * output as {@code FragCastModelZip}; rescoring with it is a separate run, told to load that zip.
 * Splitting the two is what lets FragPipe queue them as the separate steps they are.
 *
 * <p>Every model FragCast exports is used: nothing here withholds one for scoring worse than the
 * model it started from. Nothing about how the fine-tune runs is settable: every property FragCast
 * can be trained for is trained, and the only thing the caller chooses is where the result goes.
 */
public class FragCastLocalWorkflow {
    /** Flags that belong to the server workflow, refused here rather than quietly ignored. */
    private static final List<String> SERVER_FLAGS = Arrays.asList(
            "--url", "--urlTransferLearn", "--urlPredict", "--api-key");

    private static void errorMessage() {
        Print.printError("Usage: java -cp MSBooster.jar transferlearn.FragCastLocalWorkflow " +
                "--paramsList <msbooster parameters> " +
                "--library <spectral library to train on> " +
                "optional: --output-dir <output directory> --basename <model file stem> " +
                "--custom-mods <spec>");
        Print.printError("Example: java -cp MSBooster.jar transferlearn.FragCastLocalWorkflow " +
                "--paramsList msbooster_params.txt --library library.tsv " +
                "(fine-tunes on the library and prints the bundle it produced)");
        Print.printError("For server-based transfer learning, use transferlearn.SecondPassWorkflow instead.");
        System.exit(1);
    }

    public static void main(String[] args) throws Exception {
        Locale.setDefault(Locale.US);
        Thread.setDefaultUncaughtExceptionHandler((t, e) -> {
            e.printStackTrace();
            System.exit(1);
        });

        if (args.length % 2 != 0) {
            Print.printError("Malformed arguments, args of length " + args.length);
            errorMessage();
        }

        String params = "";
        String library = "";
        String outputDir = "";
        String basename = "";

        for (int i = 0; i < args.length; i += 2) {
            if (SERVER_FLAGS.contains(args[i])) {
                Print.printError(args[i] + " belongs to the server-based transfer learning in " +
                        "transferlearn.SecondPassWorkflow. This workflow runs entirely on this " +
                        "machine and contacts no server.");
                System.exit(1);
            }
            switch (args[i]) {
                case "--paramsList":
                    params = args[i + 1];
                    break;
                case "--library":
                    library = args[i + 1];
                    break;
                case "--output-dir":
                    outputDir = args[i + 1];
                    break;
                case "--basename":
                    basename = args[i + 1];
                    break;
                case "--custom-mods":
                    //FragCast resolves a peptide's modifications through these tables, so a mod that
                    //is not in them reaches it as an unrecognised delta mass and its precursor is
                    //dropped from the training library rather than learned from
                    Helpers.customModsStringToTsv(args[i + 1]);
                    PTMhandler.unimodToModMassAll = PTMhandler.makeUnimodToModMassAll(false);
                    PTMhandler.AAunimodToModMassAll = PTMhandler.makeUnimodToModMassAll(true);
                    PTMhandler.AAunimodToModMassAllKeys = PTMhandler.AAunimodToModMassAll.keySet();
                    PTMhandler.aamassToAlphapeptdeep = PTMhandler.makeModAAmassToAlphapeptdeep();
                    break;
                default:
                    Print.printError("Unknown argument " + args[i]);
                    errorMessage();
                    break;
            }
        }

        if (params.isEmpty()) {
            Print.printError("--paramsList is required.");
            errorMessage();
        }
        if (library.isEmpty()) {
            Print.printError("--library is required: it is the spectral library FragCast is " +
                    "fine-tuned on.");
            errorMessage();
        }

        printInfo(versionNumber + " FragCast Local Transfer Learning");
        fineTuneOnLibrary(params, library, outputDir, basename);
    }

    /**
     * The library path: fine-tune on a library the caller already has and stop. Nothing is rescored,
     * because a library carries no pin files - the output is the models, and their paths are printed
     * so a later run can name them.
     */
    private static void fineTuneOnLibrary(String params, String library,
                                          String outputDir, String basename) throws Exception {
        HashMap<String, String> paramsMap = new HashMap<>();
        ParameterUtils.processCommandLineInputs(
                new String[] {"--paramsList", params, "--requirePinMzml", "false"}, paramsMap);
        ParameterUtils.updateConstants(paramsMap);

        if (Constants.FragCast == null || Constants.FragCast.isEmpty()) {
            Print.printError("The FragCast parameter must point at the FragCast executable.");
            System.exit(1);
        }
        if (!new File(library).isFile()) {
            Print.printError("No such library: " + library);
            System.exit(1);
        }

        File dir = resolveOutputDirectory(outputDir, library);
        if (!dir.isDirectory() && !dir.mkdirs()) {
            Print.printError("Could not create the output directory " + dir.getPath() +
                    ". Name a writable one with --output-dir.");
            System.exit(1);
        }
        printInfo("Writing any fine-tuned models to " + dir.getPath());

        FragCastWeights before = FragCastWeights.fromConstants();
        //no pin files here, so the job's models are whatever the parameter file selects: the MS2
        //task starts from the fast weights exactly when that is the Spec model the user predicts with
        FragCastWeights accepted = FragCastTransferLearner.fineTuneAll(before, new File(library),
                dir, FragCastModels.usesFastSpecModel(Constants.spectraModel));

        if (accepted.equals(before)) {
            //every runnable task was skipped, failed, or left no usable weights, so there is no
            //new ONNX file to name in a later run's parameters
            Print.printError("No fine-tuned model was produced: every task was skipped, failed, or " +
                    "finished without usable weights. The pretrained FragCast models are unchanged.");
            System.exit(0);
        }
        //the one file this workflow exists to produce: three ONNX models are what FragCast wants,
        //one zip is what a person wants to keep, copy, or hand to FragPipe's MSBooster panel
        File bundle = FragCastTransferLearner.writeBundle(accepted, dir,
                basename.isEmpty() ? FragCastTransferLearner.DEFAULT_BUNDLE_BASENAME : basename);

        printInfo("Add this to your msbooster_params.txt to rescore with the fine-tuned models:");
        if (bundle != null) {
            printInfo("  FragCastModelZip = " + bundle.getAbsolutePath());
            //the loose ONNX files were the ingredients; the bundle is the product, and everything
            //downstream of this workflow opens the bundle. Kept only when there is no bundle to
            //open, which is the one case where they are the only copy of what was produced.
            int removed = FragCastTransferLearner.deleteLooseModels(accepted, before);
            if (removed > 0) {
                printInfo("  (" + removed + " loose ONNX file(s) removed; the bundle holds them)");
            }
        } else {
            //no bundle, so the individual models are the only way to name what was produced
            report("FragCastRtOnnx", accepted.rtOnnx, before.rtOnnx);
            report("FragCastImOnnx", accepted.imOnnx, before.imOnnx);
            report("FragCastSpecOnnx", accepted.specOnnx, before.specOnnx);
        }
        System.exit(0);
    }

    /**
     * Where this run writes. Nothing on this path establishes an output directory of its own: it
     * initialises the parameters with {@code requirePinMzml false}, which leaves
     * {@link Constants#outputDirectory} empty, and {@code new File("", name)} resolves to the root of
     * the filesystem - so the run would deposit its libraries and models in {@code C:ragcast_transfer}
     * or {@code /fragcast_transfer}, or fail to create them at all. The library is the one path the
     * caller definitely named, so its directory is the default.
     */
    static File resolveOutputDirectory(String outputDir, String library) {
        if (!outputDir.isEmpty()) {
            return new File(outputDir);
        }
        File parent = new File(library).getAbsoluteFile().getParentFile();
        return new File(parent, "fragcast_transfer");
    }

    private static void report(String param, String accepted, String before) {
        if (!accepted.equals(before)) {
            printInfo("  " + param + " = " + accepted);
        }
    }

}
