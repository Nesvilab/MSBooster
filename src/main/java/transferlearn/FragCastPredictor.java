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
import allconstants.FragCastCharges;
import allconstants.FragCastModels;
import allconstants.FragCastWeights;
import mainsteps.MainClass;
import mainsteps.ParameterUtils;
import modelcallers.FragCastModelCaller;
import modelcallers.FragCastProcess;
import peptideptmformatting.PTMhandler;
import peptideptmformatting.PeptideFormatter;
import readers.datareaders.FastaReader;
import speclib.io.ParquetToSpecLib;
import transferlearn.fragcast.FragCastModelBundle;
import utils.ProgressReporter;
import utils.ProteinLabels;

import java.io.BufferedWriter;
import java.io.File;
import java.nio.file.Files;
import java.sql.Connection;
import java.sql.DriverManager;
import java.sql.ResultSet;
import java.sql.SQLException;
import java.sql.Statement;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.List;
import java.util.Locale;

import static transferlearn.Helpers.convertParquetToLibraryTsv;
import static transferlearn.Helpers.mapProteinsListToGenes;
import static utils.Print.printError;
import static utils.Print.printInfo;

/**
 * Predict a spectral library on this machine, with FragCast.
 *
 * <p>This is {@link Predictor}'s local counterpart, and deliberately its mirror image: the same
 * flags, the same peptide sources, the same output files, the same {@code --model} zip. What differs
 * is everything on the other side of the flags. {@link Predictor} uploads peptides to a prediction
 * server, waits in its queue and downloads the result; this runs the FragCast executable in a
 * subprocess. No URL, no API key, nothing leaves the machine. Between them the two cover the same
 * job, so a caller chooses a backend by choosing a class name and, for the server, adding
 * credentials.
 *
 * <p>One difference is real rather than cosmetic, and it comes from FragCast rather than from
 * here: there is <b>no MGF</b>. FragCast writes {@code tsv} and {@code parquet} only, and the MGF
 * the server returns is written by the server itself - there is no parquet-to-MGF converter on this
 * side to stand in for it. Neither predictor takes a per-property flag: one {@code build-library}
 * run predicts retention time, ion mobility and MS2 together, and all three are kept.
 *
 * <p>Protein labels come from the input either way. A headered FragCast input may carry a
 * {@code proteins} column, so {@link #writeFragCastInput} passes the peptide list's straight through
 * and {@code PinReader.createFragCastList} carries the pin files' own, exactly as the AlphaPeptDeep
 * writer beside it does. FragCast fills {@code ProteinId}, {@code GeneName},
 * {@code AllMappedProteins}, {@code AllMappedGenes} and {@code Proteotypic} from it, and nothing here
 * labels the library afterwards.
 */
public class FragCastPredictor {
    /** Flags that belong to the server workflow, refused here rather than quietly ignored. */
    private static final List<String> SERVER_FLAGS = Arrays.asList(
            "--url", "--urlTransferLearn", "--urlPredict", "--api-key");

    private static void errorMessage() {
        printError("Usage: java -cp MSBooster.jar transferlearn.FragCastPredictor " +
                "--paramsList <msbooster parameters> " +
                "--output-format <parquet|librarytsv|speclib> " +
                "optional: --model <model.zip> --output-dir <dir> --basename <stem> --fasta <fasta> " +
                "--peptide-list-to-predict <peptides.parquet> --min-charge <int> --max-charge <int> " +
                "--decoy-tag <tag> --keep-decoys <0|1> --custom-mods <spec>");
        printError("Example: java -cp MSBooster.jar transferlearn.FragCastPredictor " +
                "--paramsList msbooster_params.txt --model FragCast-finetuned.zip " +
                "--peptide-list-to-predict peptide_list.parquet --fasta proteome.fasta " +
                "--output-format speclib --basename fragpipe-predicted-speclib");
        printError("For server-based prediction, use transferlearn.Predictor instead.");
        System.exit(1);
    }

    public static void main(String[] args) throws Exception {
        Locale.setDefault(Locale.US);
        Thread.setDefaultUncaughtExceptionHandler((t, e) -> {
            e.printStackTrace();
            System.exit(1);
        });

        if (args.length % 2 != 0) {
            printError("Malformed arguments, args of length " + args.length);
            errorMessage();
        }

        String params = "";
        String model = "";
        String peptideList = "";
        String fasta = "";
        String outputDir = "";
        String basename = "";
        String outputFormat = "parquet";
        String keepDecoys = "1"; //as Predictor defaults it, and as Constants.keepDecoys defaults
        String decoyTag = "rev_";
        int minCharge = 2;
        int maxCharge = 3;

        for (int i = 0; i < args.length; i += 2) {
            if (SERVER_FLAGS.contains(args[i])) {
                printError(args[i] + " belongs to the server-based prediction in " +
                        "transferlearn.Predictor. This workflow runs entirely on this machine and " +
                        "contacts no server.");
                System.exit(1);
            }
            switch (args[i]) {
                case "--paramsList":
                    params = args[i + 1];
                    break;
                case "--model":
                    model = args[i + 1];
                    break;
                case "--peptide-list-to-predict":
                    peptideList = args[i + 1];
                    //without this the run has no output directory of its own: requirePinMzml is
                    //false on this path, so updateOutputDirectory falls back to this field and an
                    //empty one would put spectraRT at the root of the filesystem
                    Constants.peptideListDirectory = new File(args[i + 1]).getParent();
                    break;
                case "--fasta":
                    fasta = args[i + 1];
                    break;
                case "--output-dir":
                    outputDir = args[i + 1];
                    File directory = new File(outputDir);
                    if (!directory.exists()) {
                        directory.mkdirs();
                    }
                    break;
                case "--basename":
                    basename = args[i + 1];
                    break;
                case "--output-format":
                    outputFormat = args[i + 1];
                    break;
                case "--keep-decoys":
                    keepDecoys = args[i + 1];
                    break;
                case "--decoy-tag":
                    decoyTag = args[i + 1];
                    break;
                case "--min-charge":
                    minCharge = Integer.parseInt(args[i + 1]);
                    break;
                case "--max-charge":
                    maxCharge = Integer.parseInt(args[i + 1]);
                    break;
                case "--custom-mods":
                    //the mod tables are rebuilt, not just written out: getBase below resolves a
                    //peptide's modifications through them, so a custom mod absent from them would
                    //reach FragCast as an unrecognised delta mass
                    Helpers.customModsStringToTsv(args[i + 1]);
                    PTMhandler.unimodToModMassAll = PTMhandler.makeUnimodToModMassAll(false);
                    PTMhandler.AAunimodToModMassAll = PTMhandler.makeUnimodToModMassAll(true);
                    PTMhandler.AAunimodToModMassAllKeys = PTMhandler.AAunimodToModMassAll.keySet();
                    PTMhandler.aamassToAlphapeptdeep = PTMhandler.makeModAAmassToAlphapeptdeep();
                    break;
                default:
                    printError("Unknown argument " + args[i]);
                    errorMessage();
                    break;
            }
        }

        if (params.isEmpty()) {
            printError("--paramsList is required.");
            errorMessage();
        }
        if (outputFormat.equals("mgf")) {
            printError("FragCast writes tsv and parquet only, so it cannot produce an MGF library. " +
                    "Choose parquet, librarytsv or speclib, or predict with transferlearn.Predictor " +
                    "against a prediction server.");
            System.exit(1);
        }
        final boolean annotated = outputFormat.equals("speclib") || outputFormat.equals("librarytsv");
        if (!annotated && !outputFormat.equals("parquet")) {
            printError("Unknown output format " + outputFormat + ". Must be one of parquet, " +
                    "librarytsv, or speclib.");
            System.exit(1);
        }
        if (annotated && fasta.isEmpty()) {
            printError("--fasta is required for " + outputFormat + " output.");
            System.exit(1);
        }

        printInfo(Constants.versionNumber + " FragCast Local Library Prediction");

        //Read before createInputFiles because the pin path hands MainClass explicit model overrides
        //- whatever follows --paramsList wins - so the file's own fast choice has to be captured
        //here and forwarded, not recovered from Constants afterwards.
        final boolean fastRequested = requestsFastSpecModel(params);

        //Everything the command line names is handed on as a command-line override rather than
        //written into Constants here. processCommandLineInputs applies the parameter file the moment
        //it reads --paramsList and whatever follows overrides it, so assigning a static first would
        //simply be undone: the shipped template carries an empty "FragCastModelZip =" line, and an
        //empty value is a value.
        final File inputFile = createInputFiles(peptideList, params, keepDecoys, model, decoyTag,
                minCharge, maxCharge, fastRequested);
        if (outputDir.isEmpty()) {
            outputDir = inputFile.getParent();
        }
        if (basename.isEmpty()) {
            basename = model.isEmpty() ? "fragcast_predict" : stripZip(new File(model).getName());
        }

        //the weights the bundle supplied, or the pretrained ones when no bundle was named
        final FragCastWeights weights = FragCastWeights.fromConstants();
        //what the file asked for, not Constants: the pin path just overrode Constants.spectraModel
        //with what was forwarded, and the peptide-list path applied the file's own casing - the
        //pre-read value is the one answer both paths agree on
        final boolean fast = fastRequested;
        printInfo("Predicting with " + weights);

        //Clear last run's outputs first. FragCast's COPY overwrites, but a library left from an
        //earlier run under a format this one is not writing would sit there looking current, and
        //the server path clears the same files for the same reason.
        for (String extension : new String[]{".parquet", ".tsv", ".speclib"}) {
            final File stale = new File(outputDir, basename + extension);
            if (stale.isFile() && !stale.delete()) {
                printError("Could not remove the previous " + stale + "; delete it and run again.");
                System.exit(1);
            }
        }

        //One input file, so one FragCast run. The server path splits its uploads at 500000 peptides
        //because the server takes one job per upload; a local subprocess has no such ceiling, and
        //splitting would only make the model be loaded more than once.
        final File merged = new File(outputDir, basename + ".parquet");
        printInfo("Predicting into " + merged);
        final FragCastProcess.Result result = FragCastModelCaller.predictLibrary(
                inputFile.getAbsolutePath(), merged.getAbsolutePath(), fast, weights);
        if (!result.succeeded() || !merged.isFile()) {
            printError("FragCast did not produce a library for " + inputFile +
                    " (exit code " + result.exitCode + ")");
            //everything FragCast said is already on the console, streamed as it ran - repeating the
            //tail here would only print it twice. What it never said is what the code means.
            final String diagnosis = result.diagnosis();
            if (!diagnosis.isEmpty()) {
                printError(diagnosis);
            }
            System.exit(1);
        }

        //From here on this run takes the same code the server path runs: the library is put into the
        //schema everything downstream reads, which is the server's plus the two columns only FragCast
        //can fill. The converter reads whatever its input holds, so one path serves both.
        toLibrarySchema(merged);

        if (!annotated) {
            printInfo("Total parquet file at " + merged.getAbsolutePath());
            System.exit(0);
        }

        //The gene names come from the FASTA, not from the predictor - neither backend returns them.
        printInfo("Reading FASTA file");
        final FastaReader fastaReader = new FastaReader(fasta);
        fastaReader.mapProtToGene();
        //With --peptide-list the labels are in that file; without it the peptides came from this
        //job's pin files, and the only place carrying their proteins is the input written for
        //FragCast. The server path branches the same way, on its own generated input.
        final HashMap<String, String> protMap = mapProteinsListToGenes(
                peptideList.isEmpty() ? inputFile.getAbsolutePath() : peptideList,
                fastaReader.protToGene);

        if (outputFormat.equals("librarytsv")) {
            final String tsv = new File(outputDir, basename + ".tsv").getAbsolutePath();
            try (Connection conn = DriverManager.getConnection("jdbc:duckdb:")) {
                convertParquetToLibraryTsv(merged.getAbsolutePath(), tsv, protMap, conn);
            }
            printInfo("Total file at " + tsv);
        } else {
            printInfo("Converting parquet to speclib format");
            final String speclib = new File(outputDir, basename + ".speclib").getAbsolutePath();
            new ParquetToSpecLib(merged.getAbsolutePath(), protMap, -3, true, true, true)
                    .convertAndWrite(speclib);
            printInfo("Total speclib file at " + speclib);
        }
        System.exit(0);
    }

    /**
     * The peptides to predict, written in FragCast's input dialect.
     *
     * <p>Two sources, the same two {@link Predictor} accepts. With a peptide list, every precursor is
     * already named and only needs reformatting. Without one, the peptides are whatever this run's
     * pin files hold, which is what {@code MainClass} already knows how to collect - so it is asked
     * to collect them and stop, exactly as the server path asks it to.
     */
    private static File createInputFiles(String peptideList, String params, String keepDecoys,
                                           String model, String decoyTag, int minCharge, int maxCharge,
                                           boolean fast)
            throws Exception {
        if (peptideList.isEmpty()) {
            MainClass.main(pinPathArgs(params, keepDecoys, model, decoyTag, fast));
            return new File(Constants.spectraRTPrefix + ".tsv");
        }

        final HashMap<String, String> paramsMap = new HashMap<>();
        ParameterUtils.processCommandLineInputs(peptideListArgs(params, model, decoyTag), paramsMap);
        ParameterUtils.updateConstants(paramsMap);
        FragCastModelBundle.applyToConstants();

        final File out = new File(Constants.spectraRTPrefix + ".tsv");
        writeFragCastInput(peptideList, out, minCharge, maxCharge);
        return out;
    }

    /**
     * Arguments for the run that collects this job's peptides out of its pin files.
     *
     * <p>Every one of these is also a parameter a file can name, and the file is applied the moment
     * {@code processCommandLineInputs} reads {@code --paramsList}. Whatever follows it wins, which is
     * why they are passed rather than assigned to {@link Constants} beforehand: the shipped template
     * carries an empty {@code FragCastModelZip =} line, and an empty value is a value.
     */
    static String[] pinPathArgs(String params, String keepDecoys, String model, String decoyTag,
                                boolean fast) {
        //The collection run must be told a FragCast model whatever the file names, but which Spec
        //variant only matters to the prediction that follows - so the file's own fast choice is
        //forwarded rather than flattened to the Conformer, keeping Constants.spectraModel saying
        //which weights this run actually loads. RT and IM stay under the plain name: only the MS2
        //weights come in a fast variant.
        final List<String> args = new ArrayList<>(Arrays.asList(
                "--paramsList", params,
                "--keepDecoys", keepDecoys,
                "--decoyPrefix", decoyTag,
                "--createPredFileOnly", "true",
                "--spectraModel", fast ? FragCastModels.FAST : FragCastModels.CONFORMER,
                "--rtModel", FragCastModels.CONFORMER,
                "--imModel", FragCastModels.CONFORMER));
        addModelZip(args, model);
        return args.toArray(new String[0]);
    }

    /**
     * Does the parameter file ask for FragCast's small/fast Spec model
     * ({@code spectraModel = FragCast-Fast})? That line is how a caller with no command line into
     * FragCast - FragPipe's transfer-learning step above all - picks the Spec variant, for the
     * prediction here and for {@link FragCastLocalWorkflow}'s choice of base weights alike.
     */
    static boolean requestsFastSpecModel(String paramsList) throws java.io.IOException {
        final HashMap<String, String> fileParams = new HashMap<>();
        ParameterUtils.processCommandLineInputs(new String[]{"--paramsList", paramsList}, fileParams);
        //case-insensitively, as a rescoring run reads the same line through LowercaseModelMapper -
        //the two entry points must not disagree about what a hand-written file means
        final String spectraModel = fileParams.get("spectraModel");
        return spectraModel != null && FragCastModels.FAST.equalsIgnoreCase(spectraModel.trim());
    }

    /** As above, for the run that reads a peptide list instead of this job's pin files. */
    static String[] peptideListArgs(String params, String model, String decoyTag) {
        final List<String> args = new ArrayList<>(Arrays.asList(
                "--paramsList", params,
                "--requirePinMzml", "false",
                "--decoyPrefix", decoyTag));
        addModelZip(args, model);
        return args.toArray(new String[0]);
    }

    /**
     * Name the bundle only when there is one. Passing an empty {@code --FragCastModelZip} would
     * override a path the parameter file had legitimately set with nothing.
     */
    private static void addModelZip(List<String> args, String model) {
        if (!model.isEmpty()) {
            args.add("--FragCastModelZip");
            args.add(model);
        }
    }

    /**
     * Turn a peptide list into the two columns FragCast reads.
     *
     * <p>Everything is written to one file. The 500000-peptide split {@link PredictUtils} performs is
     * there because the prediction server takes one upload per job; a local subprocess has no such
     * ceiling, and splitting would only multiply the number of times the model is loaded.
     *
     * <p>Rows whose peptide already ends in a charge keep it. Rows that do not are expanded across
     * {@code minCharge..maxCharge}, which is what those flags mean on the server side too.
     */
    static void writeFragCastInput(String peptideListParquet, File out, int minCharge, int maxCharge)
            throws SQLException, java.io.IOException {
        FragCastCharges.reset();
        final File parent = out.getAbsoluteFile().getParentFile();
        if (parent != null && !parent.isDirectory()) {
            parent.mkdirs();
        }

        try (Connection conn = DriverManager.getConnection("jdbc:duckdb:");
             Statement stmt = conn.createStatement();
             BufferedWriter writer = Files.newBufferedWriter(out.toPath())) {

            final ResultSet count = stmt.executeQuery(
                    "SELECT COUNT(*) FROM read_parquet('" + peptideListParquet.replace("'", "''") + "')");
            count.next();
            final int rows = count.getInt(1);
            printInfo(rows + " peptides to convert");
            final ProgressReporter pr = new ProgressReporter(rows);

            //The same three things the AlphaPeptDeep input carries, in FragCast's own dialect:
            //proteins is what it fills ProteinId, AllMappedProteins and Proteotypic from, and
            //is_decoy is what makes it leave the decoys out of the library, as the server does.
            writer.write("peptide" + "\t" + "charge" + "\t" + "proteins" + "\t" +
                    "is_decoy" + "\n");
            try (ResultSet rs = stmt.executeQuery("SELECT peptide, proteins, is_decoy FROM read_parquet('" +
                    peptideListParquet.replace("'", "''") + "')")) {
                while (rs.next()) {
                    String peptide = rs.getString("peptide");
                    final String proteins = rs.getString("proteins");
                    final boolean isDecoy = rs.getBoolean("is_decoy");
                    if (peptide == null || peptide.isEmpty()) {
                        pr.progress();
                        continue;
                    }
                    //the list carries the charge glued to the end of the peptide, as the server path
                    //also finds it
                    String charge = "";
                    while (!peptide.isEmpty() && Character.isDigit(peptide.charAt(peptide.length() - 1))) {
                        charge = peptide.charAt(peptide.length() - 1) + charge;
                        peptide = peptide.substring(0, peptide.length() - 1);
                    }
                    if (charge.isEmpty()) {
                        for (int z = minCharge; z <= maxCharge; z++) {
                            writeOne(writer, peptide, String.valueOf(z), proteins, isDecoy);
                        }
                    } else {
                        writeOne(writer, peptide, charge, proteins, isDecoy);
                    }
                    pr.progress();
                }
            }
        }
        FragCastCharges.reportSkipped("precursors");
    }

    private static void writeOne(BufferedWriter writer, String peptide, String charge, String proteins,
                                 boolean isDecoy) throws java.io.IOException {
        if (!FragCastCharges.canPredict(charge)) {
            return;
        }
        //getBase carries the modifications as delta masses, which is what FragCast resolves against
        //its own UniMod table
        final PeptideFormatter pf = new PeptideFormatter(peptide, charge, "apdpred");
        //normalised like every other protein column MSBooster writes: a peptide list is as free to
        //terminate its labels with a ';' as MSFragger's pin is, and FragCast reads one either way as
        //"this peptide is shared" (see ProteinLabels)
        writer.write(pf.getBase() + "\t" + charge + "\t" +
                ProteinLabels.normalize(proteins) + "\t" + isDecoy + "\n");
    }

    /**
     * Put a FragCast library into the schema everything downstream reads.
     *
     * <p>That is the schema the AlphaPeptDeep server returns, plus the two columns only FragCast can
     * fill: {@code AverageExperimentalRetentionTime} and {@code AllMappedProteins}. Keeping them is
     * what makes a predicted library carry the same protein-group information as the experimental one
     * FragSpecLib writes beside it - {@code ProteinId} names the leading protein by convention, so
     * {@code AllMappedProteins} is the only place the rest of a shared peptide's group survives.
     * Their position matches that library too: both sit after {@code FragmentLossType} and before
     * {@code Proteotypic}. The server path simply has fewer columns, and the converter reads whatever
     * its input holds rather than being told which predictor produced it.
     *
     * <p>FragCast's two GENE columns are dropped instead, because they are not gene names: it reads
     * them off the protein header, so {@code sp|Q7Z4L5|TT21B_HUMAN} yields the entry mnemonic
     * {@code TT21B} where the gene is {@code TTC21B}.
     * {@link Helpers#convertParquetToLibraryTsv} puts both back from the FASTA, which is the gene name
     * a user asked for. Narrowing the shared columns is unrelated housekeeping: FragCast writes them
     * wider than the server does.
     */
    static void toLibrarySchema(File parquet) throws SQLException {
        final File projected = new File(parquet.getAbsolutePath() + ".apd.parquet");
        try (Connection conn = DriverManager.getConnection("jdbc:duckdb:");
             Statement stmt = conn.createStatement()) {
            stmt.execute("COPY (SELECT " +
                    "CAST(PrecursorMz AS FLOAT) AS PrecursorMz, " +
                    "CAST(ProductMz AS FLOAT) AS ProductMz, " +
                    "Annotation, " +
                    "ProteinId, " +
                    "PeptideSequence, " +
                    "ModifiedPeptideSequence, " +
                    "CAST(PrecursorCharge AS SMALLINT) AS PrecursorCharge, " +
                    "CAST(LibraryIntensity AS FLOAT) AS LibraryIntensity, " +
                    "CAST(NormalizedRetentionTime AS FLOAT) AS NormalizedRetentionTime, " +
                    "CAST(PrecursorIonMobility AS VARCHAR) AS PrecursorIonMobility, " +
                    "FragmentType, " +
                    "CAST(FragmentCharge AS SMALLINT) AS FragmentCharge, " +
                    "CAST(FragmentSeriesNumber AS SMALLINT) AS FragmentSeriesNumber, " +
                    "FragmentLossType, " +
                    "CAST(AverageExperimentalRetentionTime AS VARCHAR) AS AverageExperimentalRetentionTime, " +
                    "AllMappedProteins, " +
                    "CAST(Proteotypic AS SMALLINT) AS Proteotypic " +
                    "FROM read_parquet('" + sql(parquet.getAbsolutePath()) + "')) TO '" +
                    sql(projected.getAbsolutePath()) + "' (FORMAT PARQUET)");
        } catch (SQLException e) {
            projected.delete();
            throw e;
        }
        if (!parquet.delete() || !projected.renameTo(parquet)) {
            throw new SQLException("Could not replace " + parquet + " with the converted library at " +
                    projected);
        }
    }

    private static String stripZip(String name) {
        final int i = name.lastIndexOf(".zip");
        return i == -1 ? name : name.substring(0, i);
    }

    /** Single quotes are the only thing that can break out of a DuckDB string literal. */
    private static String sql(String path) {
        return path.replace("\\", "/").replace("'", "''");
    }
}
