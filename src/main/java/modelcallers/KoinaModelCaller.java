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
import com.google.gson.Strictness;
import com.google.gson.stream.JsonReader;
import com.google.gson.stream.JsonToken;
import features.spectra.MassCalculator;
import koinaclasses.KoinaBisection;
import koinaclasses.KoinaException;
import koinaclasses.KoinaModelCheck;
import koinaclasses.KoinaModelSpec;
import koinaclasses.KoinaTask;
import org.knowm.xchart.XYChart;
import org.knowm.xchart.XYChartBuilder;
import org.knowm.xchart.style.Styler;
import peptideptmformatting.PeptideFormatter;
import predictions.PredictionEntry;
import predictions.PredictionEntryHashMap;
import readers.predictionreaders.KoinaLibReader;
import utils.ProgressReporter;

import java.awt.Font;
import java.io.File;
import java.io.IOException;
import java.io.StringReader;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.concurrent.ConcurrentHashMap;
import java.util.concurrent.ExecutionException;
import java.util.concurrent.Future;
import java.util.concurrent.ScheduledThreadPoolExecutor;
import java.util.concurrent.TimeUnit;
import java.util.concurrent.atomic.AtomicBoolean;
import java.util.concurrent.atomic.AtomicLong;
import java.util.function.Consumer;

import static allconstants.Constants.figureDirectory;
import static figures.ExtensionPlotter.plot;
import static utils.Print.printError;
import static utils.Print.printInfo;

public class KoinaModelCaller {
    private static AtomicBoolean emptyUrl = new AtomicBoolean(false);

    public KoinaModelCaller(){}

    /**
     * Sends every request file of klr's property in jsonFolder to the model, and reads the predictions into klr.
     * A request whose response is not received is sent again (KoinaTask). A request Koina refuses (HTTP 400) is not:
     * KoinaBisection sends its parts until the precursors Koina refuses are found, which klr records and which are
     * left without predictions, and the others are predicted.
     *
     * @throws KoinaException if the run cannot use the model (KoinaModelCheck), if Koina refuses every precursor
     *                        (KoinaBisection), or if a request cannot be predicted: what a task threw if it is a
     *                        KoinaException, any other failure of a task as a KoinaException naming the model. A
     *                        candidate of the best-model search ends the run the same way, since the search cannot
     *                        leave a model out
     */
    //TODO: add functions to clean this up
    public void callModel(String model, KoinaLibReader klr, String jsonFolder, ScheduledThreadPoolExecutor executorService,
                          boolean verbose, boolean makeFigure) {
        if (Constants.KoinaURL.isEmpty() && !emptyUrl.get()) {
            emptyUrl.set(true); //will only be called once
            printError("You must specify a URL to use the Koina service via the --KoinaURL parameter");
            printError("Exiting");
            System.exit(1);
        }
        //no request goes to a model before Koina was asked about it. Writing the request files asked already, so
        //this asks the server only when they were written some other way
        KoinaModelCheck.check(model);
        if (verbose) {
            printInfo("Calling " + model + " model");
        }
        long startTime = System.currentTimeMillis();

        //pass json files to http request
        File[] fileArray = new File(jsonFolder).listFiles();
        if (fileArray == null) { //no request was written
            fileArray = new File[0];
        }
        ArrayList<String> filenameArraylist = new ArrayList<>();
        for (File f : fileArray) {
            String fname = f.toString();
            if (fname.endsWith(klr.property + ".json")) {
                filenameArraylist.add(fname);
            }
        }
        int numProcesses = filenameArraylist.size();

        ProgressReporter pr = new ProgressReporter(numProcesses);
        AtomicLong waitTime;
        if (klr.spec.kind() == KoinaModelSpec.Kind.MS2) { //ms2 or ms2_aux
            waitTime = new AtomicLong(Constants.initialKoinaMillisecondsToWaitMs2);
        } else { //rt or im
            waitTime = new AtomicLong(Constants.initialKoinaMillisecondsToWaitRtIm);
        }

        //the times the request files are done, in the order they are done, from 1
        ConcurrentHashMap<Integer, Long> completionTimes = new ConcurrentHashMap<>();
        //the requests: the request files, then the parts of a refused request that KoinaBisection sends
        List<Request> requests = new ArrayList<>();
        for (String fileName : filenameArraylist) {
            requests.add(new Request(new KoinaTask(fileName, model, klr, waitTime), true));
        }
        KoinaBisection bisection = new KoinaBisection(model, klr, requests.size());
        long nextAllowedSubmissionTime = System.currentTimeMillis();
        long jobStart = System.currentTimeMillis();

        int finishedFiles = 0;
        int unfinished = requests.size();
        int inFlight = 0;
        int jobsSinceLastFailed = 0;
        int notReceived = 0; //attempts whose response was not received
        int allowedFailedAttempts = 0; //finish all tasks that have not been attempted before going to ones that have failed once, etc

        try {
            while (unfinished > 0) {
                //requests is added to while it is gone through; the ones added are gone through in the same pass
                for (int i = 0; i < requests.size(); i++) {
                    Request request = requests.get(i);
                    KoinaTask task = request.task;
                    if (request.done) {
                        continue;
                    }
                    if (request.future == null) {
                        if (task.failedAttempts <= allowedFailedAttempts) {
                            //0.02 seconds between submissions. A request sent again to confirm a refusal waits
                            //its pause on top, without holding up the others
                            long delay = Math.max(20L - (System.currentTimeMillis() - nextAllowedSubmissionTime), 0);
                            request.future = executorService.schedule(task, delay + task.pauseMs(),
                                    TimeUnit.MILLISECONDS);
                            nextAllowedSubmissionTime = System.currentTimeMillis() + delay + 20L;
                            inFlight++;
                        }
                        continue;
                    }
                    if (!request.future.isDone()) {
                        continue;
                    }

                    KoinaTask.Result result = request.future.get();
                    request.future = null;
                    inFlight--;
                    if (result == KoinaTask.Result.NOT_RECEIVED) {
                        notReceived++;
                        if (notReceived == 1) { //once: a slow or failing server would fill the log
                            printInfo("Koina model " + model + " did not answer " + task.label() + " ("
                                    + task.failure() + "). It, and any other request not answered, is sent again "
                                    + "after the others");
                        }
                        executorService.setCorePoolSize((int) Math.max(Math.ceil((double) Constants.numThreads / 5d),
                                (int) Math.floor((double) executorService.getCorePoolSize() * 0.9d)));
                        waitTime.addAndGet(5000);
                        task.failedAttempts++;
                        jobsSinceLastFailed = 0;
                        continue;
                    }

                    //predicted, refused or given up: this request is done
                    request.done = true;
                    task.completed = true;
                    unfinished--;
                    if (request.fromFile) {
                        if (verbose) {
                            pr.progress();
                        }
                        completionTimes.put(++finishedFiles, System.currentTimeMillis() - jobStart);
                    }
                    List<KoinaTask> next;
                    if (result == KoinaTask.Result.REFUSED) {
                        next = bisection.refused(task);
                    } else {
                        if (result == KoinaTask.Result.PREDICTED) {
                            jobsSinceLastFailed++;
                            if (jobsSinceLastFailed >= Constants.numThreads) {
                                executorService.setCorePoolSize(Math.min(Constants.numThreads,
                                        (int) Math.ceil(((double) executorService.getCorePoolSize() * 1.05d))));
                                jobsSinceLastFailed = 0;
                            }
                        }
                        next = bisection.finished(task, result == KoinaTask.Result.PREDICTED);
                    }
                    for (KoinaTask part : next) {
                        requests.add(new Request(part, false));
                        unfinished++;
                    }

                    if (bisection.refusesAll()) {
                        cancel(requests);
                        throw new KoinaException(bisection.refusesAllMessage());
                    }
                }
                if (inFlight == 0) {
                    allowedFailedAttempts++;
                }
                Thread.sleep(100);
            }
        } catch (ExecutionException e) { //a task threw, so the run cannot have this model's predictions
            cancel(requests);
            throw taskFailure(model, e.getCause());
        } catch (InterruptedException e) {
            cancel(requests);
            Thread.currentThread().interrupt();
            throw new KoinaException("Calling Koina model " + model + " was interrupted", e);
        } catch (RuntimeException e) {
            cancel(requests);
            throw e;
        }
        reportRefused(model, klr, bisection, verbose);
        if (notReceived > 1) {
            printInfo("Koina model " + model + " did not answer " + notReceived + " times; each of those requests "
                    + "was sent again");
        }

        long endTime = System.currentTimeMillis();
        long elapsedTime = endTime - startTime;
        if (verbose) {
            printInfo("HTTP request and parse time in milliseconds: " + elapsedTime);
        }

        //create plot
        if (makeFigure) {
            try {
                double[] xData = new double[numProcesses + 1];
                double[] yData = new double[numProcesses + 1];
                xData[0] = 0;
                yData[0] = 0;

                for (int i = 1; i < numProcesses + 1; i++) {
                    xData[i] = completionTimes.get(i) / 1000d;
                    yData[i] = i;
                }
                Arrays.sort(xData);
                XYChart chart = new XYChartBuilder().width(500).height(500).build();
                chart.setTitle("Koina timing");
                chart.setXAxisTitle("Time (sec)");
                chart.setYAxisTitle("Precursors predicted in " + KoinaModelCheck.batchSize(model) + "'s");
                chart.getStyler().setLegendPosition(Styler.LegendPosition.InsideNW);
                chart.getStyler().setYAxisDecimalPattern("0");
                chart.addSeries(model, xData, yData);

                // Increase the font sizes
                chart.getStyler().setAxisTitleFont(new Font("Helvetica", Font.PLAIN, 21));
                chart.getStyler().setLegendFont(new Font("Helvetica", Font.PLAIN, 21));
                chart.getStyler().setAxisTickLabelsFont(new Font("Helvetica", Font.PLAIN, 15));

                if (!new File(figureDirectory).exists()) {
                    new File(figureDirectory).mkdirs();
                }

                plot(chart, figureDirectory + File.separator + "Koina_timing_" + model);
            } catch (Exception e) { //the predictions are all there; only the timing plot is missing
                e.printStackTrace();
            }
        }
    }

    //a request to the model, and where it stands
    private static final class Request {
        final KoinaTask task;
        //a request file, not a part of one that KoinaBisection sends
        final boolean fromFile;
        //the attempt running or waiting to run; null when it is not
        Future<KoinaTask.Result> future;
        boolean done;

        Request(KoinaTask task, boolean fromFile) {
            this.task = task;
            this.fromFile = fromFile;
        }
    }

    //stops the requests still waiting or running once one of them has failed the run
    private static void cancel(List<Request> requests) {
        for (Request request : requests) {
            if (request.future != null) {
                request.future.cancel(true);
            }
        }
    }

    //the most precursors Koina refused that are named one per line
    private static final int REFUSED_LINES = 20;

    /**
     * Says which precursors Koina refused (KoinaBisection), when verbose, and ends the run if it refused every
     * precursor it was sent: the model then predicts nothing (D28). A candidate of the best-model search ends it
     * too, since the search cannot leave a model out.
     */
    private static void reportRefused(String model, KoinaLibReader klr, KoinaBisection bisection,
                                      boolean verbose) {
        List<KoinaBisection.Refusal> refused = bisection.isolated();
        if (refused.isEmpty()) {
            return;
        }
        if (verbose) {
            printInfo("Koina model " + model + " refused " + refused.size() + " precursor"
                    + (refused.size() == 1 ? "" : "s") + " (HTTP 400). The other precursors of their requests are "
                    + "predicted; these have no prediction and are scored as precursors the model cannot take:");
            for (int i = 0; i < Math.min(REFUSED_LINES, refused.size()); i++) {
                KoinaBisection.Refusal refusal = refused.get(i);
                printInfo("  " + refusal.peptide() + (refusal.charge() == null ? "" : " charge " + refusal.charge())
                        + ": " + KoinaBisection.quote(refusal.message()));
            }
            if (refused.size() > REFUSED_LINES) {
                printInfo("  and " + (refused.size() - REFUSED_LINES) + " more");
            }
        }
        //a request given up while the best model is looked for (KoinaTask.Result.GIVEN_UP) has no predictions
        //either; it was reported when it was given up, and the model is marked failed
        if (!klr.failed && klr.getPreds().isEmpty()) {
            throw new KoinaException("Koina model " + model + " refused "
                    + (refused.size() == 1 ? "the 1 precursor" : "all " + refused.size() + " precursors")
                    + " it was sent, with HTTP 400 " + bisection.firstRefusal() + ", so it predicts nothing");
        }
    }

    //what a task's exception means for the run: a KoinaException as it is, anything else as one naming the model
    private static KoinaException taskFailure(String model, Throwable cause) {
        KoinaException koinaException = KoinaException.in(cause);
        if (koinaException != null) {
            return koinaException;
        }
        if (cause instanceof Error) {
            throw (Error) cause;
        }
        return new KoinaException("Calling Koina model " + model + " failed: " + cause, cause);
    }

    /**
     * Reads a Koina infer response into klr's predictions. The rows of the response are the precursors of the
     * request it answers, in the order they were sent: {@code peptides} as written for the model and
     * {@code charges} as written in the request. The response to any part of a request can so be read on its own.
     * <p>
     * The outputs are found by the names in the model's spec, in whatever order Koina answers them. Their values
     * may be FP32 or FP64, and fragment annotations are BYTES. A fragment slot that is padding is left out: an
     * intensity or m/z that is not above 0 (Koina pads with -1), or an empty or "-1" annotation.
     *
     * @param charges the precursor charges; null for a model that is sent no charge
     * @throws KoinaException naming the model if the response is not one MSBooster can read. Reading the same
     *                        response again would fail the same way, so it is not worth a retry
     */
    public static void parseKoinaOutput(String koinaString, String[] peptides, String[] charges, String model,
                                        KoinaLibReader klr) {
        KoinaModelSpec spec = klr.spec;
        try {
            if (charges == null && spec.kind() != KoinaModelSpec.Kind.RT) {
                throw new IllegalArgumentException("the predictions of a " + spec.kind() + " model need the "
                        + "charges of the precursors");
            }
            if (charges != null && charges.length != peptides.length) {
                throw new IllegalArgumentException(peptides.length + " peptides were sent with " + charges.length
                        + " charges");
            }
            Map<String, Output> outputs = readOutputs(koinaString, new HashSet<>(spec.outputs()));
            switch (spec.kind()) {
                case RT:
                    assignRTs(peptides, predictions(outputs, spec.predictionOutput(), peptides.length), klr);
                    break;
                case IM:
                    assignIMs(peptides, charges, predictions(outputs, spec.predictionOutput(), peptides.length),
                            klr);
                    break;
                default: //ms2 or ms2_aux
                    parseMS2(outputs, peptides, charges, klr);
                    break;
            }
        } catch (KoinaException e) {
            throw e;
        } catch (RuntimeException e) { //the response, or an annotation in it, is not what the spec expects
            throw new KoinaException("MSBooster could not read the predictions of Koina model " + model + ": "
                    + e.getMessage(), e);
        }
    }

    private static void parseMS2(Map<String, Output> outputs, String[] peptides, String[] charges,
                                 KoinaLibReader klr) {
        KoinaModelSpec spec = klr.spec;
        int numPeptides = peptides.length;
        Output intensityOutput = output(outputs, spec.intensitiesOutput(), numPeptides, true);
        Output mzOutput = output(outputs, spec.mzOutput(), numPeptides, true);
        Output annotationOutput = output(outputs, spec.annotationOutput(), numPeptides, false);
        int vectorLength = intensityOutput.rowLength(numPeptides);
        if (mzOutput.rowLength(numPeptides) != vectorLength
                || annotationOutput.rowLength(numPeptides) != vectorLength) {
            throw new IllegalArgumentException("its outputs \"" + intensityOutput.name + "\", \"" + mzOutput.name
                    + "\" and \"" + annotationOutput.name + "\" have " + vectorLength + ", "
                    + mzOutput.rowLength(numPeptides) + " and " + annotationOutput.rowLength(numPeptides)
                    + " fragments per precursor");
        }

        float[][] allIntensities = new float[numPeptides][];
        float[][] allMZs = new float[numPeptides][];
        String[][] allFragmentIonTypes = new String[numPeptides][];
        String[][] allFullAnnotations = new String[numPeptides][];
        int[][] allFragNums = new int[numPeptides][];
        int[][] allCharges = new int[numPeptides][];

        for (int i = 0; i < numPeptides; i++) {
            ArrayList<Float> intensities = new ArrayList<>();
            ArrayList<Float> mzs = new ArrayList<>();
            ArrayList<String> fragmentIonTypes = new ArrayList<>();
            ArrayList<String> fullAnnotation = new ArrayList<>();
            ArrayList<Integer> fragNums = new ArrayList<>();
            ArrayList<Integer> fragCharges = new ArrayList<>();
            for (int j = i * vectorLength; j < (i + 1) * vectorLength; j++) {
                float intensity = intensityOutput.numbers[j];
                float mz = mzOutput.numbers[j];
                String annotation = annotationOutput.strings[j];
                if (intensity > 0 && mz > 0 && !isPadding(annotation)) { //dealing with -1 values
                    intensities.add(intensity);
                    mzs.add(mz);
                    //keep the annotation when the m/z comes from Koina. Otherwise the m/z is calculated
                    //later from the parsed fragment, in case we do not trust those from the model
                    if (spec.trustKoinaMz()) {
                        fullAnnotation.add(annotation);
                    }
                    spec.annotationParser().parse(annotation, fragCharges, fragmentIonTypes, fragNums);
                }
            }

            allIntensities[i] = toFloatArray(intensities);
            allFragmentIonTypes[i] = fragmentIonTypes.toArray(new String[0]);
            allFragNums[i] = toIntArray(fragNums);
            allCharges[i] = toIntArray(fragCharges);
            if (spec.trustKoinaMz()) {
                allMZs[i] = toFloatArray(mzs);
                allFullAnnotations[i] = fullAnnotation.toArray(new String[0]);
            }
        }

        if (spec.trustKoinaMz()) {
            assignMS2(peptides, charges, allMZs, allIntensities, allFragmentIonTypes, allFragNums, allCharges,
                    allFullAnnotations, klr);
        } else {
            assignMS2(peptides, charges, allIntensities, allFragmentIonTypes, allFragNums, allCharges, klr);
        }
    }

    //an annotation slot with no fragment in it
    private static boolean isPadding(String annotation) {
        return annotation.isEmpty() || annotation.equals("-1");
    }

    private static float[] toFloatArray(List<Float> values) {
        float[] array = new float[values.size()];
        for (int i = 0; i < array.length; i++) {
            array[i] = values.get(i);
        }
        return array;
    }

    private static int[] toIntArray(List<Integer> values) {
        int[] array = new int[values.size()];
        for (int i = 0; i < array.length; i++) {
            array[i] = values.get(i);
        }
        return array;
    }

    //one output of a Koina response. The values of an FP32 or FP64 output are in numbers, read as floats from the
    //text Koina wrote; those of a BYTES output are in strings
    private static final class Output {
        String name;
        String datatype;
        int[] shape;
        float[] numbers;
        String[] strings;
        //the values as written, while the datatype is not known yet
        List<String> raw;

        boolean isNumeric() {
            return "FP32".equals(datatype) || "FP64".equals(datatype);
        }

        boolean isBytes() {
            return "BYTES".equals(datatype);
        }

        int length() {
            return numbers != null ? numbers.length : strings != null ? strings.length : 0;
        }

        //the values per precursor
        int rowLength(int rows) {
            return rows == 0 ? 0 : length() / rows;
        }
    }

    //the prediction of an RT or IM model: one value per precursor
    private static float[] predictions(Map<String, Output> outputs, String name, int numPeptides) {
        Output output = output(outputs, name, numPeptides, true);
        if (output.length() != numPeptides) {
            throw new IllegalArgumentException("its output \"" + name + "\" has " + output.length()
                    + " values for " + numPeptides + " precursors, not one each");
        }
        return output.numbers;
    }

    //the output with this name, checked to be of the expected type and to have one row per precursor
    private static Output output(Map<String, Output> outputs, String name, int numPeptides, boolean numeric) {
        Output output = outputs.get(name);
        if (output == null) {
            throw new IllegalArgumentException("its answer has no output \"" + name + "\"; its outputs are "
                    + outputs.keySet());
        }
        if (numeric ? !output.isNumeric() : !output.isBytes()) {
            throw new IllegalArgumentException("its output \"" + name + "\" is " + output.datatype + ", not "
                    + (numeric ? "FP32 or FP64" : "BYTES"));
        }
        long size = 1;
        if (output.shape != null) {
            for (int dimension : output.shape) {
                size *= dimension;
            }
        }
        if (output.shape == null || output.shape.length == 0 || output.shape[0] != numPeptides
                || size != output.length()) {
            throw new IllegalArgumentException("its output \"" + name + "\" has shape "
                    + Arrays.toString(output.shape) + " and " + output.length() + " values, for "
                    + numPeptides + " precursors");
        }
        return output;
    }

    /**
     * The outputs of a Koina response by name. Only the values of the outputs named in {@code wanted} are kept.
     * The reader is lenient, so a NaN written without quotes reads as NaN, as it did before Gson read responses.
     */
    private static Map<String, Output> readOutputs(String koinaString, Set<String> wanted) {
        Map<String, Output> outputs = new HashMap<>();
        boolean hasOutputs = false;
        String error = null;
        try (JsonReader reader = new JsonReader(new StringReader(koinaString))) {
            reader.setStrictness(Strictness.LENIENT);
            reader.beginObject();
            while (reader.hasNext()) {
                String key = reader.nextName();
                if (key.equals("outputs")) {
                    hasOutputs = true;
                    reader.beginArray();
                    while (reader.hasNext()) {
                        Output output = readOutput(reader, wanted);
                        if (output.name != null) {
                            outputs.putIfAbsent(output.name, output);
                        }
                    }
                    reader.endArray();
                } else if (key.equals("error") && reader.peek() == JsonToken.STRING) {
                    error = reader.nextString();
                } else {
                    reader.skipValue();
                }
            }
            reader.endObject();
        } catch (IOException | IllegalStateException | NumberFormatException e) {
            throw new IllegalArgumentException("its answer is not a Koina response (" + e.getMessage() + ")", e);
        }
        if (!hasOutputs) {
            throw new IllegalArgumentException(error != null ? "Koina answered with the error \"" + error + "\""
                    : "its answer has no outputs");
        }
        return outputs;
    }

    private static Output readOutput(JsonReader reader, Set<String> wanted) throws IOException {
        Output output = new Output();
        reader.beginObject();
        while (reader.hasNext()) {
            switch (reader.nextName()) {
                case "name":
                    output.name = reader.nextString();
                    break;
                case "datatype":
                    output.datatype = reader.nextString();
                    break;
                case "shape":
                    output.shape = readShape(reader);
                    break;
                case "data":
                    if (output.name != null && !wanted.contains(output.name)) {
                        reader.skipValue();
                    } else {
                        readData(reader, output);
                    }
                    break;
                default:
                    reader.skipValue();
                    break;
            }
        }
        reader.endObject();

        if (output.raw != null) { //the data came before the datatype
            if (output.isNumeric()) {
                output.numbers = new float[output.raw.size()];
                for (int i = 0; i < output.numbers.length; i++) {
                    output.numbers[i] = Float.parseFloat(output.raw.get(i));
                }
            } else if (output.isBytes()) {
                output.strings = output.raw.toArray(new String[0]);
            }
            output.raw = null;
        }
        return output;
    }

    private static int[] readShape(JsonReader reader) throws IOException {
        List<Integer> shape = new ArrayList<>();
        reader.beginArray();
        while (reader.hasNext()) {
            shape.add(reader.nextInt());
        }
        reader.endArray();
        return toIntArray(shape);
    }

    //the values of an output's data, flat or nested by rows. A number is parsed from the text Koina wrote, so an
    //FP32 or FP64 value becomes the float nearest to it, as the text parse before Gson made it
    private static void readData(JsonReader reader, Output output) throws IOException {
        if (output.isNumeric()) {
            FloatList numbers = new FloatList(expectedLength(output));
            readValues(reader, value -> numbers.add(Float.parseFloat(value)));
            output.numbers = numbers.toArray();
        } else {
            List<String> strings = new ArrayList<>(expectedLength(output));
            readValues(reader, strings::add);
            if (output.isBytes()) {
                output.strings = strings.toArray(new String[0]);
            } else {
                output.raw = strings;
            }
        }
    }

    private static void readValues(JsonReader reader, Consumer<String> values) throws IOException {
        reader.beginArray();
        while (reader.hasNext()) {
            if (reader.peek() == JsonToken.BEGIN_ARRAY) {
                readValues(reader, values);
            } else {
                values.accept(reader.nextString());
            }
        }
        reader.endArray();
    }

    //the number of values the shape read so far promises, as the starting capacity of the list they go into
    private static int expectedLength(Output output) {
        long size = 1;
        if (output.shape != null) {
            for (int dimension : output.shape) {
                size *= Math.max(dimension, 0);
            }
        }
        return (int) Math.min(Math.max(size, 16), 1 << 24);
    }

    private static final class FloatList {
        private float[] values;
        private int size;

        FloatList(int capacity) {
            values = new float[capacity];
        }

        void add(float value) {
            if (size == values.length) {
                values = Arrays.copyOf(values, values.length * 2);
            }
            values[size++] = value;
        }

        float[] toArray() {
            return size == values.length ? values : Arrays.copyOf(values, size);
        }
    }

    private static void assignRTs(String[] peptides, float[] RTs, KoinaLibReader klr) {
        PredictionEntryHashMap preds = klr.getPreds();
        for (int i = 0; i < peptides.length; i++) {
            PeptideFormatter pf = new PeptideFormatter(peptides[i], 1, klr.spec);
            String peptide = pf.getBase() + "|";

            for (int charge = Constants.minPrecursorCharge; charge < Constants.maxPrecursorCharge + 1; charge++) {
                String peptideCharge = peptide + charge;
                PredictionEntry pe;
                if (preds.containsKey(peptideCharge)) {
                    pe = preds.get(peptideCharge);
                } else {
                    pe = new PredictionEntry(); //TODO: is this overkill?
                }
                pe.setRT(RTs[i]);
                preds.put(peptideCharge, pe);
            }
        }
    }

    private static void assignIMs(String[] peptides, String[] charges, float[] IMs, KoinaLibReader klr) {
        KoinaModelSpec.ImConversion conversion = klr.spec.imConversion();

        PredictionEntryHashMap preds = klr.getPreds();
        for (int i = 0; i < peptides.length; i++) {
            PeptideFormatter pf = new PeptideFormatter(peptides[i], charges[i], klr.spec);
            String peptide = pf.getBaseCharge();

            PredictionEntry pe;
            if (preds.containsKey(peptide)) {
                pe = preds.get(peptide);
            } else {
                pe = new PredictionEntry();
            }
            float im = IMs[i];
            if (conversion.unit() == KoinaModelSpec.ImUnit.CCS) {
                //convert from ccs to 1/K0 by the Mason Schampp equation (see KoinaModelSpec.ImConversion)
                MassCalculator mc = new MassCalculator(pf.getBase(), charges[i]);
                double reducedMass = mc.mass + mc.charge * mc.proton;
                reducedMass = reducedMass * conversion.gasMass() / (reducedMass + conversion.gasMass());
                im = (float) (im * Math.sqrt(reducedMass) / mc.charge / conversion.ccsCoefficient());
            }
            pe.setIM(im);
            preds.put(peptide, pe);
        }
    }

    private static void assignMS2(String[] peptides, String[] charges, float[][] intensities,
                                  String[][] fragmentIonTypes, int[][] fragNums, int[][] fragCharges,
                                  KoinaLibReader klr) {
        PredictionEntryHashMap preds = klr.getPreds();
        for (int i = 0; i < peptides.length; i++) {
            PeptideFormatter pf = new PeptideFormatter(peptides[i], charges[i], klr.spec);
            String peptide = pf.getBaseCharge();

            //problem with Koina not including PTM mass in m/z. Need to calculate here
            String[] pepSplit = peptide.split("\\|");
            MassCalculator mc = new MassCalculator(pepSplit[0], pepSplit[1]);
            float[] mzs = new float[intensities[i].length];
            for (int j = 0; j < intensities[i].length; j++) {
                mzs[j] = mc.calcMass(fragNums[i][j], fragmentIonTypes[i][j], fragCharges[i][j], 0);
            }

            PredictionEntry pe = new PredictionEntry(mzs, intensities[i], fragNums[i],
                    fragCharges[i], fragmentIonTypes[i]);

            if (preds.containsKey(peptide)) {
                PredictionEntry oldPe = preds.get(peptide);
                pe.setRT(oldPe.getRT());
                pe.setIM(oldPe.getIM());
            }
            preds.put(peptide, pe);
        }
    }

    private static void assignMS2(String[] peptides, String[] charges, float[][] mzs, float[][] intensities,
                                  String[][] fragmentIonTypes, int[][] fragNums, int[][] fragCharges,
                                  String[][] fullAnnotations, KoinaLibReader klr) {
        PredictionEntryHashMap preds = klr.getPreds();
        for (int i = 0; i < peptides.length; i++) {
            PeptideFormatter pf = new PeptideFormatter(peptides[i], charges[i], klr.spec);
            String peptide = pf.getBaseCharge();

            PredictionEntry pe = new PredictionEntry(mzs[i], intensities[i], fragNums[i],
                    fragCharges[i], fragmentIonTypes[i], fullAnnotations[i]);

            if (preds.containsKey(peptide)) {
                PredictionEntry oldPe = preds.get(peptide);
                pe.setRT(oldPe.getRT());
                pe.setIM(oldPe.getIM());
            }
            preds.put(peptide, pe);
        }
    }
}
