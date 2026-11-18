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

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertNotNull;
import static org.junit.jupiter.api.Assertions.assertNull;
import static org.junit.jupiter.api.Assertions.assertThrows;
import static org.junit.jupiter.api.Assertions.assertTrue;

import allconstants.Constants;
import java.lang.reflect.Field;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.Map;
import koinaclasses.KoinaBatch;
import koinaclasses.KoinaException;
import koinaclasses.KoinaModelRegistry;
import koinaclasses.KoinaModelSpec;
import org.junit.jupiter.api.Test;
import peptideptmformatting.PeptideFormatter;
import predictions.PredictionEntry;
import readers.predictionreaders.KoinaLibReader;
import testutils.KoinaGolden;

/**
 * KoinaModelCaller.parseKoinaOutput reads a Koina response by the output names of the model's spec, whatever
 * order Koina answers them in and whether the values are FP32 or FP64, against the peptides and charges it is
 * given rather than a request file. The responses here are written by hand, so each test shows the case it is
 * about; the recorded responses of every model are in KoinaGoldenTest.
 */
class KoinaModelCallerParseTest {
    //an MS2 model whose m/z MSBooster computes, and whose outputs Koina answers as annotation, mz, intensities
    private static final String MS2_MODEL = "Prosit_2019_intensity";
    private static final String RT_MODEL = "Prosit_2019_irt";
    private static final String[] PEPTIDES = {"PEPTIDEK", "LGEPTIDEK"};
    private static final String[] CHARGES = {"2", "3"};

    //two precursors with four fragment slots each. PEPTIDEK has b1+2 padded with -1; LGEPTIDEK has b2+1 at
    //intensity 0, which is not a predicted fragment either
    private static final String ANNOTATION = output("annotation", "BYTES", 2, 4,
            "\"y1+1\",\"b2+1\",\"y2+1\",\"b1+2\",\"y1+1\",\"b2+1\",\"y2+1\",\"b1+2\"");
    private static final String MZ_FP32 = output("mz", "FP32", 2, 4,
            "147.11281,227.10263,276.15509,-1.0,147.11281,227.10263,276.15509,57.02602");
    private static final String INTENSITIES_FP32 = output("intensities", "FP32", 2, 4,
            "1.0,0.25,0.123456789,-1.0,0.5,0.0,1.0,0.75");

    @Test
    void anFp32ResponseGivesEachPrecursorItsFragments() {
        KoinaLibReader klr = parse(MS2_MODEL, response(ANNOTATION, MZ_FP32, INTENSITIES_FP32), PEPTIDES, CHARGES);

        Map<String, PredictionEntry> preds = klr.getPreds();
        assertEquals(2, preds.size());
        PredictionEntry first = preds.get(key(MS2_MODEL, "PEPTIDEK", "2"));
        assertNotNull(first, "PEPTIDEK|2 in " + preds.keySet());
        assertEquals(3, first.numFragments(), "the -1 slot is padding");
        assertEquals(1.0f, intensity(first, "y", 1, 1));
        assertEquals(0.25f, intensity(first, "b", 2, 1));
        assertEquals(Float.parseFloat("0.123456789"), intensity(first, "y", 2, 1));
        assertNull(intensity(first, "b", 1, 2));

        PredictionEntry second = preds.get(key(MS2_MODEL, "LGEPTIDEK", "3"));
        assertNotNull(second, "LGEPTIDEK|3 in " + preds.keySet());
        assertEquals(3, second.numFragments(), "an intensity of 0 is no fragment");
        assertEquals(0.5f, intensity(second, "y", 1, 1));
        assertNull(intensity(second, "b", 2, 1));
        assertEquals(1.0f, intensity(second, "y", 2, 1));
        assertEquals(0.75f, intensity(second, "b", 1, 2));
    }

    @Test
    void theOutputsAreFoundByNameInAnyOrder() {
        String expected = dump(MS2_MODEL, response(INTENSITIES_FP32, MZ_FP32, ANNOTATION), PEPTIDES, CHARGES);
        List<String[]> orders = List.of(
                new String[] {INTENSITIES_FP32, ANNOTATION, MZ_FP32},
                new String[] {MZ_FP32, INTENSITIES_FP32, ANNOTATION},
                new String[] {MZ_FP32, ANNOTATION, INTENSITIES_FP32},
                new String[] {ANNOTATION, INTENSITIES_FP32, MZ_FP32},
                new String[] {ANNOTATION, MZ_FP32, INTENSITIES_FP32});
        for (String[] order : orders) {
            assertEquals(expected, dump(MS2_MODEL, response(order), PEPTIDES, CHARGES));
        }
        //an output Koina adds that the spec does not read is passed over
        String extra = output("isotopes", "INT32", 2, 1, "0,0");
        assertEquals(expected, dump(MS2_MODEL, response(extra, ANNOTATION, MZ_FP32, INTENSITIES_FP32), PEPTIDES,
                CHARGES));
    }

    @Test
    void fp64ValuesBecomeTheFloatNearestToWhatKoinaWrote() {
        String intensities = output("intensities", "FP64", 2, 4,
                "0.11706076562404633,0.25,4.612778354839975e-7,-1.0,0.5,0.0,1.0,0.000007271622962434776");
        String mz = output("mz", "FP64", 2, 4,
                "147.1128,227.1026,276.1551,-1.0,147.1128,227.1026,276.1551,57.026");
        PredictionEntry first = parse(MS2_MODEL, response(ANNOTATION, mz, intensities), PEPTIDES, CHARGES)
                .getPreds().get(key(MS2_MODEL, "PEPTIDEK", "2"));
        assertEquals(Float.parseFloat("0.11706076562404633"), intensity(first, "y", 1, 1));
        assertEquals(Float.parseFloat("4.612778354839975e-7"), intensity(first, "y", 2, 1));

        //the same values written as FP32 read the same
        String asFp32 = intensities.replace("\"FP64\"", "\"FP32\"");
        String mzAsFp32 = mz.replace("\"FP64\"", "\"FP32\"");
        assertEquals(dump(MS2_MODEL, response(ANNOTATION, mz, intensities), PEPTIDES, CHARGES),
                dump(MS2_MODEL, response(ANNOTATION, mzAsFp32, asFp32), PEPTIDES, CHARGES));

        String irt = output("irt", "FP64", 2, 1, "12.345678901234567,-3.5");
        KoinaLibReader klr = parse(RT_MODEL, response(irt), PEPTIDES, null);
        for (int charge = Constants.minPrecursorCharge; charge <= Constants.maxPrecursorCharge; charge++) {
            assertEquals(Float.parseFloat("12.345678901234567"), klr.getPreds().get("PEPTIDEK|" + charge).getRT());
            assertEquals(-3.5f, klr.getPreds().get("LGEPTIDEK|" + charge).getRT());
        }
    }

    /**
     * The IM models answer a CCS in nitrogen, which becomes 1/K0 by the Mason-Schampp equation as alphabase converts
     * it: 1/K0 = CCS * sqrt(mu) / z / 1059.62245, mu the reduced mass of the precursor ion (M + z protons) and N2
     * (28 Da). The expected values are worked out here from residue masses, not by MSBooster's MassCalculator.
     */
    @Test
    void aCcsBecomesAn1OverK0ByTheMasonSchamppEquation() {
        double peptidek = 927.45490; //PEPTIDEK, monoisotopic
        double lgeptidek = 1000.50766; //LGEPTIDEK
        for (String model : List.of("AlphaPept_ccs_generic", "IM2Deep")) {
            String ccs = output("ccs", "FP32", 2, 1, "400.0,512.5");
            KoinaLibReader klr = parse(model, response(ccs), PEPTIDES, CHARGES);
            assertEquals(2, klr.getPreds().size(), model);
            assertEquals(inverseK0(400.0, peptidek, 2), klr.getPreds().get(key(model, "PEPTIDEK", "2")).getIM(),
                    1e-5, model);
            assertEquals(inverseK0(512.5, lgeptidek, 3), klr.getPreds().get(key(model, "LGEPTIDEK", "3")).getIM(),
                    1e-5, model);
        }
        //an ion mobility, not the CCS
        assertEquals(0.98405, inverseK0(400.0, peptidek, 2), 1e-4);
    }

    /**
     * D7: the IM unit is the spec's. No Koina model MSBooster calls answers in 1/K0 today, so this one is IM2Deep's
     * spec with that unit: its predictions are used as they are, where IM2Deep's own are converted.
     */
    @Test
    void anImModelThatAnswersIn1OverK0HasItUsedAsItIs() throws Exception {
        KoinaModelSpec im2deep = KoinaModelRegistry.require("IM2Deep");
        KoinaModelSpec inverseK0 = new KoinaModelSpec(im2deep.name(), im2deep.kind(), im2deep.inputs(),
                im2deep.outputs(), KoinaModelSpec.ImConversion.INVERSE_K0, im2deep.peptideFormat(),
                im2deep.annotationParser(), im2deep.trustKoinaMz(), im2deep.modMassMap(), im2deep.modRule(),
                im2deep.encodableMods(), im2deep.limits(), im2deep.batchSize(), im2deep.minCe(), im2deep.maxCe(),
                im2deep.calibrateCe(), im2deep.instrumentMap(), im2deep.fragmentations(), im2deep.isobaricLabel(),
                im2deep.takesIsobaricLabels(), im2deep.sortBatchByLength(), im2deep.singleThreadFormatting(),
                im2deep.citationKey());
        String ims = output("ccs", "FP32", 2, 1, "0.875,1.0625");

        KoinaLibReader klr = new KoinaLibReader("IM2Deep");
        //KoinaLibReader takes its spec from the registry, which has no such model
        Field spec = KoinaLibReader.class.getField("spec");
        spec.setAccessible(true);
        spec.set(klr, inverseK0);
        KoinaModelCaller.parseKoinaOutput(response(ims), PEPTIDES, CHARGES, "IM2Deep", klr);
        assertEquals(0.875f, klr.getPreds().get(key("IM2Deep", "PEPTIDEK", "2")).getIM());
        assertEquals(1.0625f, klr.getPreds().get(key("IM2Deep", "LGEPTIDEK", "3")).getIM());

        float converted = parse("IM2Deep", response(ims), PEPTIDES, CHARGES).getPreds()
                .get(key("IM2Deep", "PEPTIDEK", "2")).getIM();
        assertEquals(inverseK0(0.875, 927.45490, 2), converted, 1e-6, "IM2Deep's own spec takes it as a CCS");
    }

    private static double inverseK0(double ccs, double neutralMass, int charge) {
        double ion = neutralMass + charge * 1.00727646688;
        double reducedMass = ion * 28.0 / (ion + 28.0);
        return ccs * Math.sqrt(reducedMass) / charge / 1059.62245;
    }

    @Test
    void aPartOfABatchIsReadOnItsOwn() {
        String whole = dump(MS2_MODEL, response(ANNOTATION, MZ_FP32, INTENSITIES_FP32), PEPTIDES, CHARGES);

        //the response Koina gives to a request of LGEPTIDEK alone: the second row of each output
        String annotation = output("annotation", "BYTES", 1, 4, "\"y1+1\",\"b2+1\",\"y2+1\",\"b1+2\"");
        String mz = output("mz", "FP32", 1, 4, "147.11281,227.10263,276.15509,57.02602");
        String intensities = output("intensities", "FP32", 1, 4, "0.5,0.0,1.0,0.75");
        KoinaLibReader klr = parse(MS2_MODEL, response(annotation, mz, intensities), new String[] {"LGEPTIDEK"},
                new String[] {"3"});

        assertEquals(1, klr.getPreds().size());
        String key = key(MS2_MODEL, "LGEPTIDEK", "3");
        assertEquals(entryLines(whole, key), entryLines(KoinaGolden.dump(klr.getPreds()), key));
    }

    @Test
    void theRowsOfAResponseMayBeNestedAndItsKeysInAnyOrder() {
        String expected = dump(MS2_MODEL, response(ANNOTATION, MZ_FP32, INTENSITIES_FP32), PEPTIDES, CHARGES);
        String annotation = "{\"data\":[[\"y1+1\",\"b2+1\",\"y2+1\",\"b1+2\"],[\"y1+1\",\"b2+1\",\"y2+1\",\"b1+2\"]],"
                + "\"shape\":[2,4],\"datatype\":\"BYTES\",\"name\":\"annotation\"}";
        String mz = "{\"data\":[[147.11281,227.10263,276.15509,-1.0],[147.11281,227.10263,276.15509,57.02602]],"
                + "\"datatype\":\"FP32\",\"name\":\"mz\",\"shape\":[2,4]}";
        assertEquals(expected, dump(MS2_MODEL, response(annotation, mz, INTENSITIES_FP32), PEPTIDES, CHARGES));
    }

    @Test
    void anEmptyOrMinusOneAnnotationIsPadding() {
        String annotation = output("annotation", "BYTES", 2, 4,
                "\"y1+1\",\"\",\"y2+1\",\"b1+2\",\"y1+1\",\"b2+1\",\"-1\",\"b1+2\"");
        String intensities = output("intensities", "FP32", 2, 4, "1.0,0.25,0.5,-1.0,0.5,0.25,1.0,0.75");
        String mz = output("mz", "FP32", 2, 4,
                "147.11281,227.10263,276.15509,-1.0,147.11281,227.10263,276.15509,57.02602");
        Map<String, PredictionEntry> preds = parse(MS2_MODEL, response(annotation, mz, intensities), PEPTIDES,
                CHARGES).getPreds();
        assertEquals(2, preds.get(key(MS2_MODEL, "PEPTIDEK", "2")).numFragments());
        assertEquals(3, preds.get(key(MS2_MODEL, "LGEPTIDEK", "3")).numFragments());
    }

    @Test
    void aResponseThatCannotBeReadEndsTheRunNamingTheModel() {
        //an annotation the model's parser does not know
        String unknown = output("annotation", "BYTES", 2, 4,
                "\"y1+1\",\"b2+1\",\"q2+1\",\"b1+2\",\"y1+1\",\"b2+1\",\"y2+1\",\"b1+2\"");
        assertFails(response(unknown, MZ_FP32, INTENSITIES_FP32), PEPTIDES, "q2+1");
        //an output the spec reads is missing
        assertFails(response(MZ_FP32, INTENSITIES_FP32), PEPTIDES, "no output \"annotation\"");
        //values of a type other than FP32 and FP64
        assertFails(response(ANNOTATION, MZ_FP32, INTENSITIES_FP32.replace("\"FP32\"", "\"INT32\"")), PEPTIDES,
                "INT32, not FP32 or FP64");
        //rows for another number of precursors than were sent
        assertFails(response(ANNOTATION, MZ_FP32, INTENSITIES_FP32), new String[] {"PEPTIDEK", "LGEPTIDEK", "ADGLR"},
                "for 3 precursors");
        //an error instead of predictions
        assertFails("{\"error\":\"unexpected inference input\"}", PEPTIDES, "unexpected inference input");
        //not JSON at all
        assertFails("<html><body>504 Gateway Time-out</body></html>", PEPTIDES, "not a Koina response");
    }

    @Test
    void theBatchOfARequestIsItsPeptidesAndCharges() {
        String request = "{\"id\":\"0\",\"inputs\":[{\"shape\":[2,1],\"data\":[[\"PEPTIDEK\"],[\"LGEPTIDEK\"]],"
                + "\"datatype\":\"BYTES\",\"name\":\"peptide_sequences\"},"
                + "{\"shape\":[2,1],\"data\":[[30],[30]],\"datatype\":\"FP32\",\"name\":\"collision_energies\"},"
                + "{\"shape\":[2,1],\"data\":[[2],[3]],\"datatype\":\"INT32\",\"name\":\"precursor_charges\"}]}";
        KoinaBatch batch = KoinaBatch.fromRequest(request);
        assertEquals(List.of(PEPTIDES), List.of(batch.peptides()));
        assertEquals(List.of(CHARGES), List.of(batch.charges()));

        KoinaBatch rt = KoinaBatch.fromRequest("{\"inputs\":[{\"name\":\"peptide_sequences\",\"datatype\":\"BYTES\","
                + "\"shape\":[1,1],\"data\":[[\"PEPTIDEK\"]]}]}");
        assertEquals(1, rt.size());
        assertNull(rt.charges());
    }

    private static void assertFails(String response, String[] peptides, String reason) {
        String[] charges = new String[peptides.length];
        Arrays.fill(charges, "2");
        KoinaException e = assertThrows(KoinaException.class, () -> parse(MS2_MODEL, response, peptides, charges));
        assertTrue(e.getMessage().contains(MS2_MODEL), e.getMessage());
        assertTrue(e.getMessage().contains(reason), e.getMessage());
    }

    private static KoinaLibReader parse(String model, String response, String[] peptides, String[] charges) {
        KoinaLibReader klr = new KoinaLibReader(model);
        KoinaModelCaller.parseKoinaOutput(response, peptides, charges, model, klr);
        return klr;
    }

    private static String dump(String model, String response, String[] peptides, String[] charges) {
        return KoinaGolden.dump(parse(model, response, peptides, charges).getPreds());
    }

    private static String key(String model, String peptide, String charge) {
        return new PeptideFormatter(peptide, charge, KoinaModelRegistry.require(model)).getBaseCharge();
    }

    //the intensity of a fragment, or null if the entry does not have it
    private static Float intensity(PredictionEntry pe, String ionType, int fragNum, int charge) {
        for (int i = 0; i < pe.numFragments(); i++) {
            if (pe.getIonTypeString(i).equals(ionType) && pe.getFragNum(i) == fragNum && pe.getCharge(i) == charge) {
                return pe.getIntensity(i);
            }
        }
        return null;
    }

    //the P line of a precursor in a dump, and the F lines after it
    private static List<String> entryLines(String dump, String key) {
        List<String> lines = new ArrayList<>();
        boolean in = false;
        for (String line : dump.split("\n")) {
            if (line.startsWith("P\t")) {
                in = line.startsWith("P\t" + key + "\t");
            }
            if (in) {
                lines.add(line);
            }
        }
        return lines;
    }

    private static String output(String name, String datatype, int rows, int columns, String data) {
        return "{\"name\":\"" + name + "\",\"datatype\":\"" + datatype + "\",\"shape\":[" + rows + "," + columns
                + "],\"data\":[" + data + "]}";
    }

    private static String response(String... outputs) {
        return "{\"id\":\"0\",\"model_name\":\"test\",\"model_version\":\"1\",\"outputs\":["
                + String.join(",", outputs) + "]}";
    }
}
