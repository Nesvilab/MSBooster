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

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertFalse;
import static org.junit.jupiter.api.Assertions.assertNotNull;
import static org.junit.jupiter.api.Assertions.assertNull;
import static org.junit.jupiter.api.Assertions.assertSame;
import static org.junit.jupiter.api.Assertions.assertThrows;
import static org.junit.jupiter.api.Assertions.assertTrue;

import allconstants.Constants;
import allconstants.LowercaseModelMapper;
import allconstants.ModelCollections;
import allconstants.NceConstants;
import java.net.URL;
import java.net.URLClassLoader;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.LinkedHashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.TreeSet;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.stream.Collectors;
import java.util.stream.Stream;
import koinaclasses.KoinaModelSpec.AnnotationParser;
import koinaclasses.KoinaModelSpec.Input;
import koinaclasses.KoinaModelSpec.Kind;
import org.json.JSONObject;
import org.junit.jupiter.api.Test;
import org.junit.jupiter.api.io.TempDir;
import peptideptmformatting.PTMhandler;
import peptideptmformatting.PeptideFormatter;
import peptideptmformatting.PeptideSkipper;
import utils.InstrumentUtils;
import writers.JSONWriter;

/**
 * The Koina model registry holds the models MSBooster had before it existed, with the rules the code had for
 * them, and the models added with it (D12), and the lists made from it are the lists there were, with those models.
 * The request and prediction goldens (KoinaGoldenTest) check the rest of today's behaviour; KoinaPrositPtmTest the
 * rules of the Prosit-PTM models.
 */
class KoinaModelRegistryTest {
    private static final Set<String> RT_MODELS = Set.of(
            "AlphaPept_rt_generic", "Prosit_2019_irt", "Prosit_2020_irt_TMT", "Prosit_2024_irt_cit",
            "Deeplc_hela_hf",
            "Prosit_2025_irt_22PTM", "Prosit_2025_irt_lac", "Chronologer_RT");
    private static final Set<String> MS2_MODELS = Set.of(
            "ms2pip_2021_HCD", "ms2pip_timsTOF2024", "ms2pip_CID_TMT", "ms2pip_TTOF5600",
            "ms2pip_Immuno_HCD", "ms2pip_iTRAQphospho",
            "AlphaPept_ms2_generic",
            "Prosit_2019_intensity", "Prosit_2020_intensity_CID", "Prosit_2020_intensity_TMT",
            "Prosit_2020_intensity_HCD", "Prosit_2023_intensity_timsTOF", "Prosit_2024_intensity_cit",
            "Prosit_2025_intensity_MultiFrag",
            "UniSpec",
            "ms2pip_timsTOF2023", "Prosit_2025_intensity_22PTM", "Prosit_2025_intensity_40PTM",
            "Prosit_2025_intensity_lac", "Altimeter_2024_intensities");
    private static final Set<String> IM_MODELS = Set.of("AlphaPept_ccs_generic", "IM2Deep");
    private static final int MODELS = RT_MODELS.size() + MS2_MODELS.size() + IM_MODELS.size();

    @Test
    void theRegistryHoldsTheModelsOfTodayAndTheListsAreMadeFromIt() {
        assertEquals(RT_MODELS, names(Kind.RT));
        assertEquals(MS2_MODELS, names(Kind.MS2));
        assertEquals(IM_MODELS, names(Kind.IM));
        assertEquals(30, KoinaModelRegistry.all().size());
        assertEquals(MODELS, KoinaModelRegistry.all().size());

        assertEquals(RT_MODELS, new HashSet<>(ModelCollections.KoinaRTmodels));
        assertEquals(MS2_MODELS, new HashSet<>(ModelCollections.KoinaMS2models));
        assertEquals(IM_MODELS, new HashSet<>(ModelCollections.KoinaIMmodels));
        assertEquals(IM_MODELS, new HashSet<>(ModelCollections.KoinaCCSmodels));
        Set<String> all = new HashSet<>(RT_MODELS);
        all.addAll(MS2_MODELS);
        all.addAll(IM_MODELS);
        assertEquals(all, new HashSet<>(ModelCollections.KoinaModels));
        assertEquals(Set.of("Prosit_2020_irt_TMT", "Prosit_2020_intensity_TMT"),
                new HashSet<>(ModelCollections.KoinaIsoLabelmodels));

        Map<String, Set<String>> fragmentations = new HashMap<>();
        fragmentations.put("Prosit_2020_intensity_TMT", Set.of("HCD", "CID"));
        fragmentations.put("Prosit_2024_intensity_cit", Set.of("HCD", "CID"));
        fragmentations.put("Prosit_2025_intensity_MultiFrag", Set.of("HCD", "ECD", "EID", "UVPD", "ETCID"));
        fragmentations.put("Prosit_2025_intensity_22PTM", Set.of("HCD", "CID"));
        fragmentations.put("Prosit_2025_intensity_40PTM", Set.of("HCD", "CID"));
        fragmentations.put("Prosit_2025_intensity_lac", Set.of("HCD", "CID"));
        assertEquals(fragmentations, new HashMap<>(ModelCollections.allowedFragmentationTypes));

        assertEquals(Set.of("Prosit", "PrositTMT", "alphapeptdeep",
                        "AlphaPept_ms2_generic", "Prosit_2019_intensity", "Prosit_2023_intensity_timsTOF",
                        "Prosit_2020_intensity_TMT", "Prosit_2020_intensity_HCD", "Prosit_2024_intensity_cit",
                        "UniSpec", "Prosit_2025_intensity_22PTM", "Prosit_2025_intensity_40PTM",
                        "Prosit_2025_intensity_lac", "Altimeter_2024_intensities"),
                new HashSet<>(NceConstants.nceModels));

        HashMap<String, String> lowercase = new LowercaseModelMapper().getLowercaseToModel();
        for (String model : all) {
            assertEquals(model, lowercase.get(model.toLowerCase()), model + " cannot be selected by name");
        }
    }

    @Test
    void findIgnoresCaseAndKnowsNoOtherModel() {
        Set<String> lowercaseNames = new HashSet<>();
        for (KoinaModelSpec spec : KoinaModelRegistry.all()) {
            assertTrue(lowercaseNames.add(spec.name().toLowerCase()), "two specs named " + spec.name());
            assertSame(spec, KoinaModelRegistry.find(spec.name()));
            assertSame(spec, KoinaModelRegistry.find(spec.name().toUpperCase()));
            assertSame(spec, KoinaModelRegistry.find(spec.name().toLowerCase()));
        }
        assertNull(KoinaModelRegistry.find("PredFull"));
        assertNull(KoinaModelRegistry.find("AlphaPeptDeep_ms2_generic")); //a citation key, not a model name
        assertNull(KoinaModelRegistry.find("DIA-NN"));
        assertNull(KoinaModelRegistry.find(""));
        assertNull(KoinaModelRegistry.find(null));
    }

    //D14: a citation key is not a model name, and neither is any alias Koina's own docs use for a model
    @Test
    void noCitationKeyOrAliasIsAModelName() {
        HashMap<String, String> lowercase = new LowercaseModelMapper().getLowercaseToModel();
        List<String> aliases = new ArrayList<>(List.of("AlphaPeptDeep_rt_generic", "AlphaPeptDeep_ms2_generic",
                "AlphaPeptDeep_ccs_generic", "ms2pip_HCD2021"));
        for (KoinaModelSpec spec : KoinaModelRegistry.all()) {
            if (!spec.citationKey().equals(spec.name())) {
                aliases.add(spec.citationKey());
            }
        }
        for (String alias : aliases) {
            assertNull(KoinaModelRegistry.find(alias), alias);
            assertFalse(lowercase.containsKey(alias.toLowerCase()), alias + " can be selected by name");
            assertFalse(ModelCollections.KoinaModels.contains(alias), alias);
        }
    }

    //the best-model search lists of the base (16df34d)
    private static final Map<String, List<String>> FIND_BEST_LISTS = Map.of(
            "RT general", List.of("DIA-NN", "FragCast", "AlphaPept_rt_generic", "Prosit_2019_irt", "Deeplc_hela_hf"),
            "RT isolabel", List.of("DIA-NN", "FragCast", "Prosit_2020_irt_TMT", "AlphaPept_rt_generic",
                    "Deeplc_hela_hf"),
            "MS2 general", List.of("DIA-NN", "FragCast", "ms2pip_2021_HCD", "ms2pip_TTOF5600", "AlphaPept_ms2_generic",
                    "Prosit_2020_intensity_CID", "Prosit_2020_intensity_HCD"),
            "MS2 isolabel", List.of("Prosit_2020_intensity_TMT", "ms2pip_CID_TMT", "ms2pip_iTRAQphospho"),
            "MS2 hla", List.of("DIA-NN", "FragCast", "ms2pip_Immuno_HCD", "AlphaPept_ms2_generic",
                    "Prosit_2020_intensity_CID", "Prosit_2020_intensity_HCD"),
            "MS2 timstof", List.of("DIA-NN", "FragCast", "ms2pip_timsTOF2024", "AlphaPept_ms2_generic",
                    "Prosit_2023_intensity_timsTOF"),
            "IM general", List.of("DIA-NN", "FragCast", "AlphaPept_ccs_generic", "IM2Deep"));

    //D18: the best-model search never tries a model added with the registry; its lists are those of the base
    @Test
    void theNewModelsStayOutOfTheBestModelSearch() {
        Set<String> newModels = Set.of("Prosit_2025_intensity_22PTM", "Prosit_2025_irt_22PTM",
                "Prosit_2025_intensity_40PTM", "Chronologer_RT", "Altimeter_2024_intensities", "ms2pip_timsTOF2023",
                "Prosit_2025_intensity_lac", "Prosit_2025_irt_lac");
        for (String model : newModels) {
            assertNotNull(KoinaModelRegistry.find(model), model);
        }

        String rt = ModelCollections.rtCollection;
        String ms2 = ModelCollections.ms2Collection;
        String im = ModelCollections.imCollection;
        Map<String, List<String>> lists = new HashMap<>();
        try {
            for (String collection : new String[] {"general", "isolabel"}) {
                ModelCollections.rtCollection = collection;
                lists.put("RT " + collection, new ArrayList<>(ModelCollections.getRtCollection()));
            }
            for (String collection : new String[] {"general", "isolabel", "hla", "timstof"}) {
                ModelCollections.ms2Collection = collection;
                lists.put("MS2 " + collection, new ArrayList<>(ModelCollections.getMs2Collection()));
            }
            ModelCollections.imCollection = "general";
            lists.put("IM general", new ArrayList<>(ModelCollections.getImCollection()));
        } finally {
            ModelCollections.rtCollection = rt;
            ModelCollections.ms2Collection = ms2;
            ModelCollections.imCollection = im;
        }

        assertEquals(FIND_BEST_LISTS, lists);
        for (Map.Entry<String, List<String>> list : lists.entrySet()) {
            for (String model : list.getValue()) {
                assertFalse(newModels.contains(model), model + " is in the " + list.getKey() + " search");
            }
        }
    }

    //what JSONWriter sent each model before the registry existed
    @Test
    void eachModelIsSentTheInputsItWasSentBefore() {
        Map<String, List<Input>> expected = new HashMap<>();
        for (String model : RT_MODELS) {
            expected.put(model, List.of());
        }
        for (String model : IM_MODELS) {
            expected.put(model, List.of(Input.PRECURSOR_CHARGES));
        }
        for (String model : MS2_MODELS) {
            expected.put(model, List.of(Input.PRECURSOR_CHARGES, Input.COLLISION_ENERGIES));
        }
        for (String model : MS2_MODELS) {
            if (model.startsWith("ms2pip")) {
                expected.put(model, List.of(Input.PRECURSOR_CHARGES));
            }
        }
        expected.put("Prosit_2020_intensity_CID", List.of(Input.PRECURSOR_CHARGES));
        expected.put("Prosit_2025_intensity_MultiFrag", List.of(Input.PRECURSOR_CHARGES, Input.FRAGMENTATION_TYPES));
        expected.put("AlphaPept_ms2_generic",
                List.of(Input.PRECURSOR_CHARGES, Input.COLLISION_ENERGIES, Input.INSTRUMENT_TYPES));
        expected.put("UniSpec", List.of(Input.PRECURSOR_CHARGES, Input.COLLISION_ENERGIES, Input.INSTRUMENT_TYPES));
        expected.put("Prosit_2020_intensity_TMT",
                List.of(Input.PRECURSOR_CHARGES, Input.COLLISION_ENERGIES, Input.FRAGMENTATION_TYPES));
        expected.put("Prosit_2024_intensity_cit",
                List.of(Input.PRECURSOR_CHARGES, Input.COLLISION_ENERGIES, Input.FRAGMENTATION_TYPES));
        //the models added with the registry: the collision energy is an NCE, the fragmentation "HCD" or "CID"
        expected.put("Prosit_2025_intensity_22PTM",
                List.of(Input.PRECURSOR_CHARGES, Input.COLLISION_ENERGIES, Input.FRAGMENTATION_TYPES));
        expected.put("Prosit_2025_intensity_40PTM",
                List.of(Input.PRECURSOR_CHARGES, Input.COLLISION_ENERGIES, Input.FRAGMENTATION_TYPES));
        //the lactylation model is sent the instrument too; Altimeter an NCE and no fragmentation
        expected.put("Prosit_2025_intensity_lac", List.of(Input.PRECURSOR_CHARGES, Input.COLLISION_ENERGIES,
                Input.INSTRUMENT_TYPES, Input.FRAGMENTATION_TYPES));
        expected.put("Altimeter_2024_intensities", List.of(Input.PRECURSOR_CHARGES, Input.COLLISION_ENERGIES));

        for (KoinaModelSpec spec : KoinaModelRegistry.all()) {
            assertEquals(expected.get(spec.name()), new ArrayList<>(spec.inputs()), spec.name());
            //a model whose collision energy is sent has it calibrated
            assertEquals(spec.sends(Input.COLLISION_ENERGIES), spec.calibrateCe(), spec.name());
        }
    }

    @Test
    void theBatchSizesAndCitationKeysAreThoseOfToday() {
        for (KoinaModelSpec spec : KoinaModelRegistry.all()) {
            int expectedBatch = 1000;
            if (spec.name().equals("UniSpec") || spec.name().equals("Prosit_2025_intensity_MultiFrag")
                    || spec.name().equals("Altimeter_2024_intensities")) {
                expectedBatch = 100;
            } else if (spec.name().equals("AlphaPept_ms2_generic")) {
                expectedBatch = 500;
            }
            assertEquals(expectedBatch, spec.batchSize(), spec.name());
            assertEquals(spec.name().startsWith("AlphaPept_"), spec.sortBatchByLength(), spec.name());
        }
        assertEquals("ms2pip_HCD2021", KoinaModelRegistry.find("ms2pip_2021_HCD").citationKey());
        assertEquals("AlphaPeptDeep_rt_generic", KoinaModelRegistry.find("AlphaPept_rt_generic").citationKey());
        assertEquals("AlphaPeptDeep_ms2_generic", KoinaModelRegistry.find("AlphaPept_ms2_generic").citationKey());
        assertEquals("AlphaPeptDeep_ccs_generic", KoinaModelRegistry.find("AlphaPept_ccs_generic").citationKey());
        assertEquals("UniSpec", KoinaModelRegistry.find("UniSpec").citationKey());
        for (String model : List.of("Prosit_2025_irt_22PTM", "Prosit_2025_intensity_22PTM",
                "Prosit_2025_intensity_40PTM", "ms2pip_timsTOF2023", "Prosit_2025_intensity_lac",
                "Prosit_2025_irt_lac", "Chronologer_RT", "Altimeter_2024_intensities")) {
            assertEquals(model, KoinaModelRegistry.find(model).citationKey());
        }
    }

    //FragPipe loads the registry from the MSBooster jar by reflection, with nothing of MSBooster run before
    @Test
    void theRegistryLoadsWithoutMSBoosterSettings() throws Exception {
        URL classes = KoinaModelRegistry.class.getProtectionDomain().getCodeSource().getLocation();
        List<String> loaded = Collections.synchronizedList(new ArrayList<>());
        //only MSBooster's own classes, without its dependencies
        try (URLClassLoader loader = new URLClassLoader(new URL[] {classes}, ClassLoader.getPlatformClassLoader()) {
            @Override
            protected Class<?> findClass(String name) throws ClassNotFoundException {
                loaded.add(name);
                return super.findClass(name);
            }
        }) {
            Class<?> registry = Class.forName("koinaclasses.KoinaModelRegistry", true, loader);
            assertEquals(MODELS, ((List<?>) registry.getMethod("all").invoke(null)).size());
            Object spec = registry.getMethod("find", String.class).invoke(null, "unispec");
            assertNotNull(spec);
            assertEquals("UniSpec", spec.getClass().getMethod("name").invoke(spec));
            assertEquals("MS2", ((Enum<?>) spec.getClass().getMethod("kind").invoke(spec)).name());
        }
        List<String> settings = loaded.stream()
                .filter(name -> name.startsWith("allconstants.") || name.startsWith("mainsteps."))
                .collect(Collectors.toList());
        assertTrue(settings.isEmpty(), "loading the registry loaded " + settings);
    }

    @Test
    void theAnnotationParsersReadTheLabelsTheyKnowAndRefuseOthers() {
        assertParsed(AnnotationParser.ION_NUMBER_CHARGE, "y12+2", 2, "y", 12);
        assertParsed(AnnotationParser.ION_NUMBER_CHARGE, "b1+1", 1, "b", 1);
        for (String label : new String[] {"c3+1", "y3", "y3+1-H2O", "a2+1", ""}) {
            assertRefused(AnnotationParser.ION_NUMBER_CHARGE, label);
        }

        assertParsed(AnnotationParser.UNISPEC, "y5-H2O^2", 2, "y-NL", 5);
        assertParsed(AnnotationParser.UNISPEC, "b3", 1, "b", 3);
        assertParsed(AnnotationParser.UNISPEC, "a2+i", 1, "a", 2);
        assertParsed(AnnotationParser.UNISPEC, "p^2", 2, "p", 0);
        assertParsed(AnnotationParser.UNISPEC, "IMA", 1, "imm", 0);
        assertParsed(AnnotationParser.UNISPEC, "Int/AE/3", 1, "int", 0);
        for (String label : new String[] {"x3", "z5^2", "c2"}) {
            assertRefused(AnnotationParser.UNISPEC, label);
        }

        assertParsed(AnnotationParser.MULTIFRAG, "Z4^2", 2, "zprime", 4);
        assertParsed(AnnotationParser.MULTIFRAG, "c3", 1, "cdot", 3);
        assertParsed(AnnotationParser.MULTIFRAG, "A7", 1, "a+1", 7);
        for (String label : new String[] {"w3", "p", "I"}) {
            assertRefused(AnnotationParser.MULTIFRAG, label);
        }
    }

    private static void assertParsed(AnnotationParser parser, String label, int charge, String ionType, int fragNum) {
        List<Integer> charges = new ArrayList<>();
        List<String> ionTypes = new ArrayList<>();
        List<Integer> fragNums = new ArrayList<>();
        parser.parse(label, charges, ionTypes, fragNums);
        assertEquals(List.of(charge), charges, label);
        assertEquals(List.of(ionType), ionTypes, label);
        assertEquals(List.of(fragNum), fragNums, label);
    }

    private static void assertRefused(AnnotationParser parser, String label) {
        assertThrows(RuntimeException.class,
                () -> parser.parse(label, new ArrayList<>(), new ArrayList<>(), new ArrayList<>()),
                parser + " must refuse " + label);
    }

    @Test
    void theInstrumentMapsAreThoseOfTheFormerInstrumentUtils() {
        String saved = Constants.instrument;
        try {
            for (KoinaModelSpec spec : KoinaModelRegistry.all()) {
                for (String instrument : new String[] {"", "QE", "qe", "QEHFX", "qehfx", "LUMOS", "lumos", "ELITE",
                        "VELOS", "NONE", "TIMSTOF", "SCIEXTOF", "THERMOTOF", "Orbitrap"}) {
                    for (String charge : new String[] {"1", "2", "3"}) {
                        Constants.instrument = instrument;
                        //D21 for the lactylation model, which is new (see KoinaPrositLacTest)
                        String expected = spec.name().equals("Prosit_2025_intensity_lac")
                                ? ("THERMOTOF".equalsIgnoreCase(instrument) ? "ASTRAL" : "LUMOS")
                                : formerInstrument(spec.name(), instrument, charge);
                        String what = spec.name() + " " + instrument + " charge " + charge;
                        assertEquals(expected, spec.instrumentMap().map(instrument, charge), what);
                        assertEquals(expected, InstrumentUtils.mapInstrumentToModelSpecific(spec, charge), what);
                    }
                }
            }
        } finally {
            Constants.instrument = saved;
        }
    }

    //InstrumentUtils.mapInstrumentToModelSpecific before the registry, as it was
    private static String formerInstrument(String model, String instrument, String charge) {
        model = model.toLowerCase();
        List<String> unispecInstruments = List.of("LUMOS", "QE", "QEHFX", "ELITE", "VELOS", "NONE");
        List<String> alphapeptdeepInstruments = List.of("QE", "LUMOS", "TIMSTOF", "SCIEXTOF");
        if (model.contains("unispec")) {
            if (!containsIgnoreCase(unispecInstruments, instrument)) {
                instrument = "NONE";
            } else if ((instrument.equals("QE") || instrument.equals("QEHFX")) && charge.equals("1")) {
                instrument = "NONE";
            }
        } else if (model.contains("alphapept")) {
            if (instrument.equals("QEHFX")) {
                instrument = "QE";
            } else if (!containsIgnoreCase(alphapeptdeepInstruments, instrument)) {
                instrument = "QE";
            }
        }
        return instrument;
    }

    private static boolean containsIgnoreCase(List<String> list, String s) {
        return list.stream().anyMatch(e -> e.equalsIgnoreCase(s));
    }

    //the models added with the registry whose limits the former PeptideSkipper had no rule for: they are checked in
    //KoinaPrositLacTest, KoinaChronologerTest and KoinaAltimeterTest, and by their goldens
    private static final Set<String> NEW_LIMITS = Set.of("Prosit_2025_intensity_lac", "Prosit_2025_irt_lac",
            "Chronologer_RT", "Altimeter_2024_intensities");

    @Test
    void theLimitsAreThoseOfTheFormerPeptideSkipper() {
        Set<String> foundUnimods = new HashSet<>(PeptideFormatter.foundUnimods);
        try {
            for (KoinaModelSpec spec : KoinaModelRegistry.all()) {
                if (NEW_LIMITS.contains(spec.name())) {
                    for (String peptide : skipperPeptides()) {
                        for (int charge = 1; charge <= 8; charge++) {
                            assertEquals(PeptideSkipper.skipPeptide(new PeptideFormatter(peptide, charge, "base"), spec),
                                    PeptideSkipper.skipPeptide(new PeptideFormatter(peptide, charge, "base"),
                                            spec.name().toUpperCase()), spec.name() + " " + peptide);
                        }
                    }
                    continue;
                }
                for (String peptide : skipperPeptides()) {
                    for (int charge = 1; charge <= 8; charge++) {
                        String what = spec.name() + " " + peptide + " charge " + charge;
                        boolean expected = formerSkip(new PeptideFormatter(peptide, charge, "base"), spec.name());

                        PeptideFormatter pf = new PeptideFormatter(peptide, charge, "base");
                        String modelPeptide = spec.limits().needsModelPeptide() ? pf.getModel(spec) : null;
                        String reason = spec.limits().skipReason(pf.getStripped(), charge, modelPeptide);
                        assertEquals(expected, reason != null, what + ": " + reason);

                        assertEquals(expected,
                                PeptideSkipper.skipPeptide(new PeptideFormatter(peptide, charge, "base"), spec), what);
                        assertEquals(expected, PeptideSkipper.skipPeptide(
                                new PeptideFormatter(peptide, charge, "base"), spec.name().toUpperCase()), what);
                    }
                }
            }
        } finally {
            PeptideFormatter.foundUnimods.clear();
            PeptideFormatter.foundUnimods.addAll(foundUnimods);
        }
    }

    //the models that are not Koina models keep their rules: the local Prosit predictions take what the Koina Prosit
    //models take, FragCast takes charges 1 to 6, the others take every precursor
    @Test
    void theModelsThatAreNotKoinaModelsKeepTheirLimits() {
        for (String model : new String[] {"Prosit", "PrositTMT", "FragCast", "FragCast-Fast", "DIA-NN",
                "alphapeptdeep", ""}) {
            for (String peptide : skipperPeptides()) {
                for (int charge = 1; charge <= 8; charge++) {
                    assertEquals(formerSkip(new PeptideFormatter(peptide, charge, "base"), model),
                            PeptideSkipper.skipPeptide(new PeptideFormatter(peptide, charge, "base"), model),
                            model + " " + peptide + " charge " + charge);
                }
            }
        }
    }

    //a PSM is predicted by a spectra, an RT and an IM model, and is skipped if any of them cannot take it
    @Test
    void aPsmIsSkippedIfAnyOfItsModelsCannotTakeIt() {
        Set<String> foundUnimods = new HashSet<>(PeptideFormatter.foundUnimods);
        try {
            String[][] modelSets = {
                    {"UniSpec", "Deeplc_hela_hf", "IM2Deep"},
                    {"Prosit_2020_intensity_HCD", "AlphaPept_rt_generic", "AlphaPept_ccs_generic"},
                    {"AlphaPept_ms2_generic", "Deeplc_hela_hf", ""},
                    {"FragCast", "Prosit_2019_irt", ""},
                    {"", "", ""},
                    {"DIA-NN", "DIA-NN", "DIA-NN"}};
            for (String[] models : modelSets) {
                for (String peptide : skipperPeptides()) {
                    for (int charge = 1; charge <= 8; charge++) {
                        boolean expected = false;
                        for (String model : models) {
                            expected |= formerSkip(new PeptideFormatter(peptide, charge, "base"), model);
                        }
                        assertEquals(expected, PeptideSkipper.skipPeptideForAnyModel(
                                        new PeptideFormatter(peptide, charge, "base"), models),
                                List.of(models) + " " + peptide + " charge " + charge);
                    }
                }
            }
        } finally {
            PeptideFormatter.foundUnimods.clear();
            PeptideFormatter.foundUnimods.addAll(foundUnimods);
        }
    }

    private static List<String> skipperPeptides() {
        List<String> peptides = new ArrayList<>(List.of(
                "PEPTIDEK", "M[15.9949]PEPTIDEK", "[42.0106]PEPTIDEK", "PEPTUDEK", "PEPTXDEK", "OPEPTIDEK",
                "PEPTBDEK", "PEPZIDEK", "PEPJIDEK", "PEPTIDEKC[57.0215]"));
        for (int length : new int[] {29, 30, 31, 39, 40, 41, 59, 60, 61}) {
            peptides.add(residues(length));
        }
        //short enough for DeepLC, but not once the mods are written as UniMod codes
        peptides.add("M[15.9949]" + residues(50) + "M[15.9949]M[15.9949]K");
        return peptides;
    }

    //PeptideSkipper.skipPeptide before the registry, as it was, for one model. Its DeepLC peptide comes from
    //getModel, as before; the request goldens (KoinaGoldenTest) check that getModel writes what it wrote
    private static boolean formerSkip(PeptideFormatter pf, String model) {
        model = model.toLowerCase();

        String stripped = pf.getStripped();
        String charge = pf.getCharge();
        String unimodName = "";
        if (model.contains("deeplc")) {
            unimodName = pf.getModel(KoinaModelRegistry.require("Deeplc_hela_hf"));
        }

        //letters
        if (model.contains("prosit") || model.contains("ms2pip") || model.contains("deeplc") ||
                model.contains("unispec") || model.contains("im2deep")) {
            for (char c : "OUBZJX".toCharArray()) {
                if (stripped.indexOf(c) != -1) {
                    return true;
                }
            }
        }
        //length
        if ((model.contains("ms2pip") || model.contains("prosit")) && stripped.length() > 30) {
            return true;
        }
        if (model.contains("unispec") && stripped.length() > 40) {
            return true;
        }
        if ((model.contains("deeplc") || model.contains("im2deep")) && stripped.length() > 60) {
            return true;
        }
        //charge
        int chargeInt = Integer.parseInt(charge);
        if (model.contains("unispec") && chargeInt > 5) {
            return true;
        }
        if (model.contains("prosit") && chargeInt > 6) {
            return true;
        }
        if (model.contains("fragcast") && (chargeInt < 1 || chargeInt > 6)) {
            return true;
        }
        //string length
        if (model.contains("deeplc")) {
            if (unimodName.length() > 60) {
                return true;
            }
        }
        return false;
    }

    //the local Prosit predictions (PinReader.createPrositList and createPrositTMTList, MspReader) write peptides
    //as the Koina Prosit models are sent them, and both write them as getProsit and getPrositTMT did before the
    //registry (base 16df34d)
    @Test
    void theLocalPrositFormatsAreThoseOfTheKoinaPrositModels() {
        Set<String> foundUnimods = new HashSet<>(PeptideFormatter.foundUnimods);
        try {
            KoinaModelSpec prosit = KoinaModelRegistry.require("Prosit_2020_intensity_HCD");
            KoinaModelSpec prositTmt = KoinaModelRegistry.require("Prosit_2020_intensity_TMT");
            //peptide, then getProsit(prositAAMods) and getPrositTMT() of the pre-registry code
            String[][] cases = {
                    {"PEPTIDEK", "PEPTIDEK", "[UNIMOD:737]-PEPTIDEK"},
                    {"PEPC[57.0215]TIDEK", "PEPC[UNIMOD:4]TIDEK", "[UNIMOD:737]-PEPC[UNIMOD:4]TIDEK"},
                    {"PEPCTIDEC", "PEPC[UNIMOD:4]TIDEC[UNIMOD:4]", "[UNIMOD:737]-PEPC[UNIMOD:4]TIDEC[UNIMOD:4]"},
                    {"M[15.9949]PEPTIDEK", "M[UNIMOD:35]PEPTIDEK", "[UNIMOD:737]-M[UNIMOD:35]PEPTIDEK"},
                    {"[42.0106]PEPTIDEK", "PEPTIDEK", "[UNIMOD:737]-PEPTIDEK"},
                    {"[229.1629]PEPTIDEK[229.1629]", "PEPTIDEK", "[UNIMOD:737]-PEPTIDEK[UNIMOD:737]"},
                    {"PEPS[229.1629]TIDEK[229.1629]", "PEPSTIDEK", "[UNIMOD:737]-PEPSTIDEK[UNIMOD:737]"},
                    {"PEPTIDEKc[-0.9840]", "PEPTIDEKc", "[UNIMOD:737]-PEPTIDEKc"},
                    {"C[39.9950]PEPTIDEK", "C[UNIMOD:4]PEPTIDEK", "[UNIMOD:737]-C[UNIMOD:4]PEPTIDEK"},
                    {"S[79.9663]PEPTIDEK", "SPEPTIDEK", "[UNIMOD:737]-SPEPTIDEK"},
                    {"PEPTIDEK[229.1629]", "PEPTIDEK", "[UNIMOD:737]-PEPTIDEK[UNIMOD:737]"}};
            for (String[] c : cases) {
                String peptide = c[0];
                assertEquals(c[1], new PeptideFormatter(peptide, 2, "base").getModel(prosit), peptide);
                assertEquals(c[1], new PeptideFormatter(peptide, 2, "base").getProsit(PTMhandler.prositAAMods),
                        peptide);
                assertEquals(c[2], new PeptideFormatter(peptide, 2, "base").getModel(prositTmt), peptide);
                assertEquals(c[2], new PeptideFormatter(peptide, 2, "base").getPrositTMT(), peptide);
            }
        } finally {
            PeptideFormatter.foundUnimods.clear();
            PeptideFormatter.foundUnimods.addAll(foundUnimods);
        }
    }

    //D27: the NCE range calibration tries is the model's own, unless minNCE or maxNCE is set
    @Test
    void theNceRangeIsTheModelsUnlessSet() {
        Integer savedMin = NceConstants.minNCE;
        Integer savedMax = NceConstants.maxNCE;
        try {
            KoinaModelSpec wide = KoinaModelSpec.builder("Wide", Kind.MS2)
                    .inputs(Input.PRECURSOR_CHARGES, Input.COLLISION_ENERGIES).calibrateCe(true).ceRange(20, 70)
                    .build();
            assertEquals(20, wide.minCe());
            assertEquals(70, wide.maxCe());

            NceConstants.minNCE = null;
            NceConstants.maxNCE = null;
            for (KoinaModelSpec spec : KoinaModelRegistry.all()) {
                assertEquals(spec.minCe(), NceConstants.getMinNCE(spec.name()), spec.name());
                assertEquals(spec.maxCe(), NceConstants.getMaxNCE(spec.name()), spec.name());
                //20-40, but the timsTOF model's collision energies, in eV
                assertEquals(20, spec.minCe(), spec.name());
                assertEquals(spec.name().equals("Prosit_2023_intensity_timsTOF") ? 70 : 40, spec.maxCe(),
                        spec.name());
            }
            assertEquals(20, NceConstants.getMinNCE("Prosit_2023_intensity_timsTOF"));
            assertEquals(70, NceConstants.getMaxNCE("Prosit_2023_intensity_timsTOF"));
            assertEquals(KoinaModelSpec.DEFAULT_MIN_CE, NceConstants.getMinNCE("alphapeptdeep"));
            assertEquals(KoinaModelSpec.DEFAULT_MAX_CE, NceConstants.getMaxNCE("alphapeptdeep"));

            NceConstants.minNCE = 25;
            assertEquals(25, NceConstants.getMinNCE("UniSpec"));
            assertEquals(40, NceConstants.getMaxNCE("UniSpec"));
            assertEquals(25, NceConstants.getMinNCE("Prosit_2023_intensity_timsTOF"));
            assertEquals(70, NceConstants.getMaxNCE("Prosit_2023_intensity_timsTOF"));
            NceConstants.maxNCE = 35;
            assertEquals(35, NceConstants.getMaxNCE("unispec"));
            assertEquals(35, NceConstants.getMaxNCE("Prosit_2023_intensity_timsTOF"),
                    "a maxNCE that is set overrides the timsTOF model's own range too");
        } finally {
            NceConstants.minNCE = savedMin;
            NceConstants.maxNCE = savedMax;
        }
    }

    private static String residues(int length) {
        String residues = "ACDEFGHIKLMNPQRSTVWY";
        StringBuilder sb = new StringBuilder();
        for (int i = 0; i < length; i++) {
            sb.append(residues.charAt(i % residues.length()));
        }
        return sb.toString();
    }

    @Test
    void aSpecRefusesInconsistentFields() {
        assertThrows(IllegalArgumentException.class,
                () -> KoinaModelSpec.builder("x", Kind.MS2).outputs("intensities", "mz").build());
        assertThrows(IllegalArgumentException.class,
                () -> KoinaModelSpec.builder("x", Kind.RT).trustKoinaMz(true).build());
        assertThrows(IllegalArgumentException.class,
                () -> KoinaModelSpec.builder("x", Kind.IM).build());
        assertThrows(IllegalArgumentException.class,
                () -> KoinaModelSpec.builder("x", Kind.RT).annotationParser(AnnotationParser.UNISPEC).build());
        assertThrows(IllegalArgumentException.class,
                () -> KoinaModelSpec.builder("x", Kind.MS2).calibrateCe(true).build());
        assertThrows(IllegalArgumentException.class,
                () -> KoinaModelSpec.builder("x", Kind.MS2).batchSize(0).build());
        assertThrows(IllegalArgumentException.class,
                () -> KoinaModelSpec.builder("x", Kind.MS2).citationKey("").build());
        assertThrows(IllegalArgumentException.class,
                () -> KoinaModelSpec.builder("x", Kind.RT).modRule(KoinaModelSpec.ModRule.ANY_KNOWN).build());

        KoinaModelSpec spec = KoinaModelSpec.builder("x", Kind.MS2)
                .inputs(Input.FRAGMENTATION_TYPES, Input.PRECURSOR_CHARGES).build();
        assertEquals(List.of(Input.PRECURSOR_CHARGES, Input.FRAGMENTATION_TYPES), new ArrayList<>(spec.inputs()),
                "inputs are sent in one order, whatever order they are given in");
        assertEquals("ms2_aux", spec.property(true));
        assertEquals("ms2", spec.property(false));
        assertEquals("rt", KoinaModelSpec.builder("y", Kind.RT).build().property(true));
    }

    @Test
    void requestsAreSplitByTheBatchSizeOfTheSpec(@TempDir Path dir) throws Exception {
        Set<String> entries = new LinkedHashSet<>();
        for (int i = 0; i < 250; i++) {
            entries.add("PEPTIDE" + residues(1 + i % 20) + "K" + i + "," + (2 + i % 2) + ",30,QE,HCD");
        }
        assertEquals(250, entries.size());
        String savedAux = Constants.auxSpectraModel;
        Constants.auxSpectraModel = "";
        ExecutorService executor = Executors.newFixedThreadPool(2);
        try {
            new JSONWriter("UniSpec", entries, false, KoinaModelRegistry.require("UniSpec").batchSize())
                    .write(true, dir.resolve("unispec").toString(), executor);
            new JSONWriter("Prosit_2020_intensity_HCD", entries, false,
                    KoinaModelRegistry.require("Prosit_2020_intensity_HCD").batchSize())
                    .write(true, dir.resolve("prosit").toString(), executor);
        } finally {
            executor.shutdownNow();
            Constants.auxSpectraModel = savedAux;
        }
        assertEquals(List.of(50, 100, 100), requestSizes(dir.resolve("unispec")));
        assertEquals(List.of(250), requestSizes(dir.resolve("prosit")));
    }

    private static List<Integer> requestSizes(Path dir) throws Exception {
        List<Integer> sizes = new ArrayList<>();
        try (Stream<Path> files = Files.list(dir)) {
            for (Path file : files.collect(Collectors.toList())) {
                JSONObject request = new JSONObject(Files.readString(file));
                sizes.add(request.getJSONArray("inputs").getJSONObject(0).getJSONArray("data").length());
            }
        }
        Collections.sort(sizes);
        return sizes;
    }

    private static Set<String> names(Kind kind) {
        return KoinaModelRegistry.all().stream()
                .filter(spec -> spec.kind() == kind)
                .map(KoinaModelSpec::name)
                .collect(Collectors.toCollection(TreeSet::new));
    }
}
