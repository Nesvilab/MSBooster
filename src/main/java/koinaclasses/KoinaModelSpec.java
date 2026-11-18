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

import peptideptmformatting.PTMhandler;
import peptideptmformatting.PeptideFormatter;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.EnumSet;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.concurrent.ConcurrentHashMap;

/**
 * What MSBooster needs to know about one Koina model: what it is sent, what it answers, how peptides are written
 * for it and read back, and which peptides it takes. {@link KoinaModelRegistry} holds one spec per model.
 * <p>
 * A spec is data and named strategies. Building one reads no MSBooster setting, so the registry loads without
 * MainClass or Constants having run (FragPipe loads it from the MSBooster jar to check its model lists). The
 * strategies that need PTMhandler tables reach them only when they are called.
 *
 * @param name                   the Koina model name, as it appears in the URL
 * @param kind                   what the model predicts
 * @param inputs                 the inputs sent besides peptide_sequences, in the order they are sent
 * @param outputs                the names of the outputs used. RT and IM: the prediction. MS2: intensities, m/z
 *                               and annotation, in that order here, whatever order Koina answers them in
 * @param imConversion           IM models: the unit Koina answers in and how it becomes 1/K0; null for others
 * @param peptideFormat          how a peptide is written for the model and read back
 * @param annotationParser       MS2 models: how a fragment annotation is read; null for others
 * @param trustKoinaMz           MS2 models: keep the m/z Koina answers with and the full annotation, instead of
 *                               computing the m/z of the parsed fragment
 * @param modMassMap             the table that turns the model's UniMod codes into masses, and masses into codes
 * @param modRule                the mods the model's peptide format writes; the others are left out of the
 *                               sequence sent, and their mass is added back to the fragment m/z afterwards
 * @param encodableMods          which of those mods the model can be sent, by what Koina reads of them; the
 *                               others are left out as those the rule does not write
 * @param limits                 the peptides and charges the model takes
 * @param batchSize              the most precursors one request carries
 * @param minCe                  the lowest collision energy the calibration tries, unless minNCE is set
 *                               (NceConstants.getMinNCE)
 * @param maxCe                  the highest collision energy the calibration tries, unless maxNCE is set
 *                               (NceConstants.getMaxNCE)
 * @param calibrateCe            whether the collision energy is calibrated for the model
 * @param instrumentMap          how the run's instrument is written for the model
 * @param fragmentations         the fragmentation types the model takes, each with the exact string sent for
 *                               it; empty if the model has no such restriction
 * @param isobaricLabel          the model is trained on isobaric-labelled peptides and labels every peptide
 * @param takesIsobaricLabels    the model can be sent peptides with an isobaric label (TMT or iTRAQ). A run whose
 *                               search has one cannot use a model that cannot (KoinaPreCheck.checkIsobaricLabel)
 * @param sortBatchByLength      the precursors are sorted by stripped peptide length before they are split
 *                               into requests
 * @param singleThreadFormatting the peptides are formatted on one thread, in pin order, because the UniMod code
 *                               a mass gets depends on the codes found before it
 *                               ({@link PeptideFormatter#foundUnimods})
 * @param citationKey            the model's path in koina_citations.yml, without the slashes and /infer
 */
public record KoinaModelSpec(
        String name,
        Kind kind,
        Set<Input> inputs,
        List<String> outputs,
        ImConversion imConversion,
        PeptideFormat peptideFormat,
        AnnotationParser annotationParser,
        boolean trustKoinaMz,
        ModMassMap modMassMap,
        ModRule modRule,
        EncodableMods encodableMods,
        Limits limits,
        int batchSize,
        int minCe,
        int maxCe,
        boolean calibrateCe,
        InstrumentMap instrumentMap,
        Map<String, String> fragmentations,
        boolean isobaricLabel,
        boolean takesIsobaricLabels,
        boolean sortBatchByLength,
        boolean singleThreadFormatting,
        String citationKey) {

    public static final String INTENSITIES = "intensities";
    public static final String MZ = "mz";
    public static final String ANNOTATION = "annotation";
    /** The collision energy range a model calibrates over unless its spec says otherwise. */
    public static final int DEFAULT_MIN_CE = 20;
    public static final int DEFAULT_MAX_CE = 40;

    public KoinaModelSpec {
        if (name == null || name.isBlank()) {
            throw new IllegalArgumentException("A Koina model spec needs a name");
        }
        requireNonNull(name, "kind", kind);
        requireNonNull(name, "peptideFormat", peptideFormat);
        requireNonNull(name, "modMassMap", modMassMap);
        requireNonNull(name, "modRule", modRule);
        requireNonNull(name, "encodableMods", encodableMods);
        requireNonNull(name, "limits", limits);
        requireNonNull(name, "instrumentMap", instrumentMap);
        if (citationKey == null || citationKey.isBlank()) {
            throw new IllegalArgumentException(name + ": no citation key");
        }

        //EnumSet iterates in declaration order, which is the order the inputs are sent in
        EnumSet<Input> inputSet = EnumSet.noneOf(Input.class);
        inputSet.addAll(inputs);
        inputs = Collections.unmodifiableSet(inputSet);
        outputs = List.copyOf(outputs);
        fragmentations = Collections.unmodifiableMap(new LinkedHashMap<>(fragmentations));

        int expectedOutputs = kind == Kind.MS2 ? 3 : 1;
        if (outputs.size() != expectedOutputs) {
            throw new IllegalArgumentException(name + ": a " + kind + " model uses " + expectedOutputs
                    + " output(s), not " + outputs);
        }
        if ((annotationParser != null) != (kind == Kind.MS2)) {
            throw new IllegalArgumentException(name + ": an annotation parser is for MS2 models only, and they need one");
        }
        if ((imConversion != null) != (kind == Kind.IM)) {
            throw new IllegalArgumentException(name + ": an IM conversion is for IM models only, and they need one");
        }
        if (trustKoinaMz && kind != Kind.MS2) {
            throw new IllegalArgumentException(name + ": only MS2 models answer with m/z");
        }
        if (batchSize < 1) {
            throw new IllegalArgumentException(name + ": batch size " + batchSize);
        }
        if (minCe > maxCe) {
            throw new IllegalArgumentException(name + ": collision energy range " + minCe + "-" + maxCe);
        }
        if (calibrateCe && !inputs.contains(Input.COLLISION_ENERGIES)) {
            throw new IllegalArgumentException(name + ": calibrates a collision energy it is not sent");
        }
        if (modRule.remembersFoundMods() && !singleThreadFormatting) {
            throw new IllegalArgumentException(name + ": its mod rule remembers the codes it finds, so its peptides "
                    + "must be formatted on one thread");
        }
    }

    private static void requireNonNull(String name, String field, Object value) {
        if (value == null) {
            throw new IllegalArgumentException(name + ": no " + field);
        }
    }

    /** A builder with the defaults of a model of this kind; see {@link Builder}. */
    public static Builder builder(String name, Kind kind) {
        return new Builder(name, kind);
    }

    public boolean sends(Input input) {
        return inputs.contains(input);
    }

    /**
     * The property name request files and KoinaLibReader use: "rt", "im", "ms2", or "ms2_aux" for an MS2 model
     * that is the auxiliary spectra model.
     */
    public String property(boolean auxSpectraModel) {
        if (kind == Kind.MS2 && auxSpectraModel) {
            return "ms2_aux";
        }
        return kind.property();
    }

    /** RT and IM models: the output that holds the prediction. */
    public String predictionOutput() {
        if (kind == Kind.MS2) {
            throw new IllegalStateException(name + " is an MS2 model");
        }
        return outputs.get(0);
    }

    /** MS2 models: the output that holds the intensities. */
    public String intensitiesOutput() {
        return ms2Output(0);
    }

    /** MS2 models: the output that holds the m/z. */
    public String mzOutput() {
        return ms2Output(1);
    }

    /** MS2 models: the output that holds the fragment annotations. */
    public String annotationOutput() {
        return ms2Output(2);
    }

    private String ms2Output(int i) {
        if (kind != Kind.MS2) {
            throw new IllegalStateException(name + " is not an MS2 model");
        }
        return outputs.get(i);
    }

    /** The fragmentation types the model takes; empty if it takes any. */
    public Set<String> allowedFragmentationTypes() {
        return fragmentations.keySet();
    }

    /**
     * The string the model is sent for this fragmentation type: the one its spec gives, or the type as it is for a
     * model without such a restriction or a type it does not list.
     */
    public String fragmentationSentAs(String type) {
        return fragmentations.getOrDefault(type, type);
    }

    /** What a model predicts. */
    public enum Kind {
        RT("rt"),
        MS2("ms2"),
        IM("im");

        private final String property;

        Kind(String property) {
            this.property = property;
        }

        /** The property name request files and KoinaLibReader use for this kind. */
        public String property() {
            return property;
        }
    }

    /** An input of a Koina request besides peptide_sequences, which every request has. */
    public enum Input {
        PRECURSOR_CHARGES("precursor_charges", "INT32"),
        COLLISION_ENERGIES("collision_energies", "FP32"),
        INSTRUMENT_TYPES("instrument_types", "BYTES"),
        FRAGMENTATION_TYPES("fragmentation_types", "BYTES");

        private final String koinaName;
        private final String datatype;

        Input(String koinaName, String datatype) {
            this.koinaName = koinaName;
            this.datatype = datatype;
        }

        public String koinaName() {
            return koinaName;
        }

        public String datatype() {
            return datatype;
        }
    }

    /** The unit an IM model answers in. */
    public enum ImUnit {
        /** collision cross section, converted to 1/K0 */
        CCS,
        /** 1/K0, used as it is */
        INVERSE_K0
    }

    /**
     * The unit an IM model answers in, and the constants that turn a CCS into 1/K0 by the Mason-Schampp equation:
     * 1/K0 = CCS * sqrt(reduced mass) / charge / ccsCoefficient, with the reduced mass of the precursor and a gas
     * of gasMass.
     */
    public record ImConversion(ImUnit unit, double ccsCoefficient, double gasMass) {
        /**
         * CCS in nitrogen, converted as alphabase does
         * (https://github.com/MannLabs/alphabase/blob/bbaecc380ae157d0f4cc87fffec097ecb7a8ceca/alphabase/peptide/mobility.py#L18)
         */
        public static final ImConversion CCS_N2 = new ImConversion(ImUnit.CCS, 1059.62245, 28.0);
        public static final ImConversion INVERSE_K0 = new ImConversion(ImUnit.INVERSE_K0, 0, 0);

        public ImConversion {
            if (unit == null) {
                throw new IllegalArgumentException("An IM conversion needs a unit");
            }
        }
    }

    /**
     * How a peptide is written for a model (PeptideFormatter.getModel, through
     * PTMhandler.formatPeptideBaseToKoina) and read back (the PeptideFormatter constructor that takes a spec).
     * Every format writes a mod as [UNIMOD:x] after its residue; they differ in how terminal mods are written and
     * in what is added to the peptide once its mods are written.
     */
    public enum PeptideFormat {
        /** N-term mods "[UNIMOD:x]-PEP", C-term mods "PEP-[UNIMOD:x]" */
        ALPHAPEPT("-", "-", false, false, false),
        /** N-term mods "[UNIMOD:x]-PEP", every bare C written as C[UNIMOD:4] */
        PROSIT("-", "", true, false, false),
        /** the Prosit format with an N-term TMT added when there is no N-term mod, and no TMT on S */
        PROSIT_TMT("-", "", true, true, false),
        /** the Prosit format, for the citrullination models */
        PROSIT_CIT("-", "", true, false, false),
        /**
         * N-term mods "[UNIMOD:x]-PEP", and C as it is: a bare C is a free cysteine, C[UNIMOD:4] a
         * carbamidomethylated one (the Prosit-PTM models, which are never sent a C-term mod)
         */
        PROSIT_PTM("-", "", false, false, false),
        /** N-term mods "[UNIMOD:x]-PEP" */
        MS2PIP("-", "", false, false, false),
        /**
         * N-term mods "[UNIMOD:x]-PEP", and C as it is: a bare C is a free cysteine, C[UNIMOD:4] a
         * carbamidomethylated one (Chronologer, which is never sent a C-term mod: Koina reads the letters of one as
         * residues)
         */
        CHRONOLOGER("-", "", false, false, false),
        /** N-term mods "[UNIMOD:x]PEP", without the hyphen */
        UNISPEC("", "", false, false, false),
        /**
         * Terminal mods on the terminal residues: N-term mods "P[UNIMOD:x]EP", C-term mods "PEP[UNIMOD:x]". Koina's
         * DeepLC preprocessing reads the atoms of one mod per residue and nothing before the first residue, so a
         * terminal residue that keeps a mod of its own keeps only that one, and the terminal mod is left out (D24)
         */
        DEEPLC("", "", false, false, true),
        /** the DeepLC format, which Koina's IM2Deep preprocessing reads the same way */
        IM2DEEP("", "", false, false, true);

        private final String ntermSeparator;
        private final String ctermSeparator;
        private final boolean camOnBareCys;
        private final boolean tmtOnNterm;
        private final boolean termModsOnResidues;

        PeptideFormat(String ntermSeparator, String ctermSeparator, boolean camOnBareCys, boolean tmtOnNterm,
                      boolean termModsOnResidues) {
            this.ntermSeparator = ntermSeparator;
            this.ctermSeparator = ctermSeparator;
            this.camOnBareCys = camOnBareCys;
            this.tmtOnNterm = tmtOnNterm;
            this.termModsOnResidues = termModsOnResidues;
        }

        /** What is written between an N-term mod and the first residue. */
        public String ntermSeparator() {
            return ntermSeparator;
        }

        /** What is written between the last residue and a C-term mod. */
        public String ctermSeparator() {
            return ctermSeparator;
        }

        /** Every C without a mod is written as carbamidomethylated, C[UNIMOD:4]. */
        public boolean camOnBareCys() {
            return camOnBareCys;
        }

        /**
         * A peptide without an N-term mod is written with an N-term TMT ([UNIMOD:737]-), and S is written without
         * TMT.
         */
        public boolean tmtOnNterm() {
            return tmtOnNterm;
        }

        /**
         * An N-term mod is written on the first residue and a C-term mod on the last, as a mod of that residue. A
         * terminal residue with a mod of its own that is written keeps it, and the terminal mod is left out.
         * Reading such a peptide back, a mod on the first residue that the mass table has as an N-term mod and not
         * as a mod of that residue is read as the N-term mod.
         */
        public boolean termModsOnResidues() {
            return termModsOnResidues;
        }
    }

    /** How a fragment annotation in a Koina response is read into ion type, fragment number and charge. */
    public enum AnnotationParser {
        /** "b3+1", "y12+2": ion type b or y, fragment number, "+", charge (Prosit, ms2pip, AlphaPeptDeep) */
        ION_NUMBER_CHARGE {
            @Override
            public void parse(String annotation, List<Integer> charges, List<String> ionTypes,
                              List<Integer> fragNums) {
                KoinaOutputParsingMethods.parseIonNumberCharge(annotation, charges, ionTypes, fragNums);
            }
        },
        /** UniSpec: b/y/a ions with neutral losses, precursor, internal and immonium ions, charge after "^" */
        UNISPEC {
            @Override
            public void parse(String annotation, List<Integer> charges, List<String> ionTypes,
                              List<Integer> fragNums) {
                KoinaOutputParsingMethods.parseUnispec(annotation, charges, ionTypes, fragNums);
            }
        },
        /** Prosit 2025 MultiFrag: a, b, c, x, y, z ions and their radical/hydrogen-shifted forms, charge after "^" */
        MULTIFRAG {
            @Override
            public void parse(String annotation, List<Integer> charges, List<String> ionTypes,
                              List<Integer> fragNums) {
                KoinaOutputParsingMethods.parsePrositMultifrag(annotation, charges, ionTypes, fragNums);
            }
        },
        /**
         * Altimeter (mzPAF): b and y ions with the charge after "^" ("y5^2"), the precursor ("p^3") and immonium
         * ions ("IY", "IC(Carbamidomethyl)", "IH-NH3+CO"). A loss is read as UniSpec's are: on a b or y ion it
         * makes b-NL or y-NL, while a precursor stays p and an immonium ion imm
         */
        ALTIMETER {
            @Override
            public void parse(String annotation, List<Integer> charges, List<String> ionTypes,
                              List<Integer> fragNums) {
                KoinaOutputParsingMethods.parseAltimeter(annotation, charges, ionTypes, fragNums);
            }
        };

        /**
         * Adds the charge, ion type and fragment number of one annotation to the lists.
         *
         * @throws IllegalArgumentException for an annotation this parser does not know
         */
        public abstract void parse(String annotation, List<Integer> charges, List<String> ionTypes,
                                   List<Integer> fragNums);
    }

    /** The table that turns a model's UniMod codes into masses (reading) and masses into codes (writing). */
    public enum ModMassMap {
        /** PTMhandler.unimodToModMassLimited and AAunimodToModMassLimited: a short list of common mods */
        LIMITED {
            @Override
            public Map<Integer, Double> unimodToMass() {
                return PTMhandler.unimodToModMassLimited;
            }

            @Override
            public Map<String, Double> aaUnimodToMass() {
                return PTMhandler.AAunimodToModMassLimited;
            }
        },
        /**
         * PTMhandler.unimodToModMassAll and AAunimodToModMassAll: alphapept_koina.csv,
         * modification_alphapeptdeep.tsv and the additional mods file
         */
        ALL {
            @Override
            public Map<Integer, Double> unimodToMass() {
                return PTMhandler.unimodToModMassAll;
            }

            @Override
            public Map<String, Double> aaUnimodToMass() {
                return PTMhandler.AAunimodToModMassAll;
            }
        },
        /**
         * {@link #ALL}, with the two pairs of the hand-written table of Koina's Prosit-PTM models that it lacks:
         * phospho on P (P21) and Gly on C (C1263) (PTMhandler.getAAunimodToModMassPrositPtm)
         */
        PROSIT_PTM {
            @Override
            public Map<Integer, Double> unimodToMass() {
                return PTMhandler.unimodToModMassAll;
            }

            @Override
            public Map<String, Double> aaUnimodToMass() {
                return PTMhandler.getAAunimodToModMassPrositPtm();
            }
        },
        /**
         * PTMhandler.unimodToModMassLac and AAunimodToModMassLac: lactyl (UNIMOD:2114) on K and the N-terminus,
         * M35 and C4. A table of their own, as the ALL table gives K the isobaric carboxyethyl (UNIMOD:378)
         */
        LAC {
            @Override
            public Map<Integer, Double> unimodToMass() {
                return PTMhandler.unimodToModMassLac;
            }

            @Override
            public Map<String, Double> aaUnimodToMass() {
                return PTMhandler.AAunimodToModMassLac;
            }
        },
        /**
         * PTMhandler.unimodToModMassChronologer and AAunimodToModMassChronologer: the pairs of Koina's Chronologer
         * table, with UNIMOD:385 on C standing for the whole pyro-carbamidomethyl cysteine (C[39.9950]), as
         * Chronologer reads it
         */
        CHRONOLOGER {
            @Override
            public Map<Integer, Double> unimodToMass() {
                return PTMhandler.unimodToModMassChronologer;
            }

            @Override
            public Map<String, Double> aaUnimodToMass() {
                return PTMhandler.AAunimodToModMassChronologer;
            }
        };

        /** UniMod code to mod mass, for reading the model's peptides. */
        public abstract Map<Integer, Double> unimodToMass();

        /** Residue (or "[" for a terminus) plus UniMod code to mod mass, for writing the model's peptides. */
        public abstract Map<String, Double> aaUnimodToMass();
    }

    /**
     * The mods a model's peptide format writes, as residue (or "[" for a terminus) plus UniMod code. A mass is
     * first matched against the preferred codes, then against the allowed ones; a mod that matches neither is left
     * out of the sequence sent.
     */
    public enum ModRule {
        /** C4 and M35 */
        PROSIT {
            @Override
            public Set<String> preferredMods() {
                return PTMhandler.prositAAMods;
            }

            @Override
            public Set<String> allowedMods() {
                return PTMhandler.prositAAMods;
            }
        },
        /** TMT (as UniMod 737) on the N-term and K first, then the Prosit mods */
        PROSIT_TMT {
            @Override
            public Set<String> preferredMods() {
                return PTMhandler.prositTmtAAMods;
            }

            @Override
            public Set<String> allowedMods() {
                return PTMhandler.prositAAMods;
            }
        },
        /** citrullination/deamidation (UniMod 7) on N, Q and R first, then the Prosit mods */
        PROSIT_CIT {
            @Override
            public Set<String> preferredMods() {
                return PTMhandler.prositCitAAMods;
            }

            @Override
            public Set<String> allowedMods() {
                return PTMhandler.prositAAMods;
            }
        },
        /** M35 and C4 */
        MS2PIP {
            @Override
            public Set<String> preferredMods() {
                return PTMhandler.ms2pipAAMods;
            }

            @Override
            public Set<String> allowedMods() {
                return PTMhandler.ms2pipAAMods;
            }
        },
        /** N-term acetyl, C4, C26, Q28, E27, M35 and phospho on S, T, Y */
        UNISPEC {
            @Override
            public Set<String> preferredMods() {
                return PTMhandler.unispecAAMods;
            }

            @Override
            public Set<String> allowedMods() {
                return PTMhandler.unispecAAMods;
            }
        },
        /**
         * every mod of the ALL mass table; the codes already found in this run first, and each code found is
         * remembered for the peptides formatted after it
         */
        ANY_KNOWN {
            @Override
            public Set<String> preferredMods() {
                return PeptideFormatter.foundUnimods;
            }

            @Override
            public Set<String> allowedMods() {
                return PTMhandler.AAunimodToModMassAllKeys;
            }

            @Override
            public boolean remembersFoundMods() {
                return true;
            }
        },
        /**
         * the pairs of the hand-written table of Koina's Prosit-PTM models first
         * ({@link KoinaPrositPtm#preferredMods}), then every mod of the {@link ModMassMap#PROSIT_PTM} table. Which of
         * them the models can be sent is up to {@link EncodableMods#PROSIT_PTM}
         */
        PROSIT_PTM {
            @Override
            public Set<String> preferredMods() {
                return KoinaPrositPtm.preferredMods();
            }

            @Override
            public Set<String> allowedMods() {
                return PTMhandler.getAAunimodToModMassPrositPtm().keySet();
            }
        },
        /** the mods of the {@link ModMassMap#LAC} table: lactyl on K and the N-terminus, M35 and C4 */
        LAC {
            @Override
            public Set<String> preferredMods() {
                return PTMhandler.AAunimodToModMassLac.keySet();
            }

            @Override
            public Set<String> allowedMods() {
                return PTMhandler.AAunimodToModMassLac.keySet();
            }
        },
        /** the mods of the {@link ModMassMap#CHRONOLOGER} table */
        CHRONOLOGER {
            @Override
            public Set<String> preferredMods() {
                return PTMhandler.AAunimodToModMassChronologer.keySet();
            }

            @Override
            public Set<String> allowedMods() {
                return PTMhandler.AAunimodToModMassChronologer.keySet();
            }
        };

        public abstract Set<String> preferredMods();

        public abstract Set<String> allowedMods();

        /**
         * Whether each code a mass is given goes into {@link PeptideFormatter#foundUnimods}, so the peptides
         * formatted after it give the same mass the same code. A model with such a rule formats its peptides on
         * one thread ({@link KoinaModelSpec#singleThreadFormatting}).
         */
        public boolean remembersFoundMods() {
            return false;
        }
    }

    /**
     * Which mods a model can be sent, by what Koina reads of them. A mod the model's {@link ModRule} writes is
     * written only if it is encodable; one that is not is left out of the peptide sent, as a mod the rule does not
     * write, and counted among the mods left out (KoinaPreCheck).
     */
    public enum EncodableMods {
        /** every mod the rule writes */
        ANY(false, true),
        /**
         * every mod the rule writes but a C-term mod. For a model whose table has N-term mods and whose peptide
         * format has no place for a C-term mod: the rule would otherwise write a C-term mod of the same mass as an
         * N-term one after the last residue, where the model reads it as a mod of that residue
         */
        NO_CTERM(false, false) {
            @Override
            public boolean encodes(String aaCode, boolean cterm) {
                return !cterm;
            }
        },
        /**
         * mods whose composition in Koina's unimod.obo ({@link KoinaUnimod}) has only C, H, N, O, S and P atoms,
         * isotopes counting as their element (DeepLC). Koina's DeepLC preprocessing reads any other symbol, such as
         * Na or Hex, as the first of these letters in it, or fails on it, and fails on a code it has no
         * composition for
         */
        KOINA_ELEMENTS(true, true),
        /**
         * {@link #KOINA_ELEMENTS}, and no isotope atoms (13C, 15N, 2H, 18O) either: Koina's IM2Deep preprocessing
         * takes the element from the first character of each token, so it fails on an isotope (IM2Deep)
         */
        KOINA_ELEMENTS_NO_ISOTOPES(true, false),
        /**
         * the mods Koina's Prosit-PTM models encode ({@link KoinaPrositPtm}): a pair of residue (or N-terminus) and
         * code of their hand-written table, or one their zero-shot encoding takes. Never a C-term mod: Koina reads
         * it as one on the N-terminus
         */
        PROSIT_PTM(false, false) {
            @Override
            public boolean encodes(String aaCode, boolean cterm) {
                return !cterm && KoinaPrositPtm.encodes(aaCode);
            }
        };

        private final boolean byComposition;
        private final boolean isotopesAllowed;
        //UniMod code to whether it is encodable, filled as codes are asked about
        private final Map<String, Boolean> known = new ConcurrentHashMap<>();

        EncodableMods(boolean byComposition, boolean isotopesAllowed) {
            this.byComposition = byComposition;
            this.isotopesAllowed = isotopesAllowed;
        }

        /**
         * Whether the model can be sent this mod on a residue or the N-terminus.
         *
         * @param aaCode the residue (or "[" for a terminus) and UniMod code, as in the mod sets of a {@link ModRule}
         */
        public boolean encodes(String aaCode) {
            return encodes(aaCode, false);
        }

        /**
         * Whether the model can be sent this mod.
         *
         * @param aaCode the residue (or "[" for a terminus) and UniMod code, as in the mod sets of a {@link ModRule}
         * @param cterm  whether it is a C-term mod
         */
        public boolean encodes(String aaCode, boolean cterm) {
            if (!byComposition) {
                return true;
            }
            return known.computeIfAbsent(aaCode.substring(1),
                    code -> KoinaUnimod.hasOnlyEncodedElements(code, isotopesAllowed));
        }
    }

    /**
     * The peptides and charges a model takes. A peptide outside them is not sent to the model.
     *
     * @param minLength             the shortest stripped peptide
     * @param maxLength             the longest stripped peptide
     * @param minCharge             the lowest precursor charge
     * @param maxCharge             the highest precursor charge
     * @param forbiddenResidues     residues the model cannot take
     * @param maxModelPeptideLength the longest peptide string, as written for the model
     * @param skipFreeCys           a peptide with a cysteine that has no mod is not sent, because the model reads
     *                              every C as carbamidomethylated or cannot take a free one
     */
    public record Limits(int minLength, int maxLength, int minCharge, int maxCharge, String forbiddenResidues,
                         int maxModelPeptideLength, boolean skipFreeCys) {
        /** no limit */
        public static final Limits NONE = new Limits(0, Integer.MAX_VALUE, 0, Integer.MAX_VALUE, "",
                Integer.MAX_VALUE, false);
        /** the residues none of the models but AlphaPeptDeep takes */
        public static final String NONSTANDARD_RESIDUES = "OUBZJX";
        /** the reason {@link #skipReason} gives for a peptide with a free cysteine */
        public static final String FREE_CYS = "free cysteine";

        public Limits {
            if (forbiddenResidues == null) {
                throw new IllegalArgumentException("Limits need forbidden residues, even if none");
            }
        }

        public Limits withLength(int minLength, int maxLength) {
            return new Limits(minLength, maxLength, minCharge, maxCharge, forbiddenResidues, maxModelPeptideLength,
                    skipFreeCys);
        }

        public Limits withCharge(int minCharge, int maxCharge) {
            return new Limits(minLength, maxLength, minCharge, maxCharge, forbiddenResidues, maxModelPeptideLength,
                    skipFreeCys);
        }

        public Limits withForbiddenResidues(String forbiddenResidues) {
            return new Limits(minLength, maxLength, minCharge, maxCharge, forbiddenResidues, maxModelPeptideLength,
                    skipFreeCys);
        }

        public Limits withMaxModelPeptideLength(int maxModelPeptideLength) {
            return new Limits(minLength, maxLength, minCharge, maxCharge, forbiddenResidues, maxModelPeptideLength,
                    skipFreeCys);
        }

        public Limits withSkipFreeCys(boolean skipFreeCys) {
            return new Limits(minLength, maxLength, minCharge, maxCharge, forbiddenResidues, maxModelPeptideLength,
                    skipFreeCys);
        }

        /** Whether {@link #skipReason} needs the peptide as written for the model. */
        public boolean needsModelPeptide() {
            return maxModelPeptideLength != Integer.MAX_VALUE;
        }

        /**
         * Why the model cannot take this precursor, or null if it can, for limits that do not skip free cysteines.
         * Limits that do need the peptide's mods: {@link #skipReason(String, String, int, String)}.
         *
         * @param stripped     the peptide without mods
         * @param charge       the precursor charge
         * @param modelPeptide the peptide as written for the model; read only if {@link #needsModelPeptide()}
         */
        public String skipReason(String stripped, int charge, String modelPeptide) {
            if (skipFreeCys) {
                throw new IllegalStateException("These limits skip free cysteines, so they need the peptide with "
                        + "its mods");
            }
            return skipReason(null, stripped, charge, modelPeptide);
        }

        /**
         * Why the model cannot take this precursor, or null if it can.
         *
         * @param base         the peptide with its mods in mass notation (PeptideFormatter.getBase), without its
         *                     C-term mod; read only if {@link #skipFreeCys()}
         * @param stripped     the peptide without mods
         * @param charge       the precursor charge
         * @param modelPeptide the peptide as written for the model; read only if {@link #needsModelPeptide()}
         */
        public String skipReason(String base, String stripped, int charge, String modelPeptide) {
            for (int i = 0; i < forbiddenResidues.length(); i++) {
                char c = forbiddenResidues.charAt(i);
                if (stripped.indexOf(c) != -1) {
                    return "residue " + c;
                }
            }
            if (stripped.length() < minLength) {
                return "shorter than " + minLength + " residues";
            }
            if (stripped.length() > maxLength) {
                return "longer than " + maxLength + " residues";
            }
            if (charge < minCharge) {
                return "charge below " + minCharge;
            }
            if (charge > maxCharge) {
                return "charge above " + maxCharge;
            }
            if (skipFreeCys && hasFreeCys(base)) {
                return FREE_CYS;
            }
            if (needsModelPeptide() && modelPeptide.length() > maxModelPeptideLength) {
                return "longer than " + maxModelPeptideLength + " characters as written for the model";
            }
            return null;
        }

        //whether a C of the peptide, in mass notation, has no mod after it
        private static boolean hasFreeCys(String base) {
            for (int i = 0; i < base.length(); i++) {
                if (base.charAt(i) == 'C' && (i + 1 == base.length() || base.charAt(i + 1) != '[')) {
                    return true;
                }
            }
            return false;
        }
    }

    /** How the run's instrument is written for a model that is sent one. */
    public enum InstrumentMap {
        /** the run's instrument */
        AS_IS {
            @Override
            public String map(String instrument, String charge) {
                return instrument;
            }
        },
        /** QE, LUMOS, TIMSTOF or SCIEXTOF; QEHFX and anything else become QE */
        ALPHAPEPT {
            @Override
            public String map(String instrument, String charge) {
                if (instrument.equals("QEHFX")) {
                    return "QE";
                } else if (!containsIgnoreCase(ALPHAPEPT_INSTRUMENTS, instrument)) {
                    return "QE";
                }
                return instrument;
            }
        },
        /** LUMOS, QE, QEHFX, ELITE, VELOS or NONE; anything else, and charge 1 on QE or QEHFX, become NONE */
        UNISPEC {
            @Override
            public String map(String instrument, String charge) {
                if (!containsIgnoreCase(UNISPEC_INSTRUMENTS, instrument)) {
                    return "NONE";
                } else if ((instrument.equals("QE") || instrument.equals("QEHFX")) && charge.equals("1")) {
                    return "NONE";
                }
                return instrument;
            }
        },
        /**
         * the Prosit lactylation model (D21): ASTRAL for an Astral (THERMOTOF), LUMOS for anything else. Koina reads
         * LUMOS, ASTRAL and ECLIPSE, upper case, and any other string as LUMOS
         */
        LAC {
            @Override
            public String map(String instrument, String charge) {
                return "THERMOTOF".equalsIgnoreCase(instrument) ? "ASTRAL" : "LUMOS";
            }
        };

        private static final List<String> ALPHAPEPT_INSTRUMENTS = List.of("QE", "LUMOS", "TIMSTOF", "SCIEXTOF");
        private static final List<String> UNISPEC_INSTRUMENTS =
                List.of("LUMOS", "QE", "QEHFX", "ELITE", "VELOS", "NONE");

        /**
         * @param instrument the run's instrument (Constants.instrument)
         * @param charge     the precursor charge
         */
        public abstract String map(String instrument, String charge);

        private static boolean containsIgnoreCase(List<String> instruments, String instrument) {
            for (String s : instruments) {
                if (s.equalsIgnoreCase(instrument)) {
                    return true;
                }
            }
            return false;
        }
    }

    /**
     * Builds a spec. The defaults are those of the plainest model of its kind: RT models are sent the peptide only
     * and answer "irt"; IM models are also sent the charge and answer "ccs" (the IM conversion has no default);
     * MS2 models are sent the charge, answer intensities, mz and annotation, and annotate fragments as "y3+1".
     * Every model defaults to the Prosit peptide format and mods (all of them encodable), no limits, 1000
     * precursors per request, a collision energy range of 20-40 that is not calibrated, the run's instrument, no
     * fragmentation restriction, isobaric labels taken, and its own name as citation key.
     */
    public static final class Builder {
        private final String name;
        private final Kind kind;
        private List<Input> inputs = new ArrayList<>();
        private List<String> outputs;
        private ImConversion imConversion;
        private PeptideFormat peptideFormat = PeptideFormat.PROSIT;
        private AnnotationParser annotationParser;
        private boolean trustKoinaMz = false;
        private ModMassMap modMassMap = ModMassMap.LIMITED;
        private ModRule modRule = ModRule.PROSIT;
        private EncodableMods encodableMods = EncodableMods.ANY;
        private Limits limits = Limits.NONE;
        private int batchSize = 1000;
        private int minCe = DEFAULT_MIN_CE;
        private int maxCe = DEFAULT_MAX_CE;
        private boolean calibrateCe = false;
        private InstrumentMap instrumentMap = InstrumentMap.AS_IS;
        private final Map<String, String> fragmentations = new LinkedHashMap<>();
        private boolean isobaricLabel = false;
        private boolean takesIsobaricLabels = true;
        private boolean sortBatchByLength = false;
        private boolean singleThreadFormatting = false;
        private String citationKey;

        private Builder(String name, Kind kind) {
            this.name = name;
            this.kind = kind;
            this.citationKey = name;
            switch (kind) {
                case RT:
                    outputs = List.of("irt");
                    break;
                case IM:
                    inputs.add(Input.PRECURSOR_CHARGES);
                    outputs = List.of("ccs");
                    break;
                case MS2:
                    inputs.add(Input.PRECURSOR_CHARGES);
                    outputs = List.of(INTENSITIES, MZ, ANNOTATION);
                    annotationParser = AnnotationParser.ION_NUMBER_CHARGE;
                    break;
            }
        }

        public Builder inputs(Input... inputs) {
            this.inputs = new ArrayList<>(Arrays.asList(inputs));
            return this;
        }

        public Builder outputs(String... outputs) {
            this.outputs = List.of(outputs);
            return this;
        }

        public Builder imConversion(ImConversion imConversion) {
            this.imConversion = imConversion;
            return this;
        }

        public Builder peptideFormat(PeptideFormat peptideFormat) {
            this.peptideFormat = peptideFormat;
            return this;
        }

        public Builder annotationParser(AnnotationParser annotationParser) {
            this.annotationParser = annotationParser;
            return this;
        }

        public Builder trustKoinaMz(boolean trustKoinaMz) {
            this.trustKoinaMz = trustKoinaMz;
            return this;
        }

        public Builder modMassMap(ModMassMap modMassMap) {
            this.modMassMap = modMassMap;
            return this;
        }

        public Builder modRule(ModRule modRule) {
            this.modRule = modRule;
            return this;
        }

        public Builder encodableMods(EncodableMods encodableMods) {
            this.encodableMods = encodableMods;
            return this;
        }

        public Builder limits(Limits limits) {
            this.limits = limits;
            return this;
        }

        public Builder batchSize(int batchSize) {
            this.batchSize = batchSize;
            return this;
        }

        public Builder ceRange(int minCe, int maxCe) {
            this.minCe = minCe;
            this.maxCe = maxCe;
            return this;
        }

        public Builder calibrateCe(boolean calibrateCe) {
            this.calibrateCe = calibrateCe;
            return this;
        }

        public Builder instrumentMap(InstrumentMap instrumentMap) {
            this.instrumentMap = instrumentMap;
            return this;
        }

        /** Fragmentation types the model takes, each sent as it is named. */
        public Builder fragmentations(String... types) {
            for (String type : types) {
                fragmentations.put(type, type);
            }
            return this;
        }

        /** A fragmentation type the model takes, and the exact string sent for it. */
        public Builder fragmentation(String type, String sentAs) {
            fragmentations.put(type, sentAs);
            return this;
        }

        public Builder isobaricLabel(boolean isobaricLabel) {
            this.isobaricLabel = isobaricLabel;
            return this;
        }

        public Builder takesIsobaricLabels(boolean takesIsobaricLabels) {
            this.takesIsobaricLabels = takesIsobaricLabels;
            return this;
        }

        public Builder sortBatchByLength(boolean sortBatchByLength) {
            this.sortBatchByLength = sortBatchByLength;
            return this;
        }

        public Builder singleThreadFormatting(boolean singleThreadFormatting) {
            this.singleThreadFormatting = singleThreadFormatting;
            return this;
        }

        public Builder citationKey(String citationKey) {
            this.citationKey = citationKey;
            return this;
        }

        public KoinaModelSpec build() {
            return new KoinaModelSpec(name, kind, inputs.isEmpty() ? EnumSet.noneOf(Input.class)
                    : EnumSet.copyOf(inputs), outputs, imConversion, peptideFormat, annotationParser,
                    trustKoinaMz, modMassMap, modRule, encodableMods, limits, batchSize, minCe, maxCe, calibrateCe,
                    instrumentMap, fragmentations, isobaricLabel, takesIsobaricLabels, sortBatchByLength,
                    singleThreadFormatting, citationKey);
        }
    }
}
