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

import koinaclasses.KoinaModelSpec.AnnotationParser;
import koinaclasses.KoinaModelSpec.Builder;
import koinaclasses.KoinaModelSpec.EncodableMods;
import koinaclasses.KoinaModelSpec.ImConversion;
import koinaclasses.KoinaModelSpec.InstrumentMap;
import koinaclasses.KoinaModelSpec.Kind;
import koinaclasses.KoinaModelSpec.Limits;
import koinaclasses.KoinaModelSpec.ModMassMap;
import koinaclasses.KoinaModelSpec.ModRule;
import koinaclasses.KoinaModelSpec.PeptideFormat;

import java.util.Collections;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Locale;
import java.util.Map;

import static koinaclasses.KoinaModelSpec.Input.COLLISION_ENERGIES;
import static koinaclasses.KoinaModelSpec.Input.FRAGMENTATION_TYPES;
import static koinaclasses.KoinaModelSpec.Input.INSTRUMENT_TYPES;
import static koinaclasses.KoinaModelSpec.Input.PRECURSOR_CHARGES;

/**
 * The Koina models MSBooster can call, one {@link KoinaModelSpec} each. The Koina model lists in
 * ModelCollections, NceConstants and LowercaseModelMapper are made from it.
 * <p>
 * Loading this class reads no MSBooster setting, so FragPipe can load it from the MSBooster jar on its own.
 */
public final class KoinaModelRegistry {
    private static final List<KoinaModelSpec> ALL = List.of(
            //RT
            alphapept("AlphaPept_rt_generic", Kind.RT)
                    .citationKey("AlphaPeptDeep_rt_generic")
                    .build(),
            prosit("Prosit_2019_irt", Kind.RT)
                    .build(),
            prosit("Prosit_2020_irt_TMT", Kind.RT)
                    .peptideFormat(PeptideFormat.PROSIT_TMT).modRule(ModRule.PROSIT_TMT)
                    .isobaricLabel(true)
                    .build(),
            prosit("Prosit_2024_irt_cit", Kind.RT)
                    .peptideFormat(PeptideFormat.PROSIT_CIT).modRule(ModRule.PROSIT_CIT)
                    .build(),
            prositPtm("Prosit_2025_irt_22PTM", Kind.RT)
                    .build(),
            prositLac("Prosit_2025_irt_lac", Kind.RT)
                    .build(),
            KoinaModelSpec.builder("Chronologer_RT", Kind.RT)
                    .outputs("rt")
                    .peptideFormat(PeptideFormat.CHRONOLOGER).modMassMap(ModMassMap.CHRONOLOGER)
                    .modRule(ModRule.CHRONOLOGER).encodableMods(EncodableMods.NO_CTERM)
                    //Koina reads the first 50 residues of a longer peptide and leaves out the rest
                    .limits(Limits.NONE.withForbiddenResidues(Limits.NONSTANDARD_RESIDUES).withLength(0, 50))
                    .build(),
            KoinaModelSpec.builder("Deeplc_hela_hf", Kind.RT)
                    .peptideFormat(PeptideFormat.DEEPLC).modMassMap(ModMassMap.ALL).modRule(ModRule.ANY_KNOWN)
                    .encodableMods(EncodableMods.KOINA_ELEMENTS)
                    .limits(Limits.NONE.withForbiddenResidues(Limits.NONSTANDARD_RESIDUES).withLength(0, 60)
                            .withMaxModelPeptideLength(60))
                    .singleThreadFormatting(true)
                    .build(),

            //MS2
            ms2pip("ms2pip_2021_HCD")
                    .citationKey("ms2pip_HCD2021")
                    .build(),
            ms2pip("ms2pip_timsTOF2024")
                    .build(),
            ms2pip("ms2pip_timsTOF2023")
                    .build(),
            ms2pip("ms2pip_CID_TMT")
                    .build(),
            ms2pip("ms2pip_TTOF5600")
                    .build(),
            ms2pip("ms2pip_Immuno_HCD")
                    .build(),
            ms2pip("ms2pip_iTRAQphospho")
                    .build(),
            alphapept("AlphaPept_ms2_generic", Kind.MS2)
                    .inputs(PRECURSOR_CHARGES, COLLISION_ENERGIES, INSTRUMENT_TYPES).calibrateCe(true)
                    .batchSize(500)
                    .citationKey("AlphaPeptDeep_ms2_generic")
                    .build(),
            prosit("Prosit_2019_intensity", Kind.MS2)
                    .inputs(PRECURSOR_CHARGES, COLLISION_ENERGIES).calibrateCe(true)
                    .build(),
            prosit("Prosit_2020_intensity_CID", Kind.MS2)
                    .inputs(PRECURSOR_CHARGES)
                    .build(),
            prosit("Prosit_2020_intensity_TMT", Kind.MS2)
                    .inputs(PRECURSOR_CHARGES, COLLISION_ENERGIES, FRAGMENTATION_TYPES).calibrateCe(true)
                    .fragmentations("HCD", "CID")
                    .peptideFormat(PeptideFormat.PROSIT_TMT).modRule(ModRule.PROSIT_TMT)
                    .isobaricLabel(true)
                    .build(),
            prosit("Prosit_2020_intensity_HCD", Kind.MS2)
                    .inputs(PRECURSOR_CHARGES, COLLISION_ENERGIES).calibrateCe(true)
                    .build(),
            prosit("Prosit_2023_intensity_timsTOF", Kind.MS2)
                    .inputs(PRECURSOR_CHARGES, COLLISION_ENERGIES).calibrateCe(true)
                    .ceRange(20, 70) //a timsTOF collision energy, in eV
                    .build(),
            prosit("Prosit_2024_intensity_cit", Kind.MS2)
                    .inputs(PRECURSOR_CHARGES, COLLISION_ENERGIES, FRAGMENTATION_TYPES).calibrateCe(true)
                    .fragmentations("HCD", "CID")
                    .peptideFormat(PeptideFormat.PROSIT_CIT).modRule(ModRule.PROSIT_CIT)
                    .build(),
            prosit("Prosit_2025_intensity_MultiFrag", Kind.MS2)
                    .inputs(PRECURSOR_CHARGES, FRAGMENTATION_TYPES)
                    .fragmentations("HCD", "ECD", "EID", "UVPD", "ETCID")
                    .annotationParser(AnnotationParser.MULTIFRAG)
                    .batchSize(100)
                    .build(),
            prositPtm("Prosit_2025_intensity_22PTM", Kind.MS2)
                    .inputs(PRECURSOR_CHARGES, COLLISION_ENERGIES, FRAGMENTATION_TYPES).calibrateCe(true)
                    .fragmentations("HCD", "CID")
                    .build(),
            prositPtm("Prosit_2025_intensity_40PTM", Kind.MS2)
                    .inputs(PRECURSOR_CHARGES, COLLISION_ENERGIES, FRAGMENTATION_TYPES).calibrateCe(true)
                    .fragmentations("HCD", "CID")
                    .build(),
            prositLac("Prosit_2025_intensity_lac", Kind.MS2)
                    .inputs(PRECURSOR_CHARGES, COLLISION_ENERGIES, INSTRUMENT_TYPES, FRAGMENTATION_TYPES)
                    .calibrateCe(true)
                    .instrumentMap(InstrumentMap.LAC)
                    .fragmentations("HCD", "CID")
                    .build(),
            //the intensities of b, y, precursor and immonium ions at an NCE of 20-40, with their m/z, 200 per
            //precursor; oxidized M and carbamidomethylated C only, every C carbamidomethylated
            KoinaModelSpec.builder("Altimeter_2024_intensities", Kind.MS2)
                    .inputs(PRECURSOR_CHARGES, COLLISION_ENERGIES).calibrateCe(true)
                    .outputs(KoinaModelSpec.INTENSITIES, KoinaModelSpec.MZ, "annotations")
                    .peptideFormat(PeptideFormat.PROSIT).modMassMap(ModMassMap.LIMITED).modRule(ModRule.PROSIT)
                    .limits(Limits.NONE.withForbiddenResidues(Limits.NONSTANDARD_RESIDUES).withLength(6, 40)
                            .withCharge(1, 7).withSkipFreeCys(true))
                    .annotationParser(AnnotationParser.ALTIMETER).trustKoinaMz(true)
                    .batchSize(100)
                    .build(),
            KoinaModelSpec.builder("UniSpec", Kind.MS2)
                    .inputs(PRECURSOR_CHARGES, COLLISION_ENERGIES, INSTRUMENT_TYPES).calibrateCe(true)
                    .instrumentMap(InstrumentMap.UNISPEC)
                    .peptideFormat(PeptideFormat.UNISPEC).modMassMap(ModMassMap.LIMITED).modRule(ModRule.UNISPEC)
                    .limits(Limits.NONE.withForbiddenResidues(Limits.NONSTANDARD_RESIDUES).withLength(0, 40)
                            .withCharge(0, 5))
                    .annotationParser(AnnotationParser.UNISPEC).trustKoinaMz(true)
                    .batchSize(100)
                    .build(),

            //IM
            alphapept("AlphaPept_ccs_generic", Kind.IM)
                    .imConversion(ImConversion.CCS_N2)
                    .citationKey("AlphaPeptDeep_ccs_generic")
                    .build(),
            KoinaModelSpec.builder("IM2Deep", Kind.IM)
                    .imConversion(ImConversion.CCS_N2)
                    .peptideFormat(PeptideFormat.IM2DEEP).modMassMap(ModMassMap.ALL).modRule(ModRule.ANY_KNOWN)
                    .encodableMods(EncodableMods.KOINA_ELEMENTS_NO_ISOTOPES)
                    .limits(Limits.NONE.withForbiddenResidues(Limits.NONSTANDARD_RESIDUES).withLength(0, 60))
                    .takesIsobaricLabels(false)
                    .singleThreadFormatting(true)
                    .build()
    );

    private static final Map<String, KoinaModelSpec> BY_LOWERCASE_NAME = indexByName(ALL);

    private KoinaModelRegistry() {}

    /** Every Koina model, RT models first, then MS2 and IM models. */
    public static List<KoinaModelSpec> all() {
        return ALL;
    }

    /** The spec of a Koina model, whatever the case of its name, or null if there is no such model. */
    public static KoinaModelSpec find(String name) {
        if (name == null) {
            return null;
        }
        return BY_LOWERCASE_NAME.get(name.toLowerCase(Locale.ROOT));
    }

    /**
     * The spec of a Koina model, whatever the case of its name, for code that is only reached with a Koina model.
     *
     * @throws IllegalArgumentException if there is no such model
     */
    public static KoinaModelSpec require(String name) {
        KoinaModelSpec spec = find(name);
        if (spec == null) {
            throw new IllegalArgumentException(name + " is not a Koina model");
        }
        return spec;
    }

    //the Prosit models: UNIMOD codes with C written as carbamidomethylated, up to 30 residues and charge 6
    private static Builder prosit(String name, Kind kind) {
        return KoinaModelSpec.builder(name, kind)
                .peptideFormat(PeptideFormat.PROSIT).modMassMap(ModMassMap.LIMITED).modRule(ModRule.PROSIT)
                .limits(Limits.NONE.withForbiddenResidues(Limits.NONSTANDARD_RESIDUES).withLength(0, 30)
                        .withCharge(0, 6));
    }

    //the Prosit-PTM models: any mod Koina's hand-written table or zero-shot encoding takes (KoinaPrositPtm), no
    //C-term mod, C as it is (a bare C is a free cysteine), up to 30 residues and charges 1 to 6. The collision energy
    //is sent as an NCE, which Koina divides by 100
    private static Builder prositPtm(String name, Kind kind) {
        return prosit(name, kind)
                .peptideFormat(PeptideFormat.PROSIT_PTM).modMassMap(ModMassMap.PROSIT_PTM).modRule(ModRule.PROSIT_PTM)
                .encodableMods(EncodableMods.PROSIT_PTM)
                .limits(Limits.NONE.withForbiddenResidues(Limits.NONSTANDARD_RESIDUES).withLength(0, 30)
                        .withCharge(1, 6));
    }

    //the Prosit lactylation models: lactyl on K and the N-terminus (UNIMOD:2114), oxidized M and carbamidomethylated
    //C, no C-term mod. Koina reads a bare C as carbamidomethylated, so a peptide with a free cysteine is not sent.
    //Up to 30 residues and charges 1 to 6; the collision energy is sent as an NCE, which Koina divides by 100
    private static Builder prositLac(String name, Kind kind) {
        return prosit(name, kind)
                .modMassMap(ModMassMap.LAC).modRule(ModRule.LAC).encodableMods(EncodableMods.NO_CTERM)
                .limits(Limits.NONE.withForbiddenResidues(Limits.NONSTANDARD_RESIDUES).withLength(0, 30)
                        .withCharge(1, 6).withSkipFreeCys(true));
    }

    //the ms2pip models: oxidized M and carbamidomethylated C, up to 30 residues, sent no collision energy
    private static Builder ms2pip(String name) {
        return KoinaModelSpec.builder(name, Kind.MS2)
                .inputs(PRECURSOR_CHARGES)
                .peptideFormat(PeptideFormat.MS2PIP).modMassMap(ModMassMap.LIMITED).modRule(ModRule.MS2PIP)
                .limits(Limits.NONE.withForbiddenResidues(Limits.NONSTANDARD_RESIDUES).withLength(0, 30));
    }

    //the AlphaPeptDeep models: any known mod, any residue, batches sorted by length
    private static Builder alphapept(String name, Kind kind) {
        return KoinaModelSpec.builder(name, kind)
                .peptideFormat(PeptideFormat.ALPHAPEPT).modMassMap(ModMassMap.ALL).modRule(ModRule.ANY_KNOWN)
                .instrumentMap(InstrumentMap.ALPHAPEPT)
                .sortBatchByLength(true)
                .singleThreadFormatting(true);
    }

    private static Map<String, KoinaModelSpec> indexByName(List<KoinaModelSpec> specs) {
        Map<String, KoinaModelSpec> byName = new LinkedHashMap<>();
        for (KoinaModelSpec spec : specs) {
            if (byName.put(spec.name().toLowerCase(Locale.ROOT), spec) != null) {
                throw new IllegalStateException("Two Koina model specs are named " + spec.name());
            }
        }
        return Collections.unmodifiableMap(byName);
    }
}
