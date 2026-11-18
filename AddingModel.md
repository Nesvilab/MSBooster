# Adding a Koina model

Every Koina model MSBooster can call is one `KoinaModelSpec` in `koinaclasses.KoinaModelRegistry`. The rest of
MSBooster reads the spec: the model lists (`ModelCollections.KoinaRTmodels`, `KoinaMS2models`, `KoinaIMmodels`,
`KoinaModels`, `KoinaIsoLabelmodels`, `allowedFragmentationTypes`, `NceConstants.nceModels`, `LowercaseModelMapper`),
the request `JSONWriter` writes, how the response is parsed, how peptides are written for the model and read back,
which PSMs and mods it is sent (`KoinaPreCheck`), the start-up check (`KoinaModelCheck`), the NCE calibration range
and the citation. Adding a model takes the steps below. A new model goes into none of the lists of the best model
search (`ModelCollections.generalRTmodels` and the others) unless it was evaluated for them; it can still be searched
by naming it in `rtSearchModelsString`, `ms2SearchModelsString` or `imSearchModelsString`.

## 1. The spec
Add one entry to `KoinaModelRegistry.ALL`, among the models of its kind (RT, MS2 or IM). Start from the builder of
its family when it has one (`prosit`, `prositPtm`, `prositLac`, `ms2pip`, `alphapept`), or from
`KoinaModelSpec.builder(name, kind)`, whose defaults are listed in `KoinaModelSpec.Builder`. Set what the model's
page on https://koina.wilhelmlab.org and its preprocessing in https://github.com/wilhelm-lab/koina (`models/<Model>/`)
say it takes:

| Builder method | What it sets |
|---|---|
| `inputs(...)` | the inputs sent besides `peptide_sequences`: `PRECURSOR_CHARGES`, `COLLISION_ENERGIES`, `INSTRUMENT_TYPES`, `FRAGMENTATION_TYPES`, sent in that order |
| `outputs(...)` | the output names read: the prediction for an RT or IM model (default `irt` or `ccs`); intensities, m/z and annotation for an MS2 model, in that order whatever order Koina answers them in |
| `imConversion(...)` | IM models: the unit Koina answers in, `ImConversion.CCS_N2` (converted to 1/K0) or `ImConversion.INVERSE_K0` |
| `peptideFormat(...)` | how terminal mods are written, whether a bare C is written as `C[UNIMOD:4]`, and the other format rules (`KoinaModelSpec.PeptideFormat`) |
| `annotationParser(...)` | MS2 models: how a fragment annotation is read (default `ION_NUMBER_CHARGE`, as in `y3+1`). An annotation the parser does not know fails the run rather than being read wrong |
| `trustKoinaMz(true)` | MS2 models: use the m/z Koina answers with, for fragments MSBooster does not calculate (UniSpec, Altimeter) |
| `modMassMap(...)` | the table that turns the model's UniMod codes into masses and back |
| `modRule(...)` | the mods the model's peptide format writes; a mass matching none of them is left out of the peptide sent |
| `encodableMods(...)` | which of those mods Koina can encode for the model; the others are left out too and counted in the log |
| `limits(...)` | the peptides and charges the model takes: residues, length, charge, length as written, free cysteines |
| `batchSize(...)` | the most precursors in one request (default 1000); the server's `max_batch_size` lowers it at run time |
| `calibrateCe(true)`, `ceRange(min, max)` | whether the collision energy is calibrated, and over which range (default 20-40); `minNCE`/`maxNCE` override the range only when set |
| `instrumentMap(...)` | how the run's instrument is written for the model |
| `fragmentations(...)`, `fragmentation(type, sentAs)` | the fragmentation types the model takes, and the exact string sent for each; with none, any type is sent as it is |
| `isobaricLabel(true)` | the model is trained on isobaric-labelled peptides and labels every peptide (the TMT models) |
| `takesIsobaricLabels(false)` | the model cannot be sent TMT or iTRAQ, so a search with one cannot use it (IM2Deep) |
| `sortBatchByLength(true)` | the precursors are sorted by length before they are split into requests (AlphaPeptDeep) |
| `singleThreadFormatting(true)` | the peptides are formatted on one thread; required by a mod rule that remembers the codes it found (`ModRule.ANY_KNOWN`) |
| `citationKey(...)` | the model's path in `koina_citations.yml`, without the slashes and `/infer` (default: the model name) |

The `KoinaModelSpec` constructor refuses a spec that contradicts itself, such as an annotation parser on a model that
is not an MS2 model, or a collision energy calibrated but not sent.

When none of the existing strategies does what the model needs, add a constant to the enum of that strategy
(`PeptideFormat`, `AnnotationParser`, `ModMassMap`, `ModRule`, `EncodableMods`, `InstrumentMap`) with what it does: a
parser in `KoinaOutputParsingMethods`, a mass table in `PTMhandler`, and so on. Building a spec must read no MSBooster
setting: FragPipe loads `KoinaModelRegistry` from the MSBooster jar on its own, without `MainClass` having run. A
strategy that needs a `PTMhandler` table reaches it only when it is called.

`KoinaModelRegistryTest` pins the registry: add the model to its model sets and to the expectations of each list,
input, batch size and citation key it checks.

## 2. The golden fixture
`koinaclasses.KoinaGoldenTest` fails for a model of the registry without a complete folder under
`src/test/resources/koina/`. Create `src/test/resources/koina/<Model>/` with:
- `case.properties`: the run settings, `instrument`, `fragmentation` and `nce`
- `psms.tsv`: synthetic PSMs only, never data of a user. One pin-style peptide (MSFragger mass notation, flanking
residues included), charge and note per line. Cover the model's rules: a mod it is sent, one it leaves out, a peptide
it does not take

Then record the model's answer from the live Koina server:

    mvn -B test -Dtest=KoinaGoldenTest -Dkoina.golden.record=<Model>

This sends the request to Koina (`-Dkoina.golden.url=...` picks another server) and writes `request.json`,
`response.json`, `metadata.json`, `config.json`, `parsed.tsv` and `transferred.tsv`. Review them before committing.
The folder layout and the update mode are described in `src/test/resources/koina/README.md`. Rules of the model that
the goldens do not show, such as the reason a peptide is not sent, get unit tests of their own, as in
`KoinaAltimeterTest`, `KoinaChronologerTest`, `KoinaPrositLacTest` and `KoinaPrositPtmTest`.

## 3. The citation
`src/main/resources/koina_citations.yml` is Koina's https://koina.wilhelmlab.org/openapi.yml. MSBooster prints the
`### Citation` section of the model's path, `/<citationKey>/infer`, and writes it into the parameter file
(`citations.KoinaYamlParser`). If the path is missing, download the file again and change the date in its header.
If Koina gives no citation to print, write one into the section and list the edit in the header comment of the file,
as is done for the models there. `KoinaYamlParserTest.everyKoinaModelHasACitation` fails for a model whose key has
no citation.

## 4. The documentation
Add the model to the table in [Koina.md](Koina.md), with what it is sent, its mods and the peptides it takes.

## 5. FragPipe
FragPipe offers the Koina models in the MSBooster panel of the Validation tab. In
`FragPipe-GUI/src/main/java/org/nesvilab/fragpipe/`:
- `tools/msbooster/MSBoosterPanel.java`: append the display name and the MSBooster name to `RT_MODELS`,
`SPECTRA_MODELS` or `IM_MODELS`, and add the model's text to `KOINA_MODEL_MODS`: the mods it is sent and the
peptides it takes, written by hand from the spec. The text is shown as an HTML tooltip, so it holds no `<`, `>` or `&`
- `ui.translations`: append `msbooster.<rt|spectra|im>-model::<display name> <=> <MSBooster name>`, in the same
order

`KoinaModelModsTest` checks that every model offered has a `KOINA_MODEL_MODS` text, `MSBoosterTranslationsTest` that
every model offered has a translation, and `MSBoosterJarRegistryTest` that `KoinaModelRegistry.find` in the MSBooster
jar of `tools/` returns a spec of the matching kind for every model offered. Until an MSBooster with the model is
pinned in `components.json`, check that last one against the jar built from this repository:

    gradlew test -Dmsbooster.jar=<path to MSBooster jar> --tests *MSBoosterJarRegistryTest
