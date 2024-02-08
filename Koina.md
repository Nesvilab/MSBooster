# Using Koina with MSBooster

## Why use Koina?
[Koina](https://koina.proteomicsdb.org/) (described in [Picciani et al.,_Proteomics_, 2023](https://pubmed.ncbi.nlm.nih.gov/37672792/))
is an online service that allows us to leverage a remote GPU server for predictions
from multiple MS2 and RT models, thereby democratizing deep learning. The following
models are currently supported by MSBooster:

| RT models            | 
|----------------------|
| Deeplc_hela_hf       |
| AlphaPept_rt_generic |
| Prosit_2019_irt      |
| Prosit_2020_irt_TMT  |

| MS2 models                    |
|-------------------------------|
| ms2pip_2021_HCD               |
| AlphaPept_ms2_generic         |
| Prosit_2019_intensity         |
| Prosit_2023_intensity_timsTOF |
| Prosit_2020_intensity_CID     |
| Prosit_2020_intensity_TMT     |
| Prosit_2020_intensity_HCD     |

## Using Koina
Using Koina is currently only supported in the command line interface. You can invoke
it with a command like this
```
java -Xmx14G -jar MSBooster-1.1.27.jar --paramsList msbooster_params.txt
```
An example parameter file for "--paramsList" can be found at
[msbooster_params.txt](msbooster_params.txt). Add the lines below to the
parameter file:
```
useKoina = true
rtModel = _____
spectraModel = _____
```
Substitute the blanks with one of the models listed in the tables above.
You can mix and match, or even use "DIA-NN" for one of them if you want
to call the DIA-NN prediction model specified by the "DiaNN" parameter
in paramsList.

If you want to use an edited.pin from Koina in FragPipe, please load your workflow
in FragPipe. Uncheck all boxes from MSBooster onwards. Run MSBooster via the command
line. Move the original .pin file somewhere else, and rename the edited.pin file
to the original name (ex. move "search.pin" to a different folder, and rename
"search_edited.pin" to "search.pin"). Then, in FragPipe, reload the workflow, but
uncheck the MSFragger and MSBooster checkboxes.

## Figure ideas for manuscript
1. Workflow figure
    1. How we work with Koina
    2. Models we support
    3. Table 1: Types of supported peptides for each (at least on Koina)
    4. Supplementary (timing of analysis of each dataset)
    5. Incorporation into FragPipe (plus screenshot)
2. Number of peptides and proteins identified in each dataset (refer to
   [the datasets section](#datasets-to-test-on))
    1. Table 2: datasets each is trained on. Upset plot to show intersection
    2. 10 iterations of Percolator for each
    3. Show performance for RT and MS2 models separately
    4. Then, pick the best combo and compare it to using same model for RT and MS2
    5. RT calibration for unsupported PTMs (phospho, TMT, etc.)
3. Future directions
    1. Hosting other instances of server
    2. APD for transfer learning

## Datasets to test on
- DDA tryptic (proteome maps)
    - [mouse](https://www.nature.com/articles/s41592-022-01526-y)
        - PDAC
        - phospho and whole proteome
        - "micro-flow LC–MS/MS system using a modified Vanquish pump (Thermo Fisher
          Scientific) coupled to a Q Exactive Orbitrap HF-X (Thermo Fisher Scientific)"
        - /storage/yangkl/fasta/2023-12-19-decoys-reviewed-isoforms-contam-UP000000589.fas
        - 66 PDAC lines
        - 8 phospho files
    - [Arabidopsis](https://www.nature.com/articles/s41586-020-2094-2)
        - 24 fractions
        - "Nanoflow LC-MS/MS was performed by coupling a Dionex 3000 (Thermo Fisher
          Scientific) to a QExactive Orbitrap HF (Thermo Fisher Scientific)"
        - 110 min gradient
        - 12th fraction of 30 tissues
- HLA
    - [DDA Class I and II](https://jitc.bmj.com/content/9/4/e002071)
        - "90,428 HLA-I and 142,625 HLA-II peptides"
        - "Orbitrap Fusion Lumos mass spectrometer (Thermo Fisher Scientific, San Jose,
          California, USA) equipped with a Nanospray Flex Ion Source (Thermo Fisher
          Scientific) coupled to an Ultimate 3000 RSLC Nano UHPLC System (Thermo Fisher
          Scientific)"
        - Use mass shifts for AA substitutions?
        - no C+57
    - [DDA Class I and II timsTOF](https://www.mcponline.org/article/S1535-9476(23)00073-7/fulltext)
    - [DDA Class I stepped HCD](https://www.sciencedirect.com/science/article/pii/S2666166721000927)
    - [DIA](https://www.mcponline.org/article/S1535-9476(21)00053-0/fulltext)
- [timsTOF-Pro machine](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6283298/)
    - Prosit already has a model for comparison
    - 1 tryptic DDA
    - "four replicate injections each of HeLa lysate acquired at five different
      TIMS ramp (accumulation) times"

## Workflow for getting results
- Use most recent FragPipe version
- Run in FragPipe GUI
    - Run MSFragger search first
- Run MSBooster via bash
    - Copy pin and pepxml files into sep folders
    - Rename MSBooster pin file
    - Run RT and MS2 models separately
- Run rest of workflow in FragPipe headless
    - 10 Percolator iterations, changing "--seed"
- Repeat the process, but picking the best combo of individual features
    - Compare against paired models

## Timing

## General recommendations for models to use by type of experiment

## Notes
- DeepLC, APD_rt and ms2 did change predictions with added phospho on T
    - Not ms2pip, prosit
    - need to check DIANN
- ms2pip could not predict longer than 30
  - Prosit is fine, but maybe not for TMT?
- ms2pip has other models on web server, but one on Koina is only HCD
- for TOF, use APD and Prosit
  - for comparison, use APD with Lumos and Prosit HCD
- for TMT, use DIA-NN, Prosit
  - Prosit assumes nterm must be labeled
  - less than 31 long for prosit
- TMT phospho, DIANN
- C57 changes with deeplc+ms2pip, APD
  - not prosit
  - check DIANN

- figure for percent of peptides in testing data that overlap with model 
training?

- supported peptides
  - HCD vs CID (fragmentation)
- table for info needed by each model
- subsection for TMT
  - does Prosit do better for peptides with n-term TMT?
- TMT analysis uses 1 group of 11 samples