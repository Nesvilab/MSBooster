# Using Koina with MSBooster

## Why use Koina?
[Koina](https://koina.proteomicsdb.org/) (described in [Picciani et al.,
_Proteomics_, 2023](https://pubmed.ncbi.nlm.nih.gov/37672792/)) 
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
3. Future directions
   1. Hosting other instances of server
   2. APD for transfer learning

## Datasets to test on
- DDA tryptic (proteome maps)
  - [mouse](https://www.nature.com/articles/s41592-022-01526-y)
    - 32 fractions
    - 90 min gradient
    - "micro-flow LC–MS/MS system using a modified Vanquish pump (Thermo Fisher 
    Scientific) coupled to a Q Exactive Orbitrap HF-X (Thermo Fisher Scientific)"
  - [Arabidopsis](https://www.nature.com/articles/s41586-020-2094-2)
    - 24 fractions
    - "Nanoflow LC-MS/MS was performed by coupling a Dionex 3000 (Thermo Fisher 
    Scientific) to a QExactive Orbitrap HF (Thermo Fisher Scientific)"
    - 110 min gradient
- HLA 
  - [DDA Class I and II](https://jitc.bmj.com/content/9/4/e002071)
    - "90,428 HLA-I and 142,625 HLA-II peptides"
    - "Orbitrap Fusion Lumos mass spectrometer (Thermo Fisher Scientific, San Jose, 
    California, USA) equipped with a Nanospray Flex Ion Source (Thermo Fisher 
    Scientific) coupled to an Ultimate 3000 RSLC Nano UHPLC System (Thermo Fisher 
    Scientific)"
  - [DDA Class I and II timsTOF](https://www.mcponline.org/article/S1535-9476(23)00073-7/fulltext)
  - [DIA](https://www.mcponline.org/article/S1535-9476(21)00053-0/fulltext)
- [timsTOF-Pro machine](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6283298/)
  - Prosit already has a model for comparison
  - 1 tryptic DDA
  - "four replicate injections each of HeLa lysate acquired at five different 
  TIMS ramp (accumulation) times"

## Workflow for getting results

## Timing

## General recommendations for models to use by type of experiment
