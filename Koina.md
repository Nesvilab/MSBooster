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

Using Koina is currently only supported in the command line interface. You can invoke
it with a command like this
```
java -Xmx14G -jar MSBooster-1.1.27.jar --paramsList msbooster_params.txt
```
An example parameter file for --paramsList