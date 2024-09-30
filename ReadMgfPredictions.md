# Reading in predictions from any model via MGF files

If you wish to use predictions produced outside of MSBooster for PSM rescoring (e.g. from your own prediction model), 
you can upload them as an mgf file. It should be formatted as below:
```
BEGIN IONS
TITLE=FFAPWC[57.0215]GHC[57.0215]K
CHARGE=2
RT=30.12642
1/K0=0.97269547
295.14407	0.4740689
307.1435	0.11716866
366.18118	0.16534789
424.17056	0.04255961
444.2024	0.117110334
463.23395	0.02398945
472.69696	0.72797287
501.22385	0.2157223
508.2155	0.3965
581.7497	0.060824234
582.23114	0.019092904
649.3133	0.020508597
661.2546	0.46449134
809.34393	0.0181899
847.33386	0.33339772
866.3654	0.017972104
944.3866	1.0
1003.4243	0.01978566
1015.4237	0.3828275
1162.4922	0.032336686
1163.4551	0.013988824
END IONS
```

The mgf file can include any of RT, 1/K0, or MS/MS information. The number of digits in the mass shift of the peptide
should match what is in the MSFragger pin file. Furthermore, there should be an mgf entry for each unique precursor in
the pin file. This means that if your model does not support a mass shift included in the MSFragger search, you must
 add an mgf entry with those shifted MS/MS fragments. If you are unable to create such a proxy prediction (e.g. pin file
includes a peptide with amino acid U that your model cannot predict, you can set RT and/or 1/K0 to 0 and exclude any 
fragments from the mgf entry).

Provide these files with the parameters <code>spectraPredFile</code>, <code>RTPredFile</code>, and 
<code>IMPredFile</code>. Even if all your predictions are in the same file, this parameter must be included for each
predicted peptide property you wish to use. You do not need to specify the name of the model you are using in
<code>spectraModel</code> etc.