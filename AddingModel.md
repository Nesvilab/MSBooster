1. Add it to the Koina models hashsets (NCE included) and model search strings
2. Pass inputs to it in JSONWriter
3. Add PTM mass, unimod code, AAunimodToModMass, and AAmods in PTMhandler

In some cases if it's a very new model (i.e. not another Prosit or ms2pip model, etc) it will require more changes
4. Update PeptideSkipper with restrictions
5. PeptideFormatter changes to write out correct sequence string for new model