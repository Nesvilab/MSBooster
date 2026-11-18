# Koina golden tests

One folder per Koina model. `koinaclasses.KoinaGoldenTest` runs each folder's synthetic PSMs through
today's request writer, serves the folder's recorded Koina response from a local stub server, and
compares the request and the parsed predictions with the files here. The test never calls Koina.

| File | What it is |
|---|---|
| `case.properties` | Run settings: instrument, fragmentation, nce |
| `psms.tsv` | Synthetic PSMs: pin-style peptide in MSFragger mass notation, charge, and a note on what each one tests |
| `request.json` | The request JSONWriter writes for these PSMs |
| `response.json` | Koina's answer to that request, byte for byte |
| `metadata.json`, `config.json` | Koina's answers to `GET <model>` and `GET <model>/config` at recording time. The stub serves them to the start-up check (`koinaclasses.KoinaModelCheck`), which needs every output the spec reads in the metadata and cuts the request at the config's `max_batch_size` when it is below the spec's batch size |
| `parsed.tsv` | Predictions as KoinaModelCaller parses them, keyed by the model's own peptide |
| `transferred.tsv` | The same predictions after `transferKoinaPreds` moved them onto the PSMs' peptides |

The goldens change only on purpose (see `testutils.KoinaGolden`):

```
# rewrite request.json, parsed.tsv and transferred.tsv from the current code and the stored response
mvn -B test -Dtest=KoinaGoldenTest -Dkoina.golden.update=<Model>[,<Model>...]|all

# also send the request to the live Koina server and store its response, metadata and config
# (a new model folder, or a request that changed); -Dkoina.golden.url=... picks another server
mvn -B test -Dtest=KoinaGoldenTest -Dkoina.golden.record=<Model>[,<Model>...]|all
```

Review the diff of a folder before committing it. Only synthetic peptides go into `psms.tsv`.
