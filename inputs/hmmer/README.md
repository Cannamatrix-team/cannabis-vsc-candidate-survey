# Archived HMMER search outputs

This directory contains the 42 HMMER 3.4 `--tblout` files used by the 23 target bins in `queries/query_panel.tsv`. They were copied unchanged from CannamatrixAI commit `ca6a37a5240e8e0c85b29912d9136fc087cd7d1d`. The original round and family directory structure is retained.

The complete GMO v1 high-confidence proteome is not duplicated here. Download the exact [`GMO.v1.primary_high_confidence.proteins.fasta.gz`](https://resources.michael.salk.edu/jbrowsers/data/csat/assemblies/GMO/genes_v1/GMO.v1.primary_high_confidence.proteins.fasta.gz) file listed in the Genes V1 Primary table on the [Michael Lab Cannabis pangenome portal](https://resources.michael.salk.edu/resources/cannabis_genomes/index.html), then decompress it. The verified input contains 55,790 FASTA records and has SHA-256 `953fbdcb7f10cc02af10e39a8e4d23108b6854b4ea90d700bf2e24c5735d1da7`. Run:

```bash
make verify-search VSC_PROTEOME=/path/to/GMO.v1.primary_high_confidence.proteins.fasta
```

`workflow/rebuild_search_assignments.py` applies the paper thresholds of bit score at least 50 and E-value at most 1e-5, records whether each protein passed a Pfam-profile search, a full-length-reference search, or both, and reproduces `data/search_assignments_1005.tsv` exactly.

The exact query models are packaged under `inputs/search_queries/`. To rerun all 42 searches with HMMER 3.4 and verify the resulting assignments, create the environment described in the main README and run `make verify-hmmer` with the same `VSC_PROTEOME` argument.
