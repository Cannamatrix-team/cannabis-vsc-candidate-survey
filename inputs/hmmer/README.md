# Archived HMMER search outputs

This directory contains the 42 HMMER 3.4 `--tblout` files used by the 23 target bins in `queries/query_panel.tsv`. They were copied unchanged from CannamatrixAI commit `ca6a37a5240e8e0c85b29912d9136fc087cd7d1d`. The original round and family directory structure is retained.

The complete GMO v1 high-confidence proteome is not duplicated here. The verified input contains 55,790 FASTA records and has SHA-256 `953fbdcb7f10cc02af10e39a8e4d23108b6854b4ea90d700bf2e24c5735d1da7`. Obtain it from the GMO v1 source dataset identified in the main README, stage it as an uncompressed FASTA, and run:

```bash
make verify-search VSC_PROTEOME=/path/to/GMO.v1.primary_high_confidence.proteins.fasta
```

`workflow/rebuild_search_assignments.py` applies the paper thresholds of bit score at least 50 and E-value at most 1e-5, records whether each protein passed a Pfam-profile search, a full-length-reference search, or both, and reproduces `data/search_assignments_1005.tsv` exactly.

These files support candidate extraction from completed HMMER searches. They do not themselves rerun `hmmsearch`; packaging the query models and a frozen HMMER execution command is a separate upstream step.
