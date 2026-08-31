# Manuscript figures and source tables

This directory freezes the presentation layer used by the reconciled PLOS ONE manuscript. The four figures and four manuscript-only tables are regenerated entirely from the public 941-candidate evidence package; they do not read the private `results_v2` tree.

The comparison authority is `vsc-manuscript` commit `c45f33f`. `figure_checksums.tsv` records the eight PNG/PDF hashes from that frozen manuscript package. `tables/` contains the four exact source tables not already present elsewhere in this repository.

Create the pinned Linux plotting environment and verify the artifacts:

```bash
conda env create -f manuscript/environment.yml
conda activate cannabis-vsc-manuscript
make verify-manuscript
```

Generated files are written under `build/manuscript/`. The verifier requires byte-identical PNG and PDF figures and byte-identical TSV tables.

The generator uses:

- the search assignments, stage ledger, category counts, and motif rules for Figures 1 and 2 and the reporting tables;
- the packaged alignment/tree, assignments, and reference metadata for Figure 3;
- the deposited expression sensitivity, exact mappings, and expression bridge for Figure 4;
- the actual BLASTp/post-analysis parameter boundary documented by the expression runner.

The plotting environment is pinned because Matplotlib, Biopython tree rendering, FreeType, and PDF metadata can change output bytes even when the scientific values are unchanged.
