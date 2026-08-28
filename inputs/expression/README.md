# Expression-mapping inputs

The public runner requires the Cannabis Expression Atlas v1.1 protein FASTA as an external input. The full Atlas protein database is not redistributed here.

## Atlas protein FASTA

The validated source is the official Docker image `barbosaxavierkevelin/cannabisexpressionatlas:v1.1`, pinned to manifest digest:

```text
sha256:4df81a1ae4a755278c8c22f40877c7859ecd03f2dd445dde3ecead7405e4eaef
```

The required file is `/build_zone/blast_db/CEA_protein_sequences_all.faa` inside that image. It contains 27,893 proteins and must have SHA-256:

```text
5b998b6bd3a3c59a1538a7cd3dc4e9f006e41b99855e9a5c352c52524a134e4e
```

With Docker installed, extract it without running the application:

```bash
image=barbosaxavierkevelin/cannabisexpressionatlas@sha256:4df81a1ae4a755278c8c22f40877c7859ecd03f2dd445dde3ecead7405e4eaef
container=$(docker create "$image")
docker cp "$container":/build_zone/blast_db/CEA_protein_sequences_all.faa ./CEA_protein_sequences_all.faa
docker rm "$container"
sha256sum CEA_protein_sequences_all.faa
```

The application layer containing the FASTA is independently addressable as `sha256:65f6472e8cecab9ef7360fe00db80aca177f47008ecd667281128b230991a400`.

The runner pins the Linux Bioconda build `blast=2.5.0=hc0b0e79_3`. Its 1,062 manuscript-scope HSP rows, sorted canonically to remove output-order differences, have SHA-256 `54c9dd60577c546e74ea5f0e00ac86c86639249a018ce1a436f955a828752acd`.

## Minimal gene metadata

`atlas_gene_metadata_423.tsv` contains only the classification, tissue label, and Tau fields for the 423 Atlas genes reached at the least stringent deposited sensitivity threshold. It was imported from `VSC_ATLAS_gene_metadata.csv` at CannamatrixAI commit `ca6a37a5240e8e0c85b29912d9136fc087cd7d1d`; the complete 482-row historical source has SHA-256 `1c04882e3f518d781121e2b384edc39e10a29a42560be5b38209acea92b735c1`.

The extraction is recorded in `workflow/import_expression_metadata.py`. Missing values were normalized from `NA` to empty fields; no classification, tissue, or Tau values were otherwise changed.

## Historical command semantics

The Atlas application invoked protein BLAST as:

```text
blastp -query QUERY -db prot_atlas_blast_db -evalue 0.001 -outfmt 6 -matrix BLOSUM62 -qcov_hsp_perc 80 -max_target_seqs 1
```

Although the saved UI options also named `perc_identity: 100`, the Atlas v1.1 protein-search branch did not pass that option to `blastp`. The manuscript package therefore applies 100% identity after selecting the highest-bitscore HSP per query. The same downstream operation produces the deposited threshold-sensitivity series.
