# Cannabis VSC candidate survey

Data and supporting analysis artifacts for the manuscript:

> A genome-scale, hypothesis-guided survey of protein candidates relevant to sulfur metabolism and volatile chemistry in *Cannabis sativa*

Authors: Emanuel Maminakis, Logan Geffen, Kevelin Barbosa-Xavier, and Suliman Sharif.

## Status

This is a pre-release working repository. No GitHub release or Zenodo DOI has been created yet. Contents may change before the manuscript-supporting release is frozen.

No reuse license has been selected yet.

## Current contents

- `data/candidate_catalog_941.tsv`: the paper-specific catalog of 941 retained GMO v1 proteins across 20 reporting categories.
- `data/candidate_sequences_941.fasta`: full-length amino-acid sequences for those 941 proteins.
- `data/candidate_counts_by_category.tsv`: counts and hypothesis-group crosswalk for the 20 categories.
- `data/expression/`: exact-identity expression-atlas mappings and identity-threshold sensitivity results.
- `phylogeny/pf01053/`: alignment, tree, reference metadata, and nearest-reference assignments for the six PF01053 candidates.

The initial repository deliberately excludes the complete private analysis repository, the full 55,790-protein GMO v1 input proteome, Pfam libraries, large GWAS intermediates, and manuscript audit/session artifacts.

## Candidate catalog construction

The canonical analysis manifest contained 1,046 retained rows across 21 categories. The final manuscript does not report the 105 AAT rows, leaving 941 proteins across its 20 reporting categories. `candidate_catalog_941.tsv` is that exact paper-specific subset joined to the corresponding final-pass validation fields by `gene_id`.

The catalog contains:

- discovery fields: candidate identifier, reporting category, search source, discovery round, sequence scope, matched protein, HMMER score and E-value, and sequence length;
- PF01053 resolution fields: original family pool, resolution group and method, status, nearest reference, and nearest-reference family;
- validation fields: MEME/FIMO status, motif support, expected-domain status, annotation flags, active-gate flags, and final-pass status.

The two MGL-nearest proteins have blank `pfam_status` values because the generated MGL-nearest reporting category had no configured expected-domain rule. Both sequences carry PF01053 and were retained by the gate default, as described in the manuscript.

## Expression tables

The stringent expression cross-reference used a highest-bitscore BLASTp HSP at 100% amino-acid identity across at least 80% of the GMO query. It produced 168 candidate mappings to 128 distinct Cannabis Expression Atlas genes.

`expression_identity_sensitivity.tsv` reports how the mapping counts change at lower identity thresholds. The source atlas is described by Barbosa-Xavier et al. (2024), DOI: [10.1111/ppl.70010](https://doi.org/10.1111/ppl.70010).

## Provenance

The underlying search and annotation artifacts came from CannamatrixAI commit `ca6a37a5240e8e0c85b29912d9136fc087cd7d1d` on branch `KB-VSC_pipeline_validation`. GMO v1 protein models originate from the Cannabis pangenome resource described by Lynch et al. (2025), DOI: [10.1038/s41586-025-09065-0](https://doi.org/10.1038/s41586-025-09065-0).

Checksums for the current files are recorded in `SHA256SUMS`.
