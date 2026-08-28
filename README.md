# Cannabis VSC candidate survey

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.21685901.svg)](https://doi.org/10.5281/zenodo.21685901)

Supporting data for the manuscript:

> Computational Identification of Candidate Gene Families for Volatile Sulfur Compound Biosynthesis in *Cannabis sativa* Using Profile Hidden Markov Models

Authors: Emanuel Maminakis, Logan Geffen, Kevelin Barbosa-Xavier, and Suliman Sharif.

## Status

Version 0.1.0 is the first manuscript-supporting release and is archived on Zenodo at [10.5281/zenodo.21685901](https://doi.org/10.5281/zenodo.21685901).

## Evidence package

The files preserve the path from the stated query panel to the paper-specific catalog:

1. `queries/query_panel.tsv` documents the 23 target bins, Pfam profiles, full-length reference queries, biochemical contexts, admission rule, and any post-search resolution step.
2. `data/search_assignments_1005.tsv` contains the 1,005 initial target-bin assignments, representing 975 unique GMO v1 proteins. A protein may occur in more than one target bin. The support columns distinguish Pfam-only, reference-only, and dual-route hits.
3. `data/candidate_stage_ledger_975.tsv` provides one row per unique protein after family resolution and carries the motif, expected-domain, gate, and final outcome fields through the validation stages.
4. `data/candidate_catalog_941.tsv` is the final paper-specific catalog of 941 retained proteins across 20 reporting categories.
5. `data/candidate_sequences_941.fasta` contains the corresponding full-length amino-acid sequences.
6. `data/candidate_counts_by_category.tsv` gives category counts and the manuscript hypothesis-group crosswalk.

The stage ledger records 34 removals for incomplete expected-domain architecture: 22 ADH, 5 alliinase, 2 FMO, and 5 LOX candidates. Two MGL-nearest proteins had no configured expected-domain rule and were retained by the gate default; this is explicitly marked in `final_outcome`.

## Search and validation details

Searches were run with HMMER 3.4. Candidate extraction retained hits with a bit score of at least 50 and an E-value of at most 1e-5 from any Pfam-profile or full-length-reference search listed for that target bin. `search_assignments_1005.tsv` reconstructs the passing search route from the original HMMER tabular outputs at those thresholds.

Overlapping PF01053 search bins were resolved with the joint maximum-likelihood phylogeny in `phylogeny/pf01053/`. The stage ledger then joins the resolved assignments to the MEME/FIMO motif results, expected-domain validation, and final gate by `gene_id`. It contains 975 motif records, including one explicitly flagged singleton that was not motif-evaluable, and 973 expected-domain records; the two absent records are the MGL-nearest cases described above.

## Expression tables

`data/expression/` contains exact-identity expression-atlas mappings and identity-threshold sensitivity results. The stringent cross-reference used the highest-bitscore BLASTp HSP at 100% amino-acid identity across at least 80% of the GMO query, producing 168 candidate mappings to 128 distinct Cannabis Expression Atlas genes.

The source atlas is described by Barbosa-Xavier et al. (2024), DOI: [10.1111/ppl.70010](https://doi.org/10.1111/ppl.70010).

## Provenance and upstream data

The underlying search and annotation artifacts came from CannamatrixAI commit `ca6a37a5240e8e0c85b29912d9136fc087cd7d1d` on branch `KB-VSC_pipeline_validation`.

GMO v1 protein models originate from the Cannabis pangenome resource described by Lynch et al. (2025), DOI: [10.1038/s41586-025-09065-0](https://doi.org/10.1038/s41586-025-09065-0). The source annotation dataset is available from Figshare at DOI: [10.25452/figshare.plus.25909024.v1](https://doi.org/10.25452/figshare.plus.25909024.v1). The complete 55,790-protein input proteome is not duplicated here; this package includes only the 941 paper-specific sequences.

The 24 exact Pfam profiles and 13 UniProt reference sequences used by the search panel are preserved under `inputs/search_queries/`. These are only the search-specific assets, not complete copies of either database. See `THIRD_PARTY_NOTICE.md` for source-specific attribution and licensing.

The complete private analysis repository, large GWAS intermediates, and manuscript audit/session artifacts are not part of this evidence package.

## Reproducibility workflow

The public workflow starts from the deposited tables and provides a fast check of the released result set:

```bash
make verify-release
```

This verifies the archived checksums; reconciles the search, stage-ledger, catalog, FASTA, expression, and PF01053 identifiers and counts; and regenerates the two deposited reporting-count tables for comparison. To write those regenerated tables under `build/reporting/`, run:

```bash
make reporting
```

Only Python 3 and GNU `sha256sum` are required for this release-level workflow. [Workflow provenance and remaining gaps](workflow/PROVENANCE.md) distinguishes these checks from the upstream analysis stages.

The first upstream reconstruction rebuilds all 1,005 search assignments from the archived HMMER tabular outputs. The complete GMO v1 high-confidence protein FASTA remains an external input:

```bash
make verify-search VSC_PROTEOME=/path/to/GMO.v1.primary_high_confidence.proteins.fasta
```

This command regenerates `build/search_assignments_1005.tsv` and requires an exact match with the deposited table. See [the HMMER input notes](inputs/hmmer/README.md) for scope and provenance.

The next upstream check reruns all 42 searches from the packaged query assets. Create and activate the HMMER 3.4 environment, then provide the same external proteome:

```bash
conda env create -f environment.yml
conda activate cannabis-vsc-survey
make verify-hmmer VSC_PROTEOME=/path/to/GMO.v1.primary_high_confidence.proteins.fasta
```

This verifies every regenerated HMMER hit identifier and score against the archived tables, then confirms that they produce the same 1,005 search assignments. HMMER 3.4 on Linux can round an E-value slightly differently from the archived macOS run or report values below `1e-300` as zero. The comparison permits only those numeric representation differences; candidate membership and all other assignment fields must match.

The PF01053 family-resolution workflow is self-contained:

```bash
make verify-pf01053
```

This runs a deterministic MAFFT 7.520, trimAl 1.5.0, and IQ-TREE 3.1.1 reconstruction from the packaged 41 references and six candidates. It requires the same CBL, CGS, and MGL nearest-reference assignments as the deposit. It also reruns IQ-TREE from the deposited trimmed alignment with the archived seed and requires the same unrooted topology, selected model, and log-likelihood.

The historical MAFFT command used multithreaded iterative refinement, whose scheduling can change alignment columns and reference-only branches across reruns. The public reconstruction therefore uses one MAFFT thread and treats the six deposited nearest-reference assignments as the full-run invariant. The deposited trimmed alignment provides the separate reproducibility boundary for the archived tree.

## Citation and Zenodo metadata

`CITATION.cff` supplies GitHub-readable citation metadata and the versioned Zenodo DOI. `.zenodo.json` supplies the metadata Zenodo uses when tagged GitHub releases are archived.

## License and integrity

Original tables, documentation, and selection/arrangement in this repository are licensed under [CC BY 4.0](LICENSE). Third-party source materials retain their original terms as described in `THIRD_PARTY_NOTICE.md`.

`SHA256SUMS` records SHA-256 checksums for every distributed file.
