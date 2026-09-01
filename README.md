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
5. `data/candidate_sequences_975.fasta` contains the exact full-length input to the manuscript-scope motif screen; `data/candidate_sequences_941.fasta` contains the final retained subset.
6. `motifs/` freezes the 19 evaluable family-specific MEME models and all family classification rules. The remaining GSH2 family is the explicitly flagged singleton.
7. `domains/family_rules.tsv` freezes the expected Pfam architecture for every reporting category and reuses the 24 exact profiles under `inputs/search_queries/pfam/`.
8. `gate/` documents the two manuscript-active final gates and the inactive historical annotations excluded from the public runner.
9. `data/candidate_counts_by_category.tsv` gives category counts and the manuscript hypothesis-group crosswalk.
10. `inputs/expression/atlas_gene_metadata_423.tsv` is the minimal Atlas annotation subset needed to regenerate the exact-identity bridge and all six identity-sensitivity rows from the external Atlas protein FASTA.
11. `manuscript/` freezes the four submitted figure checksums and four manuscript-only source tables, plus a public-input-only generator for reproducing them.

The stage ledger records 34 removals for incomplete expected-domain architecture: 22 ADH, 5 alliinase, 2 FMO, and 5 LOX candidates. Two MGL-nearest proteins had no configured expected-domain rule and were retained by the gate default; this is explicitly marked in `final_outcome`.

## Search and validation details

Searches were run with HMMER 3.4. Candidate extraction retained hits with a bit score of at least 50 and an E-value of at most 1e-5 from any Pfam-profile or full-length-reference search listed for that target bin. `search_assignments_1005.tsv` reconstructs the passing search route from the original HMMER tabular outputs at those thresholds.

Overlapping PF01053 search bins were resolved with the joint maximum-likelihood phylogeny in `phylogeny/pf01053/`. The stage ledger then joins the resolved assignments to the MEME/FIMO motif results, expected-domain validation, and final gate by `gene_id`. It contains 975 motif records, including one explicitly flagged singleton that was not motif-evaluable, and 973 expected-domain records; the two absent records are the MGL-nearest cases described above.

## Expression tables

`data/expression/` contains exact-identity expression-atlas mappings and identity-threshold sensitivity results. The stringent cross-reference used the highest-bitscore BLASTp HSP at 100% amino-acid identity across at least 80% of the GMO query, producing 168 candidate mappings to 128 distinct Cannabis Expression Atlas genes.

The source atlas is described by Barbosa-Xavier et al. (2024), DOI: [10.1111/ppl.70010](https://doi.org/10.1111/ppl.70010).

## Provenance and upstream data

The underlying search and annotation artifacts came from CannamatrixAI commit `ca6a37a5240e8e0c85b29912d9136fc087cd7d1d` on branch `KB-VSC_pipeline_validation`.

GMO v1 protein models originate from the Cannabis pangenome resource described by Lynch et al. (2025), DOI: [10.1038/s41586-025-09065-0](https://doi.org/10.1038/s41586-025-09065-0). The authors' [Michael Lab Cannabis pangenome portal](https://resources.michael.salk.edu/resources/cannabis_genomes/index.html) provides the genome and annotation files, including the exact [GMO v1 primary high-confidence protein FASTA](https://resources.michael.salk.edu/jbrowsers/data/csat/assemblies/GMO/genes_v1/GMO.v1.primary_high_confidence.proteins.fasta.gz). The complete 55,790-protein input proteome is not duplicated here; this package includes the 975 motif-stage sequences and the 941 final retained sequences.

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

The motif-stage workflow is self-contained and pins MEME Suite 5.5.9. The fast verification rescans the exact archived de novo models with FIMO and requires all 975 deposited motif records:

```bash
make verify-motifs
```

The exhaustive verification also reruns de novo MEME discovery for all 19 evaluable families, requires every motif identity, E-value, and probability matrix to match the archived models, then requires the same 975 records after FIMO classification:

```bash
make verify-motif-discovery
```

This full audit took about 46 minutes with four MEME threads during package validation. MEME can emit harmless warnings when optional EPS logos cannot be converted to PNG; motif XML, FIMO output, and classification are unaffected.

Expected-domain validation reruns HMMER 3.4 against the same 24 packaged Pfam profiles and requires exact agreement for all 973 deposited records:

```bash
make verify-domains
```

The other two proteins are the explicitly recorded MGL no-rule cases. The check requires the deposited outcome distribution of 939 complete, 32 partial, and 2 absent domain architectures.

The final gate can then be rerun from the reconstructed motif and domain outputs:

```bash
make verify-gate
```

This requires all 975 gate dispositions, the same 34 removals, and a byte-for-byte match to the deposited 941-protein FASTA. The historical stability and contamination switches were off, and selection context was annotation-only, so they are not part of the manuscript-active runner.

Expression mapping requires the external Cannabis Expression Atlas v1.1 protein FASTA. Extract and checksum the pinned file as described in [the expression input notes](inputs/expression/README.md), then run:

```bash
make verify-expression ATLAS_PROTEINS=/path/to/CEA_protein_sequences_all.faa
```

This rebuilds the protein database with BLAST 2.5.0, runs the historical BLASTp command against the 941 retained proteins, selects the highest-bitscore HSP per query, and requires exact agreement with the deposited 168 mappings, 128 Atlas genes, and six identity-sensitivity rows. Numeric comparison permits only equivalent floating-point renderings of BLAST e-values.

The four manuscript figures and four remaining presentation tables can be regenerated without either external proteome:

```bash
conda env create -f manuscript/environment.yml
conda activate cannabis-vsc-manuscript
make verify-manuscript
```

This requires byte-identical agreement with the frozen reconciled-manuscript PNG/PDF figures and TSV tables. See [the manuscript artifact notes](manuscript/README.md) for the input mapping and plotting-environment boundary.

## Citation and Zenodo metadata

`CITATION.cff` supplies GitHub-readable citation metadata and the versioned Zenodo DOI. `.zenodo.json` supplies the metadata Zenodo uses when tagged GitHub releases are archived.

## License and integrity

Original tables, documentation, and selection/arrangement in this repository are licensed under [CC BY 4.0](LICENSE). Third-party source materials retain their original terms as described in `THIRD_PARTY_NOTICE.md`.

`SHA256SUMS` records SHA-256 checksums for every distributed file.
