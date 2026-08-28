# Workflow provenance

This repository's deposited data release is the authority for the manuscript result set: 1,005 search assignments, 975 unique proteins, and 941 retained candidates in 20 reporting categories. The historical private source is CannamatrixAI commit `ca6a37a5240e8e0c85b29912d9136fc087cd7d1d` on `KB-VSC_pipeline_validation`.

The historical commit is a provenance reference, not yet a frozen end-to-end workflow. It contains older result states, unrelated GWAS and structural-analysis stages, automatic dependency installation, and machine-specific paths. Importing that tree wholesale could silently replace the released 941-candidate result with a stale result set.

## Current reproducible scope

| Stage | Public implementation | What it proves |
|---|---|---|
| Release integrity and reconciliation | `workflow/verify_release.py` | Archived files are unchanged; identifier sets and reported counts agree across the search assignments, stage ledger, catalog, FASTA, expression tables, and PF01053 assignments. |
| HMMER search execution | `workflow/run_hmmer_searches.py`; `inputs/search_queries/` | HMMER 3.4 reruns all 42 Pfam-profile and full-length-reference searches from the packaged query assets and external GMO v1 proteome. Hit identifiers and scores match the archived searches; the comparator permits only documented E-value rounding or underflow across platforms. |
| Candidate extraction from HMMER outputs | `workflow/rebuild_search_assignments.py`; `inputs/hmmer/` | The 42 archived HMMER result tables and external GMO v1 high-confidence proteome reproduce all 1,005 deposited search assignments exactly. |
| Shared PF01053 resolution | `workflow/run_pf01053_resolution.py`; `workflow/compare_pf01053_outputs.py`; `phylogeny/pf01053/` | A deterministic fresh alignment and tree reproduce all six CBL, CGS, and MGL nearest-reference assignments. A second inference from the deposited trimmed alignment reproduces the archived unrooted topology, LG+R2 model, and log-likelihood. |
| Publication reporting counts | `workflow/generate_reporting_tables.py` | Candidate counts are recalculated from the catalog and expression mappings are recalculated from the exact-identity bridge; both outputs match the deposited tables. |

Together, these checks reconstruct the 1,005 initial search assignments from the full GMO v1 proteome and rerun the shared PF01053 family resolution. They do not yet rerun the downstream motif screen, expected-domain validation, final gate, or expression mapping.

The historical MAFFT stage used multithreaded L-INS-i iterative refinement. Repeated runs can differ in alignment columns and reference-only branches because thread scheduling changes the refinement path. The public full rerun fixes MAFFT to one thread and verifies the scientific invariant used downstream: the same nearest reference and family for every candidate. The archived trimmed alignment remains deposited so the IQ-TREE inference itself can be checked independently against the archived topology, model, and likelihood.

## Upstream stages to import and reconcile

Paths below are relative to the historical CannamatrixAI commit.

| Stage | Historical source | Required reconciliation before public use |
|---|---|---|
| MEME/FIMO motif screen | `synthase_features/pipeline/2_meme_validation/`; `synthase_features/pipeline/phase6_config.yaml` | Separate source from the bundled tool installation, freeze motif parameters and versions, and reproduce the 975 motif records including the GSH2 singleton disposition. |
| Expected-domain validation | `synthase_features/pipeline/2b_domain_validation/`; `synthase_features/pipeline/shared/domain_validation.py` | Freeze the Pfam models and HMMER version, remove automatic installs, and reproduce 973 domain records plus the two explicitly absent MGL-nearest records. |
| Final validation gate | `synthase_features/pipeline/2e_validation_gate/`; `synthase_features/pipeline/shared/validation_gate.py` | Limit the gate to manuscript-active rules and reproduce the 34 removals and 941 retained proteins exactly. |
| Expression-atlas mapping | `synthase_features/pipeline/2_plots/vsc_expression_overview.py`; `synthase_features/pipeline/results_v2/2h_expression_overview/` | Replace hard-coded paths with inputs, capture the BLAST database build and command, freeze BLAST version/options, and reproduce the 168 mappings to 128 atlas genes. |
| Manuscript figures and remaining source tables | Maintained in the separate `vsc-manuscript` repository | Port only figures supported by public inputs, declare plotting dependencies, and compare regenerated files with the submitted figures. |

The GWAS, structure-prediction, Foldseek, structural-QC, geometry, and feature-export stages are outside this manuscript's analysis and should not be added to this workflow.

## External inputs

The complete GMO v1 proteome is not duplicated in this repository. Its source dataset is identified in the main README, and its expected SHA-256 checksum is recorded in the HMMER input notes. The exact Pfam models and UniProt reference queries needed for this search panel are packaged under `inputs/search_queries/` with their original terms documented in `THIRD_PARTY_NOTICE.md`.
