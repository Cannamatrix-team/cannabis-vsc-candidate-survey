# Workflow provenance

This repository's deposited data release is the authority for the manuscript result set: 1,005 search assignments, 975 unique proteins, and 941 retained candidates in 20 reporting categories. The historical private source is CannamatrixAI commit `ca6a37a5240e8e0c85b29912d9136fc087cd7d1d` on `KB-VSC_pipeline_validation`.

The historical commit is a provenance reference, not yet a frozen end-to-end workflow. It contains older result states, unrelated GWAS and structural-analysis stages, automatic dependency installation, and machine-specific paths. Importing that tree wholesale could silently replace the released 941-candidate result with a stale result set.

## Current reproducible scope

| Stage | Public implementation | What it proves |
|---|---|---|
| Release integrity and reconciliation | `workflow/verify_release.py` | Archived files are unchanged; identifier sets and reported counts agree across the search assignments, stage ledger, catalog, FASTA, expression tables, and PF01053 assignments. |
| Publication reporting counts | `workflow/generate_reporting_tables.py` | Candidate counts are recalculated from the catalog and expression mappings are recalculated from the exact-identity bridge; both outputs match the deposited tables. |

These are release-level checks. They do not reconstruct the candidates from the full GMO v1 proteome.

## Upstream stages to import and reconcile

Paths below are relative to the historical CannamatrixAI commit.

| Stage | Historical source | Required reconciliation before public use |
|---|---|---|
| HMMER search and extraction | `synthase_features/hmmer_results/VSC/scripts/`; `synthase_features/hmmer_results/VSC/round2_results/`; `synthase_features/pipeline/1_hmmer_search/run_phase.sh` | Freeze both search rounds, the extractor, HMMER version, 55,790-protein GMO v1 input, query files, and thresholds; reproduce `search_assignments_1005.tsv`. |
| Shared PF01053 resolution | `synthase_features/pipeline/1b_family_resolution/` | Remove automatic installs, freeze MAFFT/trimAl/IQ-TREE versions and seeds/options, and reproduce the six deposited assignments and archived tree. |
| MEME/FIMO motif screen | `synthase_features/pipeline/2_meme_validation/`; `synthase_features/pipeline/phase6_config.yaml` | Separate source from the bundled tool installation, freeze motif parameters and versions, and reproduce the 975 motif records including the GSH2 singleton disposition. |
| Expected-domain validation | `synthase_features/pipeline/2b_domain_validation/`; `synthase_features/pipeline/shared/domain_validation.py` | Freeze the Pfam models and HMMER version, remove automatic installs, and reproduce 973 domain records plus the two explicitly absent MGL-nearest records. |
| Final validation gate | `synthase_features/pipeline/2e_validation_gate/`; `synthase_features/pipeline/shared/validation_gate.py` | Limit the gate to manuscript-active rules and reproduce the 34 removals and 941 retained proteins exactly. |
| Expression-atlas mapping | `synthase_features/pipeline/2_plots/vsc_expression_overview.py`; `synthase_features/pipeline/results_v2/2h_expression_overview/` | Replace hard-coded paths with inputs, capture the BLAST database build and command, freeze BLAST version/options, and reproduce the 168 mappings to 128 atlas genes. |
| Manuscript figures and remaining source tables | Maintained in the separate `vsc-manuscript` repository | Port only figures supported by public inputs, declare plotting dependencies, and compare regenerated files with the submitted figures. |

The GWAS, structure-prediction, Foldseek, structural-QC, geometry, and feature-export stages are outside this manuscript's analysis and should not be added to this workflow.

## External inputs

The complete GMO v1 proteome is not duplicated in this repository. Its source dataset is identified in the main README. Pfam libraries are also not redistributed; the exact profile accessions and versions are in `queries/query_panel.tsv`. Full-length reference queries are identified there by UniProt accession. A full rerun must fetch or stage those inputs under their original terms and verify their checksums before analysis.
