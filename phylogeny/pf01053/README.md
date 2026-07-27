# PF01053 nearest-reference analysis

The shared PF01053 search recovered six unique GMO v1 proteins. Candidate and reference sequences were aligned with MAFFT, trimmed with trimAl, and analyzed by maximum likelihood with IQ-TREE 3.1.1. The LG+R2 model was selected by Bayesian information criterion, with 1,000 ultrafast bootstrap replicates.

The 41-reference panel and shortest-patristic-distance assignment produced two CBL-nearest, two CGS-nearest, and two MGL-nearest candidates.

- `assignment_summary.tsv`: candidate assignments and nearest-reference distances.
- `reference_metadata.tsv`: reference identifiers, families, organisms, and review status.
- `group_candidates.tsv` and `.fasta`: the six GMO v1 candidates.
- `combined_for_phylogeny.fasta`: candidates and reference sequences used for alignment.
- `alignment/`: full and trimmed alignments.
- `tree/group_tree.treefile`: maximum-likelihood tree.
- `tree/group_tree.contree`: consensus tree with branch support.
- `tree/group_tree.iqtree`: IQ-TREE model-selection and run report.
