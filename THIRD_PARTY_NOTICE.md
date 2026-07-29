# Third-party sources and attribution

This repository combines Cannamatrix-created annotations and tables with selected material derived from public biological resources. The repository-level CC BY 4.0 license applies to Cannamatrix's original contributions and does not replace the terms attached to upstream material.

## Cannabis GMO v1 protein models

`data/candidate_sequences_941.fasta` is a filtered extract of the GMO v1 primary high-confidence protein models. Candidate identifiers and derived tables also refer to those models.

- Source publication: Lynch et al. (2025), *Nature*, DOI: [10.1038/s41586-025-09065-0](https://doi.org/10.1038/s41586-025-09065-0)
- Source annotation data: *Cannabis Pangenome Annotation Data*, DOI: [10.25452/figshare.plus.25909024.v1](https://doi.org/10.25452/figshare.plus.25909024.v1)
- Source collection: *Cannabis Pangenome*, DOI: [10.25452/figshare.plus.c.7248427.v1](https://doi.org/10.25452/figshare.plus.c.7248427.v1)
- Source-data license: CC0 1.0 as stated on the Figshare record

## UniProt reference proteins

The PF01053 phylogeny includes reference protein sequences obtained from UniProt. `queries/query_panel.tsv` identifies the full-length reference queries by UniProt accession and entry name; it does not redistribute the remaining query FASTA files.

- Resource: [UniProt](https://www.uniprot.org/)
- License: [CC BY 4.0](https://www.uniprot.org/help/license)

## Pfam profiles

`queries/query_panel.tsv` records Pfam accessions and model versions used in the searches. The Pfam HMM libraries are not redistributed in this repository.

- Resource: [Pfam through InterPro](https://www.ebi.ac.uk/interpro/entry/pfam/)
- License information: [InterPro data licensing](https://www.ebi.ac.uk/interpro/about/license/)

## Cannabis Expression Atlas

The full Cannabis Expression Atlas database is not redistributed. Files under `data/expression/` are derived cross-reference and sensitivity tables produced by matching GMO v1 candidates to the atlas proteins.

- Source publication: Barbosa-Xavier et al. (2024), *Physiologia Plantarum*, DOI: [10.1111/ppl.70010](https://doi.org/10.1111/ppl.70010)

Users should cite the relevant upstream resource together with this repository when reusing source-derived sequences or annotations.
