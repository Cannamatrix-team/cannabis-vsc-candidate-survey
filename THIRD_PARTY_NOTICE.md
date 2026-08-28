# Third-party sources and attribution

This repository combines Cannamatrix-created annotations and tables with selected material derived from public biological resources. The repository-level CC BY 4.0 license applies to Cannamatrix's original contributions and does not replace the terms attached to upstream material.

## Cannabis GMO v1 protein models

`data/candidate_sequences_975.fasta` and `data/candidate_sequences_941.fasta` are filtered extracts of the GMO v1 primary high-confidence protein models. Candidate identifiers and derived tables also refer to those models.

- Source publication: Lynch et al. (2025), *Nature*, DOI: [10.1038/s41586-025-09065-0](https://doi.org/10.1038/s41586-025-09065-0)
- Source annotation data: *Cannabis Pangenome Annotation Data*, DOI: [10.25452/figshare.plus.25909024.v1](https://doi.org/10.25452/figshare.plus.25909024.v1)
- Source collection: *Cannabis Pangenome*, DOI: [10.25452/figshare.plus.c.7248427.v1](https://doi.org/10.25452/figshare.plus.c.7248427.v1)
- Source-data license: CC0 1.0 as stated on the Figshare record

## UniProt reference proteins

`inputs/search_queries/uniprot/` contains the 13 UniProt protein sequences used as full-length search queries. The PF01053 phylogeny also includes UniProt reference sequences. `queries/query_panel.tsv` records each search query's accession, entry name, and sequence version.

- Resource: [UniProt](https://www.uniprot.org/)
- License: [CC BY 4.0](https://www.uniprot.org/help/license)

## Pfam profiles

`inputs/search_queries/pfam/` contains the 24 exact Pfam models used by the search panel, rather than a complete Pfam library. `queries/query_panel.tsv` records their accessions and model versions.

- Resource: [Pfam through InterPro](https://www.ebi.ac.uk/interpro/entry/pfam/)
- License: [CC0 1.0](https://www.ebi.ac.uk/interpro/about/license/)

## Cannabis Expression Atlas

The full Cannabis Expression Atlas database is not redistributed. Files under `data/expression/` are derived cross-reference and sensitivity tables produced by matching GMO v1 candidates to the atlas proteins.

- Source publication: Barbosa-Xavier et al. (2024), *Physiologia Plantarum*, DOI: [10.1111/ppl.70010](https://doi.org/10.1111/ppl.70010)

## MEME Suite

The compressed XML files under `motifs/meme/` are Cannamatrix analysis outputs generated with MEME Suite 5.5.9; they do not redistribute the MEME software. The reproducibility environment obtains MEME Suite from Bioconda under the software's upstream custom license.

- Software: [MEME Suite](https://meme-suite.org/)
- Version: 5.5.9
- License: [MEME Suite copyright and license](https://meme-suite.org/meme/doc/copyright.html)

Users should cite the relevant upstream resource together with this repository when reusing source-derived sequences or annotations.
