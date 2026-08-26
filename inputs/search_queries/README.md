# HMMER search queries

This directory contains only the query assets used by the 42 searches in `queries/query_panel.tsv`:

- `pfam/`: 24 Pfam profile HMMs. Each filename includes the exact model accession and version.
- `uniprot/`: 13 full-length UniProt reference proteins. Each filename is the UniProt accession; the panel records the entry name and sequence version.

These are the exact files preserved in CannamatrixAI commit `ca6a37a5240e8e0c85b29912d9136fc087cd7d1d`, the historical source for the deposited analysis. They are packaged here so the searches can be rerun without relying on the current contents of changing upstream databases.

The full GMO v1 high-confidence proteome remains an external input. See `inputs/hmmer/README.md` for its source and checksum.

Pfam data are available under [CC0 1.0](https://www.ebi.ac.uk/interpro/about/license/). UniProt data are available under [CC BY 4.0](https://www.uniprot.org/help/license). See `THIRD_PARTY_NOTICE.md` for attribution.
