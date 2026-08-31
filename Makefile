.PHONY: check-atlas check-blast-tools check-domain-tools check-fimo-tools check-hmmer-tools check-manuscript-tools check-meme-tools check-phylogeny-tools check-proteome rebuild-search reporting run-domains run-expression run-gate run-hmmer run-manuscript run-motif-discovery run-motifs run-pf01053 verify-domains verify-expression verify-gate verify-hmmer verify-manuscript verify-motif-discovery verify-motifs verify-pf01053 verify-release verify-search

DOMAIN_THREADS ?= 4
MOTIF_THREADS ?= 4

check-proteome:
	@test -n "$(VSC_PROTEOME)" || { echo "Set VSC_PROTEOME=/path/to/GMO.v1.primary_high_confidence.proteins.fasta"; exit 2; }
	@test -f "$(VSC_PROTEOME)" || { echo "Proteome not found: $(VSC_PROTEOME)"; exit 2; }

check-atlas:
	@test -n "$(ATLAS_PROTEINS)" || { echo "Set ATLAS_PROTEINS=/path/to/CEA_protein_sequences_all.faa"; exit 2; }
	@test -f "$(ATLAS_PROTEINS)" || { echo "Atlas protein FASTA not found: $(ATLAS_PROTEINS)"; exit 2; }

check-blast-tools:
	@command -v makeblastdb >/dev/null || { echo "makeblastdb not found; create and activate environment.yml"; exit 2; }
	@command -v blastp >/dev/null || { echo "blastp not found; create and activate environment.yml"; exit 2; }

check-hmmer-tools:
	@command -v hmmsearch >/dev/null || { echo "hmmsearch not found; create and activate environment.yml"; exit 2; }
	@command -v phmmer >/dev/null || { echo "phmmer not found; create and activate environment.yml"; exit 2; }

check-domain-tools:
	@command -v hmmscan >/dev/null || { echo "hmmscan not found; create and activate environment.yml"; exit 2; }
	@command -v hmmpress >/dev/null || { echo "hmmpress not found; create and activate environment.yml"; exit 2; }

check-phylogeny-tools:
	@command -v mafft >/dev/null || { echo "mafft not found; create and activate environment.yml"; exit 2; }
	@command -v trimal >/dev/null || { echo "trimal not found; create and activate environment.yml"; exit 2; }
	@command -v iqtree >/dev/null || { echo "iqtree not found; create and activate environment.yml"; exit 2; }

check-fimo-tools:
	@command -v fimo >/dev/null || { echo "fimo not found; create and activate environment.yml"; exit 2; }

check-meme-tools: check-fimo-tools
	@command -v meme >/dev/null || { echo "meme not found; create and activate environment.yml"; exit 2; }

check-manuscript-tools:
	@python3 -c "import Bio, matplotlib, numpy" || { echo "Create and activate manuscript/environment.yml"; exit 2; }

rebuild-search: check-proteome
	python3 workflow/rebuild_search_assignments.py \
		--hmmer-root inputs/hmmer \
		--proteome "$(VSC_PROTEOME)"

run-hmmer: check-proteome check-hmmer-tools
	python3 workflow/run_hmmer_searches.py --proteome "$(VSC_PROTEOME)"

run-pf01053: check-phylogeny-tools
	python3 workflow/run_pf01053_resolution.py

run-domains: check-domain-tools
	python3 workflow/run_domain_validation.py --threads "$(DOMAIN_THREADS)"

run-expression: check-atlas check-blast-tools
	python3 workflow/run_expression_mapping.py \
		--atlas-proteins "$(ATLAS_PROTEINS)"

run-manuscript: check-manuscript-tools
	python3 manuscript/generate_paper_figures.py

run-gate: run-motifs run-domains
	python3 workflow/run_validation_gate.py

run-motifs: check-fimo-tools
	python3 workflow/run_motif_screen.py

run-motif-discovery: check-meme-tools
	python3 workflow/run_motif_screen.py \
		--discover \
		--meme-threads "$(MOTIF_THREADS)" \
		--output-root build/motif-discovery

reporting:
	python3 workflow/generate_reporting_tables.py

verify-search: rebuild-search
	diff -u data/search_assignments_1005.tsv build/search_assignments_1005.tsv

verify-hmmer: run-hmmer
	python3 workflow/compare_hmmer_outputs.py
	python3 workflow/rebuild_search_assignments.py \
		--hmmer-root build/hmmer \
		--proteome "$(VSC_PROTEOME)" \
		--output build/search_assignments_fresh_hmmer.tsv
	python3 workflow/compare_search_assignments.py \
		data/search_assignments_1005.tsv \
		build/search_assignments_fresh_hmmer.tsv

verify-pf01053: run-pf01053
	python3 workflow/compare_pf01053_outputs.py build/pf01053
	python3 workflow/run_pf01053_resolution.py \
		--alignment phylogeny/pf01053/alignment/trimmed.fasta \
		--output-root build/pf01053-archive-tree
	python3 workflow/compare_pf01053_outputs.py \
		build/pf01053-archive-tree \
		--require-archived-topology

verify-domains: run-domains
	python3 workflow/compare_domain_outputs.py

verify-expression: run-expression
	python3 workflow/compare_expression_outputs.py

verify-manuscript: run-manuscript
	python3 workflow/compare_manuscript_artifacts.py

verify-gate: run-gate
	python3 workflow/compare_gate_outputs.py

verify-motifs: run-motifs
	python3 workflow/compare_motif_outputs.py

verify-motif-discovery: run-motif-discovery
	python3 workflow/compare_motif_models.py
	python3 workflow/compare_motif_outputs.py build/motif-discovery/motif_screen.tsv

verify-release:
	sha256sum -c SHA256SUMS
	python3 workflow/verify_release.py
	python3 workflow/generate_reporting_tables.py --check
