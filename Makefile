.PHONY: check-hmmer-tools check-phylogeny-tools check-proteome rebuild-search reporting run-hmmer run-pf01053 verify-hmmer verify-pf01053 verify-release verify-search

check-proteome:
	@test -n "$(VSC_PROTEOME)" || { echo "Set VSC_PROTEOME=/path/to/GMO.v1.primary_high_confidence.proteins.fasta"; exit 2; }
	@test -f "$(VSC_PROTEOME)" || { echo "Proteome not found: $(VSC_PROTEOME)"; exit 2; }

check-hmmer-tools:
	@command -v hmmsearch >/dev/null || { echo "hmmsearch not found; create and activate environment.yml"; exit 2; }
	@command -v phmmer >/dev/null || { echo "phmmer not found; create and activate environment.yml"; exit 2; }

check-phylogeny-tools:
	@command -v mafft >/dev/null || { echo "mafft not found; create and activate environment.yml"; exit 2; }
	@command -v trimal >/dev/null || { echo "trimal not found; create and activate environment.yml"; exit 2; }
	@command -v iqtree >/dev/null || { echo "iqtree not found; create and activate environment.yml"; exit 2; }

rebuild-search: check-proteome
	python3 workflow/rebuild_search_assignments.py \
		--hmmer-root inputs/hmmer \
		--proteome "$(VSC_PROTEOME)"

run-hmmer: check-proteome check-hmmer-tools
	python3 workflow/run_hmmer_searches.py --proteome "$(VSC_PROTEOME)"

run-pf01053: check-phylogeny-tools
	python3 workflow/run_pf01053_resolution.py

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

verify-release:
	sha256sum -c SHA256SUMS
	python3 workflow/verify_release.py
	python3 workflow/generate_reporting_tables.py --check
