.PHONY: check-hmmer-tools check-proteome rebuild-search reporting run-hmmer verify-hmmer verify-release verify-search

check-proteome:
	@test -n "$(VSC_PROTEOME)" || { echo "Set VSC_PROTEOME=/path/to/GMO.v1.primary_high_confidence.proteins.fasta"; exit 2; }
	@test -f "$(VSC_PROTEOME)" || { echo "Proteome not found: $(VSC_PROTEOME)"; exit 2; }

check-hmmer-tools:
	@command -v hmmsearch >/dev/null || { echo "hmmsearch not found; create and activate environment.yml"; exit 2; }
	@command -v phmmer >/dev/null || { echo "phmmer not found; create and activate environment.yml"; exit 2; }

rebuild-search: check-proteome
	python3 workflow/rebuild_search_assignments.py \
		--hmmer-root inputs/hmmer \
		--proteome "$(VSC_PROTEOME)"

run-hmmer: check-proteome check-hmmer-tools
	python3 workflow/run_hmmer_searches.py --proteome "$(VSC_PROTEOME)"

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

verify-release:
	sha256sum -c SHA256SUMS
	python3 workflow/verify_release.py
	python3 workflow/generate_reporting_tables.py --check
