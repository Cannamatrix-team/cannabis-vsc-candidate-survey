.PHONY: rebuild-search reporting verify-release verify-search

rebuild-search:
	@test -n "$(VSC_PROTEOME)" || { echo "Set VSC_PROTEOME=/path/to/GMO.v1.primary_high_confidence.proteins.fasta"; exit 2; }
	python3 workflow/rebuild_search_assignments.py \
		--hmmer-root inputs/hmmer \
		--proteome "$(VSC_PROTEOME)"

reporting:
	python3 workflow/generate_reporting_tables.py

verify-search: rebuild-search
	diff -u data/search_assignments_1005.tsv build/search_assignments_1005.tsv

verify-release:
	sha256sum -c SHA256SUMS
	python3 workflow/verify_release.py
	python3 workflow/generate_reporting_tables.py --check
