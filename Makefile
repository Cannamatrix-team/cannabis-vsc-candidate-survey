.PHONY: reporting verify-release

reporting:
	python3 workflow/generate_reporting_tables.py

verify-release:
	sha256sum -c SHA256SUMS
	python3 workflow/verify_release.py
	python3 workflow/generate_reporting_tables.py --check
