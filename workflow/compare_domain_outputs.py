#!/usr/bin/env python3

import argparse
import csv
import math
from collections import Counter
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
TEXT_FIELDS = [
    "pfam_n_domains",
    "pfam_domain_architecture",
    "pfam_confirmation_rule",
    "pfam_required_ids",
    "pfam_required_n_domains",
    "pfam_required_n_present",
    "domain_confirmation_status",
    "domain_validated",
]
FLOAT_FIELDS = ["pfam_best_evalue", "pfam_total_score"]


def read_tsv(path):
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def equivalent_float(left, right):
    if left == right:
        return True
    if left == "" or right == "":
        return False
    return math.isclose(float(left), float(right), rel_tol=1e-12, abs_tol=0.0)


def main():
    parser = argparse.ArgumentParser(
        description="Compare regenerated expected-domain records with the deposit."
    )
    parser.add_argument(
        "generated",
        type=Path,
        nargs="?",
        default=ROOT / "build/domains/domain_validation.tsv",
    )
    parser.add_argument(
        "--ledger", type=Path, default=ROOT / "data/candidate_stage_ledger_975.tsv"
    )
    args = parser.parse_args()

    ledger = read_tsv(args.ledger)
    expected = {
        row["gene_id"]: row
        for row in ledger
        if row["expected_domain_record_present"] == "true"
    }
    observed_rows = read_tsv(args.generated)
    observed = {row["gene_id"]: row for row in observed_rows}
    if len(observed) != len(observed_rows):
        raise SystemExit("Generated expected-domain output has duplicate identifiers")
    if set(observed) != set(expected):
        missing = sorted(set(expected) - set(observed))
        extra = sorted(set(observed) - set(expected))
        raise SystemExit(
            f"Domain identifier mismatch: missing={missing}, extra={extra}"
        )

    mismatches = []
    for gene_id, expected_row in expected.items():
        observed_row = observed[gene_id]
        for field in TEXT_FIELDS:
            if observed_row[field] != expected_row[field]:
                mismatches.append(
                    (gene_id, field, expected_row[field], observed_row[field])
                )
        for field in FLOAT_FIELDS:
            if not equivalent_float(observed_row[field], expected_row[field]):
                mismatches.append(
                    (gene_id, field, expected_row[field], observed_row[field])
                )

    absent_records = [
        row for row in ledger if row["expected_domain_record_present"] == "false"
    ]
    if len(absent_records) != 2 or {
        row["reporting_category"] for row in absent_records
    } != {"mgl"}:
        raise SystemExit("Expected exactly two explicit MGL no-record dispositions")
    if mismatches:
        for mismatch in mismatches[:20]:
            print("\t".join(mismatch))
        raise SystemExit(f"Domain comparison failed with {len(mismatches)} mismatches")

    statuses = Counter(row["domain_confirmation_status"] for row in observed.values())
    print(
        f"Expected-domain records match: {len(observed)} proteins; "
        f"statuses={dict(statuses)}; explicit_no_record=2"
    )


if __name__ == "__main__":
    main()
