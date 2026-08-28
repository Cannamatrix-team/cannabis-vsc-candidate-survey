#!/usr/bin/env python3

import argparse
import csv
import math
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
TEXT_FIELDS = [
    "subject_acc",
    "family_slug",
    "Classification",
    "Specific tissues",
]
INTEGER_FIELDS = [
    "alignment_length",
    "mismatches",
    "gap_opens",
    "q_start",
    "q_end",
    "s_start",
    "s_end",
    "query_length",
]
FLOAT_FIELDS = ["p_identity", "evalue", "bit_score", "Tau", "query_coverage"]


def read_tsv(path):
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def equivalent_float(left, right):
    if left == right:
        return True
    if left == "" or right == "":
        return False
    return math.isclose(float(left), float(right), rel_tol=1e-12, abs_tol=0.0)


def compare_bridge(generated, deposited):
    observed_rows = read_tsv(generated)
    expected_rows = read_tsv(deposited)
    observed = {row["query_acc"]: row for row in observed_rows}
    expected = {row["query_acc"]: row for row in expected_rows}
    if len(observed) != len(observed_rows):
        raise SystemExit("Generated expression bridge has duplicate queries")
    if set(observed) != set(expected):
        missing = sorted(set(expected) - set(observed))
        extra = sorted(set(observed) - set(expected))
        raise SystemExit(
            f"Expression query mismatch: missing={missing}, extra={extra}"
        )

    mismatches = []
    for query, expected_row in expected.items():
        observed_row = observed[query]
        for field in TEXT_FIELDS + INTEGER_FIELDS:
            if observed_row[field] != expected_row[field]:
                mismatches.append(
                    (query, field, expected_row[field], observed_row[field])
                )
        for field in FLOAT_FIELDS:
            if not equivalent_float(observed_row[field], expected_row[field]):
                mismatches.append(
                    (query, field, expected_row[field], observed_row[field])
                )
    if mismatches:
        for mismatch in mismatches[:20]:
            print("\t".join(mismatch))
        raise SystemExit(
            f"Expression bridge comparison failed with {len(mismatches)} mismatches"
        )
    print(
        f"Expression bridge matches: {len(observed)} candidate mappings; "
        f"{len({row['subject_acc'] for row in observed_rows})} Atlas genes"
    )


def compare_sensitivity(generated, deposited):
    observed = read_tsv(generated)
    expected = read_tsv(deposited)
    if observed != expected:
        raise SystemExit("Expression identity-sensitivity table differs from the deposit")
    print(f"Expression sensitivity matches: {len(observed)} identity thresholds")


def main():
    parser = argparse.ArgumentParser(
        description="Compare regenerated expression mappings with the deposit."
    )
    parser.add_argument(
        "--generated-root", type=Path, default=ROOT / "build/expression"
    )
    parser.add_argument(
        "--deposited-root", type=Path, default=ROOT / "data/expression"
    )
    args = parser.parse_args()
    compare_bridge(
        args.generated_root / "expression_bridge_exact100_qcov80.tsv",
        args.deposited_root / "expression_bridge_exact100_qcov80.tsv",
    )
    compare_sensitivity(
        args.generated_root / "expression_identity_sensitivity.tsv",
        args.deposited_root / "expression_identity_sensitivity.tsv",
    )


if __name__ == "__main__":
    main()
