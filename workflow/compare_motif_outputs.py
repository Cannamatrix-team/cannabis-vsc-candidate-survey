#!/usr/bin/env python3

import argparse
import csv
import math
from collections import Counter
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
TEXT_FIELDS = [
    "meme_status",
    "meme_status_quality_aware",
    "motif_support_quality",
    "meme_evaluable",
    "singleton_passthrough",
    "n_motifs_found",
    "n_motifs_expected",
    "motif_ids",
]
FLOAT_FIELDS = ["motif_best_evalue", "motif_coverage"]


def read_tsv(path):
    with path.open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    return {row["gene_id"]: row for row in rows}


def equivalent_float(left, right):
    if left == right:
        return True
    if left == "" or right == "":
        return False
    return math.isclose(float(left), float(right), rel_tol=1e-12, abs_tol=0.0)


def main():
    parser = argparse.ArgumentParser(
        description="Compare reconstructed motif records with the deposited ledger."
    )
    parser.add_argument(
        "generated",
        type=Path,
        nargs="?",
        default=ROOT / "build/motifs/motif_screen.tsv",
    )
    parser.add_argument(
        "--ledger", type=Path, default=ROOT / "data/candidate_stage_ledger_975.tsv"
    )
    args = parser.parse_args()

    expected = read_tsv(args.ledger)
    observed = read_tsv(args.generated)
    if set(observed) != set(expected):
        missing = sorted(set(expected) - set(observed))
        extra = sorted(set(observed) - set(expected))
        raise SystemExit(f"Motif identifier mismatch: missing={missing}, extra={extra}")

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

    if mismatches:
        for mismatch in mismatches[:20]:
            print("\t".join(mismatch))
        raise SystemExit(f"Motif comparison failed with {len(mismatches)} mismatches")

    statuses = Counter(row["meme_status"] for row in observed.values())
    evaluable = Counter(row["meme_evaluable"] for row in observed.values())
    passthrough = Counter(row["singleton_passthrough"] for row in observed.values())
    print(
        f"Motif records match: {len(observed)} proteins; "
        f"statuses={dict(statuses)}; evaluable={dict(evaluable)}; "
        f"singleton_passthrough={dict(passthrough)}"
    )


if __name__ == "__main__":
    main()
