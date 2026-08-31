#!/usr/bin/env python3

import argparse
import csv
import hashlib
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]


def sha256(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    parser = argparse.ArgumentParser(
        description="Compare regenerated manuscript artifacts with the frozen submission."
    )
    parser.add_argument(
        "--generated-root", type=Path, default=ROOT / "build/manuscript"
    )
    args = parser.parse_args()

    with (ROOT / "manuscript/figure_checksums.tsv").open(
        newline="", encoding="utf-8"
    ) as handle:
        figures = list(csv.DictReader(handle, delimiter="\t"))
    for row in figures:
        path = args.generated_root / "figures" / row["artifact"]
        observed = sha256(path)
        if observed != row["sha256"]:
            raise SystemExit(
                f"Figure mismatch for {row['artifact']}: "
                f"expected {row['sha256']}, found {observed}"
            )
    print(f"Manuscript figures match: {len(figures)} PNG/PDF artifacts")

    expected_root = ROOT / "manuscript/tables"
    generated_root = args.generated_root / "tables"
    expected = sorted(path.name for path in expected_root.glob("*.tsv"))
    generated = sorted(path.name for path in generated_root.glob("*.tsv"))
    if generated != expected:
        raise SystemExit(
            f"Manuscript table set mismatch: expected={expected}, generated={generated}"
        )
    for name in expected:
        if (expected_root / name).read_bytes() != (generated_root / name).read_bytes():
            raise SystemExit(f"Manuscript table mismatch: {name}")
    print(f"Manuscript source tables match: {len(expected)} artifacts")


if __name__ == "__main__":
    main()
