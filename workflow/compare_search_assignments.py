#!/usr/bin/env python3

import argparse
import csv
from pathlib import Path


def read_table(path):
    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        return reader.fieldnames, list(reader)


def same_evalue(left, right):
    left = float(left)
    right = float(right)
    if left == right:
        return True
    return min(left, right) == 0 and max(left, right) < 1e-300


def main():
    parser = argparse.ArgumentParser(
        description="Compare fresh-HMMER assignments with the deposited table."
    )
    parser.add_argument("expected", type=Path)
    parser.add_argument("generated", type=Path)
    args = parser.parse_args()

    expected_fields, expected_rows = read_table(args.expected)
    generated_fields, generated_rows = read_table(args.generated)
    if generated_fields != expected_fields:
        raise SystemExit("FAIL: assignment columns differ")
    if len(generated_rows) != len(expected_rows):
        raise SystemExit(
            f"FAIL: expected {len(expected_rows)} assignments, "
            f"found {len(generated_rows)}"
        )

    for line_number, (expected, generated) in enumerate(
        zip(expected_rows, generated_rows), start=2
    ):
        for field in expected_fields:
            if field == "evalue":
                matches = same_evalue(expected[field], generated[field])
            else:
                matches = expected[field] == generated[field]
            if not matches:
                raise SystemExit(
                    f"FAIL: assignment line {line_number} differs in {field}: "
                    f"{expected[field]!r} != {generated[field]!r}"
                )

    print(f"Verified {len(expected_rows)} fresh-HMMER search assignments.")


if __name__ == "__main__":
    main()
