#!/usr/bin/env python3

import argparse
import math
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]


def contents(path):
    program = None
    version = None
    records = []
    with path.open(encoding="utf-8") as handle:
        for line in handle:
            if line.startswith("# Program:"):
                program = line.split(":", 1)[1].strip()
            elif line.startswith("# Version:"):
                version = line.split(":", 1)[1].strip()
            elif line.strip() and not line.startswith("#"):
                records.append(line.split())
    return program, version, records


def same_evalue(left, right):
    left = float(left)
    right = float(right)
    if left == right:
        return True
    if min(left, right) == 0:
        return max(left, right) < 1e-300
    return math.isclose(left, right, rel_tol=0.1)


def same_contents(generated, archived):
    generated_program, generated_version, generated_records = generated
    archived_program, archived_version, archived_records = archived
    if (generated_program, generated_version) != (archived_program, archived_version):
        return False
    if len(generated_records) != len(archived_records):
        return False
    for generated_record, archived_record in zip(generated_records, archived_records):
        if len(generated_record) != len(archived_record):
            return False
        for index, (generated_value, archived_value) in enumerate(
            zip(generated_record, archived_record)
        ):
            if index in (4, 7):
                if not same_evalue(generated_value, archived_value):
                    return False
            elif generated_value != archived_value:
                return False
    return True


def main():
    parser = argparse.ArgumentParser(
        description="Compare regenerated HMMER tblout records with the archived searches."
    )
    parser.add_argument(
        "--archived-root", type=Path, default=ROOT / "inputs" / "hmmer"
    )
    parser.add_argument(
        "--generated-root", type=Path, default=ROOT / "build" / "hmmer"
    )
    args = parser.parse_args()

    archived = sorted(args.archived_root.rglob("*.tbl"))
    for expected in archived:
        relative = expected.relative_to(args.archived_root)
        generated = args.generated_root / relative
        if not same_contents(contents(generated), contents(expected)):
            raise SystemExit(f"FAIL: regenerated HMMER output differs: {relative}")
        print(f"PASS: {relative}")
    print(f"Verified {len(archived)} regenerated HMMER searches.")


if __name__ == "__main__":
    main()
