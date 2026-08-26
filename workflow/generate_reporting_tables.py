#!/usr/bin/env python3

import argparse
import csv
from collections import Counter
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]


def read_tsv(relative_path):
    with (ROOT / relative_path).open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        return reader.fieldnames, list(reader)


def write_tsv(path, fieldnames, rows):
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=fieldnames, delimiter="\t", lineterminator="\n"
        )
        writer.writeheader()
        writer.writerows(rows)


def category_counts():
    _, catalog = read_tsv("data/candidate_catalog_941.tsv")
    fieldnames, template = read_tsv("data/candidate_counts_by_category.tsv")
    counts = Counter(row["family"].lower() for row in catalog)
    slugs = {row["original_pipeline_slug"] for row in template}
    if set(counts) != slugs:
        raise SystemExit("Catalog categories do not match the reporting template")
    rows = []
    for row in template:
        generated = dict(row)
        generated["candidate_count"] = str(counts[row["original_pipeline_slug"]])
        rows.append(generated)
    return fieldnames, rows


def expression_counts():
    _, bridge = read_tsv("data/expression/expression_bridge_exact100_qcov80.tsv")
    fieldnames, template = read_tsv(
        "data/expression/expression_exact100_candidate_mappings_by_bin.tsv"
    )
    queries_by_family = {}
    for row in bridge:
        queries_by_family.setdefault(row["family_slug"], set()).add(row["query_acc"])
    if set(queries_by_family) != {row["family_slug"] for row in template}:
        raise SystemExit("Expression categories do not match the reporting template")
    rows = []
    for row in template:
        generated = dict(row)
        generated["candidate_mappings"] = str(
            len(queries_by_family[row["family_slug"]])
        )
        rows.append(generated)
    return fieldnames, rows


def require_matches(relative_path, generated):
    _, deposited = read_tsv(relative_path)
    if generated != deposited:
        raise SystemExit(f"FAIL: regenerated {relative_path} differs from the deposit")
    print(f"PASS: regenerated {relative_path}")


def parse_args():
    parser = argparse.ArgumentParser(
        description="Regenerate publication reporting tables from the deposited candidate data."
    )
    parser.add_argument("--check", action="store_true", help="Compare with deposited tables")
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=ROOT / "build" / "reporting",
        help="Output directory used without --check",
    )
    return parser.parse_args()


def main():
    args = parse_args()
    category_fields, categories = category_counts()
    expression_fields, expression = expression_counts()
    outputs = [
        (
            "data/candidate_counts_by_category.tsv",
            category_fields,
            categories,
        ),
        (
            "data/expression/expression_exact100_candidate_mappings_by_bin.tsv",
            expression_fields,
            expression,
        ),
    ]
    if args.check:
        for relative_path, _, rows in outputs:
            require_matches(relative_path, rows)
        return
    for relative_path, fieldnames, rows in outputs:
        destination = args.output_dir / Path(relative_path).name
        write_tsv(destination, fieldnames, rows)
        print(f"Wrote {destination}")


if __name__ == "__main__":
    main()
