#!/usr/bin/env python3
import argparse
import csv
import math
import re
from pathlib import Path

from Bio import Phylo


ROOT = Path(__file__).resolve().parents[1]
EXPECTED = ROOT / "phylogeny" / "pf01053"
ASSIGNMENT_FIELDS = (
    "assigned_family",
    "nearest_reference_id",
    "nearest_reference_family",
    "family_resolution_method",
    "family_resolution_status",
)


def read_assignments(path):
    with path.open(newline="") as handle:
        return {row["gene_id"]: row for row in csv.DictReader(handle, delimiter="\t")}


def compare_assignments(generated_root):
    expected = read_assignments(EXPECTED / "assignment_summary.tsv")
    generated = read_assignments(generated_root / "assignment_summary.tsv")
    if generated.keys() != expected.keys():
        raise SystemExit("FAIL: generated candidate identifiers differ from the deposit")
    for gene_id, expected_row in expected.items():
        generated_row = generated[gene_id]
        for field in ASSIGNMENT_FIELDS:
            if generated_row[field] != expected_row[field]:
                raise SystemExit(
                    f"FAIL: {gene_id} {field}: expected {expected_row[field]}, got {generated_row[field]}"
                )
    counts = {}
    for row in generated.values():
        counts[row["assigned_family"]] = counts.get(row["assigned_family"], 0) + 1
    if counts != {"CBL": 2, "CGS": 2, "MGL": 2}:
        raise SystemExit(f"FAIL: unexpected family counts: {counts}")
    print("PASS: six deposited nearest-reference assignments (CBL 2, CGS 2, MGL 2)")


def unrooted_splits(tree):
    terminals = frozenset(node.name for node in tree.get_terminals())
    splits = set()
    for clade in tree.get_nonterminals():
        side = frozenset(node.name for node in clade.get_terminals())
        other = terminals - side
        if len(side) < 2 or len(other) < 2:
            continue
        canonical = min((side, other), key=lambda item: (len(item), sorted(item)))
        splits.add(canonical)
    return terminals, splits


def iqtree_summary(path):
    text = path.read_text()
    model = re.search(r"Best-fit model according to BIC: (\S+)", text)
    likelihood = re.search(r"Log-likelihood of the tree: (-?\d+(?:\.\d+)?)", text)
    if not model or not likelihood:
        raise SystemExit(f"FAIL: could not parse {path}")
    return model.group(1), float(likelihood.group(1))


def compare_archived_topology(generated_root):
    expected_tree = Phylo.read(EXPECTED / "tree" / "group_tree.treefile", "newick")
    generated_tree = Phylo.read(generated_root / "tree" / "group_tree.treefile", "newick")
    expected_terminals, expected_splits = unrooted_splits(expected_tree)
    generated_terminals, generated_splits = unrooted_splits(generated_tree)
    if generated_terminals != expected_terminals:
        raise SystemExit("FAIL: generated tree terminals differ from the archived tree")
    if generated_splits != expected_splits:
        raise SystemExit("FAIL: generated unrooted topology differs from the archived tree")

    expected_model, expected_likelihood = iqtree_summary(EXPECTED / "tree" / "group_tree.iqtree")
    generated_model, generated_likelihood = iqtree_summary(generated_root / "tree" / "group_tree.iqtree")
    if generated_model != expected_model:
        raise SystemExit(f"FAIL: expected model {expected_model}, got {generated_model}")
    if not math.isclose(generated_likelihood, expected_likelihood, abs_tol=1e-4):
        raise SystemExit(
            f"FAIL: expected log-likelihood {expected_likelihood}, got {generated_likelihood}"
        )
    print(
        f"PASS: archived unrooted topology ({len(generated_splits)} splits), "
        f"{generated_model} model, and log-likelihood"
    )


def main():
    parser = argparse.ArgumentParser(description="Compare a PF01053 reconstruction with the deposit")
    parser.add_argument("generated_root", type=Path)
    parser.add_argument("--require-archived-topology", action="store_true")
    args = parser.parse_args()
    compare_assignments(args.generated_root)
    if args.require_archived_topology:
        compare_archived_topology(args.generated_root)


if __name__ == "__main__":
    main()
