#!/usr/bin/env python3
import argparse
import csv
import shutil
import subprocess
from pathlib import Path

from Bio import Phylo


ROOT = Path(__file__).resolve().parents[1]
PF01053 = ROOT / "phylogeny" / "pf01053"


def require_tool(name):
    path = shutil.which(name)
    if not path:
        raise SystemExit(f"{name} not found; create and activate environment.yml")
    return path


def run_alignment(input_fasta, output_root, mafft_threads):
    alignment_dir = output_root / "alignment"
    alignment_dir.mkdir(parents=True, exist_ok=True)
    aligned = alignment_dir / "aligned.fasta"
    trimmed = alignment_dir / "trimmed.fasta"
    logs = output_root / "logs"
    logs.mkdir(parents=True, exist_ok=True)

    with aligned.open("wb") as stdout, (logs / "mafft.log").open("wb") as stderr:
        subprocess.run(
            [require_tool("mafft"), "--auto", "--thread", str(mafft_threads), "--reorder", str(input_fasta)],
            check=True,
            stdout=stdout,
            stderr=stderr,
        )
    with (logs / "trimal.log").open("wb") as log:
        subprocess.run(
            [require_tool("trimal"), "-in", str(aligned), "-out", str(trimmed), "-automated1"],
            check=True,
            stdout=log,
            stderr=subprocess.STDOUT,
        )
    return trimmed


def run_tree(alignment, output_root, iqtree_threads, bootstraps, seed):
    tree_dir = output_root / "tree"
    logs = output_root / "logs"
    tree_dir.mkdir(parents=True, exist_ok=True)
    logs.mkdir(parents=True, exist_ok=True)
    prefix = tree_dir / "group_tree"
    with (logs / "iqtree.log").open("wb") as log:
        subprocess.run(
            [
                require_tool("iqtree"),
                "-s",
                str(alignment),
                "-B",
                str(bootstraps),
                "-T",
                str(iqtree_threads),
                "-seed",
                str(seed),
                "--prefix",
                str(prefix),
                "-redo",
            ],
            check=True,
            stdout=log,
            stderr=subprocess.STDOUT,
        )
    return Path(f"{prefix}.treefile")


def load_references():
    with (PF01053 / "reference_metadata.tsv").open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def load_candidates():
    with (PF01053 / "group_candidates.tsv").open(newline="") as handle:
        return [row["gene_id"] for row in csv.DictReader(handle, delimiter="\t")]


def write_assignments(tree_path, output_root):
    tree = Phylo.read(tree_path, "newick")
    references = load_references()
    rows = []
    for gene_id in load_candidates():
        nearest = min(references, key=lambda row: tree.distance(gene_id, row["leaf_id"]))
        rows.append(
            {
                "gene_id": gene_id,
                "assigned_family": nearest["family"],
                "nearest_reference_id": nearest["reference_id"],
                "nearest_reference_family": nearest["family"],
                "nearest_reference_distance": tree.distance(gene_id, nearest["leaf_id"]),
                "family_resolution_method": "phylogeny_nearest_reference",
                "family_resolution_status": "resolved",
            }
        )

    output = output_root / "assignment_summary.tsv"
    with output.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=rows[0], delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
    print(f"Wrote {len(rows)} assignments to {output}")


def main():
    parser = argparse.ArgumentParser(description="Reproduce the deposited PF01053 family resolution")
    parser.add_argument("--input-fasta", type=Path, default=PF01053 / "combined_for_phylogeny.fasta")
    parser.add_argument("--alignment", type=Path)
    parser.add_argument("--output-root", type=Path, default=ROOT / "build" / "pf01053")
    parser.add_argument("--mafft-threads", type=int, default=1)
    parser.add_argument("--iqtree-threads", type=int, default=8)
    parser.add_argument("--bootstraps", type=int, default=1000)
    parser.add_argument("--seed", type=int, default=793803)
    args = parser.parse_args()

    args.output_root.mkdir(parents=True, exist_ok=True)
    alignment = args.alignment
    if alignment is None:
        alignment = run_alignment(args.input_fasta, args.output_root, args.mafft_threads)
    tree_path = run_tree(alignment, args.output_root, args.iqtree_threads, args.bootstraps, args.seed)
    write_assignments(tree_path, args.output_root)


if __name__ == "__main__":
    main()
