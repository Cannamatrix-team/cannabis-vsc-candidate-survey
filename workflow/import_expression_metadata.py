#!/usr/bin/env python3

import argparse
import csv
import hashlib
import io
import subprocess
import sys
from pathlib import Path

from run_expression_mapping import best_hits, eligible_hits, fasta_lengths, read_blast


ROOT = Path(__file__).resolve().parents[1]
HISTORICAL_COMMIT = "ca6a37a5240e8e0c85b29912d9136fc087cd7d1d"
HISTORICAL_PATH = (
    "synthase_features/pipeline/results_v2/2h_expression_overview/"
    "VSC_ATLAS_gene_metadata.csv"
)
HISTORICAL_SHA256 = "1c04882e3f518d781121e2b384edc39e10a29a42560be5b38209acea92b735c1"
FIELDS = ["GeneID", "Classification", "Specific tissues", "Tau"]


def main():
    parser = argparse.ArgumentParser(
        description="Import the minimal Atlas metadata used by expression sensitivity checks."
    )
    parser.add_argument("--historical-repo", type=Path, required=True)
    parser.add_argument("--commit", default=HISTORICAL_COMMIT)
    parser.add_argument(
        "--blast", type=Path, default=ROOT / "build/expression/blast.tsv"
    )
    parser.add_argument(
        "--candidates",
        type=Path,
        default=ROOT / "data/candidate_sequences_941.fasta",
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=ROOT / "inputs/expression/atlas_gene_metadata_423.tsv",
    )
    args = parser.parse_args()

    data = subprocess.check_output(
        [
            "git",
            "-C",
            str(args.historical_repo),
            "show",
            f"{args.commit}:{HISTORICAL_PATH}",
        ]
    )
    digest = hashlib.sha256(data).hexdigest()
    if digest != HISTORICAL_SHA256:
        raise SystemExit(
            f"Historical metadata checksum mismatch: expected {HISTORICAL_SHA256}, "
            f"found {digest}"
        )

    lengths = fasta_lengths(args.candidates)
    hits = best_hits(read_blast(args.blast))
    subjects = {
        row["subject_acc"] for row in eligible_hits(hits, lengths, 80)
    }
    csv.field_size_limit(sys.maxsize)
    selected = []
    for row in csv.DictReader(io.StringIO(data.decode())):
        if row["GeneID"] in subjects:
            selected.append(
                {
                    field: "" if row[field] == "NA" else row[field]
                    for field in FIELDS
                }
            )
    selected.sort(key=lambda row: row["GeneID"])
    if {row["GeneID"] for row in selected} != subjects or len(selected) != 423:
        raise SystemExit("Historical metadata does not cover the 423 Atlas genes")

    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=FIELDS, delimiter="\t", lineterminator="\n"
        )
        writer.writeheader()
        writer.writerows(selected)
    print(f"Wrote {len(selected)} Atlas metadata rows to {args.output}")


if __name__ == "__main__":
    main()
