#!/usr/bin/env python3

import argparse
import csv
import subprocess
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
PROTEOME_SUBDIR = "GMO.v1.primary_high_confidence.proteins"


def read_panel(path):
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def output_directory(root, panel_row):
    if panel_row["discovery_round"] == "round2":
        root /= "round2_results"
    return root / panel_row["pipeline_family"] / PROTEOME_SUBDIR


def run(command):
    print("Running", command[0], Path(command[-2]).name)
    subprocess.run(command, check=True)


def main():
    parser = argparse.ArgumentParser(
        description="Rerun the 42 HMMER searches used by the VSC candidate survey."
    )
    parser.add_argument("--proteome", type=Path, required=True)
    parser.add_argument(
        "--query-panel", type=Path, default=ROOT / "queries" / "query_panel.tsv"
    )
    parser.add_argument(
        "--queries-root",
        type=Path,
        default=ROOT / "inputs" / "search_queries",
    )
    parser.add_argument(
        "--output-root", type=Path, default=ROOT / "build" / "hmmer"
    )
    args = parser.parse_args()

    for panel_row in read_panel(args.query_panel):
        destination = output_directory(args.output_root, panel_row)
        destination.mkdir(parents=True, exist_ok=True)

        for profile in filter(None, panel_row["pfam_profiles"].split(";")):
            query = args.queries_root / "pfam" / f"{profile}.hmm"
            name = profile.split(".", 1)[0]
            run(
                [
                    "hmmsearch",
                    "-Z",
                    "567483",
                    "--incE",
                    "0.001",
                    "--cpu",
                    "1",
                    "-A",
                    str(destination / f"{name}.sto"),
                    "--tblout",
                    str(destination / f"{name}.tbl"),
                    "-o",
                    str(destination / f"{name}.log"),
                    str(query),
                    str(args.proteome),
                ]
            )

        references = filter(
            None, panel_row["full_length_reference_queries"].split(";")
        )
        for reference in references:
            accession, entry_name, _ = reference.split("/")
            query = args.queries_root / "uniprot" / f"{accession}.fasta"
            name = f"{entry_name}.fasta"
            run(
                [
                    "phmmer",
                    "--incE",
                    "0.001",
                    "--cpu",
                    "1",
                    "-A",
                    str(destination / f"{name}.sto"),
                    "--tblout",
                    str(destination / f"{name}.tbl"),
                    "-o",
                    str(destination / f"{name}.log"),
                    str(query),
                    str(args.proteome),
                ]
            )


if __name__ == "__main__":
    main()
