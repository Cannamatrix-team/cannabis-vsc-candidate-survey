#!/usr/bin/env python3

import argparse
import csv
import gzip
import hashlib
import subprocess
import xml.etree.ElementTree as ET
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
HISTORICAL_COMMIT = "ca6a37a5240e8e0c85b29912d9136fc087cd7d1d"
HISTORICAL_ROOT = "synthase_features/pipeline/results_v2"


def git_show(repo, commit, path):
    return subprocess.check_output(["git", "-C", str(repo), "show", f"{commit}:{path}"])


def read_ledger(path):
    with path.open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    return rows, {row["gene_id"] for row in rows}


def filter_fasta(data, selected_ids):
    records = []
    current = []
    current_id = None

    for line in data.decode().splitlines(keepends=True):
        if line.startswith(">"):
            if current_id in selected_ids:
                records.extend(current)
            current = [line]
            current_id = line[1:].split()[0]
        else:
            current.append(line)
    if current_id in selected_ids:
        records.extend(current)

    found = {line[1:].split()[0] for line in records if line.startswith(">")}
    if found != selected_ids:
        missing = sorted(selected_ids - found)
        extra = sorted(found - selected_ids)
        raise ValueError(f"FASTA identifier mismatch: missing={missing}, extra={extra}")
    return "".join(records).encode()


def xml_summary(data, expected_ids):
    root = ET.fromstring(data)
    training_set = root.find("training_set")
    model = root.find("model")
    sequence_ids = {node.get("name") for node in training_set.findall("sequence")}
    if sequence_ids != expected_ids:
        missing = sorted(expected_ids - sequence_ids)
        extra = sorted(sequence_ids - expected_ids)
        raise ValueError(
            f"MEME training-set mismatch: missing={missing}, extra={extra}"
        )
    return {
        "meme_version": root.get("version", ""),
        "n_sequences": len(sequence_ids),
        "n_motifs": len(list(root.iter("motif"))),
        "command": model.findtext("command_line", "").strip(),
    }


def main():
    parser = argparse.ArgumentParser(
        description="Import the exact manuscript-scope MEME inputs and models."
    )
    parser.add_argument("--historical-repo", type=Path, required=True)
    parser.add_argument("--commit", default=HISTORICAL_COMMIT)
    parser.add_argument(
        "--ledger", type=Path, default=ROOT / "data/candidate_stage_ledger_975.tsv"
    )
    parser.add_argument(
        "--fasta-output", type=Path, default=ROOT / "data/candidate_sequences_975.fasta"
    )
    parser.add_argument("--models-output", type=Path, default=ROOT / "motifs/meme")
    args = parser.parse_args()

    ledger_rows, selected_ids = read_ledger(args.ledger)
    family_ids = {}
    for row in ledger_rows:
        family_ids.setdefault(row["reporting_category"], set()).add(row["gene_id"])

    combined_path = f"{HISTORICAL_ROOT}/resolved_input/all_vsc_candidates.fasta"
    fasta = filter_fasta(
        git_show(args.historical_repo, args.commit, combined_path), selected_ids
    )
    args.fasta_output.parent.mkdir(parents=True, exist_ok=True)
    args.fasta_output.write_bytes(fasta)

    args.models_output.mkdir(parents=True, exist_ok=True)
    manifest_rows = []
    for family, gene_ids in sorted(family_ids.items()):
        if len(gene_ids) == 1:
            manifest_rows.append(
                {
                    "family": family,
                    "meme_evaluable": "false",
                    "n_sequences": 1,
                    "n_motifs": 0,
                    "meme_version": "5.5.9",
                    "model_file": "",
                    "model_xml_sha256": "",
                    "historical_command": "singleton passthrough; MEME/FIMO not run",
                }
            )
            continue

        model_path = (
            f"{HISTORICAL_ROOT}/2_meme_validation/output/{family}/"
            "meme_discovery/meme.xml"
        )
        data = git_show(args.historical_repo, args.commit, model_path)
        summary = xml_summary(data, gene_ids)
        if summary["meme_version"] != "5.5.9" or summary["n_motifs"] != 10:
            raise ValueError(f"Unexpected MEME model metadata for {family}: {summary}")

        destination = args.models_output / f"{family}.xml.gz"
        destination.write_bytes(gzip.compress(data, compresslevel=9, mtime=0))
        manifest_rows.append(
            {
                "family": family,
                "meme_evaluable": "true",
                "n_sequences": summary["n_sequences"],
                "n_motifs": summary["n_motifs"],
                "meme_version": summary["meme_version"],
                "model_file": destination.relative_to(ROOT).as_posix(),
                "model_xml_sha256": hashlib.sha256(data).hexdigest(),
                "historical_command": summary["command"],
            }
        )

    manifest = args.models_output.parent / "model_manifest.tsv"
    with manifest.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(
            handle,
            delimiter="\t",
            fieldnames=manifest_rows[0],
            lineterminator="\n",
        )
        writer.writeheader()
        writer.writerows(manifest_rows)

    print(f"Wrote {len(selected_ids)} sequences to {args.fasta_output}")
    print(
        f"Wrote {len(manifest_rows) - 1} archived MEME models to {args.models_output}"
    )
    print(f"Wrote model manifest to {manifest}")


if __name__ == "__main__":
    main()
