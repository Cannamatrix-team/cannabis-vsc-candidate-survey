#!/usr/bin/env python3

import argparse
import csv
import json
from collections import Counter
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
MEME_PASS_STATUSES = {"confirmed", "partial"}


def read_fasta(path):
    records = []
    header = None
    sequence = []
    with path.open(encoding="utf-8") as handle:
        for line in handle:
            line = line.rstrip("\n")
            if line.startswith(">"):
                if header is not None:
                    records.append((header, "".join(sequence)))
                header = line[1:]
                sequence = []
            else:
                sequence.append(line.strip())
    if header is not None:
        records.append((header, "".join(sequence)))
    return records


def metadata(header):
    values = {"gene_id": header.split()[0]}
    for token in header.split()[1:]:
        if "=" in token:
            key, value = token.split("=", 1)
            values[key] = value
    return values


def read_tsv(path):
    with path.open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    indexed = {row["gene_id"]: row for row in rows}
    if len(indexed) != len(rows):
        raise ValueError(f"Duplicate gene identifiers in {path}")
    return indexed


def as_bool(value):
    return str(value).lower() == "true"


def write_fasta(records, path):
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8") as handle:
        for header, sequence in records:
            handle.write(f">{header}\n")
            for start in range(0, len(sequence), 60):
                handle.write(sequence[start : start + 60] + "\n")


def gate_result(header, motif, domain):
    values = metadata(header)
    gene_id = values["gene_id"]
    biological_pass = motif["meme_status"] in MEME_PASS_STATUSES
    singleton = as_bool(motif["singleton_passthrough"])
    operational_pass = biological_pass or singleton
    if domain:
        pfam_status = domain["domain_confirmation_status"]
        pfam_pass = pfam_status == "complete"
        domain_record = "true"
    else:
        pfam_status = ""
        pfam_pass = True
        domain_record = "false"

    final_pass = operational_pass and pfam_pass
    if not operational_pass:
        outcome = "removed_without_operational_motif_support"
    elif not domain:
        outcome = "retained_by_default_no_configured_domain_rule"
    elif not pfam_pass:
        outcome = "removed_incomplete_expected_domain_architecture"
    else:
        outcome = "retained"

    return {
        "gene_id": gene_id,
        "family": values["family"].lower(),
        "meme_biological_pass": biological_pass,
        "meme_operational_pass": operational_pass,
        "meme_pass": biological_pass,
        "expected_domain_record_present": domain_record,
        "pfam_pass": pfam_pass,
        "pfam_status": pfam_status,
        "pfam_gate_active": True,
        "final_pass": final_pass,
        "final_outcome": outcome,
    }


def main():
    parser = argparse.ArgumentParser(
        description="Apply the manuscript-active motif and expected-domain gate."
    )
    parser.add_argument(
        "--fasta", type=Path, default=ROOT / "data/candidate_sequences_975.fasta"
    )
    parser.add_argument(
        "--motifs", type=Path, default=ROOT / "build/motifs/motif_screen.tsv"
    )
    parser.add_argument(
        "--domains",
        type=Path,
        default=ROOT / "build/domains/domain_validation.tsv",
    )
    parser.add_argument("--output-root", type=Path, default=ROOT / "build/gate")
    args = parser.parse_args()

    records = read_fasta(args.fasta)
    gene_ids = [metadata(header)["gene_id"] for header, _ in records]
    if len(gene_ids) != len(set(gene_ids)):
        raise ValueError("Candidate FASTA has duplicate identifiers")
    motifs = read_tsv(args.motifs)
    domains = read_tsv(args.domains)
    if set(motifs) != set(gene_ids):
        raise ValueError("Motif and candidate FASTA identifiers differ")
    if not set(domains) < set(gene_ids):
        raise ValueError(
            "Expected-domain identifiers are not a proper candidate subset"
        )

    results = [
        gate_result(
            header,
            motifs[metadata(header)["gene_id"]],
            domains.get(metadata(header)["gene_id"]),
        )
        for header, _ in records
    ]
    passing_ids = {row["gene_id"] for row in results if row["final_pass"]}
    retained = [
        record for record in records if metadata(record[0])["gene_id"] in passing_ids
    ]

    args.output_root.mkdir(parents=True, exist_ok=True)
    summary_path = args.output_root / "validation_summary.tsv"
    with summary_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(
            handle, delimiter="\t", fieldnames=results[0], lineterminator="\n"
        )
        writer.writeheader()
        writer.writerows(results)
    fasta_path = args.output_root / "validated_candidates.fasta"
    write_fasta(retained, fasta_path)

    summary = {
        "sequences": len(results),
        "retained": len(retained),
        "final_outcome": Counter(row["final_outcome"] for row in results),
    }
    (args.output_root / "summary.json").write_text(
        json.dumps(summary, indent=2) + "\n", encoding="utf-8"
    )
    print(f"Wrote {len(results)} gate records to {summary_path}")
    print(f"Wrote {len(retained)} retained sequences to {fasta_path}")


if __name__ == "__main__":
    main()
