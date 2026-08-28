#!/usr/bin/env python3

import argparse
import csv
from collections import Counter
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
FIELDS = [
    "meme_biological_pass",
    "meme_operational_pass",
    "meme_pass",
    "expected_domain_record_present",
    "pfam_pass",
    "pfam_status",
    "pfam_gate_active",
    "final_pass",
    "final_outcome",
]


def read_tsv(path):
    with path.open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    return {row["gene_id"]: row for row in rows}


def main():
    parser = argparse.ArgumentParser(
        description="Compare regenerated validation-gate outputs with the deposit."
    )
    parser.add_argument(
        "--generated-summary",
        type=Path,
        default=ROOT / "build/gate/validation_summary.tsv",
    )
    parser.add_argument(
        "--generated-fasta",
        type=Path,
        default=ROOT / "build/gate/validated_candidates.fasta",
    )
    parser.add_argument(
        "--ledger", type=Path, default=ROOT / "data/candidate_stage_ledger_975.tsv"
    )
    parser.add_argument(
        "--deposited-fasta",
        type=Path,
        default=ROOT / "data/candidate_sequences_941.fasta",
    )
    args = parser.parse_args()

    expected = read_tsv(args.ledger)
    observed = read_tsv(args.generated_summary)
    if set(observed) != set(expected):
        missing = sorted(set(expected) - set(observed))
        extra = sorted(set(observed) - set(expected))
        raise SystemExit(f"Gate identifier mismatch: missing={missing}, extra={extra}")

    mismatches = []
    for gene_id, expected_row in expected.items():
        observed_row = observed[gene_id]
        for field in FIELDS:
            if observed_row[field] != expected_row[field]:
                mismatches.append(
                    (gene_id, field, expected_row[field], observed_row[field])
                )
    if mismatches:
        for mismatch in mismatches[:20]:
            print("\t".join(mismatch))
        raise SystemExit(f"Gate comparison failed with {len(mismatches)} mismatches")
    if args.generated_fasta.read_bytes() != args.deposited_fasta.read_bytes():
        raise SystemExit("Generated retained FASTA differs from the deposited FASTA")

    outcomes = Counter(row["final_outcome"] for row in observed.values())
    print(
        f"Validation gate matches: {len(observed)} proteins; "
        f"outcomes={dict(outcomes)}; retained_FASTA=941"
    )


if __name__ == "__main__":
    main()
