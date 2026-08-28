#!/usr/bin/env python3

import argparse
import csv
import json
import re
import shutil
import subprocess
from collections import Counter, defaultdict
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
HMMER_VERSION = "3.4"
SEQUENCE_EVALUE = 1e-5
DOMAIN_EVALUE = 1e-3


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


def write_fasta(records, path):
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8") as handle:
        for header, sequence in records:
            handle.write(f">{header}\n")
            for start in range(0, len(sequence), 80):
                handle.write(sequence[start : start + 80] + "\n")


def read_rules(path):
    with path.open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    return {row["family"]: row for row in rows}


def require_tools():
    tools = {name: shutil.which(name) for name in ("hmmscan", "hmmpress")}
    missing = [name for name, path in tools.items() if not path]
    if missing:
        raise SystemExit(
            f"{', '.join(missing)} not found; create and activate environment.yml"
        )
    help_text = subprocess.run(
        [tools["hmmscan"], "-h"], check=True, capture_output=True, text=True
    ).stdout
    match = re.search(r"HMMER\s+(\S+)", help_text)
    if not match or match.group(1) != HMMER_VERSION:
        found = match.group(1) if match else "unknown"
        raise SystemExit(f"HMMER {HMMER_VERSION} is required; found {found}")
    return tools


def profile_path(root, pfam_id):
    matches = list(root.glob(f"{pfam_id}.*.hmm"))
    if len(matches) != 1:
        raise ValueError(
            f"Expected one packaged profile for {pfam_id}, found {matches}"
        )
    return matches[0]


def build_database(hmmpress, profile_root, pfam_ids, destination):
    destination.parent.mkdir(parents=True, exist_ok=True)
    with destination.open("wb") as handle:
        for pfam_id in pfam_ids:
            data = profile_path(profile_root, pfam_id).read_bytes()
            handle.write(data)
            if not data.endswith(b"\n"):
                handle.write(b"\n")
    subprocess.run([hmmpress, "-f", str(destination)], check=True, capture_output=True)


def run_scan(hmmscan, database, fasta, output, threads):
    command = [
        hmmscan,
        "--domtblout",
        str(output),
        "--noali",
        "-E",
        str(SEQUENCE_EVALUE),
        "--domE",
        str(DOMAIN_EVALUE),
        "--cpu",
        str(threads),
        str(database),
        str(fasta),
    ]
    print("Running", " ".join(command))
    subprocess.run(command, check=True, capture_output=True)


def parse_hits(path, pfam_ids):
    expected = set(pfam_ids)
    hits = defaultdict(list)
    with path.open(encoding="utf-8") as handle:
        for line in handle:
            if line.startswith("#"):
                continue
            fields = line.split()
            if len(fields) < 22:
                continue
            target_name, accession, gene_id = fields[0], fields[1], fields[3]
            candidates = (target_name.split(".", 1)[0], accession.split(".", 1)[0])
            pfam_id = next((value for value in candidates if value in expected), None)
            if pfam_id:
                hits[gene_id].append(
                    {
                        "pfam_id": pfam_id,
                        "dom_evalue": float(fields[12]),
                        "dom_score": float(fields[13]),
                    }
                )
    return hits


def domain_features(gene_id, family, hits, rule):
    pfam_ids = rule["pfam_ids"].split(",")
    required_ids = rule["required_pfam_ids"].split(",")
    best_hits = {}
    for pfam_id in pfam_ids:
        candidates = [hit for hit in hits if hit["pfam_id"] == pfam_id]
        if candidates:
            best_hits[pfam_id] = max(candidates, key=lambda hit: hit["dom_score"])

    present_ids = [pfam_id for pfam_id in pfam_ids if pfam_id in best_hits]
    n_required = sum(pfam_id in best_hits for pfam_id in required_ids)
    if len(present_ids) == len(pfam_ids):
        architecture = "both" if len(pfam_ids) == 2 else "all"
    elif not present_ids:
        architecture = "none"
    else:
        architecture = "+".join(present_ids)

    if rule["confirmation_rule"] == "any_required":
        status = "complete" if n_required else "absent"
    elif not n_required:
        status = "absent"
    elif n_required == len(required_ids):
        status = "complete"
    else:
        status = "partial"

    selected = list(best_hits.values())
    return {
        "gene_id": gene_id,
        "family": family,
        "pfam_n_domains": len(present_ids),
        "pfam_domain_architecture": architecture,
        "pfam_best_evalue": min((hit["dom_evalue"] for hit in selected), default=""),
        "pfam_total_score": round(sum(hit["dom_score"] for hit in selected), 1),
        "pfam_confirmation_rule": rule["confirmation_rule"],
        "pfam_required_ids": rule["required_pfam_ids"],
        "pfam_required_n_domains": len(required_ids),
        "pfam_required_n_present": n_required,
        "domain_confirmation_status": status,
        "domain_validated": status == "complete",
    }


def main():
    parser = argparse.ArgumentParser(
        description="Reproduce the manuscript-scope expected-domain validation."
    )
    parser.add_argument(
        "--fasta", type=Path, default=ROOT / "data/candidate_sequences_975.fasta"
    )
    parser.add_argument("--rules", type=Path, default=ROOT / "domains/family_rules.tsv")
    parser.add_argument(
        "--profiles", type=Path, default=ROOT / "inputs/search_queries/pfam"
    )
    parser.add_argument("--output-root", type=Path, default=ROOT / "build/domains")
    parser.add_argument("--threads", type=int, default=4)
    parser.add_argument("--family", action="append", default=[])
    args = parser.parse_args()

    tools = require_tools()
    rules = read_rules(args.rules)
    records = read_fasta(args.fasta)
    by_family = defaultdict(list)
    for header, sequence in records:
        values = metadata(header)
        family = values.get("family", "").lower()
        if family not in rules:
            raise ValueError(
                f"No expected-domain rule for {family}: {values['gene_id']}"
            )
        by_family[family].append((header, sequence))

    if set(rules) != set(by_family):
        raise ValueError(
            f"Rule and FASTA families differ: rules={sorted(rules)}, "
            f"FASTA={sorted(by_family)}"
        )
    selected = set(args.family) if args.family else set(by_family)
    unknown = selected - set(by_family)
    if unknown:
        raise ValueError(f"Unknown families: {sorted(unknown)}")

    args.output_root.mkdir(parents=True, exist_ok=True)
    results = []
    for family in sorted(selected):
        rule = rules[family]
        family_records = by_family[family]
        if rule["record_expected"] != "true":
            print(f"{family}: no historical expected-domain record")
            continue

        family_root = args.output_root / family
        fasta = family_root / f"{family}_candidates.fasta"
        database = family_root / "pfam_combined.hmm"
        domtbl = family_root / "hmmscan_domtbl.out"
        write_fasta(family_records, fasta)
        pfam_ids = rule["pfam_ids"].split(",")
        build_database(tools["hmmpress"], args.profiles, pfam_ids, database)
        run_scan(tools["hmmscan"], database, fasta, domtbl, args.threads)
        hits = parse_hits(domtbl, pfam_ids)
        for header, _ in family_records:
            gene_id = metadata(header)["gene_id"]
            results.append(domain_features(gene_id, family, hits[gene_id], rule))

    order = {
        metadata(header)["gene_id"]: index for index, (header, _) in enumerate(records)
    }
    results.sort(key=lambda row: order[row["gene_id"]])
    output = args.output_root / "domain_validation.tsv"
    with output.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(
            handle, delimiter="\t", fieldnames=results[0], lineterminator="\n"
        )
        writer.writeheader()
        writer.writerows(results)

    summary = {
        "hmmer_version": HMMER_VERSION,
        "families_scanned": sum(
            rules[family]["record_expected"] == "true" for family in selected
        ),
        "sequences": len(results),
        "domain_confirmation_status": Counter(
            row["domain_confirmation_status"] for row in results
        ),
    }
    (args.output_root / "summary.json").write_text(
        json.dumps(summary, indent=2) + "\n", encoding="utf-8"
    )
    print(f"Wrote {len(results)} expected-domain records to {output}")


if __name__ == "__main__":
    main()
