#!/usr/bin/env python3

import argparse
import csv
import gzip
import json
import shutil
import subprocess
import xml.etree.ElementTree as ET
from collections import Counter, defaultdict
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
MEME_VERSION = "5.5.9"
FIMO_PVALUE = 1e-4
STRONG_EVALUE = 1e-3
WEAK_EVALUE = 1.0


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


def header_metadata(header):
    tokens = header.split()
    metadata = {"gene_id": tokens[0]}
    for token in tokens[1:]:
        if "=" in token:
            key, value = token.split("=", 1)
            metadata[key] = value
    return metadata


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


def tool_version(command, flag):
    result = subprocess.run([command, flag], check=True, capture_output=True, text=True)
    return (result.stdout + result.stderr).strip().splitlines()[-1]


def require_tools(discover):
    fimo = shutil.which("fimo")
    if not fimo:
        raise SystemExit("fimo not found; create and activate environment.yml")
    if tool_version(fimo, "--version") != MEME_VERSION:
        raise SystemExit(f"FIMO {MEME_VERSION} is required")

    meme = None
    if discover:
        meme = shutil.which("meme")
        if not meme:
            raise SystemExit("meme not found; create and activate environment.yml")
        if tool_version(meme, "-version") != MEME_VERSION:
            raise SystemExit(f"MEME {MEME_VERSION} is required")
    return meme, fimo


def parse_model(path):
    root = ET.parse(path).getroot()
    motifs = []
    for node in root.iter("motif"):
        evalue = float(node.get("e_value", 1.0))
        if evalue <= STRONG_EVALUE:
            quality = "strong"
        elif evalue <= WEAK_EVALUE:
            quality = "weak"
        else:
            quality = "very_weak"
        motifs.append(
            {
                "id": node.get("id", ""),
                "alt": node.get("alt", ""),
                "name": node.get("name", ""),
                "evalue": evalue,
                "quality": quality,
            }
        )
    return motifs


def parse_fimo(path):
    hits = defaultdict(list)
    with path.open(newline="", encoding="utf-8") as handle:
        rows = csv.DictReader(
            (line for line in handle if not line.startswith("#")), delimiter="\t"
        )
        for row in rows:
            if not row.get("sequence_name") or not row.get("p-value"):
                continue
            hits[row["sequence_name"].split()[0]].append(
                {
                    "motif_id": row.get("motif_id", ""),
                    "motif_alt_id": row.get("motif_alt_id", ""),
                    "pvalue": float(row["p-value"]),
                }
            )
    return hits


def classify(gene_id, hits, motifs, rule):
    aliases = {}
    quality_by_id = {}
    for motif in motifs:
        quality_by_id[motif["id"]] = motif["quality"]
        for alias in (motif["id"], motif["alt"], motif["name"]):
            if alias:
                aliases[alias] = motif["id"]

    significant = [hit for hit in hits if hit["pvalue"] < FIMO_PVALUE]
    found = set()
    hit_quality = Counter()
    for hit in significant:
        motif_id = aliases.get(hit["motif_alt_id"]) or aliases.get(hit["motif_id"])
        if not motif_id:
            raise ValueError(f"Unrecognized motif in FIMO output for {gene_id}: {hit}")
        found.add(motif_id)
        hit_quality[quality_by_id[motif_id]] += 1

    found_by_quality = {
        quality: sorted(
            motif_id for motif_id in found if quality_by_id[motif_id] == quality
        )
        for quality in ("strong", "weak", "very_weak")
    }
    expected_by_quality = Counter(motif["quality"] for motif in motifs)
    coverage = len(found) / len(motifs)
    strong_coverage = len(found_by_quality["strong"]) / max(
        expected_by_quality["strong"], 1
    )
    min_confirmed = int(rule["min_motifs_confirmed"])
    min_partial = int(rule["min_motifs_partial"])
    confirmed_coverage = float(rule["min_coverage_confirmed"])
    partial_coverage = float(rule["min_coverage_partial"])

    if len(found) >= min_confirmed and coverage >= confirmed_coverage:
        status = "confirmed"
    elif len(found) >= min_partial and coverage >= partial_coverage:
        status = "partial"
    else:
        status = "absent"

    n_strong = len(found_by_quality["strong"])
    n_weak = len(found_by_quality["weak"])
    n_very_weak = len(found_by_quality["very_weak"])
    if n_strong >= min_confirmed and strong_coverage >= confirmed_coverage:
        quality_status = "confirmed"
    elif (n_strong >= min_partial and strong_coverage >= partial_coverage) or n_weak:
        quality_status = "partial"
    else:
        quality_status = "absent"

    if n_strong and not n_weak and not n_very_weak:
        support_quality = "strong_only"
    elif n_strong:
        support_quality = "strong_plus_weaker"
    elif n_weak:
        support_quality = "weak_only"
    elif n_very_weak:
        support_quality = "very_weak_only"
    else:
        support_quality = "no_significant_motifs"

    return {
        "gene_id": gene_id,
        "meme_status": status,
        "meme_status_quality_aware": quality_status,
        "motif_support_quality": support_quality,
        "meme_evaluable": "True",
        "singleton_passthrough": "False",
        "n_motifs_found": len(found),
        "n_motifs_expected": len(motifs),
        "motif_ids": ",".join(sorted(found)),
        "motif_best_evalue": min((hit["pvalue"] for hit in significant), default=1.0),
        "motif_coverage": round(coverage, 4),
        "fimo_hits": len(significant),
    }


def singleton_result(gene_id):
    return {
        "gene_id": gene_id,
        "meme_status": "absent",
        "meme_status_quality_aware": "absent",
        "motif_support_quality": "single_sequence",
        "meme_evaluable": "False",
        "singleton_passthrough": "True",
        "n_motifs_found": 0,
        "n_motifs_expected": 0,
        "motif_ids": "",
        "motif_best_evalue": "",
        "motif_coverage": 0.0,
        "fimo_hits": 0,
    }


def run_meme(meme, fasta, output, threads):
    size = sum(len(sequence) for _, sequence in read_fasta(fasta))
    command = [
        meme,
        str(fasta),
        "-protein",
        "-oc",
        str(output),
        "-nmotifs",
        "10",
        "-minw",
        "6",
        "-maxw",
        "50",
        "-mod",
        "zoops",
        "-nostatus",
        "-maxsize",
        str(max(100000, size + 10000)),
        "-p",
        str(threads),
    ]
    print("Running", " ".join(command))
    subprocess.run(command, check=True)
    return output / "meme.xml"


def main():
    parser = argparse.ArgumentParser(
        description="Reproduce the manuscript-scope MEME/FIMO motif screen."
    )
    parser.add_argument(
        "--fasta", type=Path, default=ROOT / "data/candidate_sequences_975.fasta"
    )
    parser.add_argument("--rules", type=Path, default=ROOT / "motifs/family_rules.tsv")
    parser.add_argument("--models", type=Path, default=ROOT / "motifs/meme")
    parser.add_argument("--output-root", type=Path, default=ROOT / "build/motifs")
    parser.add_argument(
        "--discover",
        action="store_true",
        help="rerun de novo MEME discovery instead of using the archived models",
    )
    parser.add_argument("--meme-threads", type=int, default=4)
    parser.add_argument("--family", action="append", default=[])
    args = parser.parse_args()

    meme, fimo = require_tools(args.discover)
    rules = read_rules(args.rules)
    records = read_fasta(args.fasta)
    by_family = defaultdict(list)
    for header, sequence in records:
        metadata = header_metadata(header)
        family = metadata.get("family", "").lower()
        if family not in rules:
            raise ValueError(f"No motif rule for {family}: {metadata['gene_id']}")
        by_family[family].append((header, sequence))

    selected = set(args.family) if args.family else set(by_family)
    unknown = selected - set(by_family)
    if unknown:
        raise ValueError(f"Unknown families: {sorted(unknown)}")

    args.output_root.mkdir(parents=True, exist_ok=True)
    results = []
    for family in sorted(selected):
        family_records = by_family[family]
        family_dir = args.output_root / family
        fasta = family_dir / f"{family}_candidates.fasta"
        write_fasta(family_records, fasta)
        print(f"{family}: {len(family_records)} sequences")

        if len(family_records) == 1:
            gene_id = header_metadata(family_records[0][0])["gene_id"]
            results.append(singleton_result(gene_id))
            continue

        if args.discover:
            model = run_meme(
                meme, fasta, family_dir / "meme_discovery", args.meme_threads
            )
        else:
            archived = args.models / f"{family}.xml.gz"
            model = family_dir / "meme.xml"
            model.write_bytes(gzip.decompress(archived.read_bytes()))

        fimo_dir = family_dir / "fimo_scan"
        command = [
            fimo,
            "--oc",
            str(fimo_dir),
            "--thresh",
            str(FIMO_PVALUE),
            "--verbosity",
            "1",
            str(model),
            str(fasta),
        ]
        print("Running", " ".join(command))
        subprocess.run(command, check=True)

        motifs = parse_model(model)
        if len(motifs) != 10:
            raise ValueError(f"Expected 10 motifs for {family}, found {len(motifs)}")
        hits = parse_fimo(fimo_dir / "fimo.tsv")
        for header, _ in family_records:
            gene_id = header_metadata(header)["gene_id"]
            results.append(
                classify(gene_id, hits.get(gene_id, []), motifs, rules[family])
            )

    order = {
        header_metadata(header)["gene_id"]: i for i, (header, _) in enumerate(records)
    }
    results.sort(key=lambda row: order[row["gene_id"]])
    output = args.output_root / "motif_screen.tsv"
    with output.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(
            handle, delimiter="\t", fieldnames=results[0], lineterminator="\n"
        )
        writer.writeheader()
        writer.writerows(results)

    summary = {
        "mode": "fresh_discovery" if args.discover else "archived_model_rescan",
        "meme_suite_version": MEME_VERSION,
        "families": len(selected),
        "sequences": len(results),
        "meme_status": Counter(row["meme_status"] for row in results),
        "meme_evaluable": Counter(row["meme_evaluable"] for row in results),
        "singleton_passthrough": Counter(
            row["singleton_passthrough"] for row in results
        ),
    }
    (args.output_root / "summary.json").write_text(
        json.dumps(summary, indent=2) + "\n", encoding="utf-8"
    )
    print(f"Wrote {len(results)} motif records to {output}")


if __name__ == "__main__":
    main()
