#!/usr/bin/env python3

import argparse
import csv
import hashlib
import json
import re
import shutil
import subprocess
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
BLAST_VERSION = "2.5.0"
ATLAS_FASTA_SHA256 = "5b998b6bd3a3c59a1538a7cd3dc4e9f006e41b99855e9a5c352c52524a134e4e"
ATLAS_SEQUENCE_COUNT = 27893
EXPECTED_HSP_COUNT = 1062
EXPECTED_QUERY_COUNT = 900
EXPECTED_CANONICAL_HSP_SHA256 = "54c9dd60577c546e74ea5f0e00ac86c86639249a018ce1a436f955a828752acd"
IDENTITY_THRESHOLDS = (100, 99, 98, 95, 90, 80)
BLAST_FIELDS = [
    "query_acc",
    "subject_acc",
    "p_identity",
    "alignment_length",
    "mismatches",
    "gap_opens",
    "q_start",
    "q_end",
    "s_start",
    "s_end",
    "evalue",
    "bit_score",
]
BRIDGE_FIELDS = BLAST_FIELDS + [
    "family_slug",
    "Classification",
    "Specific tissues",
    "Tau",
    "query_length",
    "query_coverage",
]
SENSITIVITY_FIELDS = [
    "minimum_percent_identity",
    "candidate_mappings",
    "unique_atlas_genes",
    "tissue_specific_candidate_mappings",
    "root_specific_candidate_mappings",
    "trichome_specific_candidate_mappings",
    "candidate_mappings_with_any_trichome_label",
]


def sha256(path):
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def fasta_lengths(path):
    lengths = {}
    identifier = None
    length = 0
    with path.open(encoding="utf-8") as handle:
        for line in handle:
            line = line.strip()
            if line.startswith(">"):
                if identifier is not None:
                    lengths[identifier] = length
                identifier = line[1:].split(maxsplit=1)[0]
                if identifier in lengths:
                    raise ValueError(f"Duplicate FASTA identifier: {identifier}")
                length = 0
            else:
                length += len(line)
    if identifier is not None:
        lengths[identifier] = length
    return lengths


def read_tsv(path):
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def read_blast(path):
    rows = []
    with path.open(newline="", encoding="utf-8") as handle:
        for values in csv.reader(handle, delimiter="\t"):
            if len(values) != len(BLAST_FIELDS):
                raise ValueError(f"Unexpected BLAST row with {len(values)} fields")
            rows.append(dict(zip(BLAST_FIELDS, values)))
    return rows


def canonical_blast_sha256(rows):
    lines = [
        "\t".join(row[field] for field in BLAST_FIELDS) + "\n" for row in rows
    ]
    return hashlib.sha256("".join(sorted(lines)).encode()).hexdigest()


def best_hits(rows):
    selected = {}
    for row in sorted(rows, key=lambda value: float(value["bit_score"]), reverse=True):
        selected.setdefault(row["query_acc"], row)
    return list(selected.values())


def eligible_hits(rows, lengths, minimum_identity):
    return [
        row
        for row in rows
        if float(row["p_identity"]) >= minimum_identity
        and int(row["alignment_length"]) / lengths[row["query_acc"]] >= 0.8
    ]


def write_tsv(path, fieldnames, rows):
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=fieldnames, delimiter="\t", lineterminator="\n"
        )
        writer.writeheader()
        writer.writerows(rows)


def require_tools():
    tools = {name: shutil.which(name) for name in ("makeblastdb", "blastp")}
    missing = [name for name, path in tools.items() if not path]
    if missing:
        raise SystemExit(
            f"{', '.join(missing)} not found; create and activate environment.yml"
        )
    output = subprocess.run(
        [tools["blastp"], "-version"], check=True, capture_output=True, text=True
    ).stdout
    match = re.search(r"blastp:\s+(\d+\.\d+\.\d+)\+?", output)
    found = match.group(1) if match else "unknown"
    if found != BLAST_VERSION:
        raise SystemExit(f"BLAST {BLAST_VERSION} is required; found {found}")
    return tools


def validate_inputs(atlas, candidates, catalog, metadata):
    atlas_digest = sha256(atlas)
    if atlas_digest != ATLAS_FASTA_SHA256:
        raise SystemExit(
            f"Atlas FASTA checksum mismatch: expected {ATLAS_FASTA_SHA256}, "
            f"found {atlas_digest}"
        )
    atlas_count = len(fasta_lengths(atlas))
    if atlas_count != ATLAS_SEQUENCE_COUNT:
        raise SystemExit(
            f"Atlas FASTA record mismatch: expected {ATLAS_SEQUENCE_COUNT}, "
            f"found {atlas_count}"
        )

    lengths = fasta_lengths(candidates)
    families = {row["gene_id"]: row["family"].lower() for row in read_tsv(catalog)}
    if set(lengths) != set(families):
        raise SystemExit("Candidate FASTA and catalog identifiers differ")

    metadata_rows = read_tsv(metadata)
    metadata_by_gene = {row["GeneID"]: row for row in metadata_rows}
    if len(metadata_by_gene) != len(metadata_rows):
        raise SystemExit("Atlas metadata has duplicate GeneID values")
    return atlas_digest, lengths, families, metadata_by_gene


def run_blast(tools, atlas, candidates, output_root):
    database = output_root / "atlas-db" / "atlas"
    database.parent.mkdir(parents=True, exist_ok=True)
    make_command = [
        tools["makeblastdb"],
        "-in",
        str(atlas),
        "-dbtype",
        "prot",
        "-out",
        str(database),
    ]
    print("Running", " ".join(make_command))
    subprocess.run(make_command, check=True, capture_output=True)

    raw_output = output_root / "blast.tsv"
    blast_command = [
        tools["blastp"],
        "-query",
        str(candidates),
        "-db",
        str(database),
        "-evalue",
        "0.001",
        "-outfmt",
        "6",
        "-matrix",
        "BLOSUM62",
        "-qcov_hsp_perc",
        "80",
        "-max_target_seqs",
        "1",
        "-out",
        str(raw_output),
    ]
    print("Running", " ".join(blast_command))
    subprocess.run(blast_command, check=True)
    return raw_output, make_command, blast_command


def bridge_rows(hits, lengths, families, metadata):
    rows = []
    for hit in eligible_hits(hits, lengths, 100):
        gene_metadata = metadata[hit["subject_acc"]]
        row = dict(hit)
        row["p_identity"] = str(float(hit["p_identity"]))
        row["family_slug"] = families[hit["query_acc"]]
        row["Classification"] = gene_metadata["Classification"]
        row["Specific tissues"] = gene_metadata["Specific tissues"]
        row["Tau"] = (
            str(float(gene_metadata["Tau"])) if gene_metadata["Tau"] else ""
        )
        row["query_length"] = str(lengths[hit["query_acc"]])
        row["query_coverage"] = str(
            int(hit["alignment_length"]) / lengths[hit["query_acc"]]
        )
        rows.append(row)
    return sorted(rows, key=lambda row: row["query_acc"])


def sensitivity_rows(hits, lengths, metadata):
    rows = []
    for threshold in IDENTITY_THRESHOLDS:
        selected = eligible_hits(hits, lengths, threshold)
        joined = [(row, metadata[row["subject_acc"]]) for row in selected]
        rows.append(
            {
                "minimum_percent_identity": threshold,
                "candidate_mappings": len(selected),
                "unique_atlas_genes": len(
                    {row["subject_acc"] for row in selected}
                ),
                "tissue_specific_candidate_mappings": sum(
                    gene["Classification"] == "Tissue-Specific"
                    for _, gene in joined
                ),
                "root_specific_candidate_mappings": sum(
                    gene["Classification"] == "Tissue-Specific"
                    and "Root" in gene["Specific tissues"]
                    for _, gene in joined
                ),
                "trichome_specific_candidate_mappings": sum(
                    gene["Classification"] == "Tissue-Specific"
                    and "Trichome" in gene["Specific tissues"]
                    for _, gene in joined
                ),
                "candidate_mappings_with_any_trichome_label": sum(
                    "Trichome" in gene["Specific tissues"] for _, gene in joined
                ),
            }
        )
    return rows


def main():
    parser = argparse.ArgumentParser(
        description="Reproduce Cannabis Expression Atlas mappings for the 941 candidates."
    )
    parser.add_argument("--atlas-proteins", type=Path, required=True)
    parser.add_argument(
        "--candidates",
        type=Path,
        default=ROOT / "data/candidate_sequences_941.fasta",
    )
    parser.add_argument(
        "--catalog", type=Path, default=ROOT / "data/candidate_catalog_941.tsv"
    )
    parser.add_argument(
        "--metadata",
        type=Path,
        default=ROOT / "inputs/expression/atlas_gene_metadata_423.tsv",
    )
    parser.add_argument(
        "--output-root", type=Path, default=ROOT / "build/expression"
    )
    args = parser.parse_args()

    tools = require_tools()
    atlas_digest, lengths, families, metadata = validate_inputs(
        args.atlas_proteins, args.candidates, args.catalog, args.metadata
    )
    args.output_root.mkdir(parents=True, exist_ok=True)
    raw_output, make_command, blast_command = run_blast(
        tools, args.atlas_proteins, args.candidates, args.output_root
    )
    raw_rows = read_blast(raw_output)
    hits = best_hits(raw_rows)
    canonical_hsp_digest = canonical_blast_sha256(raw_rows)
    if (
        len(raw_rows) != EXPECTED_HSP_COUNT
        or len(hits) != EXPECTED_QUERY_COUNT
        or canonical_hsp_digest != EXPECTED_CANONICAL_HSP_SHA256
    ):
        raise SystemExit(
            "Raw BLAST output differs from the historical 941-query projection: "
            f"hsps={len(raw_rows)}, queries={len(hits)}, "
            f"canonical_sha256={canonical_hsp_digest}"
        )
    required_metadata = {
        row["subject_acc"] for row in eligible_hits(hits, lengths, 80)
    }
    if set(metadata) != required_metadata:
        missing = sorted(required_metadata - set(metadata))
        extra = sorted(set(metadata) - required_metadata)
        raise SystemExit(
            f"Atlas metadata boundary mismatch: missing={missing}, extra={extra}"
        )

    bridge = bridge_rows(hits, lengths, families, metadata)
    sensitivity = sensitivity_rows(hits, lengths, metadata)
    bridge_output = args.output_root / "expression_bridge_exact100_qcov80.tsv"
    sensitivity_output = args.output_root / "expression_identity_sensitivity.tsv"
    write_tsv(bridge_output, BRIDGE_FIELDS, bridge)
    write_tsv(sensitivity_output, SENSITIVITY_FIELDS, sensitivity)

    summary = {
        "atlas_fasta_sha256": atlas_digest,
        "atlas_sequences": ATLAS_SEQUENCE_COUNT,
        "blast_version": BLAST_VERSION,
        "makeblastdb_command": make_command,
        "blastp_command": blast_command,
        "raw_hsps": len(raw_rows),
        "canonical_hsp_sha256": canonical_hsp_digest,
        "queries_with_hits": len(hits),
        "exact100_candidate_mappings": len(bridge),
        "exact100_unique_atlas_genes": len(
            {row["subject_acc"] for row in bridge}
        ),
    }
    (args.output_root / "summary.json").write_text(
        json.dumps(summary, indent=2) + "\n", encoding="utf-8"
    )
    print(
        f"Wrote {len(bridge)} exact-identity mappings and "
        f"{len(sensitivity)} sensitivity rows to {args.output_root}"
    )


if __name__ == "__main__":
    main()
