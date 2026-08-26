#!/usr/bin/env python3

import argparse
import csv
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
PROTEOME_SUBDIR = "GMO.v1.primary_high_confidence.proteins"
TARGET_ORDER = [
    "ADH",
    "Alliinase",
    "GGT",
    "GST",
    "CSCBL-like",
    "Lachrymatory-factor synthase",
    "LOX",
    "MGL-Cannabis",
    "MGL-Durio",
    "MGL-Artocarpus",
    "Myrosinase",
    "AOC",
    "AOS/HPL-motivated P450 search",
    "CGS",
    "CBL",
    "FMO",
    "GSH1",
    "GSH2",
    "OASTL",
    "Rhodanese",
    "SAMS",
    "SULT",
    "UGT",
]
FIELDS = [
    "gene_id",
    "family",
    "biochemical_context",
    "target_bin",
    "length",
    "score",
    "evalue",
    "source",
    "sequence_scope",
    "discovery_round",
    "matched_protein_id",
    "pfam_hit",
    "full_length_reference_hit",
    "passing_pfam_profiles",
    "passing_full_length_reference_queries",
    "search_route",
]


def read_tsv(path):
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def fasta_lengths(path):
    lengths = {}
    identifier = None
    length = 0
    with path.open(encoding="utf-8") as handle:
        for line in handle:
            if line.startswith(">"):
                if identifier is not None:
                    lengths[identifier] = length
                identifier = line[1:].split(maxsplit=1)[0]
                length = 0
            else:
                length += len(line.strip())
    if identifier is not None:
        lengths[identifier] = length
    return lengths


def passing_hits(path, score_min, evalue_max):
    hits = {}
    with path.open(encoding="utf-8") as handle:
        for line in handle:
            if line.startswith("#"):
                continue
            fields = line.split()
            if not fields:
                continue
            evalue = float(fields[4])
            score = float(fields[5])
            if score >= score_min and evalue <= evalue_max:
                gene_id = fields[0]
                previous = hits.get(gene_id)
                if previous is None or score > previous[0]:
                    hits[gene_id] = (score, evalue)
    return hits


def source_tables(panel_row, hmmer_root):
    family_root = hmmer_root
    if panel_row["discovery_round"] == "round2":
        family_root /= "round2_results"
    family_root /= panel_row["pipeline_family"]
    table_root = family_root / PROTEOME_SUBDIR

    profiles = [item for item in panel_row["pfam_profiles"].split(";") if item]
    references = [
        item for item in panel_row["full_length_reference_queries"].split(";") if item
    ]
    profile_tables = [
        (profile, table_root / f"{profile.split('.', 1)[0]}.tbl")
        for profile in profiles
    ]
    reference_tables = [
        (reference, table_root / f"{reference.split('/')[1]}.fasta.tbl")
        for reference in references
    ]
    return profile_tables, reference_tables


def build_rows(panel, hmmer_root, lengths, score_min, evalue_max):
    panel_by_target = {row["target_bin"]: row for row in panel}
    rows = []
    for target in TARGET_ORDER:
        panel_row = panel_by_target[target]
        profile_tables, reference_tables = source_tables(panel_row, hmmer_root)
        profiles_by_gene = {}
        references_by_gene = {}
        best_by_gene = {}

        for label, path in profile_tables:
            for gene_id, result in passing_hits(path, score_min, evalue_max).items():
                profiles_by_gene.setdefault(gene_id, []).append(label)
                if gene_id not in best_by_gene or result[0] > best_by_gene[gene_id][0]:
                    best_by_gene[gene_id] = result
        for label, path in reference_tables:
            for gene_id, result in passing_hits(path, score_min, evalue_max).items():
                references_by_gene.setdefault(gene_id, []).append(label)
                if gene_id not in best_by_gene or result[0] > best_by_gene[gene_id][0]:
                    best_by_gene[gene_id] = result

        ordered_hits = sorted(best_by_gene.items(), key=lambda item: (-item[1][0], item[0]))
        for gene_id, (score, evalue) in ordered_hits:
            profiles = profiles_by_gene.get(gene_id, [])
            references = references_by_gene.get(gene_id, [])
            route = "pfam_and_reference" if profiles and references else "pfam_only"
            if references and not profiles:
                route = "reference_only"
            rows.append(
                {
                    "gene_id": gene_id,
                    "family": panel_row["pipeline_family"],
                    "biochemical_context": panel_row["biochemical_context"],
                    "target_bin": target,
                    "length": str(lengths[gene_id]),
                    "score": str(score),
                    "evalue": str(evalue),
                    "source": "proteome_v1_exact",
                    "sequence_scope": "full_length",
                    "discovery_round": panel_row["discovery_round"],
                    "matched_protein_id": gene_id,
                    "pfam_hit": str(bool(profiles)).lower(),
                    "full_length_reference_hit": str(bool(references)).lower(),
                    "passing_pfam_profiles": ";".join(profiles),
                    "passing_full_length_reference_queries": ";".join(references),
                    "search_route": route,
                }
            )
    return rows


def write_tsv(path, rows):
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=FIELDS, delimiter="\t", lineterminator="\n"
        )
        writer.writeheader()
        writer.writerows(rows)


def parse_args():
    parser = argparse.ArgumentParser(
        description="Rebuild the deposited VSC search-assignment table from HMMER tblout files."
    )
    parser.add_argument("--hmmer-root", type=Path, required=True)
    parser.add_argument("--proteome", type=Path, required=True)
    parser.add_argument(
        "--query-panel", type=Path, default=ROOT / "queries" / "query_panel.tsv"
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=ROOT / "build" / "search_assignments_1005.tsv",
    )
    parser.add_argument("--score-min", type=float, default=50.0)
    parser.add_argument("--evalue-max", type=float, default=1e-5)
    return parser.parse_args()


def main():
    args = parse_args()
    panel = read_tsv(args.query_panel)
    lengths = fasta_lengths(args.proteome)
    rows = build_rows(panel, args.hmmer_root, lengths, args.score_min, args.evalue_max)
    write_tsv(args.output, rows)
    print(f"Wrote {len(rows)} assignments to {args.output}")


if __name__ == "__main__":
    main()
