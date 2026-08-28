#!/usr/bin/env python3

import argparse
import gzip
import math
import xml.etree.ElementTree as ET
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


def motif_records(data):
    root = ET.fromstring(data)
    records = []
    for motif in root.iter("motif"):
        probabilities = []
        matrix = motif.find("probabilities/alphabet_matrix")
        for row in matrix.findall("alphabet_array"):
            probabilities.append([float(value.text) for value in row.findall("value")])
        records.append(
            {
                "id": motif.get("id"),
                "alt": motif.get("alt"),
                "name": motif.get("name"),
                "width": int(motif.get("width")),
                "sites": int(motif.get("sites")),
                "evalue": float(motif.get("e_value")),
                "probabilities": probabilities,
            }
        )
    return records


def main():
    parser = argparse.ArgumentParser(
        description="Compare fresh de novo MEME models with the archived models."
    )
    parser.add_argument(
        "--fresh-root", type=Path, default=ROOT / "build/motif-discovery"
    )
    parser.add_argument("--archived-root", type=Path, default=ROOT / "motifs/meme")
    args = parser.parse_args()

    failures = []
    for archived in sorted(args.archived_root.glob("*.xml.gz")):
        family = archived.name.removesuffix(".xml.gz")
        fresh = args.fresh_root / family / "meme_discovery/meme.xml"
        if not fresh.exists():
            failures.append(f"{family}: fresh model missing")
            continue
        expected = motif_records(gzip.decompress(archived.read_bytes()))
        observed = motif_records(fresh.read_bytes())
        if len(expected) != len(observed):
            failures.append(f"{family}: motif count {len(observed)} != {len(expected)}")
            continue
        for expected_motif, observed_motif in zip(expected, observed):
            for field in ("id", "alt", "name", "width", "sites"):
                if expected_motif[field] != observed_motif[field]:
                    failures.append(f"{family} {expected_motif['id']}: {field} differs")
            if not math.isclose(
                expected_motif["evalue"],
                observed_motif["evalue"],
                rel_tol=1e-12,
                abs_tol=0.0,
            ):
                failures.append(f"{family} {expected_motif['id']}: e-value differs")
            for expected_row, observed_row in zip(
                expected_motif["probabilities"], observed_motif["probabilities"]
            ):
                if len(expected_row) != len(observed_row) or any(
                    not math.isclose(left, right, rel_tol=1e-12, abs_tol=1e-15)
                    for left, right in zip(expected_row, observed_row)
                ):
                    failures.append(
                        f"{family} {expected_motif['id']}: probability matrix differs"
                    )
                    break

    if failures:
        for failure in failures[:20]:
            print(failure)
        raise SystemExit(
            f"MEME model comparison failed with {len(failures)} differences"
        )
    print("Fresh MEME discovery matches all 19 archived family models")


if __name__ == "__main__":
    main()
