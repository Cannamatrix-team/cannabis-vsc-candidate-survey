#!/usr/bin/env python3

import csv
import gzip
import hashlib
from collections import Counter
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


def read_tsv(relative_path):
    with (ROOT / relative_path).open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def require(actual, expected, label):
    if actual != expected:
        raise SystemExit(f"FAIL: {label}: expected {expected!r}, found {actual!r}")
    print(f"PASS: {label}")


def fasta_ids(relative_path):
    identifiers = []
    with (ROOT / relative_path).open(encoding="utf-8") as handle:
        for line in handle:
            if line.startswith(">"):
                identifiers.append(line[1:].split(maxsplit=1)[0])
    return identifiers


def main():
    search = read_tsv("data/search_assignments_1005.tsv")
    ledger = read_tsv("data/candidate_stage_ledger_975.tsv")
    catalog = read_tsv("data/candidate_catalog_941.tsv")
    category_counts = read_tsv("data/candidate_counts_by_category.tsv")
    motif_models = read_tsv("motifs/model_manifest.tsv")
    domain_rules = read_tsv("domains/family_rules.tsv")
    expression = read_tsv("data/expression/expression_bridge_exact100_qcov80.tsv")
    expression_by_bin = read_tsv(
        "data/expression/expression_exact100_candidate_mappings_by_bin.tsv"
    )
    expression_sensitivity = read_tsv(
        "data/expression/expression_identity_sensitivity.tsv"
    )
    pf01053 = read_tsv("phylogeny/pf01053/assignment_summary.tsv")

    search_ids = {row["gene_id"] for row in search}
    ledger_ids = {row["gene_id"] for row in ledger}
    require(len(search), 1005, "initial target-bin assignments")
    require(len(search_ids), 975, "unique searched proteins")
    require(len(ledger), 975, "stage-ledger rows")
    require(len(ledger_ids), 975, "unique stage-ledger proteins")
    require(ledger_ids, search_ids, "search and stage-ledger protein sets")

    stage_sequence_ids = fasta_ids("data/candidate_sequences_975.fasta")
    require(len(stage_sequence_ids), 975, "motif-stage FASTA records")
    require(len(set(stage_sequence_ids)), 975, "unique motif-stage FASTA identifiers")
    require(
        set(stage_sequence_ids),
        ledger_ids,
        "motif-stage FASTA and stage-ledger protein sets",
    )

    require(len(motif_models), 20, "motif model manifest families")
    require(
        Counter(row["meme_evaluable"] for row in motif_models),
        Counter({"true": 19, "false": 1}),
        "motif model evaluability",
    )
    require(
        sum(int(row["n_sequences"]) for row in motif_models),
        975,
        "motif model sequence total",
    )
    for row in motif_models:
        if not row["model_file"]:
            require(row["family"], "gsh2", "motif-model singleton family")
            continue
        model = gzip.decompress((ROOT / row["model_file"]).read_bytes())
        require(
            hashlib.sha256(model).hexdigest(),
            row["model_xml_sha256"],
            f"{row['family']} motif model checksum",
        )

    ledger_families = {row["reporting_category"] for row in ledger}
    require(len(domain_rules), 20, "expected-domain rule families")
    require(
        {row["family"] for row in domain_rules},
        ledger_families,
        "expected-domain rule and ledger families",
    )
    require(
        Counter(row["record_expected"] for row in domain_rules),
        Counter({"true": 19, "false": 1}),
        "expected-domain record rules",
    )
    no_record_rules = [row for row in domain_rules if row["record_expected"] == "false"]
    require(len(no_record_rules), 1, "expected-domain no-record rule count")
    require(no_record_rules[0]["family"], "mgl", "expected-domain no-record family")

    configured_profiles = {
        pfam_id
        for row in domain_rules
        for pfam_id in row["pfam_ids"].split(",")
        if pfam_id
    }
    packaged_profiles = {
        path.name.split(".", 1)[0]
        for path in (ROOT / "inputs/search_queries/pfam").glob("*.hmm")
    }
    require(len(configured_profiles), 24, "configured expected-domain profiles")
    require(
        configured_profiles,
        packaged_profiles,
        "configured and packaged expected-domain profiles",
    )

    require(
        Counter(row["expected_domain_record_present"] for row in ledger),
        Counter({"true": 973, "false": 2}),
        "expected-domain record presence",
    )
    no_domain_records = [
        row for row in ledger if row["expected_domain_record_present"] == "false"
    ]
    require(
        {row["reporting_category"] for row in no_domain_records},
        {"mgl"},
        "expected-domain no-record ledger family",
    )
    require(
        all(
            not row["pfam_confirmation_rule"] and not row["pfam_required_ids"]
            for row in no_domain_records
        ),
        True,
        "expected-domain no-record fields",
    )
    domain_statuses = Counter(
        row["domain_confirmation_status"]
        for row in ledger
        if row["expected_domain_record_present"] == "true"
    )
    require(
        domain_statuses,
        Counter({"complete": 939, "partial": 32, "absent": 2}),
        "expected-domain confirmation statuses",
    )
    rules_by_family = {row["family"]: row for row in domain_rules}
    for family in sorted(ledger_families - {"mgl"}):
        family_rows = [row for row in ledger if row["reporting_category"] == family]
        observed = {
            (row["pfam_required_ids"], row["pfam_confirmation_rule"])
            for row in family_rows
        }
        rule = rules_by_family[family]
        require(
            observed,
            {(rule["required_pfam_ids"], rule["confirmation_rule"])},
            f"{family} expected-domain rule",
        )

    def bool_text(value):
        return "True" if value else "False"

    biological_pass = {
        row["gene_id"]: row["meme_status"] in {"confirmed", "partial"} for row in ledger
    }
    operational_pass = {
        row["gene_id"]: biological_pass[row["gene_id"]]
        or row["singleton_passthrough"] == "True"
        for row in ledger
    }
    expected_pfam_pass = {
        row["gene_id"]: row["expected_domain_record_present"] == "false"
        or row["domain_confirmation_status"] == "complete"
        for row in ledger
    }
    expected_final_pass = {
        row["gene_id"]: operational_pass[row["gene_id"]]
        and expected_pfam_pass[row["gene_id"]]
        for row in ledger
    }
    require(
        all(
            row["meme_biological_pass"] == bool_text(biological_pass[row["gene_id"]])
            and row["meme_operational_pass"]
            == bool_text(operational_pass[row["gene_id"]])
            and row["meme_pass"] == bool_text(biological_pass[row["gene_id"]])
            for row in ledger
        ),
        True,
        "motif gate semantics",
    )
    require(
        all(
            row["pfam_pass"] == bool_text(expected_pfam_pass[row["gene_id"]])
            and row["pfam_gate_active"] == "True"
            for row in ledger
        ),
        True,
        "expected-domain gate semantics",
    )
    require(
        all(
            row["final_pass"] == bool_text(expected_final_pass[row["gene_id"]])
            for row in ledger
        ),
        True,
        "final validation-gate semantics",
    )

    expected_outcomes = Counter(
        {
            "retained": 939,
            "retained_by_default_no_configured_domain_rule": 2,
            "removed_incomplete_expected_domain_architecture": 34,
        }
    )
    outcomes = Counter(row["final_outcome"] for row in ledger)
    retained_ids = {row["gene_id"] for row in ledger if row["final_pass"] == "True"}
    require(outcomes, expected_outcomes, "stage-ledger final outcomes")
    require(len(retained_ids), 941, "retained stage-ledger proteins")

    catalog_ids = [row["gene_id"] for row in catalog]
    require(len(catalog), 941, "catalog rows")
    require(len(set(catalog_ids)), 941, "unique catalog proteins")
    require(
        set(catalog_ids), retained_ids, "catalog and retained stage-ledger protein sets"
    )
    require(
        {row["final_pass"] for row in catalog}, {"True"}, "catalog final-pass values"
    )

    expected_categories = Counter(
        {
            row["original_pipeline_slug"]: int(row["candidate_count"])
            for row in category_counts
        }
    )
    catalog_categories = Counter(row["family"].lower() for row in catalog)
    require(len(category_counts), 20, "reporting categories")
    require(sum(expected_categories.values()), 941, "category-count total")
    require(catalog_categories, expected_categories, "catalog category counts")

    sequence_ids = fasta_ids("data/candidate_sequences_941.fasta")
    require(len(sequence_ids), 941, "FASTA records")
    require(len(set(sequence_ids)), 941, "unique FASTA identifiers")
    require(set(sequence_ids), set(catalog_ids), "FASTA and catalog protein sets")

    expression_queries = {row["query_acc"] for row in expression}
    tissue_specific = [
        row for row in expression if row["Classification"] == "Tissue-Specific"
    ]
    require(len(expression), 168, "exact-identity expression mappings")
    require(
        len({row["subject_acc"] for row in expression}),
        128,
        "distinct mapped atlas genes",
    )
    require(
        expression_queries <= set(catalog_ids),
        True,
        "expression queries belong to catalog",
    )
    require(
        all(
            float(row["p_identity"]) == 100.0 and float(row["query_coverage"]) >= 0.8
            for row in expression
        ),
        True,
        "expression identity and coverage thresholds",
    )
    require(len(tissue_specific), 36, "tissue-specific candidate mappings")
    require(
        len({row["subject_acc"] for row in tissue_specific}),
        29,
        "distinct tissue-specific atlas genes",
    )
    require(
        sum(int(row["candidate_mappings"]) for row in expression_by_bin),
        168,
        "expression mappings by-bin total",
    )
    identity_100 = next(
        row
        for row in expression_sensitivity
        if row["minimum_percent_identity"] == "100"
    )
    require(
        identity_100["candidate_mappings"], "168", "100%-identity sensitivity mappings"
    )
    require(
        identity_100["unique_atlas_genes"], "128", "100%-identity sensitivity genes"
    )

    require(len(pf01053), 6, "PF01053 resolved candidates")
    require(
        Counter(row["assigned_family"] for row in pf01053),
        Counter({"CBL": 2, "CGS": 2, "MGL": 2}),
        "PF01053 family assignments",
    )

    print("Verified the deposited 975-protein stage ledger and 941-candidate release.")


if __name__ == "__main__":
    main()
