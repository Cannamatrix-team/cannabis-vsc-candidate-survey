#!/usr/bin/env python3

import argparse
import csv
from collections import Counter
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from Bio import Phylo
from matplotlib import rcParams
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch


ROOT = Path(__file__).resolve().parents[1]
FINAL_ORDER = [
    "adh",
    "alliinase",
    "aoc",
    "aos_hpl",
    "cbl",
    "cgs",
    "fmo",
    "ggt",
    "gsh1",
    "gsh2",
    "gst",
    "lf_synthase",
    "lox",
    "mgl",
    "myrosinase",
    "oastl",
    "rhodanese",
    "sams",
    "sult",
    "ugt",
]

TEAL = "#087F73"
BLUE = "#386CB0"
ORANGE = "#E69F00"
PURPLE = "#7A5195"
RED = "#B2473E"
GREY = "#6B7479"
LIGHT_GREY = "#E8ECEE"
INK = "#243238"

rcParams.update(
    {
        "font.family": "DejaVu Sans",
        "font.size": 10,
        "axes.titlesize": 12,
        "axes.labelsize": 10,
        "xtick.labelsize": 9,
        "ytick.labelsize": 9,
        "figure.facecolor": "white",
        "axes.facecolor": "white",
        "savefig.facecolor": "white",
        "savefig.dpi": 300,
    }
)

SHORT_DISPLAY = {
    "adh": "ADH architecture",
    "alliinase": "Alliinase/ACC-synthase-like",
    "aoc": "AOC domain",
    "aos_hpl": "P450 domain",
    "cbl": "PF01053 nearest-reference CBL",
    "cgs": "PF01053 nearest-reference CGS",
    "fmo": "FMO-like domain",
    "ggt": "GGT domain",
    "gsh1": "Glutamate-cysteine ligase",
    "gsh2": "Glutathione synthetase",
    "gst": "GST_N domain",
    "lf_synthase": "Polyketide_cyc2 domain",
    "lox": "LOX architecture",
    "mgl": "PF01053 nearest-reference MGL",
    "myrosinase": "GH1 beta-glucosidase",
    "oastl": "PF00291 PALP domain",
    "rhodanese": "Rhodanese domain",
    "sams": "S-adenosylmethionine synthase",
    "sult": "SULT domain",
    "ugt": "UDPGT domain",
}

CATEGORY_ROWS = [
    ("Sulfur amino-acid and glutathione supply", "OASTL", "PF00291 PALP domain", "PF00291", "oastl"),
    ("Sulfur amino-acid and glutathione supply", "SAMS", "S-adenosylmethionine synthase", "PF00438 + PF02772 + PF02773", "sams"),
    ("Sulfur amino-acid and glutathione supply", "GSH1", "Glutamate-cysteine ligase", "PF04107", "gsh1"),
    ("Sulfur amino-acid and glutathione supply", "CBL", "PF01053 nearest-reference CBL", "PF01053", "cbl"),
    ("Sulfur amino-acid and glutathione supply", "CGS", "PF01053 nearest-reference CGS", "PF01053", "cgs"),
    ("Sulfur amino-acid and glutathione supply", "GSH2", "Glutathione synthetase", "PF03199 + PF03917", "gsh2"),
    ("Conjugation, transfer, and redox chemistry", "UGT", "UDPGT domain", "PF00201", "ugt"),
    ("Conjugation, transfer, and redox chemistry", "GST", "GST_N domain", "PF02798", "gst"),
    ("Conjugation, transfer, and redox chemistry", "FMO", "FMO-like domain", "PF00743", "fmo"),
    ("Conjugation, transfer, and redox chemistry", "SULT", "SULT domain", "PF00685", "sult"),
    ("Conjugation, transfer, and redox chemistry", "Rhodanese", "Rhodanese domain", "PF00581", "rhodanese"),
    ("Conjugation, transfer, and redox chemistry", "GGT", "GGT domain", "PF01019", "ggt"),
    ("Sulfur release and comparative reference architectures", "Myrosinase", "GH1 beta-glucosidase", "PF00232", "myrosinase"),
    ("Sulfur release and comparative reference architectures", "Lachrymatory-factor synthase", "Polyketide_cyc2 domain", "PF10604", "lf_synthase"),
    ("Sulfur release and comparative reference architectures", "Alliinase", "Alliinase/ACC-synthase-like", "PF04863 + PF04864", "alliinase"),
    ("Sulfur release and comparative reference architectures", "MGL-Artocarpus, MGL-Cannabis, MGL-Durio", "PF01053 nearest-reference MGL", "PF01053", "mgl"),
    ("Sulfur release and comparative reference architectures", "CSCBL-like", "no reporting category", "PF01053", "cscbl_like"),
    ("Oxylipin and volatile-product chemistry", "AOS/HPL-motivated P450 search", "P450 domain", "PF00067", "aos_hpl"),
    ("Oxylipin and volatile-product chemistry", "ADH", "ADH architecture", "PF00107 + PF08240", "adh"),
    ("Oxylipin and volatile-product chemistry", "LOX", "LOX architecture", "PF00305 + PF01477", "lox"),
    ("Oxylipin and volatile-product chemistry", "AOC", "AOC domain", "PF06351", "aoc"),
]

TARGET_CROSSWALK = [
    ("SAMS", "", "S-adenosylmethionine synthase", "one-to-one"),
    ("CBL", "PF01053", "PF01053 nearest-reference CBL", "shared-pool phylogenetic assignment"),
    ("CGS", "PF01053", "PF01053 nearest-reference CGS", "shared-pool phylogenetic assignment"),
    ("GSH1", "", "Glutamate-cysteine ligase", "one-to-one"),
    ("GSH2", "", "Glutathione synthetase", "one-to-one"),
    ("OASTL", "", "PF00291 PALP domain", "one-to-one"),
    ("GST", "", "GST_N domain", "one-to-one"),
    ("GGT", "", "GGT domain", "one-to-one"),
    ("UGT", "", "UDPGT domain", "one-to-one"),
    ("SULT", "", "SULT domain", "one-to-one"),
    ("FMO", "", "FMO-like domain", "one-to-one"),
    ("Rhodanese", "", "Rhodanese domain", "one-to-one"),
    ("CSCBL-like", "PF01053", "", "shared-pool reference group; no nearest assignments"),
    ("MGL-Artocarpus", "PF01053", "PF01053 nearest-reference MGL", "three MGL targets collapse to one reporting category"),
    ("MGL-Cannabis", "PF01053", "PF01053 nearest-reference MGL", "three MGL targets collapse to one reporting category"),
    ("MGL-Durio", "PF01053", "PF01053 nearest-reference MGL", "three MGL targets collapse to one reporting category"),
    ("Alliinase", "", "Alliinase/ACC-synthase-like", "one-to-one"),
    ("Lachrymatory-factor synthase", "", "Polyketide_cyc2 domain", "one-to-one"),
    ("Myrosinase", "", "GH1 beta-glucosidase", "one-to-one"),
    ("LOX", "", "LOX architecture", "one-to-one"),
    ("AOS/HPL-motivated P450 search", "", "P450 domain", "one-to-one"),
    ("AOC", "", "AOC domain", "one-to-one"),
    ("ADH", "", "ADH architecture", "one-to-one"),
]

GROUPS = [
    ("Sulfur amino-acid and\nglutathione-supply context", ["sams", "cgs", "cbl", "gsh1", "gsh2", "oastl"], TEAL),
    ("Conjugation, transfer, and\nredox context", ["gst", "ggt", "ugt", "sult", "fmo", "rhodanese"], BLUE),
    ("Sulfur-release and comparative\nreference architectures", ["mgl", "alliinase", "lf_synthase", "myrosinase"], ORANGE),
    ("Oxylipin and\nvolatile-product context", ["lox", "aos_hpl", "aoc", "adh"], PURPLE),
]


def read_tsv(path):
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path, fieldnames, rows):
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=fieldnames, delimiter="\t", lineterminator="\n"
        )
        writer.writeheader()
        writer.writerows(rows)


def load_counts():
    rows = read_tsv(ROOT / "data/candidate_counts_by_category.tsv")
    observed = {
        row["original_pipeline_slug"]: int(row["candidate_count"]) for row in rows
    }
    if set(observed) != set(FINAL_ORDER) or sum(observed.values()) != 941:
        raise ValueError("Candidate reporting categories do not reconcile to 941")
    return {family: observed[family] for family in FINAL_ORDER}


def domain_scope():
    rows = read_tsv(ROOT / "data/candidate_stage_ledger_975.tsv")
    retained = [row for row in rows if row["final_pass"] == "True"]
    recorded_complete = sum(
        row["expected_domain_record_present"] == "true"
        and row["domain_confirmation_status"] == "complete"
        for row in retained
    )
    defaults = {
        row["gene_id"]
        for row in retained
        if row["expected_domain_record_present"] == "false"
    }
    expected_defaults = {
        "GMO.00001695.v1.g225230.t1",
        "GMO.00002556.v1.g618800.t1",
    }
    if recorded_complete != 939 or defaults != expected_defaults or len(retained) != 941:
        raise ValueError("Expected-domain scope does not reconcile to 939 plus 2")
    return {
        "recorded_complete": recorded_complete,
        "mgl_default_passes": len(defaults),
        "final_total": len(retained),
    }


def save(fig, output, name):
    output.mkdir(parents=True, exist_ok=True)
    fig.savefig(
        output / f"{name}.png",
        bbox_inches="tight",
        pad_inches=0.12,
        metadata={"Software": "generate_paper_figures.py"},
    )
    fig.savefig(
        output / f"{name}.pdf",
        bbox_inches="tight",
        pad_inches=0.12,
        metadata={"Creator": "generate_paper_figures.py", "CreationDate": None},
    )
    plt.close(fig)


def write_source_tables(output, counts):
    reporting = [
        {
            "biochemical_context": context,
            "table1_pre_search_target": target,
            "reporting_category": category,
            "pfam_profile_or_architecture": pfam,
            "pipeline_slug": slug,
            "candidate_count": counts.get(slug, 0),
        }
        for context, target, category, pfam, slug in CATEGORY_ROWS
    ]
    occupied = [row for row in reporting if row["reporting_category"] != "no reporting category"]
    if len(occupied) != 20 or sum(row["candidate_count"] for row in reporting) != 941:
        raise ValueError("Manuscript reporting table does not reconcile")
    write_tsv(output / "reporting_category_counts.tsv", list(reporting[0]), reporting)

    crosswalk = [
        {
            "table1_pre_search_target": target,
            "shared_candidate_pool": shared,
            "reporting_category": category,
            "disposition": disposition,
        }
        for target, shared, category, disposition in TARGET_CROSSWALK
    ]
    if len(crosswalk) != 23 or len({row["reporting_category"] for row in crosswalk if row["reporting_category"]}) != 20:
        raise ValueError("Pre-search target crosswalk does not reconcile")
    write_tsv(output / "target_to_reporting_category_crosswalk.tsv", list(crosswalk[0]), crosswalk)

    rules = {row["family"]: row for row in read_tsv(ROOT / "motifs/family_rules.tsv")}
    motif_rows = []
    for label, scope, family in [
        ("Default", "All evaluable categories except ADH, GST, UGT, and Alliinase", "aoc"),
        ("ADH override", "adh", "adh"),
        ("GST override", "gst", "gst"),
        ("UGT override", "ugt", "ugt"),
        ("Alliinase/ACC-synthase-like override", "alliinase", "alliinase"),
    ]:
        rule = rules[family]
        motif_rows.append(
            {
                "rule_scope": label,
                "configured_categories": scope,
                "confirmed_minimum_distinct_motifs": rule["min_motifs_confirmed"],
                "confirmed_minimum_detected_motif_fraction": rule["min_coverage_confirmed"],
                "partial_minimum_distinct_motifs": rule["min_motifs_partial"],
                "partial_minimum_detected_motif_fraction": rule["min_coverage_partial"],
                "disposition": "Confirmed and partial pass",
            }
        )
    motif_rows.append(
        {
            "rule_scope": "Singleton",
            "configured_categories": "gsh2",
            "confirmed_minimum_distinct_motifs": "NA",
            "confirmed_minimum_detected_motif_fraction": "NA",
            "partial_minimum_distinct_motifs": "NA",
            "partial_minimum_detected_motif_fraction": "NA",
            "disposition": "One sequence; motif discovery unevaluable, proceeded directly to the domain gate",
        }
    )
    write_tsv(output / "supplementary_table_s1_motif_screen_rules.tsv", list(motif_rows[0]), motif_rows)

    parameters = [
        ("search", "program", "blastp"),
        ("search", "similarity_matrix", "BLOSUM62"),
        ("search", "e_value_maximum", "0.001"),
        ("search", "query_coverage_per_hsp_minimum_percent", "80"),
        ("search", "maximum_target_sequences_per_query", "1"),
        ("analysis", "primary_hsp_identity_percent", "100"),
        ("analysis", "primary_query_coverage_minimum_percent", "80"),
    ]
    write_tsv(
        output / "expression_mapping_parameters.tsv",
        ["stage", "parameter", "value"],
        [dict(zip(["stage", "parameter", "value"], row)) for row in parameters],
    )


def figure_search_framework(output, scope):
    fig, ax = plt.subplots(figsize=(9.0, 11.2))
    ax.set_xlim(0, 9)
    ax.set_ylim(0, 11.2)
    ax.axis("off")
    ax.text(4.5, 10.95, "Genome-scale candidate discovery workflow", ha="center", va="top", fontsize=16, fontweight="bold", color=INK)
    ax.text(4.5, 10.55, "Comparative biochemistry guided a sequence-based survey of the GMO v1 proteome", ha="center", va="top", fontsize=9.5, color=GREY)
    steps = [
        ("Define the candidate-search panel", "Allium sulfur metabolism, grape thiol precursors, Brassicaceae glucosinolate\nturnover, and oxylipin and volatile-product chemistry", TEAL),
        ("Search the GMO v1 high-confidence proteome", "HMMER full-sequence score ≥ 50 and E-value ≤ 1 × 10⁻⁵\n1,005 search-level assignments", BLUE),
        ("Recover and harmonize full-length candidates", "Exact GMO v1 identifiers  |  975 unique proteins", BLUE),
        ("Resolve the shared PF01053 group", "MAFFT → trimAl → IQ-TREE  |  2 CBL-nearest, 2 CGS-nearest, 2 MGL-nearest", ORANGE),
        ("Screen sequence consistency and expected domains", "MEME/FIMO motif screening + expected Pfam architecture\n" f"{scope['final_total']} proteins in 20 reporting categories\n" f"{scope['recorded_complete']} recorded complete, {scope['mgl_default_passes']} MGL-nearest by rule default", PURPLE),
        ("Cross-reference the Cannabis Expression Atlas", "100% amino-acid identity and ≥ 80% GMO-query coverage\n168 mappings to 128 atlas genes", TEAL),
    ]
    x, width, height = 0.75, 7.5, 1.12
    for index, ((title, detail, color), y) in enumerate(zip(steps, [9.05, 7.55, 6.05, 4.55, 3.05, 1.55]), start=1):
        ax.add_patch(FancyBboxPatch((x, y), width, height, boxstyle="round,pad=0.025,rounding_size=0.10", facecolor="white", edgecolor=color, linewidth=1.8))
        ax.add_patch(FancyBboxPatch((x + 0.22, y + 0.27), 0.58, 0.58, boxstyle="round,pad=0.02,rounding_size=0.12", facecolor=color, edgecolor=color, linewidth=1.0))
        ax.text(x + 0.51, y + 0.56, str(index), ha="center", va="center", fontsize=11, fontweight="bold", color="white")
        ax.text(x + 1.0, y + 0.80, title, ha="left", va="center", fontsize=11.2, fontweight="bold", color=INK)
        ax.text(x + 1.0, y + 0.32, detail, ha="left", va="center", fontsize=8.4, color=GREY, linespacing=1.22)
        if index < len(steps):
            ax.add_patch(FancyArrowPatch((4.5, y - 0.04), (4.5, y - 0.34), arrowstyle="-|>", mutation_scale=13, linewidth=1.2, color=GREY))
    ax.text(4.5, 0.72, "Output: a prioritized candidate catalog for full-length phylogeny and biochemical testing", ha="center", va="center", fontsize=10.2, fontweight="bold", color=INK, bbox=dict(boxstyle="round,pad=0.40", facecolor="#F5F8F8", edgecolor=TEAL, linewidth=1.4))
    save(fig, output, "figure1_search_framework")


def figure_candidate_funnel(output, counts, scope):
    fig = plt.figure(figsize=(12.5, 7.0))
    gs = fig.add_gridspec(1, 2, width_ratios=[0.82, 1.9], wspace=0.58)
    ax0 = fig.add_subplot(gs[0, 0])
    ax1 = fig.add_subplot(gs[0, 1])
    funnel = [
        ("Search-level candidate\nassignments", 1005, BLUE),
        ("Unique proteins after shared-\nprofile resolution", 975, TEAL),
        ("Proteins entering expected-\ndomain screening", 975, GREY),
        ("Retained after expected-\ndomain screening", scope["final_total"], ORANGE),
    ]
    y = np.arange(len(funnel))[::-1]
    widths = [value / 1005 * 0.95 for _, value, _ in funnel]
    for yi, width, (label, value, color) in zip(y, widths, funnel):
        ax0.barh(yi, width, height=0.62, color=color, edgecolor=INK, linewidth=0.7)
        ax0.text(width / 2, yi, f"{value:,}", ha="center", va="center", color="white", fontweight="bold", fontsize=12)
        ax0.text(-0.03, yi, label, ha="right", va="center", color=INK, fontsize=9.5)
    ax0.set_xlim(-0.55, 1.0)
    ax0.set_ylim(-0.7, 3.7)
    ax0.axis("off")
    ax0.set_title("A  Candidate filtering and retention", loc="left", fontweight="bold", color=INK)
    fig.text(0.02, 0.012, "Motif criteria were met by 974 evaluable proteins; the single GSH2 candidate entered expected-domain screening directly. The Pfam check removed 34 proteins.\n" f"Of the {scope['final_total']} retained, {scope['recorded_complete']} have a recorded complete expected architecture; the {scope['mgl_default_passes']} MGL-nearest proteins were retained by default in a category with no expected-domain rule.", fontsize=8.5, color=GREY, ha="left", va="bottom", linespacing=1.3)

    items = sorted(counts.items(), key=lambda item: item[1])
    labels = [SHORT_DISPLAY[family] for family, _ in items]
    values = [value for _, value in items]
    family_to_color = {family: color for _, families, color in GROUPS for family in families}
    colors = [family_to_color[family] for family, _ in items]
    yy = np.arange(len(items))
    ax1.barh(yy, values, color=colors, edgecolor=INK, linewidth=0.5, height=0.72)
    ax1.set_xscale("log")
    ax1.set_xlim(0.8, 700)
    ax1.set_yticks(yy)
    ax1.set_yticklabels(labels, fontsize=8.0)
    ax1.set_xlabel("Candidates (log scale)")
    ax1.grid(axis="x", color=LIGHT_GREY, linewidth=0.8)
    ax1.set_axisbelow(True)
    ax1.spines[["top", "right"]].set_visible(False)
    for yi, value in zip(yy, values):
        ax1.text(value * 1.08, yi, str(value), va="center", fontsize=8.5, color=INK)
    ax1.set_title("B  Final candidates by reporting category", loc="left", fontweight="bold", color=INK)
    save(fig, output, "figure2_candidate_flow_and_counts")


def figure_pf01053_tree(output):
    base = ROOT / "phylogeny/pf01053"
    tree = Phylo.read(base / "tree/group_tree.treefile", "newick")
    assignments = {row["gene_id"]: row for row in read_tsv(base / "assignment_summary.tsv")}
    references = {row["leaf_id"]: row for row in read_tsv(base / "reference_metadata.tsv")}

    def label_func(clade):
        if not clade.is_terminal():
            return None
        name = clade.name
        if name.startswith("GMO."):
            return f"{name}  [{assignments[name]['assigned_family']}-nearest]"
        if name in references:
            row = references[name]
            return f"{row['reference_id']}  [{row['family']}]"
        return name

    def label_color(text):
        return RED if text and text.startswith("GMO.") else GREY

    tree.ladderize()
    fig, ax = plt.subplots(figsize=(10.5, 11.5))
    Phylo.draw(tree, axes=ax, do_show=False, label_func=label_func, label_colors=label_color, show_confidence=False, branch_labels=None)
    ax.set_title("Nearest-reference partition of PF01053 candidates", loc="left", fontweight="bold", color=INK, pad=12)
    ax.set_xlabel("Branch length (substitutions per site)")
    ax.set_ylabel("")
    ax.set_yticks([])
    ax.spines[["top", "right", "left"]].set_visible(False)
    for line in ax.get_lines():
        line.set_linewidth(0.8)
    for text in ax.texts:
        text.set_fontsize(7.2)
        if text.get_text().startswith("GMO."):
            text.set_fontweight("bold")
    ax.text(0.995, 0.015, "Shortest patristic distance: 2 CBL-nearest, 2 CGS-nearest, 2 MGL-nearest", transform=ax.transAxes, ha="right", va="bottom", fontsize=9.0, color=INK, bbox=dict(boxstyle="round,pad=0.3", facecolor="white", edgecolor=LIGHT_GREY))
    ax.text(0.995, 0.115, "Candidate-nearest-reference common-ancestor support\nCBL-nearest: 100, 62   |   CGS-nearest: 100, unavailable\nMGL-nearest: 99, unavailable", transform=ax.transAxes, ha="right", va="bottom", fontsize=8.5, color=INK, bbox=dict(boxstyle="round,pad=0.35", facecolor="#F7F9FA", edgecolor=LIGHT_GREY))
    save(fig, output, "figure3_pf01053_phylogeny")


def figure_expression_mapping(output):
    sensitivity = read_tsv(ROOT / "data/expression/expression_identity_sensitivity.tsv")
    mappings = read_tsv(ROOT / "data/expression/expression_exact100_candidate_mappings_by_bin.tsv")
    bridge = read_tsv(ROOT / "data/expression/expression_bridge_exact100_qcov80.tsv")
    x = [int(row["minimum_percent_identity"]) for row in sensitivity]

    fig = plt.figure(figsize=(11.5, 6.3))
    gs = fig.add_gridspec(1, 2, width_ratios=[1.0, 1.35], wspace=0.35)
    ax0 = fig.add_subplot(gs[0, 0])
    ax1 = fig.add_subplot(gs[0, 1])
    for field, color, label in [
        ("candidate_mappings", BLUE, "GMO candidate mappings"),
        ("unique_atlas_genes", TEAL, "Unique atlas genes"),
        ("tissue_specific_candidate_mappings", ORANGE, "Tissue-specific candidate mappings"),
    ]:
        ax0.plot(x, [int(row[field]) for row in sensitivity], marker="o", color=color, linewidth=2.0, label=label)
    ax0.invert_xaxis()
    ax0.set_xlabel("Minimum amino-acid identity (%)")
    ax0.set_ylabel("Count")
    ax0.grid(color=LIGHT_GREY, linewidth=0.8)
    ax0.set_axisbelow(True)
    ax0.spines[["top", "right"]].set_visible(False)
    ax0.legend(frameon=False, fontsize=8.5, loc="upper right")
    ax0.set_title("A  Identity threshold changes mapping counts", loc="left", fontweight="bold", color=INK, fontsize=10.5)

    families = [row["family_slug"] for row in mappings]
    values = [int(row["candidate_mappings"]) for row in mappings]
    labels = [SHORT_DISPLAY[family] for family in families]
    family_to_color = {family: color for _, group, color in GROUPS for family in group}
    colors = [family_to_color[family] for family in families]
    yy = np.arange(len(families))
    ax1.barh(yy, values, color=colors, edgecolor=INK, linewidth=0.5, height=0.72)
    ax1.set_yticks(yy)
    ax1.set_yticklabels(labels)
    ax1.set_xlabel("100%-identity mappings (query coverage ≥ 80%)")
    ax1.grid(axis="x", color=LIGHT_GREY, linewidth=0.8)
    ax1.set_axisbelow(True)
    ax1.spines[["top", "right"]].set_visible(False)
    for yi, value in zip(yy, values):
        ax1.text(value + 0.7, yi, str(value), va="center", fontsize=8.5, color=INK)
    ax1.set_title("B  Stringent mappings by candidate category", loc="left", fontweight="bold", color=INK, fontsize=10.5)

    best_by_subject = {}
    for row in sorted(bridge, key=lambda value: (value["subject_acc"], -float(value["bit_score"]))):
        best_by_subject.setdefault(row["subject_acc"], row)
    unique = list(best_by_subject.values())
    class_counts = Counter(row["Classification"] or "Unclassified" for row in unique)
    tissue = [row for row in unique if row["Classification"] == "Tissue-Specific"]
    root_n = sum("Root" in row["Specific tissues"] for row in tissue)
    trichome_n = sum("Trichome" in row["Specific tissues"] for row in tissue)
    fig.text(0.01, 0.01, f"At 100% aligned-region identity and query coverage ≥ 80%: {len(bridge):,} candidate mappings to {len(best_by_subject):,} atlas genes; " f"{class_counts['Tissue-Specific']} unique atlas genes were tissue-specific ({root_n} root; {trichome_n} trichome).", ha="left", va="bottom", fontsize=9.0, color=INK)
    save(fig, output, "figure4_expression_mapping")


def main():
    parser = argparse.ArgumentParser(description="Regenerate VSC manuscript figures and source tables from public inputs.")
    parser.add_argument("--output-root", type=Path, default=ROOT / "build/manuscript")
    args = parser.parse_args()
    counts = load_counts()
    scope = domain_scope()
    write_source_tables(args.output_root / "tables", counts)
    figure_search_framework(args.output_root / "figures", scope)
    figure_candidate_funnel(args.output_root / "figures", counts, scope)
    figure_pf01053_tree(args.output_root / "figures")
    figure_expression_mapping(args.output_root / "figures")
    print(f"Wrote manuscript figures and source tables to {args.output_root}")


if __name__ == "__main__":
    main()
