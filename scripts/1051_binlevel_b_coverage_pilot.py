#!/usr/bin/env python3
"""Build an auditable Bin50 B-lineage RNA coverage pilot.

This script reuses the frozen immune candidate pool and the previously
extracted raw-count-derived center and center-excluded 50 um neighborhood
marker scores. It does not rerun RCTD, rebuild the graph, recluster bins, or
use disease labels to choose thresholds.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.backends.backend_pdf import PdfPages


TARGET_PANELS = {"B__core", "B__extended_only"}
MIN_EXTENDED_GENE_N = 2
MIN_CORE_CORROBORATION_GENE_N = 1
MIN_CORE_ONLY_GENE_N = 2


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--score-file", required=True, type=Path)
    parser.add_argument("--assignments", required=True, type=Path)
    parser.add_argument("--output-dir", required=True, type=Path)
    parser.add_argument("--chunksize", type=int, default=250_000)
    return parser.parse_args()


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def evidence_tier(
    extended_gene_n: pd.Series,
    core_gene_n: pd.Series,
) -> pd.Series:
    extended_supported = extended_gene_n >= MIN_EXTENDED_GENE_N
    core_corroborated = core_gene_n >= MIN_CORE_CORROBORATION_GENE_N
    core_only_supported = core_gene_n >= MIN_CORE_ONLY_GENE_N
    any_detected = (extended_gene_n + core_gene_n) > 0
    return pd.Series(
        np.select(
            [
                extended_supported & core_corroborated,
                extended_supported,
                ~extended_supported & core_only_supported,
                any_detected,
            ],
            [
                "extended_supported_core_corroborated",
                "extended_supported_no_core",
                "core_supported_no_extended",
                "low_information_single_gene",
            ],
            default="no_detectable_b_coverage",
        ),
        index=extended_gene_n.index,
    )


def read_b_gene_scores(path: Path, chunksize: int) -> pd.DataFrame:
    usecols = [
        "chip",
        "internal_id",
        "Bin_ID",
        "patient",
        "ROI",
        "x_um",
        "y_um",
        "original_coarse_niche",
        "nUMI",
        "nGene",
        "neighbor_count",
        "raw_center_umi",
        "raw_center_detected_genes",
        "neighborhood_total_umi",
        "neighborhood_mean_detected_genes",
        "feature_type",
        "panel",
        "gene",
        "compartment",
        "score",
        "detection_fraction",
    ]
    retained: list[pd.DataFrame] = []
    for chunk in pd.read_csv(
        path,
        sep="\t",
        compression="gzip",
        usecols=usecols,
        chunksize=chunksize,
        low_memory=False,
    ):
        keep = chunk["feature_type"].eq("gene") & chunk["panel"].isin(
            TARGET_PANELS
        )
        if keep.any():
            retained.append(chunk.loc[keep].copy())
    if not retained:
        raise RuntimeError("No B core/extended gene rows were found")
    return pd.concat(retained, ignore_index=True)


def build_bin_table(scores: pd.DataFrame) -> pd.DataFrame:
    metadata_columns = [
        "chip",
        "internal_id",
        "Bin_ID",
        "patient",
        "ROI",
        "x_um",
        "y_um",
        "original_coarse_niche",
        "nUMI",
        "nGene",
        "neighbor_count",
        "raw_center_umi",
        "raw_center_detected_genes",
        "neighborhood_total_umi",
        "neighborhood_mean_detected_genes",
    ]
    metadata = scores[metadata_columns].drop_duplicates()
    if metadata["internal_id"].duplicated().any():
        raise RuntimeError("Bin metadata are not unique by internal_id")

    grouped = (
        scores.assign(detected=scores["detection_fraction"].gt(0))
        .groupby(
            ["internal_id", "compartment", "panel"],
            as_index=False,
        )
        .agg(
            available_gene_n=("gene", "nunique"),
            detected_gene_n=("detected", "sum"),
            mean_gene_score=("score", "mean"),
            mean_detection_fraction=("detection_fraction", "mean"),
        )
    )
    wide = grouped.pivot(
        index="internal_id",
        columns=["compartment", "panel"],
        values=[
            "available_gene_n",
            "detected_gene_n",
            "mean_gene_score",
            "mean_detection_fraction",
        ],
    )
    wide.columns = [
        f"{compartment}_{panel.replace('B__', '')}_{metric}"
        for metric, compartment, panel in wide.columns
    ]
    wide = wide.reset_index()
    result = metadata.merge(wide, on="internal_id", validate="one_to_one")

    for compartment in ["center", "neighborhood"]:
        extended = result[
            f"{compartment}_extended_only_detected_gene_n"
        ].fillna(0)
        core = result[f"{compartment}_core_detected_gene_n"].fillna(0)
        result[f"{compartment}_b_extended_supported"] = (
            extended >= MIN_EXTENDED_GENE_N
        )
        result[f"{compartment}_b_core_corroborated"] = (
            core >= MIN_CORE_CORROBORATION_GENE_N
        )
        result[f"{compartment}_b_core_only_supported"] = (
            core >= MIN_CORE_ONLY_GENE_N
        )
        result[f"{compartment}_b_evidence_tier"] = evidence_tier(
            extended,
            core,
        )
        result[f"{compartment}_b_extended_gain"] = (
            result[f"{compartment}_b_extended_supported"]
            & ~result[f"{compartment}_b_core_only_supported"]
        )

    result["primary_b_coverage"] = result[
        "neighborhood_b_extended_supported"
    ]
    result["primary_b_coverage_core_corroborated"] = (
        result["primary_b_coverage"]
        & result["neighborhood_b_core_corroborated"]
    )
    result["center_neighbor_concordance"] = np.select(
        [
            result["center_b_extended_supported"]
            & result["neighborhood_b_extended_supported"],
            result["center_b_extended_supported"]
            & ~result["neighborhood_b_extended_supported"],
            ~result["center_b_extended_supported"]
            & result["neighborhood_b_extended_supported"],
        ],
        ["both", "center_only", "neighborhood_only"],
        default="neither",
    )
    return result


def add_assignments(bin_table: pd.DataFrame, path: Path) -> pd.DataFrame:
    assignments = pd.read_csv(
        path,
        sep="\t",
        compression="gzip",
        usecols=["chip", "internal_id", "k", "cluster"],
    )
    pivot = assignments.pivot(
        index=["chip", "internal_id"],
        columns="k",
        values="cluster",
    ).reset_index()
    pivot = pivot.rename(columns={3: "immune_k3", 4: "immune_k4", 5: "immune_k5"})
    return bin_table.merge(
        pivot,
        on=["chip", "internal_id"],
        how="left",
        validate="one_to_one",
    )


def summarize_patient(bin_table: pd.DataFrame) -> pd.DataFrame:
    summary = (
        bin_table.groupby(["chip", "patient", "ROI"], as_index=False)
        .agg(
            candidate_bin_n=("internal_id", "size"),
            median_n_umi=("nUMI", "median"),
            center_extended_supported_n=("center_b_extended_supported", "sum"),
            center_core_only_supported_n=("center_b_core_only_supported", "sum"),
            center_extended_gain_n=("center_b_extended_gain", "sum"),
            neighborhood_extended_supported_n=(
                "neighborhood_b_extended_supported",
                "sum",
            ),
            neighborhood_core_only_supported_n=(
                "neighborhood_b_core_only_supported",
                "sum",
            ),
            neighborhood_extended_gain_n=("neighborhood_b_extended_gain", "sum"),
            neighborhood_core_corroborated_n=(
                "neighborhood_b_core_corroborated",
                "sum",
            ),
        )
    )
    count_columns = [
        column
        for column in summary.columns
        if column.endswith("_n") and column != "candidate_bin_n"
    ]
    for column in count_columns:
        summary[column.replace("_n", "_fraction")] = (
            summary[column] / summary["candidate_bin_n"]
        )
    return summary


def summarize_chip(bin_table: pd.DataFrame) -> pd.DataFrame:
    return (
        bin_table.groupby("chip", as_index=False)
        .agg(
            patient_n=("patient", "nunique"),
            candidate_bin_n=("internal_id", "size"),
            median_n_umi=("nUMI", "median"),
            neighborhood_extended_supported_n=(
                "neighborhood_b_extended_supported",
                "sum",
            ),
            neighborhood_core_only_supported_n=(
                "neighborhood_b_core_only_supported",
                "sum",
            ),
            neighborhood_extended_gain_n=("neighborhood_b_extended_gain", "sum"),
            neighborhood_core_corroborated_n=(
                "primary_b_coverage_core_corroborated",
                "sum",
            ),
        )
        .assign(
            neighborhood_extended_supported_fraction=lambda x: (
                x["neighborhood_extended_supported_n"] / x["candidate_bin_n"]
            ),
            neighborhood_core_only_supported_fraction=lambda x: (
                x["neighborhood_core_only_supported_n"] / x["candidate_bin_n"]
            ),
            neighborhood_extended_gain_fraction=lambda x: (
                x["neighborhood_extended_gain_n"] / x["candidate_bin_n"]
            ),
        )
    )


def plot_patient_summary(summary: pd.DataFrame, path: Path) -> None:
    plot = summary.sort_values(["chip", "patient"]).copy()
    labels = plot["chip"] + " | " + plot["patient"]
    core = plot["neighborhood_core_only_supported_fraction"]
    extended = plot["neighborhood_extended_supported_fraction"]
    x = np.arange(len(plot))
    fig, ax = plt.subplots(figsize=(13, 6.5))
    ax.bar(x - 0.2, core, width=0.4, label="Core-only rule", color="#7F8C8D")
    ax.bar(
        x + 0.2,
        extended,
        width=0.4,
        label="Extended coverage rule",
        color="#17827B",
    )
    ax.set_xticks(x, labels, rotation=65, ha="right")
    ax.set_ylabel("Fraction of immune candidate Bin50")
    ax.set_title("Patient-level B RNA coverage in center-excluded 50 um neighborhoods")
    ax.legend(frameon=False)
    ax.grid(axis="y", color="#DDE5E8", linewidth=0.7)
    fig.tight_layout()
    fig.savefig(path.with_suffix(".png"), dpi=180)
    fig.savefig(path.with_suffix(".pdf"))
    plt.close(fig)


def plot_spatial_maps(bin_table: pd.DataFrame, path: Path) -> None:
    colors = {
        "extended_supported_core_corroborated": "#0B6E69",
        "extended_supported_no_core": "#E67E22",
        "core_supported_no_extended": "#3973B7",
        "low_information_single_gene": "#B8C4CC",
        "no_detectable_b_coverage": "#ECEFF1",
    }
    order = list(colors)
    with PdfPages(path) as pdf:
        for (chip, patient), local in bin_table.groupby(
            ["chip", "patient"],
            sort=True,
        ):
            fig, axes = plt.subplots(1, 2, figsize=(12, 5.5))
            for axis, compartment in zip(axes, ["center", "neighborhood"]):
                tier_column = f"{compartment}_b_evidence_tier"
                for tier in order:
                    selected = local[tier_column].eq(tier)
                    if not selected.any():
                        continue
                    axis.scatter(
                        local.loc[selected, "x_um"],
                        local.loc[selected, "y_um"],
                        s=3.5,
                        c=colors[tier],
                        label=tier,
                        linewidths=0,
                        rasterized=True,
                    )
                axis.set_title(f"{compartment.capitalize()} B coverage")
                axis.set_aspect("equal")
                axis.invert_yaxis()
                axis.set_xlabel("x (um)")
                axis.set_ylabel("y (um)")
            handles, labels = axes[1].get_legend_handles_labels()
            fig.legend(
                handles,
                labels,
                loc="lower center",
                ncol=3,
                frameon=False,
                fontsize=8,
            )
            fig.suptitle(f"{chip} | {patient}")
            fig.tight_layout(rect=[0, 0.08, 1, 0.94])
            pdf.savefig(fig, dpi=180)
            plt.close(fig)


def write_report(
    output_dir: Path,
    bin_table: pd.DataFrame,
    patient_summary: pd.DataFrame,
    chip_summary: pd.DataFrame,
) -> None:
    total = len(bin_table)
    primary_n = int(bin_table["primary_b_coverage"].sum())
    gain_n = int(bin_table["neighborhood_b_extended_gain"].sum())
    corroborated_n = int(
        bin_table["primary_b_coverage_core_corroborated"].sum()
    )
    all_patient_supported = int(
        (patient_summary["neighborhood_extended_supported_n"] > 0).sum()
    )
    chip_summary_text = chip_summary.to_csv(sep="\t", index=False)
    report = f"""# Bin-level B-lineage RNA coverage pilot

## Analysis boundary

This pilot reused the frozen 29,645-bin immune candidate pool and previously
extracted raw-count-derived marker scores. It did not rerun RCTD, rebuild the
50 um graph, recluster bins, use disease labels, or alter the extended marker
library.

## Fixed rule

- Primary regional endpoint: center-excluded 50 um neighborhood B coverage.
- Extended support: at least {MIN_EXTENDED_GENE_N} detected extended-only B genes.
- Core corroboration: at least {MIN_CORE_CORROBORATION_GENE_N} detected core B gene.
- Core does not veto extended support.
- Shared/functional broad genes are not used.
- Single-gene evidence remains low-information rather than being forced negative.

## Main counts

- Candidate Bin50: {total:,}
- Neighborhood extended-supported Bin50: {primary_n:,} ({primary_n / total:.2%})
- Extended-supported with core corroboration: {corroborated_n:,}
- Extended gain beyond a >=2-core-gene rule: {gain_n:,}
- Patients with at least one extended-supported neighborhood Bin50: {all_patient_supported}/{len(patient_summary)}

## Interpretation boundary

The output is a B-lineage RNA coverage mask, not a B-cell identity call, B-cell
niche, aggregate, or TLS annotation. Center and neighborhood evidence are
stored separately. Disease comparison is not part of this pilot.

## Chip summary

```text
{chip_summary_text}```
"""
    (output_dir / "reports" / "08_B_COVERAGE_PILOT_KEY_FINDINGS.md").write_text(
        report,
        encoding="utf-8",
    )


def main() -> None:
    args = parse_args()
    output_dir = args.output_dir.resolve()
    for directory in [
        output_dir,
        output_dir / "figures",
        output_dir / "reports",
        output_dir / "qc",
        output_dir / "logs",
    ]:
        directory.mkdir(parents=True, exist_ok=True)

    scores = read_b_gene_scores(args.score_file, args.chunksize)
    bin_table = build_bin_table(scores)
    bin_table = add_assignments(bin_table, args.assignments)
    if len(bin_table) != 29_645:
        raise RuntimeError(f"Expected 29,645 candidate bins, observed {len(bin_table):,}")
    if bin_table["internal_id"].duplicated().any():
        raise RuntimeError("Final Bin-level table contains duplicated internal_id")

    output_table = output_dir / "01_BIN_LEVEL_B_COVERAGE_MASK.tsv.gz"
    bin_table.to_csv(output_table, sep="\t", index=False, compression="gzip")

    patient_summary = summarize_patient(bin_table)
    patient_summary.to_csv(
        output_dir / "02_PATIENT_LEVEL_B_COVERAGE_SUMMARY.tsv",
        sep="\t",
        index=False,
    )
    chip_summary = summarize_chip(bin_table)
    chip_summary.to_csv(
        output_dir / "03_CHIP_LEVEL_B_COVERAGE_SUMMARY.tsv",
        sep="\t",
        index=False,
    )

    coarse = (
        bin_table.groupby(
            ["chip", "original_coarse_niche"],
            as_index=False,
        )
        .agg(
            candidate_bin_n=("internal_id", "size"),
            neighborhood_extended_supported_n=("primary_b_coverage", "sum"),
            neighborhood_extended_gain_n=("neighborhood_b_extended_gain", "sum"),
        )
    )
    coarse["neighborhood_extended_supported_fraction"] = (
        coarse["neighborhood_extended_supported_n"] / coarse["candidate_bin_n"]
    )
    coarse.to_csv(
        output_dir / "04_COARSE_NICHE_B_COVERAGE.tsv",
        sep="\t",
        index=False,
    )

    immune_rows: list[pd.DataFrame] = []
    for k in [3, 4, 5]:
        column = f"immune_k{k}"
        summary = (
            bin_table.groupby(["chip", column], as_index=False)
            .agg(
                candidate_bin_n=("internal_id", "size"),
                neighborhood_extended_supported_n=("primary_b_coverage", "sum"),
                neighborhood_extended_gain_n=("neighborhood_b_extended_gain", "sum"),
            )
            .rename(columns={column: "cluster"})
        )
        summary.insert(1, "k", k)
        summary["neighborhood_extended_supported_fraction"] = (
            summary["neighborhood_extended_supported_n"]
            / summary["candidate_bin_n"]
        )
        immune_rows.append(summary)
    pd.concat(immune_rows, ignore_index=True).to_csv(
        output_dir / "05_IMMUNE_K345_B_COVERAGE.tsv",
        sep="\t",
        index=False,
    )

    bin_table["patient_umi_quintile"] = bin_table.groupby("patient")[
        "nUMI"
    ].transform(lambda x: pd.qcut(x.rank(method="first"), 5, labels=False) + 1)
    umi_summary = (
        bin_table.groupby(
            ["chip", "patient", "patient_umi_quintile"],
            as_index=False,
        )
        .agg(
            candidate_bin_n=("internal_id", "size"),
            median_n_umi=("nUMI", "median"),
            neighborhood_extended_supported_n=("primary_b_coverage", "sum"),
            neighborhood_extended_gain_n=("neighborhood_b_extended_gain", "sum"),
        )
    )
    umi_summary["neighborhood_extended_supported_fraction"] = (
        umi_summary["neighborhood_extended_supported_n"]
        / umi_summary["candidate_bin_n"]
    )
    umi_summary.to_csv(
        output_dir / "06_PATIENT_UMI_QUINTILE_B_COVERAGE.tsv",
        sep="\t",
        index=False,
    )

    audit = pd.DataFrame(
        {
            "audit_item": [
                "candidate_bin_n",
                "unique_internal_id_n",
                "patient_n",
                "chip_n",
                "center_and_neighborhood_stored_separately",
                "extended_core_separated",
                "shared_functional_genes_used",
                "disease_labels_used",
                "graph_rebuilt",
                "reclustering_performed",
            ],
            "observed": [
                len(bin_table),
                bin_table["internal_id"].nunique(),
                bin_table["patient"].nunique(),
                bin_table["chip"].nunique(),
                True,
                True,
                False,
                False,
                False,
                False,
            ],
            "expected": [29_645, 29_645, 21, 3, True, True, False, False, False, False],
        }
    )
    audit["passed"] = audit["observed"].astype(str).eq(
        audit["expected"].astype(str)
    )
    audit.to_csv(
        output_dir / "qc" / "07_B_COVERAGE_PILOT_AUDIT.tsv",
        sep="\t",
        index=False,
    )
    if not audit["passed"].all():
        raise RuntimeError("B coverage pilot audit failed")

    plot_patient_summary(
        patient_summary,
        output_dir / "figures" / "01_PATIENT_B_COVERAGE_CORE_VS_EXTENDED_EN",
    )
    plot_spatial_maps(
        bin_table,
        output_dir / "figures" / "02_ALL_PATIENT_B_COVERAGE_SPATIAL_MAPS_EN.pdf",
    )
    write_report(output_dir, bin_table, patient_summary, chip_summary)

    manifest = pd.DataFrame(
        [
            {
                "purpose": "marker_score_input",
                "path": str(args.score_file.resolve()),
                "size_bytes": args.score_file.stat().st_size,
                "sha256": sha256(args.score_file),
            },
            {
                "purpose": "frozen_assignment_input",
                "path": str(args.assignments.resolve()),
                "size_bytes": args.assignments.stat().st_size,
                "sha256": sha256(args.assignments),
            },
        ]
    )
    manifest.to_csv(output_dir / "input_manifest.tsv", sep="\t", index=False)
    session = {
        "script": str(Path(__file__).resolve()),
        "min_extended_gene_n": MIN_EXTENDED_GENE_N,
        "min_core_corroboration_gene_n": MIN_CORE_CORROBORATION_GENE_N,
        "min_core_only_gene_n": MIN_CORE_ONLY_GENE_N,
        "primary_endpoint": "center-excluded 50 um neighborhood",
        "random_steps": "none",
    }
    (output_dir / "SESSION_INFO.json").write_text(
        json.dumps(session, indent=2),
        encoding="utf-8",
    )
    (output_dir / "ANALYSIS_COMPLETE_PENDING_USER_APPROVAL.ok").write_text(
        "ANALYSIS_COMPLETE_PENDING_USER_APPROVAL\n",
        encoding="utf-8",
    )
    print(f"B coverage pilot complete: {output_dir}")


if __name__ == "__main__":
    main()
