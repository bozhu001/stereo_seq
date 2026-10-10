#!/usr/bin/env python3
"""Calibrate Bin50 B-lineage coverage against frozen multi-resolution anchors.

The analysis is deliberately label-free with respect to disease.  Positive
anchors are bins assigned to an independently validated B-enriched cluster at
at least two of K=3/4/5; bins assigned at zero resolutions form the calibration
background and one-resolution bins are left ambiguous.  All score thresholds
are estimated within patient and neighborhood-UMI quintile.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


PANELS = {"B__core", "B__extended_only"}
QUANTILES = (0.85, 0.90, 0.95)
MIN_GENE_COUNTS = (1, 2, 3, 4, 5)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--score-file", required=True, type=Path)
    parser.add_argument("--assignments", required=True, type=Path)
    parser.add_argument("--cluster-validation", required=True, type=Path)
    parser.add_argument("--output-dir", required=True, type=Path)
    parser.add_argument("--chunksize", type=int, default=250_000)
    return parser.parse_args()


def read_scores(path: Path, chunksize: int) -> tuple[pd.DataFrame, pd.DataFrame]:
    columns = [
        "chip", "internal_id", "patient", "ROI", "x_um", "y_um",
        "original_coarse_niche", "nUMI", "neighbor_count",
        "neighborhood_total_umi", "feature_type", "panel", "gene",
        "compartment", "score", "detection_fraction",
    ]
    kept: list[pd.DataFrame] = []
    for chunk in pd.read_csv(
        path,
        sep="\t",
        compression="gzip",
        usecols=columns,
        chunksize=chunksize,
        low_memory=False,
    ):
        mask = (
            chunk["feature_type"].eq("gene")
            & chunk["panel"].isin(PANELS)
            & chunk["compartment"].eq("neighborhood")
        )
        if mask.any():
            kept.append(chunk.loc[mask].copy())
    if not kept:
        raise RuntimeError("No neighborhood B-panel gene rows found")
    scores = pd.concat(kept, ignore_index=True)
    metadata_columns = [
        "chip", "internal_id", "patient", "ROI", "x_um", "y_um",
        "original_coarse_niche", "nUMI", "neighbor_count",
        "neighborhood_total_umi",
    ]
    metadata = scores[metadata_columns].drop_duplicates()
    if metadata["internal_id"].duplicated().any():
        raise RuntimeError("Metadata are not unique by internal_id")
    return scores, metadata


def build_anchor_table(
    metadata: pd.DataFrame,
    assignments_path: Path,
    validation_path: Path,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    validation = pd.read_csv(validation_path, sep="\t")
    supported = validation[
        validation["major_lineage"].eq("B")
        & validation["tier"].eq("extended_only")
        & validation["tier_supported"].eq(True)
    ][["chip", "k", "cluster"]].drop_duplicates()
    if supported.empty:
        raise RuntimeError("No validated B extended-only clusters found")

    assignments = pd.read_csv(
        assignments_path,
        sep="\t",
        compression="gzip",
        usecols=["chip", "internal_id", "k", "cluster"],
    )
    assignments = assignments.merge(
        supported.assign(b_supported_resolution=True),
        on=["chip", "k", "cluster"],
        how="left",
        validate="many_to_one",
    )
    assignments["b_supported_resolution"] = assignments[
        "b_supported_resolution"
    ].eq(True)
    anchor = (
        assignments.groupby(["chip", "internal_id"], as_index=False)
        .agg(
            evaluated_resolution_n=("k", "nunique"),
            b_supported_resolution_n=("b_supported_resolution", "sum"),
        )
    )
    if not anchor["evaluated_resolution_n"].eq(3).all():
        raise RuntimeError("Every bin must have K=3/4/5 assignments")
    anchor["anchor_class"] = np.select(
        [
            anchor["b_supported_resolution_n"].ge(2),
            anchor["b_supported_resolution_n"].eq(0),
        ],
        ["positive_consensus", "background_consensus"],
        default="ambiguous_one_resolution",
    )
    result = metadata.merge(
        anchor,
        on=["chip", "internal_id"],
        validate="one_to_one",
    )
    result["umi_quintile"] = (
        result.groupby("patient", group_keys=False)["neighborhood_total_umi"]
        .transform(
            lambda values: pd.qcut(
                values.rank(method="first"),
                5,
                labels=False,
            )
            + 1
        )
        .astype(int)
    )
    return result, supported


def patient_matched_gene_metrics(
    scores: pd.DataFrame,
    anchors: pd.DataFrame,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    local = scores.merge(
        anchors[["internal_id", "anchor_class", "umi_quintile"]],
        on="internal_id",
        validate="many_to_one",
    )
    local = local[local["anchor_class"].isin(
        ["positive_consensus", "background_consensus"]
    )].copy()
    local["detected"] = local["detection_fraction"].gt(0)
    stratum = (
        local.groupby(
            [
                "chip", "patient", "umi_quintile", "panel", "gene",
                "anchor_class",
            ],
            as_index=False,
        )
        .agg(
            bin_n=("internal_id", "size"),
            mean_score=("score", "mean"),
            detection_fraction=("detected", "mean"),
        )
    )
    wide = stratum.pivot(
        index=["chip", "patient", "umi_quintile", "panel", "gene"],
        columns="anchor_class",
        values=["bin_n", "mean_score", "detection_fraction"],
    ).reset_index()
    wide.columns = [
        "__".join(str(value) for value in column if str(value) != "")
        if isinstance(column, tuple)
        else column
        for column in wide.columns
    ]
    required = [
        "bin_n__positive_consensus", "bin_n__background_consensus",
        "mean_score__positive_consensus", "mean_score__background_consensus",
        "detection_fraction__positive_consensus",
        "detection_fraction__background_consensus",
    ]
    wide = wide.dropna(subset=required).copy()
    wide["score_delta"] = (
        wide["mean_score__positive_consensus"]
        - wide["mean_score__background_consensus"]
    )
    wide["detection_delta"] = (
        wide["detection_fraction__positive_consensus"]
        - wide["detection_fraction__background_consensus"]
    )

    rows: list[dict[str, object]] = []
    for keys, frame in wide.groupby(["chip", "patient", "panel", "gene"]):
        weights = frame["bin_n__positive_consensus"].to_numpy(float)
        rows.append(
            {
                "chip": keys[0],
                "patient": keys[1],
                "panel": keys[2],
                "gene": keys[3],
                "evaluable_umi_stratum_n": len(frame),
                "matched_score_delta": np.average(frame["score_delta"], weights=weights),
                "matched_detection_delta": np.average(
                    frame["detection_delta"], weights=weights
                ),
            }
        )
    patient = pd.DataFrame(rows)
    summary = (
        patient.groupby(["panel", "gene"], as_index=False)
        .agg(
            evaluable_patient_n=("patient", "nunique"),
            positive_score_delta_patient_n=(
                "matched_score_delta", lambda values: int((values > 0).sum())
            ),
            positive_detection_delta_patient_n=(
                "matched_detection_delta", lambda values: int((values > 0).sum())
            ),
            median_matched_score_delta=("matched_score_delta", "median"),
            median_matched_detection_delta=("matched_detection_delta", "median"),
        )
    )
    chip = (
        patient.groupby(["panel", "gene", "chip"], as_index=False)
        .agg(chip_median_score_delta=("matched_score_delta", "median"))
    )
    chip_support = (
        chip.groupby(["panel", "gene"], as_index=False)
        .agg(positive_chip_n=("chip_median_score_delta", lambda x: int((x > 0).sum())))
    )
    summary = summary.merge(chip_support, on=["panel", "gene"], validate="one_to_one")
    summary["positive_patient_fraction"] = (
        summary["positive_score_delta_patient_n"]
        / summary["evaluable_patient_n"]
    )
    summary["informative_for_bin_rule"] = (
        summary["panel"].eq("B__extended_only")
        & summary["positive_patient_fraction"].ge(2 / 3)
        & summary["median_matched_score_delta"].gt(0.001)
        & summary["median_matched_detection_delta"].gt(0)
        & summary["positive_chip_n"].eq(3)
    )
    return patient, summary


def build_panel_table(
    scores: pd.DataFrame,
    anchors: pd.DataFrame,
    informative_genes: list[str],
) -> pd.DataFrame:
    if len(informative_genes) < 2:
        raise RuntimeError("Fewer than two informative extended genes survived")
    extended = scores[
        scores["panel"].eq("B__extended_only")
        & scores["gene"].isin(informative_genes)
    ].copy()
    extended["detected"] = extended["detection_fraction"].gt(0)
    panel = (
        extended.groupby("internal_id", as_index=False)
        .agg(
            informative_gene_n=("gene", "nunique"),
            informative_detected_gene_n=("detected", "sum"),
            informative_panel_score=("score", "mean"),
        )
    )
    result = anchors.merge(panel, on="internal_id", validate="one_to_one")
    expected_n = len(informative_genes)
    if not result["informative_gene_n"].eq(expected_n).all():
        raise RuntimeError("Informative gene availability differs across bins")
    return result


def evaluate_rules(panel: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    background = panel[panel["anchor_class"].eq("background_consensus")]
    target_strata = panel[
        ["chip", "patient", "umi_quintile"]
    ].drop_duplicates()
    thresholds: list[pd.DataFrame] = []
    for quantile in QUANTILES:
        patient_quintile = (
            background.groupby(["patient", "umi_quintile"], as_index=False)
            ["informative_panel_score"]
            .quantile(quantile)
            .rename(columns={"informative_panel_score": "background_score_threshold"})
        )
        patient_quintile["threshold_source"] = "patient_umi_quintile"
        patient_all = (
            background.groupby("patient", as_index=False)["informative_panel_score"]
            .quantile(quantile)
            .rename(columns={"informative_panel_score": "patient_all_threshold"})
        )
        chip_quintile = (
            background.groupby(["chip", "umi_quintile"], as_index=False)
            ["informative_panel_score"]
            .quantile(quantile)
            .rename(columns={"informative_panel_score": "chip_quintile_threshold"})
        )
        chip_all = (
            background.groupby("chip", as_index=False)["informative_panel_score"]
            .quantile(quantile)
            .rename(columns={"informative_panel_score": "chip_all_threshold"})
        )
        current = target_strata.merge(
            patient_quintile,
            on=["patient", "umi_quintile"],
            how="left",
            validate="one_to_one",
        ).merge(
            patient_all,
            on="patient",
            how="left",
            validate="many_to_one",
        ).merge(
            chip_quintile,
            on=["chip", "umi_quintile"],
            how="left",
            validate="many_to_one",
        ).merge(
            chip_all,
            on="chip",
            how="left",
            validate="many_to_one",
        )
        missing = current["background_score_threshold"].isna()
        current.loc[missing, "background_score_threshold"] = current.loc[
            missing, "patient_all_threshold"
        ]
        current.loc[missing & current["patient_all_threshold"].notna(), "threshold_source"] = "patient_all_umi"
        missing = current["background_score_threshold"].isna()
        current.loc[missing, "background_score_threshold"] = current.loc[
            missing, "chip_quintile_threshold"
        ]
        current.loc[missing & current["chip_quintile_threshold"].notna(), "threshold_source"] = "chip_umi_quintile"
        missing = current["background_score_threshold"].isna()
        current.loc[missing, "background_score_threshold"] = current.loc[
            missing, "chip_all_threshold"
        ]
        current.loc[missing & current["chip_all_threshold"].notna(), "threshold_source"] = "chip_all_umi"
        if current["background_score_threshold"].isna().any():
            raise RuntimeError("Hierarchical background threshold fallback failed")
        current = current.drop(
            columns=[
                "patient_all_threshold", "chip_quintile_threshold",
                "chip_all_threshold",
            ]
        )
        current["background_quantile"] = quantile
        thresholds.append(current)
    threshold_table = pd.concat(thresholds, ignore_index=True)

    rows: list[dict[str, object]] = []
    for quantile in QUANTILES:
        local_thresholds = threshold_table[
            threshold_table["background_quantile"].eq(quantile)
        ].drop(columns="background_quantile")
        local = panel.merge(
            local_thresholds,
            on=["chip", "patient", "umi_quintile"],
            how="left",
            validate="many_to_one",
        )
        if local["background_score_threshold"].isna().any():
            raise RuntimeError("Missing patient-depth background threshold")
        for min_gene_n in MIN_GENE_COUNTS:
            called = (
                local["informative_detected_gene_n"].ge(min_gene_n)
                & local["informative_panel_score"].gt(
                    local["background_score_threshold"]
                )
            )
            for keys, frame in local.assign(called=called).groupby(
                ["chip", "anchor_class"], sort=False
            ):
                rows.append(
                    {
                        "background_quantile": quantile,
                        "min_detected_gene_n": min_gene_n,
                        "chip": keys[0],
                        "anchor_class": keys[1],
                        "bin_n": len(frame),
                        "called_n": int(frame["called"].sum()),
                        "called_fraction": float(frame["called"].mean()),
                    }
                )
    return pd.DataFrame(rows), threshold_table


def plot_gene_deltas(summary: pd.DataFrame, path: Path) -> None:
    data = summary[summary["panel"].eq("B__extended_only")].sort_values(
        "median_matched_score_delta"
    )
    colors = np.where(data["informative_for_bin_rule"], "#1479B8", "#B8B8B8")
    fig, ax = plt.subplots(figsize=(8.5, 6.5))
    ax.barh(data["gene"], data["median_matched_score_delta"], color=colors)
    ax.axvline(0, color="black", linewidth=0.8)
    ax.set_xlabel("Median within-patient, UMI-matched score delta")
    ax.set_ylabel("Extended-only B gene")
    ax.set_title("B-anchor enrichment identifies informative extended genes")
    fig.tight_layout()
    fig.savefig(path.with_suffix(".png"), dpi=220)
    fig.savefig(path.with_suffix(".pdf"))
    plt.close(fig)


def plot_rule_grid(rule_table: pd.DataFrame, path: Path) -> None:
    aggregate = (
        rule_table.groupby(
            ["background_quantile", "min_detected_gene_n", "anchor_class"],
            as_index=False,
        )
        .agg(bin_n=("bin_n", "sum"), called_n=("called_n", "sum"))
    )
    aggregate["called_fraction"] = aggregate["called_n"] / aggregate["bin_n"]
    fig, axes = plt.subplots(1, 3, figsize=(12, 4), sharey=True)
    for ax, quantile in zip(axes, QUANTILES):
        subset = aggregate[aggregate["background_quantile"].eq(quantile)]
        for anchor_class, color in [
            ("positive_consensus", "#1479B8"),
            ("ambiguous_one_resolution", "#E69F00"),
            ("background_consensus", "#777777"),
        ]:
            line = subset[subset["anchor_class"].eq(anchor_class)]
            ax.plot(
                line["min_detected_gene_n"],
                line["called_fraction"],
                marker="o",
                label=anchor_class.replace("_", " "),
                color=color,
            )
        ax.set_title(f"Background q={quantile:.2f}")
        ax.set_xlabel("Minimum informative genes detected")
        ax.set_ylim(0, 1)
    axes[0].set_ylabel("Called fraction")
    axes[-1].legend(frameon=False, fontsize=8, loc="upper right")
    fig.tight_layout()
    fig.savefig(path.with_suffix(".png"), dpi=220)
    fig.savefig(path.with_suffix(".pdf"))
    plt.close(fig)


def main() -> None:
    args = parse_args()
    output = args.output_dir.resolve()
    for directory in [output, output / "figures", output / "qc", output / "logs"]:
        directory.mkdir(parents=True, exist_ok=True)

    scores, metadata = read_scores(args.score_file, args.chunksize)
    anchors, supported_clusters = build_anchor_table(
        metadata,
        args.assignments,
        args.cluster_validation,
    )
    patient_metrics, gene_summary = patient_matched_gene_metrics(scores, anchors)
    informative = gene_summary.loc[
        gene_summary["informative_for_bin_rule"], "gene"
    ].sort_values().tolist()
    panel = build_panel_table(scores, anchors, informative)
    rule_table, thresholds = evaluate_rules(panel)

    anchors.to_csv(output / "01_MULTIRESOLUTION_B_ANCHORS.tsv.gz", sep="\t", index=False, compression="gzip")
    supported_clusters.to_csv(output / "02_SUPPORTED_B_CLUSTERS.tsv", sep="\t", index=False)
    patient_metrics.to_csv(output / "03_PATIENT_MATCHED_GENE_METRICS.tsv", sep="\t", index=False)
    gene_summary.to_csv(output / "04_GENE_INFORMATIVENESS_AUDIT.tsv", sep="\t", index=False)
    panel.to_csv(output / "05_BIN_INFORMATIVE_PANEL_SCORES.tsv.gz", sep="\t", index=False, compression="gzip")
    thresholds.to_csv(output / "06_PATIENT_UMI_BACKGROUND_THRESHOLDS.tsv", sep="\t", index=False)
    rule_table.to_csv(output / "07_RULE_SENSITIVITY_GRID.tsv", sep="\t", index=False)

    anchor_summary = (
        anchors.groupby(["chip", "anchor_class"], as_index=False)
        .agg(bin_n=("internal_id", "size"), patient_n=("patient", "nunique"))
    )
    anchor_summary.to_csv(output / "08_ANCHOR_SUMMARY.tsv", sep="\t", index=False)

    plot_gene_deltas(gene_summary, output / "figures" / "01_EXTENDED_GENE_MATCHED_DELTAS_EN")
    plot_rule_grid(rule_table, output / "figures" / "02_RULE_SENSITIVITY_GRID_EN")

    audit = pd.DataFrame(
        {
            "audit_item": [
                "unique_bin_n", "patient_n", "chip_n", "resolution_n_per_bin",
                "disease_labels_unused", "positive_anchor_requires_ge2_resolutions",
                "background_anchor_requires_zero_resolutions", "informative_gene_n",
            ],
            "observed": [
                anchors["internal_id"].nunique(), anchors["patient"].nunique(),
                anchors["chip"].nunique(),
                int(anchors["evaluated_resolution_n"].min()), True,
                bool(anchors.loc[anchors["anchor_class"].eq("positive_consensus"), "b_supported_resolution_n"].ge(2).all()),
                bool(anchors.loc[anchors["anchor_class"].eq("background_consensus"), "b_supported_resolution_n"].eq(0).all()),
                len(informative),
            ],
            "expected": [29645, 21, 3, 3, True, True, True, ">=2"],
        }
    )
    audit["passed"] = [
        audit.loc[index, "observed"] >= 2
        if audit.loc[index, "audit_item"] == "informative_gene_n"
        else str(audit.loc[index, "observed"]) == str(audit.loc[index, "expected"])
        for index in audit.index
    ]
    audit.to_csv(output / "qc" / "09_RULE_CALIBRATION_AUDIT.tsv", sep="\t", index=False)
    if not audit["passed"].all():
        raise RuntimeError("Rule-calibration audit failed")

    session = {
        "analysis": "Bin-level B coverage rule calibration",
        "disease_labels_used": False,
        "stochastic_steps": False,
        "positive_anchor": "B-supported at >=2 of frozen K=3/4/5",
        "background_anchor": "B-supported at 0 of frozen K=3/4/5",
        "informative_genes": informative,
        "quantiles_tested": QUANTILES,
        "minimum_gene_counts_tested": MIN_GENE_COUNTS,
    }
    (output / "SESSION_INFO.json").write_text(
        json.dumps(session, indent=2), encoding="utf-8"
    )
    print(f"Rule calibration complete: {output}")
    print(f"Informative extended genes ({len(informative)}): {', '.join(informative)}")


if __name__ == "__main__":
    main()
