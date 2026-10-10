#!/usr/bin/env python3
"""Spatially validate q90/q95 calibrated B-coverage rules and freeze a mask."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.backends.backend_pdf import PdfPages


QUANTILES = (0.90, 0.95)
MIN_DETECTED_GENE_N = 2
PRIMARY_QUANTILE = 0.90
PRIMARY_MIN_REGION_COMPONENT_N = 5
SENSITIVITY_MIN_REGION_COMPONENT_N = 3


class UnionFind:
    def __init__(self, items: list[str]) -> None:
        self.parent = {item: item for item in items}
        self.size = {item: 1 for item in items}

    def find(self, item: str) -> str:
        root = item
        while self.parent[root] != root:
            root = self.parent[root]
        while self.parent[item] != item:
            next_item = self.parent[item]
            self.parent[item] = root
            item = next_item
        return root

    def union(self, first: str, second: str) -> None:
        root_first = self.find(first)
        root_second = self.find(second)
        if root_first == root_second:
            return
        if self.size[root_first] < self.size[root_second]:
            root_first, root_second = root_second, root_first
        self.parent[root_second] = root_first
        self.size[root_first] += self.size[root_second]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--panel-scores", required=True, type=Path)
    parser.add_argument("--thresholds", required=True, type=Path)
    parser.add_argument("--pilot-mask", required=True, type=Path)
    parser.add_argument("--edges", required=True, type=Path)
    parser.add_argument("--output-dir", required=True, type=Path)
    return parser.parse_args()


def add_rule_calls(panel: pd.DataFrame, thresholds: pd.DataFrame) -> pd.DataFrame:
    result = panel.copy()
    for quantile in QUANTILES:
        suffix = f"q{int(100 * quantile)}"
        local_thresholds = thresholds[
            thresholds["background_quantile"].eq(quantile)
        ][
            [
                "chip", "patient", "umi_quintile",
                "background_score_threshold", "threshold_source",
            ]
        ].rename(
            columns={
                "background_score_threshold": f"{suffix}_score_threshold",
                "threshold_source": f"{suffix}_threshold_source",
            }
        )
        result = result.merge(
            local_thresholds,
            on=["chip", "patient", "umi_quintile"],
            validate="many_to_one",
        )
        result[f"b_bin_support_{suffix}"] = (
            result["informative_detected_gene_n"].ge(MIN_DETECTED_GENE_N)
            & result["informative_panel_score"].gt(
                result[f"{suffix}_score_threshold"]
            )
        )
    return result


def prepare_edges(
    edges_path: Path,
    metadata: pd.DataFrame,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    mapping = metadata[["chip", "Bin_ID", "internal_id"]].drop_duplicates()
    if mapping.duplicated(["chip", "Bin_ID"]).any():
        raise RuntimeError("chip + Bin_ID does not uniquely identify a candidate bin")
    edges = pd.read_csv(edges_path, sep="\t", compression="gzip")
    edges = edges.merge(
        mapping.rename(
            columns={
                "Bin_ID": "neighbor_bin_id",
                "internal_id": "neighbor_internal_id",
            }
        ),
        on=["chip", "neighbor_bin_id"],
        how="left",
        validate="many_to_one",
    )
    total_edge_n = len(edges)
    outside_candidate_edge_n = int(edges["neighbor_internal_id"].isna().sum())
    edges = edges.dropna(subset=["neighbor_internal_id"]).copy()
    known = set(metadata["internal_id"])
    if not set(edges["center_internal_id"]).issubset(known):
        raise RuntimeError("Unknown center_internal_id in edge table")
    audit = pd.DataFrame(
        {
            "edge_class": ["all_center_to_neighbor", "candidate_to_candidate", "neighbor_outside_candidate_pool"],
            "edge_n": [total_edge_n, len(edges), outside_candidate_edge_n],
        }
    )
    return edges[["chip", "center_internal_id", "neighbor_internal_id"]], audit


def add_components(
    table: pd.DataFrame,
    edges: pd.DataFrame,
    call_column: str,
    suffix: str,
) -> pd.DataFrame:
    result = table.copy()
    component_ids = pd.Series(pd.NA, index=result.index, dtype="string")
    component_sizes = pd.Series(0, index=result.index, dtype="int64")
    for chip, local in result.groupby("chip", sort=False):
        called_ids = local.loc[local[call_column], "internal_id"].tolist()
        called_set = set(called_ids)
        union_find = UnionFind(called_ids)
        local_edges = edges[
            edges["chip"].eq(chip)
            & edges["center_internal_id"].isin(called_set)
            & edges["neighbor_internal_id"].isin(called_set)
        ]
        for edge in local_edges.itertuples(index=False):
            union_find.union(edge.center_internal_id, edge.neighbor_internal_id)
        roots = {item: union_find.find(item) for item in called_ids}
        root_sizes = pd.Series(list(roots.values())).value_counts().to_dict()
        local_index = result.index[result["chip"].eq(chip) & result[call_column]]
        local_ids = result.loc[local_index, "internal_id"]
        component_ids.loc[local_index] = local_ids.map(roots).astype("string")
        component_sizes.loc[local_index] = local_ids.map(
            lambda item: root_sizes[roots[item]]
        ).astype(int)
    result[f"component_id_{suffix}"] = component_ids
    result[f"component_size_{suffix}"] = component_sizes
    result[f"b_region_support_{suffix}_ge3"] = (
        result[call_column]
        & result[f"component_size_{suffix}"].ge(
            SENSITIVITY_MIN_REGION_COMPONENT_N
        )
    )
    result[f"b_region_support_{suffix}_ge5"] = (
        result[call_column]
        & result[f"component_size_{suffix}"].ge(
            PRIMARY_MIN_REGION_COMPONENT_N
        )
    )
    return result


def summarize(table: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    rows: list[dict[str, object]] = []
    patient_rows: list[dict[str, object]] = []
    for quantile in QUANTILES:
        suffix = f"q{int(100 * quantile)}"
        call = f"b_bin_support_{suffix}"
        region_ge3 = f"b_region_support_{suffix}_ge3"
        region_ge5 = f"b_region_support_{suffix}_ge5"
        for keys, frame in table.groupby(["chip", "anchor_class"], sort=False):
            rows.append(
                {
                    "background_quantile": quantile,
                    "chip": keys[0],
                    "anchor_class": keys[1],
                    "bin_n": len(frame),
                    "called_n": int(frame[call].sum()),
                    "called_fraction": float(frame[call].mean()),
                    "region_supported_ge3_n": int(frame[region_ge3].sum()),
                    "region_supported_ge3_fraction": float(frame[region_ge3].mean()),
                    "region_supported_ge5_n": int(frame[region_ge5].sum()),
                    "region_supported_ge5_fraction": float(frame[region_ge5].mean()),
                    "called_in_ge5_region_fraction": float(
                        frame.loc[frame[call], region_ge5].mean()
                    ) if frame[call].any() else np.nan,
                }
            )
        for keys, frame in table.groupby(["chip", "patient"], sort=False):
            patient_rows.append(
                {
                    "background_quantile": quantile,
                    "chip": keys[0],
                    "patient": keys[1],
                    "candidate_bin_n": len(frame),
                    "called_n": int(frame[call].sum()),
                    "called_fraction": float(frame[call].mean()),
                    "region_supported_ge3_n": int(frame[region_ge3].sum()),
                    "region_supported_ge3_fraction": float(frame[region_ge3].mean()),
                    "region_supported_ge5_n": int(frame[region_ge5].sum()),
                    "region_supported_ge5_fraction": float(frame[region_ge5].mean()),
                }
            )
    return pd.DataFrame(rows), pd.DataFrame(patient_rows)


def plot_spatial(table: pd.DataFrame, path: Path) -> None:
    with PdfPages(path) as pdf:
        for (chip, patient), frame in table.groupby(["chip", "patient"], sort=True):
            fig, axes = plt.subplots(1, 2, figsize=(11, 5))
            for ax, quantile in zip(axes, QUANTILES):
                suffix = f"q{int(100 * quantile)}"
                ax.scatter(frame["x_um"], frame["y_um"], s=2, c="#D9D9D9", rasterized=True)
                called = frame[frame[f"b_bin_support_{suffix}"]]
                ax.scatter(called["x_um"], called["y_um"], s=4, c="#F2A900", rasterized=True, label="bin support")
                region = frame[frame[f"b_region_support_{suffix}_ge5"]]
                ax.scatter(region["x_um"], region["y_um"], s=5, c="#1479B8", rasterized=True, label="component >=3")
                ax.set_title(f"{suffix}: {len(region):,} bins in components >=5")
                ax.set_aspect("equal")
                ax.invert_yaxis()
                ax.axis("off")
            axes[0].legend(frameon=False, fontsize=8, loc="best")
            fig.suptitle(f"{chip} | {patient}: calibrated B-lineage RNA coverage")
            fig.tight_layout()
            pdf.savefig(fig)
            plt.close(fig)


def plot_chip_overviews(table: pd.DataFrame, output_dir: Path) -> None:
    for chip, chip_frame in table.groupby("chip", sort=True):
        patients = sorted(chip_frame["patient"].unique())
        fig, axes = plt.subplots(
            len(patients),
            2,
            figsize=(10, 3.2 * len(patients)),
            squeeze=False,
        )
        for row_index, patient in enumerate(patients):
            frame = chip_frame[chip_frame["patient"].eq(patient)]
            for column_index, quantile in enumerate(QUANTILES):
                suffix = f"q{int(100 * quantile)}"
                ax = axes[row_index, column_index]
                ax.scatter(
                    frame["x_um"], frame["y_um"], s=1.5,
                    c="#D9D9D9", rasterized=True,
                )
                called = frame[frame[f"b_bin_support_{suffix}"]]
                ax.scatter(
                    called["x_um"], called["y_um"], s=3,
                    c="#F2A900", rasterized=True,
                )
                region = frame[frame[f"b_region_support_{suffix}_ge5"]]
                ax.scatter(
                    region["x_um"], region["y_um"], s=4,
                    c="#1479B8", rasterized=True,
                )
                ax.set_title(
                    f"{patient} | {suffix} | region n={len(region):,}",
                    fontsize=9,
                )
                ax.set_aspect("equal")
                ax.invert_yaxis()
                ax.axis("off")
        fig.suptitle(
            f"{chip}: calibrated B-lineage RNA coverage (blue = component >=5)",
            fontsize=13,
            y=0.995,
        )
        fig.tight_layout(rect=[0, 0, 1, 0.98])
        fig.savefig(
            output_dir / f"02_{chip}_Q90_Q95_SPATIAL_OVERVIEW_EN.png",
            dpi=180,
        )
        plt.close(fig)


def main() -> None:
    args = parse_args()
    output = args.output_dir.resolve()
    for directory in [output, output / "figures", output / "qc", output / "logs"]:
        directory.mkdir(parents=True, exist_ok=True)

    panel = pd.read_csv(args.panel_scores, sep="\t", compression="gzip")
    thresholds = pd.read_csv(args.thresholds, sep="\t")
    pilot = pd.read_csv(args.pilot_mask, sep="\t", compression="gzip")
    keep = [
        "internal_id", "Bin_ID", "center_core_detected_gene_n",
        "neighborhood_core_detected_gene_n", "neighborhood_b_core_corroborated",
        "neighborhood_b_core_only_supported",
    ]
    table = add_rule_calls(panel, thresholds).merge(
        pilot[keep], on="internal_id", validate="one_to_one"
    )
    edges, edge_audit = prepare_edges(args.edges, table)
    for quantile in QUANTILES:
        suffix = f"q{int(100 * quantile)}"
        table = add_components(table, edges, f"b_bin_support_{suffix}", suffix)

    table["primary_b_coverage"] = table["b_region_support_q90_ge5"]
    table["primary_b_coverage_core_corroborated"] = (
        table["primary_b_coverage"]
        & table["neighborhood_b_core_corroborated"]
    )
    table["primary_extended_gain_beyond_core_only"] = (
        table["primary_b_coverage"]
        & ~table["neighborhood_b_core_only_supported"]
    )

    spatial_summary, patient_summary = summarize(table)
    table.to_csv(output / "01_FINAL_CALIBRATED_B_COVERAGE_MASK.tsv.gz", sep="\t", index=False, compression="gzip")
    spatial_summary.to_csv(output / "02_SPATIAL_RULE_COMPARISON.tsv", sep="\t", index=False)
    patient_summary.to_csv(output / "03_PATIENT_COVERAGE_SUMMARY.tsv", sep="\t", index=False)
    edge_audit.to_csv(output / "03B_EDGE_SCOPE_AUDIT.tsv", sep="\t", index=False)

    primary_summary = (
        table.groupby("chip", as_index=False)
        .agg(
            candidate_bin_n=("internal_id", "size"),
            calibrated_bin_support_n=("b_bin_support_q90", "sum"),
            primary_region_support_n=("primary_b_coverage", "sum"),
            core_only_support_n=("neighborhood_b_core_only_supported", "sum"),
            extended_gain_n=("primary_extended_gain_beyond_core_only", "sum"),
            core_corroborated_n=("primary_b_coverage_core_corroborated", "sum"),
        )
    )
    for column in [
        "calibrated_bin_support_n", "primary_region_support_n",
        "core_only_support_n", "extended_gain_n", "core_corroborated_n",
    ]:
        primary_summary[column.replace("_n", "_fraction")] = (
            primary_summary[column] / primary_summary["candidate_bin_n"]
        )
    primary_summary.to_csv(output / "04_PRIMARY_MASK_CHIP_SUMMARY.tsv", sep="\t", index=False)
    plot_spatial(table, output / "figures" / "01_Q90_Q95_SPATIAL_COMPARISON_EN.pdf")
    plot_chip_overviews(table, output / "figures")

    audit = pd.DataFrame(
        {
            "audit_item": [
                "unique_bin_n", "candidate_to_candidate_edge_n_positive",
                "outside_candidate_neighbors_excluded_explicitly",
                "primary_quantile", "minimum_informative_gene_n",
                "primary_minimum_component_n",
                "sensitivity_minimum_component_n", "disease_labels_unused",
                "ge5_is_subset_of_ge3",
                "q90_is_superset_of_q95_bin_calls",
            ],
            "observed": [
                table["internal_id"].nunique(), bool(len(edges) > 0), True,
                PRIMARY_QUANTILE,
                MIN_DETECTED_GENE_N, PRIMARY_MIN_REGION_COMPONENT_N,
                SENSITIVITY_MIN_REGION_COMPONENT_N, True,
                bool((~table["b_region_support_q90_ge5"] | table["b_region_support_q90_ge3"]).all()),
                bool((~table["b_bin_support_q95"] | table["b_bin_support_q90"]).all()),
            ],
            "expected": [29645, True, True, 0.90, 2, 5, 3, True, True, True],
        }
    )
    audit["passed"] = audit["observed"].astype(str).eq(audit["expected"].astype(str))
    audit.to_csv(output / "qc" / "05_FINAL_MASK_AUDIT.tsv", sep="\t", index=False)
    if not audit["passed"].all():
        raise RuntimeError("Final calibrated mask audit failed")

    session = {
        "primary_definition": "q90 patient/UMI-matched score, >=2 informative genes, connected component >=5",
        "sensitivity_definitions": [
            "q90 patient/UMI-matched score, >=2 informative genes, connected component >=3",
            "q95 patient/UMI-matched score, >=2 informative genes, connected component >=5",
        ],
        "interpretation": "B-lineage RNA coverage region; not a B-cell identity or TLS call",
        "disease_labels_used": False,
        "stochastic_steps": False,
    }
    (output / "SESSION_INFO.json").write_text(json.dumps(session, indent=2), encoding="utf-8")
    print(f"Final calibrated B coverage mask complete: {output}")


if __name__ == "__main__":
    main()
