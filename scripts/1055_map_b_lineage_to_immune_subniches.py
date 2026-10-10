#!/usr/bin/env python3
"""Map the frozen B-lineage coverage mask onto frozen K3/K4/K5 subniches."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


EXPECTED_CANDIDATES = {"L3": 14839, "J2": 11893, "Y40102K8": 2913}
EXPECTED_PRIMARY_B = {"L3": 1386, "J2": 545, "Y40102K8": 738}
K_VALUES = (3, 4, 5)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--mask", required=True, type=Path)
    parser.add_argument("--assignments", required=True, type=Path)
    parser.add_argument("--output-dir", required=True, type=Path)
    return parser.parse_args()


def save_table(table: pd.DataFrame, path: Path) -> None:
    table.to_csv(path, sep="\t", index=False)


def validate_inputs(
    mask: pd.DataFrame,
    assignments: pd.DataFrame,
) -> pd.DataFrame:
    required_mask = {
        "chip", "internal_id", "patient", "ROI", "primary_b_coverage",
        "primary_b_coverage_core_corroborated",
        "primary_extended_gain_beyond_core_only", "component_id_q90",
        "component_size_q90", "informative_panel_score",
        "informative_detected_gene_n",
    }
    required_assignments = {
        "chip", "internal_id", "patient", "ROI", "k", "cluster",
        "original_coarse_niche",
    }
    if missing := required_mask.difference(mask.columns):
        raise RuntimeError(f"Mask columns missing: {sorted(missing)}")
    if missing := required_assignments.difference(assignments.columns):
        raise RuntimeError(f"Assignment columns missing: {sorted(missing)}")
    if mask["internal_id"].duplicated().any():
        raise RuntimeError("Mask internal_id is not unique")
    if assignments.duplicated(["internal_id", "k"]).any():
        raise RuntimeError("Assignment internal_id + k is not unique")
    if set(assignments["k"].unique()) != set(K_VALUES):
        raise RuntimeError("Assignments do not contain exactly K=3/4/5")

    audit_rows: list[dict[str, object]] = []
    for chip, expected_n in EXPECTED_CANDIDATES.items():
        mask_chip = mask[mask["chip"].eq(chip)]
        assign_chip = assignments[assignments["chip"].eq(chip)]
        expected_b = EXPECTED_PRIMARY_B[chip]
        checks = {
            "candidate_count": (len(mask_chip), expected_n),
            "primary_b_count": (int(mask_chip["primary_b_coverage"].sum()), expected_b),
            "patient_count": (mask_chip["patient"].nunique(), 7),
            "assignment_count": (len(assign_chip), expected_n * 3),
            "assignment_id_count": (assign_chip["internal_id"].nunique(), expected_n),
        }
        for item, (observed, expected) in checks.items():
            audit_rows.append(
                {
                    "chip": chip,
                    "audit_item": item,
                    "observed": observed,
                    "expected": expected,
                    "passed": observed == expected,
                }
            )
    audit = pd.DataFrame(audit_rows)
    if not audit["passed"].all():
        raise RuntimeError("Input audit failed")
    if set(mask["internal_id"]) != set(assignments["internal_id"]):
        raise RuntimeError("Mask and assignment IDs differ")
    return audit


def build_mapping(mask: pd.DataFrame, assignments: pd.DataFrame) -> pd.DataFrame:
    mask_columns = [
        "chip", "internal_id", "primary_b_coverage",
        "primary_b_coverage_core_corroborated",
        "primary_extended_gain_beyond_core_only", "component_id_q90",
        "component_size_q90", "informative_panel_score",
        "informative_detected_gene_n",
    ]
    mapping = assignments.merge(
        mask[mask_columns],
        on=["chip", "internal_id"],
        how="left",
        validate="many_to_one",
    )
    if mapping["primary_b_coverage"].isna().any():
        raise RuntimeError("Missing B-mask values after join")
    return mapping


def summarize_patient(mapping: pd.DataFrame) -> pd.DataFrame:
    rows: list[dict[str, object]] = []
    group_columns = ["chip", "k", "patient"]
    for (chip, k, patient), patient_frame in mapping.groupby(
        group_columns, sort=True
    ):
        for cluster, frame in patient_frame.groupby("cluster", sort=True):
            other = patient_frame[~patient_frame["cluster"].eq(cluster)]
            b_n = int(frame["primary_b_coverage"].sum())
            other_b_n = int(other["primary_b_coverage"].sum())
            fraction = b_n / len(frame)
            other_fraction = other_b_n / len(other) if len(other) else np.nan
            smoothed_fraction = (b_n + 0.5) / (len(frame) + 1)
            smoothed_other = (other_b_n + 0.5) / (len(other) + 1)
            b_frame = frame[frame["primary_b_coverage"]]
            rows.append(
                {
                    "chip": chip,
                    "k": k,
                    "patient": patient,
                    "cluster": cluster,
                    "candidate_bin_n": len(frame),
                    "b_coverage_bin_n": b_n,
                    "b_coverage_fraction": fraction,
                    "other_cluster_bin_n": len(other),
                    "other_cluster_b_bin_n": other_b_n,
                    "other_cluster_b_fraction": other_fraction,
                    "b_fraction_delta_vs_other": fraction - other_fraction,
                    "b_enrichment_ratio_vs_other": smoothed_fraction / smoothed_other,
                    "direction_enriched_vs_other": fraction > other_fraction,
                    "core_corroborated_b_n": int(
                        b_frame["primary_b_coverage_core_corroborated"].sum()
                    ),
                    "extended_gain_b_n": int(
                        b_frame["primary_extended_gain_beyond_core_only"].sum()
                    ),
                }
            )
    return pd.DataFrame(rows)


def summarize_clusters(
    mapping: pd.DataFrame,
    patient_summary: pd.DataFrame,
) -> pd.DataFrame:
    rows: list[dict[str, object]] = []
    for (chip, k, cluster), frame in mapping.groupby(
        ["chip", "k", "cluster"], sort=True
    ):
        patient = patient_summary[
            patient_summary["chip"].eq(chip)
            & patient_summary["k"].eq(k)
            & patient_summary["cluster"].eq(cluster)
        ]
        b_frame = frame[frame["primary_b_coverage"]].copy()
        b_n = len(b_frame)
        chip_b_n = int(
            mapping[
                mapping["chip"].eq(chip) & mapping["k"].eq(k)
            ]["primary_b_coverage"].sum()
        )
        patient_b_counts = b_frame.groupby("patient").size()
        largest_patient_share = (
            float(patient_b_counts.max() / b_n) if b_n else np.nan
        )
        component_counts = b_frame.dropna(subset=["component_id_q90"]).groupby(
            ["patient", "component_id_q90"]
        ).size()
        largest_component_share = (
            float(component_counts.max() / b_n) if b_n and len(component_counts) else np.nan
        )
        recurrent_n = int(patient["direction_enriched_vs_other"].sum())
        median_ratio = float(patient["b_enrichment_ratio_vs_other"].median())
        median_delta = float(patient["b_fraction_delta_vs_other"].median())
        recurrent_supported = (
            recurrent_n >= 5
            and median_ratio >= 1.25
            and median_delta > 0
            and b_n >= 25
        )
        if recurrent_supported and (
            not np.isfinite(largest_patient_share) or largest_patient_share < 0.50
        ):
            mapping_class = "recurrent_B_enriched_patient_balanced"
        elif recurrent_supported:
            mapping_class = "recurrent_B_enriched_patient_imbalanced"
        elif recurrent_n >= 5 and median_ratio > 1 and median_delta > 0:
            mapping_class = "recurrent_weak_B_enrichment"
        elif recurrent_n <= 2 and median_ratio < 1:
            mapping_class = "B_depleted_or_not_supported"
        else:
            mapping_class = "heterogeneous_or_inconclusive"
        rows.append(
            {
                "chip": chip,
                "k": k,
                "cluster": cluster,
                "candidate_bin_n": len(frame),
                "b_coverage_bin_n": b_n,
                "b_coverage_fraction": b_n / len(frame),
                "fraction_of_chip_b_coverage": b_n / chip_b_n,
                "patient_n": frame["patient"].nunique(),
                "patient_with_b_n": b_frame["patient"].nunique(),
                "patient_direction_enriched_n": recurrent_n,
                "median_patient_b_fraction": float(
                    patient["b_coverage_fraction"].median()
                ),
                "median_patient_other_b_fraction": float(
                    patient["other_cluster_b_fraction"].median()
                ),
                "median_patient_b_fraction_delta": median_delta,
                "median_patient_enrichment_ratio": median_ratio,
                "core_corroborated_b_n": int(
                    b_frame["primary_b_coverage_core_corroborated"].sum()
                ),
                "core_corroborated_fraction_of_b": float(
                    b_frame["primary_b_coverage_core_corroborated"].mean()
                ) if b_n else np.nan,
                "extended_gain_b_n": int(
                    b_frame["primary_extended_gain_beyond_core_only"].sum()
                ),
                "extended_gain_fraction_of_b": float(
                    b_frame["primary_extended_gain_beyond_core_only"].mean()
                ) if b_n else np.nan,
                "largest_patient_contribution": largest_patient_share,
                "component_n": int(len(component_counts)),
                "largest_component_contribution": largest_component_share,
                "mapping_class": mapping_class,
            }
        )
    result = pd.DataFrame(rows)
    return result.sort_values(
        ["chip", "k", "median_patient_enrichment_ratio", "b_coverage_fraction"],
        ascending=[True, True, False, False],
    )


def summarize_transitions(mapping: pd.DataFrame) -> pd.DataFrame:
    wide = mapping.pivot(
        index=["chip", "internal_id", "patient", "primary_b_coverage"],
        columns="k",
        values="cluster",
    ).reset_index()
    rows: list[dict[str, object]] = []
    for parent_k, child_k in ((3, 4), (4, 5)):
        for (chip, parent, child), frame in wide.groupby(
            ["chip", parent_k, child_k], sort=True
        ):
            parent_all = wide[
                wide["chip"].eq(chip) & wide[parent_k].eq(parent)
            ]
            child_all = wide[
                wide["chip"].eq(chip) & wide[child_k].eq(child)
            ]
            rows.append(
                {
                    "chip": chip,
                    "parent_k": parent_k,
                    "parent_cluster": parent,
                    "child_k": child_k,
                    "child_cluster": child,
                    "overlap_bin_n": len(frame),
                    "parent_fraction_to_child": len(frame) / len(parent_all),
                    "child_fraction_from_parent": len(frame) / len(child_all),
                    "b_coverage_bin_n": int(frame["primary_b_coverage"].sum()),
                    "b_coverage_fraction": float(frame["primary_b_coverage"].mean()),
                    "patient_n": frame["patient"].nunique(),
                }
            )
    return pd.DataFrame(rows).sort_values(
        ["chip", "parent_k", "parent_cluster", "overlap_bin_n"],
        ascending=[True, True, True, False],
    )


def summarize_components(mapping: pd.DataFrame) -> pd.DataFrame:
    b_mapping = mapping[mapping["primary_b_coverage"]].copy()
    rows: list[dict[str, object]] = []
    for keys, frame in b_mapping.groupby(
        ["chip", "k", "cluster", "patient", "component_id_q90"],
        dropna=False,
        sort=True,
    ):
        rows.append(
            {
                "chip": keys[0],
                "k": keys[1],
                "cluster": keys[2],
                "patient": keys[3],
                "component_id_q90": keys[4],
                "component_total_size": int(frame["component_size_q90"].max()),
                "component_bin_n_in_subniche": len(frame),
                "component_fraction_in_subniche": (
                    len(frame) / int(frame["component_size_q90"].max())
                ),
                "core_corroborated_n": int(
                    frame["primary_b_coverage_core_corroborated"].sum()
                ),
                "extended_gain_n": int(
                    frame["primary_extended_gain_beyond_core_only"].sum()
                ),
            }
        )
    return pd.DataFrame(rows)


def plot_cluster_heatmap(summary: pd.DataFrame, output_dir: Path) -> None:
    chips = ["L3", "J2", "Y40102K8"]
    fig, axes = plt.subplots(1, 3, figsize=(14, 6), constrained_layout=True)
    image = None
    for ax, chip in zip(axes, chips):
        local = summary[summary["chip"].eq(chip)]
        clusters = sorted(local["cluster"].unique())
        matrix = np.full((len(clusters), len(K_VALUES)), np.nan)
        for row_index, cluster in enumerate(clusters):
            for column_index, k in enumerate(K_VALUES):
                match = local[local["cluster"].eq(cluster) & local["k"].eq(k)]
                if len(match):
                    matrix[row_index, column_index] = match.iloc[0][
                        "median_patient_enrichment_ratio"
                    ]
        log2_matrix = np.log2(matrix)
        image = ax.imshow(log2_matrix, cmap="RdBu_r", vmin=-2, vmax=2, aspect="auto")
        ax.set_xticks(range(len(K_VALUES)), [f"K={k}" for k in K_VALUES])
        ax.set_yticks(range(len(clusters)), clusters)
        ax.set_title(chip)
        for row_index in range(matrix.shape[0]):
            for column_index in range(matrix.shape[1]):
                value = matrix[row_index, column_index]
                if np.isfinite(value):
                    ax.text(column_index, row_index, f"{value:.2f}", ha="center", va="center", fontsize=8)
    if image is not None:
        fig.colorbar(image, ax=axes, shrink=0.75, label="log2 median patient enrichment ratio")
    fig.suptitle("B-lineage coverage enrichment across frozen immune subniches")
    fig.savefig(output_dir / "figures/01_B_ENRICHMENT_BY_SUBNICHE.pdf")
    fig.savefig(output_dir / "figures/01_B_ENRICHMENT_BY_SUBNICHE.png", dpi=240)
    plt.close(fig)


def plot_patient_fractions(patient: pd.DataFrame, output_dir: Path) -> None:
    chips = ["L3", "J2", "Y40102K8"]
    fig, axes = plt.subplots(3, 3, figsize=(14, 12), constrained_layout=True)
    for row_index, chip in enumerate(chips):
        for column_index, k in enumerate(K_VALUES):
            ax = axes[row_index, column_index]
            local = patient[patient["chip"].eq(chip) & patient["k"].eq(k)].copy()
            clusters = sorted(local["cluster"].unique())
            for cluster_index, cluster in enumerate(clusters):
                values = local[local["cluster"].eq(cluster)]["b_coverage_fraction"]
                jitter = np.linspace(-0.08, 0.08, len(values)) if len(values) > 1 else np.array([0.0])
                ax.scatter(
                    np.full(len(values), cluster_index) + jitter,
                    values,
                    s=22,
                    alpha=0.8,
                )
                ax.plot(
                    [cluster_index - 0.18, cluster_index + 0.18],
                    [values.median(), values.median()],
                    color="black",
                    linewidth=2,
                )
            ax.set_xticks(range(len(clusters)), clusters)
            ax.set_ylim(bottom=0)
            ax.set_title(f"{chip}, K={k}")
            if column_index == 0:
                ax.set_ylabel("Patient B-coverage fraction")
    fig.suptitle("Patient-level B-lineage coverage by immune subniche")
    fig.savefig(output_dir / "figures/02_PATIENT_B_COVERAGE_BY_SUBNICHE.pdf")
    fig.savefig(output_dir / "figures/02_PATIENT_B_COVERAGE_BY_SUBNICHE.png", dpi=240)
    plt.close(fig)


def write_findings(
    summary: pd.DataFrame,
    transitions: pd.DataFrame,
    output_dir: Path,
) -> None:
    lines = [
        "# B-lineage覆盖与immune subniche映射：关键结果",
        "",
        "## 分析边界",
        "",
        "本轮只交叉映射冻结的29,645-bin Strategy C免疫候选池、K=3/4/5 assignments和最终2,669-bin B-lineage覆盖mask。未重新运行RCTD、邻域构建、聚类、K选择、marker筛选或疾病比较。患者是重复单位；Bin不作为独立生物学重复。",
        "",
        "## 患者重复的B富集候选",
        "",
    ]
    selected = summary[summary["mapping_class"].str.startswith("recurrent_B")]
    if selected.empty:
        lines.append("没有subniche达到预设的跨患者B富集描述规则。")
    else:
        for row in selected.itertuples(index=False):
            lines.append(
                f"- {row.chip} K={row.k} {row.cluster}: "
                f"B覆盖={row.b_coverage_bin_n}/{row.candidate_bin_n} "
                f"({row.b_coverage_fraction:.1%}); "
                f"捕获该芯片B-mask={row.fraction_of_chip_b_coverage:.1%}；"
                f"患者同向={row.patient_direction_enriched_n}/7；"
                f"患者中位富集比={row.median_patient_enrichment_ratio:.2f}；"
                f"core佐证={row.core_corroborated_fraction_of_b:.1%}；"
                f"extended gain={row.extended_gain_fraction_of_b:.1%}；"
                f"分类={row.mapping_class}。"
            )
    lines.extend(["", "## K3→K4→K5主分支", ""])
    for chip in ("L3", "J2", "Y40102K8"):
        chain_parts: list[str] = []
        for parent_k in (3, 4):
            local = transitions[
                transitions["chip"].eq(chip)
                & transitions["parent_k"].eq(parent_k)
            ]
            if local.empty:
                continue
            top = local.sort_values("b_coverage_bin_n", ascending=False).iloc[0]
            chain_parts.append(
                f"K{int(top.parent_k)} {top.parent_cluster}→"
                f"K{int(top.child_k)} {top.child_cluster} "
                f"(B={int(top.b_coverage_bin_n)}, overlap={int(top.overlap_bin_n)})"
            )
        lines.append(f"- {chip}: " + "; ".join(chain_parts) + "。")
    lines.extend(
        [
            "",
            "## 解释限制",
            "",
            "这里的分类是映射审计标签，不是最终immune subniche命名。B-lineage mask表示RNA覆盖而非纯B细胞判定；K3/K4/K5来自同一候选池且相互嵌套。最终全免疫注释仍需在T、NK、Plasma、Myeloid和Mast/Basophil覆盖扩展后联合冻结。",
        ]
    )
    (output_dir / "07_B_SUBNICHE_MAPPING_KEY_FINDINGS_ZH.md").write_text(
        "\n".join(lines), encoding="utf-8"
    )


def main() -> None:
    args = parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)
    (args.output_dir / "figures").mkdir(exist_ok=True)

    mask = pd.read_csv(args.mask, sep="\t", compression="gzip", low_memory=False)
    assignments = pd.read_csv(
        args.assignments, sep="\t", compression="gzip", low_memory=False
    )
    audit = validate_inputs(mask, assignments)
    mapping = build_mapping(mask, assignments)
    patient_summary = summarize_patient(mapping)
    cluster_summary = summarize_clusters(mapping, patient_summary)
    transitions = summarize_transitions(mapping)
    components = summarize_components(mapping)

    save_table(audit, args.output_dir / "00_INPUT_AND_ID_AUDIT.tsv")
    mapping.to_csv(
        args.output_dir / "01_BIN_K345_B_COVERAGE_MAPPING.tsv.gz",
        sep="\t",
        index=False,
        compression="gzip",
    )
    save_table(patient_summary, args.output_dir / "02_PATIENT_SUBNICHE_B_COVERAGE.tsv")
    save_table(cluster_summary, args.output_dir / "03_SUBNICHE_B_ENRICHMENT_SUMMARY.tsv")
    save_table(transitions, args.output_dir / "04_K345_LINEAGE_TRANSITIONS.tsv")
    save_table(components, args.output_dir / "05_B_COMPONENT_SUBNICHE_INTERSECTIONS.tsv")
    prioritized = cluster_summary.sort_values(
        ["mapping_class", "patient_direction_enriched_n", "median_patient_enrichment_ratio"],
        ascending=[True, False, False],
    )
    save_table(prioritized, args.output_dir / "06_B_SUBNICHE_MAPPING_PRIORITY.tsv")

    plot_cluster_heatmap(cluster_summary, args.output_dir)
    plot_patient_fractions(patient_summary, args.output_dir)
    write_findings(cluster_summary, transitions, args.output_dir)

    final_checks = pd.DataFrame(
        [
            {"audit_item": "input checks passed", "value": bool(audit["passed"].all()), "passed": bool(audit["passed"].all())},
            {"audit_item": "mapping row count", "value": len(mapping), "passed": len(mapping) == sum(EXPECTED_CANDIDATES.values()) * 3},
            {"audit_item": "no missing cluster labels", "value": int(mapping["cluster"].isna().sum()), "passed": not mapping["cluster"].isna().any()},
            {"audit_item": "primary B count reproduced at each K", "value": int(mapping["primary_b_coverage"].sum()), "passed": int(mapping["primary_b_coverage"].sum()) == sum(EXPECTED_PRIMARY_B.values()) * 3},
            {"audit_item": "patient-level summaries complete", "value": len(patient_summary), "passed": not patient_summary.empty},
            {"audit_item": "disease comparison performed", "value": False, "passed": True},
            {"audit_item": "clustering or K selection performed", "value": False, "passed": True},
        ]
    )
    save_table(final_checks, args.output_dir / "08_FINAL_AUDIT.tsv")
    if not final_checks["passed"].all():
        raise RuntimeError("Final audit failed")

    session = {
        "mask": str(args.mask),
        "assignments": str(args.assignments),
        "candidate_n": int(mask.shape[0]),
        "primary_b_n": int(mask["primary_b_coverage"].sum()),
        "mapping_row_n": int(mapping.shape[0]),
        "recurrent_b_enriched_n": int(
            cluster_summary["mapping_class"].str.startswith("recurrent_B").sum()
        ),
    }
    (args.output_dir / "SESSION_INFO.json").write_text(
        json.dumps(session, indent=2), encoding="utf-8"
    )
    print(json.dumps(session, indent=2))


if __name__ == "__main__":
    main()
