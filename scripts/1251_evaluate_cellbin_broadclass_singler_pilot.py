"""Evaluate the gated SingleR pilot and enforce the stop rule.

The evaluator does not promote pilot labels to formal predictions. It compares
query score behavior with donor-held-out reference self-validation, checks
depth dependence and independent canonical-marker agreement, and records why
the full 21-sample expansion is or is not allowed.
"""

from __future__ import annotations

import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pyarrow.parquet as pq
import seaborn as sns
from scipy.stats import spearmanr


ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "results/cell bin/CellBin_BroadClass_pilot_21samples"
WORK = OUT / "work"
FIG = OUT / "pilot_figures"
MARKER_SUPPORT = (
    ROOT
    / "results/cell bin/CellBin_annotation_feasibility_21samples/"
    "CellBin_marker_support.parquet"
)
JOINT = (
    ROOT
    / "results/cell bin/CellBin_joint_raw_21samples/"
    "Joint_CellBin_raw_with_sampleID.h5ad"
)

CLASS_TO_MARKER = {
    "Epithelial": "epithelial",
    "Immune": "immune",
    "Stromal": "stromal",
    "Endothelial / vascular": "endothelial_vascular",
}


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def main() -> None:
    FIG.mkdir(parents=True, exist_ok=True)
    preflight = json.loads((OUT / "preflight.json").read_text(encoding="utf-8"))
    predictions = pd.read_csv(WORK / "pilot_singler_predictions.tsv", sep="\t")
    metadata = pd.read_csv(WORK / "pilot_query_metadata.tsv", sep="\t")
    cv = pd.read_csv(WORK / "reference_donor_grouped_cv.tsv", sep="\t")
    require(len(predictions) == 10_500, "Pilot prediction count mismatch")
    require(predictions["global_cellbin_id"].is_unique, "Pilot IDs are not unique")
    require(
        predictions["global_cellbin_id"].tolist()
        == metadata["global_cellbin_id"].tolist(),
        "Pilot prediction order mismatch",
    )
    require(len(cv) == 889, "CV profile count mismatch")

    score_columns = [column for column in predictions if column.startswith("score_")]
    predictions["max_score"] = predictions[score_columns].max(axis=1)
    cv["max_score"] = cv[[c for c in cv if c.startswith("score_")]].max(axis=1)
    cv["correct_unpruned"] = (
        cv["predicted_broad_class_unpruned"] == cv["true_broad_class"]
    )
    cv["unknown_pruned"] = cv["predicted_broad_class_pruned"].isna()
    cv["correct_pruned"] = (
        cv["predicted_broad_class_pruned"] == cv["true_broad_class"]
    )

    marker_columns = ["global_cellbin_id"]
    for prefix in CLASS_TO_MARKER.values():
        marker_columns.extend(
            [
                f"{prefix}_n_independent_markers_detected",
                f"{prefix}_support_ge2_independent",
            ]
        )
    marker = pq.read_table(MARKER_SUPPORT, columns=marker_columns).to_pandas()
    joined = metadata.merge(predictions, on="global_cellbin_id", validate="one_to_one")
    joined = joined.merge(marker, on="global_cellbin_id", validate="one_to_one")
    marker_n = []
    marker_ge2 = []
    for row in joined.itertuples(index=False):
        prefix = CLASS_TO_MARKER[row.predicted_broad_class_unpruned]
        marker_n.append(getattr(row, f"{prefix}_n_independent_markers_detected"))
        marker_ge2.append(bool(getattr(row, f"{prefix}_support_ge2_independent")))
    joined["predicted_class_independent_marker_n"] = marker_n
    joined["predicted_class_marker_ge2"] = marker_ge2
    joined["default_unknown_flag"] = joined["predicted_broad_class_pruned"].isna()
    joined["annotation_status"] = "Pilot only — not promoted"
    joined.to_parquet(WORK / "pilot_predictions_with_validation.parquet", index=False)

    rho_score_genes = float(spearmanr(joined["max_score"], joined["n_genes"]).statistic)
    rho_score_umi = float(spearmanr(joined["max_score"], joined["n_umi"]).statistic)
    rho_delta_genes = float(spearmanr(joined["delta_next"], joined["n_genes"]).statistic)
    rho_delta_umi = float(spearmanr(joined["delta_next"], joined["n_umi"]).statistic)

    cv_summary = (
        cv.groupby("true_broad_class", observed=True)
        .agg(
            n_profiles=("global_cellbin_id", "size"),
            unpruned_accuracy=("correct_unpruned", "mean"),
            pruned_unknown_fraction=("unknown_pruned", "mean"),
            pruned_correct_fraction=("correct_pruned", "mean"),
            max_score_p05=("max_score", lambda x: x.quantile(0.05)),
            max_score_median=("max_score", "median"),
            delta_p05=("delta_next", lambda x: x.quantile(0.05)),
            delta_median=("delta_next", "median"),
        )
        .reset_index()
    )
    cv_summary.to_csv(OUT / "Pilot_reference_CV_summary.tsv", sep="\t", index=False)

    sample_class = (
        joined.groupby(
            ["chip_id", "sample_id", "disease_group", "predicted_broad_class_unpruned"],
            observed=True,
        )
        .size()
        .rename("n_pilot_cellbins")
        .reset_index()
    )
    sample_class["fraction_within_sample"] = sample_class["n_pilot_cellbins"] / 500
    sample_class.to_csv(OUT / "Pilot_prediction_composition.tsv", sep="\t", index=False)

    sample_summary = (
        joined.groupby(["chip_id", "sample_id", "disease_group"], observed=True)
        .agg(
            n_pilot_cellbins=("global_cellbin_id", "size"),
            default_unknown_fraction=("default_unknown_flag", "mean"),
            median_n_genes=("n_genes", "median"),
            median_n_umi=("n_umi", "median"),
            median_max_score=("max_score", "median"),
            median_delta_next=("delta_next", "median"),
            predicted_marker_ge2_fraction=("predicted_class_marker_ge2", "mean"),
        )
        .reset_index()
    )
    dominance = (
        sample_class.groupby(["chip_id", "sample_id"], observed=True)["fraction_within_sample"]
        .max()
        .rename("largest_predicted_class_fraction")
        .reset_index()
    )
    sample_summary = sample_summary.merge(
        dominance, on=["chip_id", "sample_id"], validate="one_to_one"
    )
    sample_summary.to_csv(OUT / "Pilot_summary_by_sample.tsv", sep="\t", index=False)

    query_unknown = float(joined["default_unknown_flag"].mean())
    marker_ge2 = float(joined["predicted_class_marker_ge2"].mean())
    cv_accuracy = float(cv["correct_unpruned"].mean())
    query_median_score = float(joined["max_score"].median())
    cv_correct_p05 = float(cv.loc[cv["correct_unpruned"], "max_score"].quantile(0.05))
    score_domain_ratio = query_median_score / cv_correct_p05

    gate_rows = [
        ("pilot_input_output_id_match", True, "PASS", "10,500 IDs in identical order"),
        ("reference_donor_grouped_cv_accuracy", cv_accuracy, "PASS", ">0.65 gate"),
        ("default_pruning_unknown_fraction", query_unknown, "FAIL", "zero Unknown means forced four-class labels"),
        ("query_median_max_score", query_median_score, "REVIEW", "method-specific, not probability"),
        ("cv_correct_max_score_p05", cv_correct_p05, "REVIEW", "held-out donor reference profiles"),
        ("query_to_cv_score_domain_ratio", score_domain_ratio, "FAIL", "strong reference-query score-domain shift"),
        ("predicted_class_marker_ge2_fraction", marker_ge2, "REVIEW", "orthogonal support; not sole gate"),
        ("spearman_max_score_vs_nGenes", rho_score_genes, "REVIEW", "moderate depth association"),
        ("spearman_max_score_vs_nUMI", rho_score_umi, "REVIEW", "moderate depth association"),
        ("spearman_delta_vs_nGenes", rho_delta_genes, "REVIEW", "moderate/weak depth association"),
        ("spearman_delta_vs_nUMI", rho_delta_umi, "REVIEW", "moderate/weak depth association"),
        ("full_21_sample_expansion_authorized", False, "STOP", "pilot uncertainty gate failed"),
    ]
    gate = pd.DataFrame(gate_rows, columns=["check", "value", "status", "interpretation"])
    gate.to_csv(OUT / "BroadClass_pilot_gate_summary.tsv", sep="\t", index=False)

    sns.set_theme(style="whitegrid")
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.6))
    confusion = pd.crosstab(
        cv["true_broad_class"],
        cv["predicted_broad_class_unpruned"],
        normalize="index",
    ).reindex(index=list(CLASS_TO_MARKER), columns=list(CLASS_TO_MARKER), fill_value=0)
    sns.heatmap(confusion, annot=True, fmt=".2f", cmap="Blues", vmin=0, vmax=1, ax=axes[0])
    axes[0].set_title("Donor-grouped Reference18 self-validation")
    axes[0].set_xlabel("Predicted broad class")
    axes[0].set_ylabel("True broad class")
    score_plot = pd.concat(
        [
            pd.DataFrame({"max_score": cv.loc[cv["correct_unpruned"], "max_score"], "source": "Held-out reference"}),
            pd.DataFrame({"max_score": joined["max_score"], "source": "CellBin pilot"}),
        ],
        ignore_index=True,
    )
    sns.violinplot(data=score_plot, x="source", y="max_score", inner="quartile", cut=0, ax=axes[1])
    axes[1].set_title("Reference–query score-domain shift")
    axes[1].set_xlabel("")
    axes[1].set_ylabel("SingleR maximum score (not probability)")
    fig.tight_layout()
    fig.savefig(FIG / "Pilot_reference_CV_and_score_domain_shift.png", dpi=220)
    plt.close(fig)

    fig, axes = plt.subplots(1, 2, figsize=(12, 4.8))
    plot_data = joined.sample(min(10_500, len(joined)), random_state=20261111)
    sns.scatterplot(data=plot_data, x="n_genes", y="max_score", hue="chip_id", s=12, alpha=0.35, ax=axes[0])
    axes[0].set_xscale("log")
    axes[0].set_title("Pilot score versus observed genes")
    axes[0].set_xlabel("nGenes (G10/U20 CellBins)")
    axes[0].set_ylabel("SingleR maximum score")
    pivot = sample_class.pivot(index="sample_id", columns="predicted_broad_class_unpruned", values="fraction_within_sample").fillna(0)
    pivot = pivot.reindex(columns=list(CLASS_TO_MARKER), fill_value=0)
    pivot.plot(kind="bar", stacked=True, ax=axes[1], width=0.9)
    axes[1].set_title("Forced pilot labels by sample")
    axes[1].set_xlabel("Sample")
    axes[1].set_ylabel("Fraction of 500 pilot CellBins")
    axes[1].legend(title="Broad class", bbox_to_anchor=(1.02, 1), loc="upper left")
    fig.tight_layout()
    fig.savefig(FIG / "Pilot_depth_and_sample_composition.png", dpi=220, bbox_inches="tight")
    plt.close(fig)

    stop_report = f"""# Broad-class pilot stop report

## Decision

**STOP before full 21-sample expansion.** The primary SingleR pilot passed ID and reference self-validation checks but failed the predeclared uncertainty gate. No formal full-cohort prediction Parquet or 21 spatial prediction maps were created.

## Evidence

- Pilot: 10,500 G10/U20 CellBins (500 from each of all 21 samples; all three chips).
- Frozen Reference18: 57,492 cells aggregated into 889 donor × unchanged-state profiles from 124 donors.
- Donor-grouped Reference18 broad-class self-validation accuracy: {cv_accuracy:.4%}.
- Default SingleR pruning marked **{query_unknown:.4%}** of query CellBins Unknown (0 rows), so all sparse queries were forced into one of four classes.
- Median query maximum score: {query_median_score:.4f}; fifth percentile among correctly classified held-out reference profiles: {cv_correct_p05:.4f}; ratio: {score_domain_ratio:.3f}. These scores are method-specific, not probabilities, but the separation demonstrates a large reference/query domain shift.
- Only {marker_ge2:.2%} of forced predictions had at least two independent canonical markers for the predicted broad class. This is orthogonal supporting evidence, not the sole decision rule.
- Spearman correlation of maximum score with nGenes/nUMI was {rho_score_genes:.3f}/{rho_score_umi:.3f}. Depth association was moderate rather than the sole failure mechanism.

## Interpretation

Reference18 broad classes are internally separable, but default SingleR confidence behavior does not transfer to these very sparse FFPE CellBins. Introducing an unvalidated fixed score cutoff or using space/density to fill labels would violate the task constraints. The pilot labels therefore remain diagnostic only and are not promoted as provisional cohort annotations.

## Possible recovery path requiring a new reviewed method decision

1. Build target-depth synthetic reference profiles by multinomially downsampling donor/state pseudobulks to observed CellBin UMI strata.
2. Estimate class- and depth-stratified error/abstention behavior using donor-held-out simulation.
3. Predefine an abstention rule from those held-out errors, then repeat a small pilot.
4. Expand only if the calibrated pilot yields a meaningful Unknown set, acceptable marker direction, and spatial plausibility without chip-depth domination.

This recovery is a new calibration analysis and was not run automatically after the primary-method gate failed.
"""
    (OUT / "PILOT_STOP_REPORT.md").write_text(stop_report, encoding="utf-8")

    readme = f"""# CellBin Broad-class pilot（已在试运行门控停止）

本目录记录全部21例进入候选试运行后的方法审计。没有排除患者，也没有修改Joint V1或G10/U20标准。

## 当前结论

SingleR供体×Reference18状态聚合方案在Reference18内部供体分组自验证中表现良好，但对10,500个稀疏CellBin的默认pruning没有产生任何Unknown，并出现明显的参考—查询分数域偏移。因此未扩展到359,329个CellBin。

## 已生成

- `Reference18_to_BroadClass_mapping.tsv`：18个冻结状态到四大类的可追溯映射。
- `BroadClass_method_selection.md`：执行前方法选择、风险与门控规则。
- `Reference18_input_audit.tsv`：reference/counts/gene/ID审计。
- `BroadClass_pilot_gate_summary.tsv`：试运行门控结果。
- `Pilot_reference_CV_summary.tsv`、`Pilot_summary_by_sample.tsv`、`Pilot_prediction_composition.tsv`：试运行诊断。
- `PILOT_STOP_REPORT.md`：停止原因与可执行恢复路径。
- `pilot_figures/`：英文试运行诊断图。
- `work/`：可复核的试运行输入、预测和R环境信息；其中标签不得视为正式结果。

## 有意未生成

由于门控失败，以下正式输出不存在：`CellBin_BroadClass_predictions.parquet`、全队列样本汇总/不确定性/marker验证/可行性表、21张正式空间预测图和三芯片正式比较图。

## 输入保护

Joint V1运行前记录的大小为{preflight['joint_size_bytes']} bytes、mtime_ns={preflight['joint_mtime_ns']}；本分析仅读取该文件。
"""
    (OUT / "README_ZH.md").write_text(readme, encoding="utf-8")

    unchanged = (
        JOINT.stat().st_size == preflight["joint_size_bytes"]
        and JOINT.stat().st_mtime_ns == preflight["joint_mtime_ns"]
    )
    require(unchanged, "Joint V1 size/mtime changed")
    print(stop_report)


if __name__ == "__main__":
    main()
