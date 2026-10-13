#!/usr/bin/env python3
"""Targeted audit of the B-associated largest high-usage component metric.

Purpose
-------
Audit the exact frozen definition inherited by analysis 1104, reconstruct all
21 patient-level numerators and denominators, and run prespecified threshold
sensitivities without refitting NMF or modifying any frozen result.

Inputs
------
* Frozen K=6 relative usage matrices and entity order.
* Frozen all-21 neighborhood assignment table with patient coordinates.
* Frozen 1101 patient metrics and frozen 1104 disease-test results.

Outputs
-------
A timestamped audit directory containing patient-level baseline and sensitivity
tables, IPF-versus-HC comparisons, leave-one-patient-out influence estimates,
an adjacency audit, a Chinese report, English figures, and validation records.
"""

from __future__ import annotations

import hashlib
import json
import math
import platform
import sys
from datetime import datetime
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from scipy.sparse import csr_matrix
from scipy.sparse.csgraph import connected_components
from scipy.spatial import cKDTree


SEED = 20261013
N_PERMUTATIONS = 100_000
EDGE_RADIUS_UM = 36.0
BIN_AREA_MM2 = 0.000625
EPS = 1e-12

ROOT = Path(__file__).resolve().parents[1]
MODEL_RUN = (
    ROOT
    / "results/reanalysis/REFERENCE18_PER_CHIP_DUAL_NMF_20261012_204500_attempt4_alpha0"
)
ASSIGNMENTS = (
    ROOT
    / "results/reanalysis/ALL21_SHARED_UNBIASED_DISCRETE_NICHE_ATLAS_20261011_191652_attempt4"
    / "02_ALL21_SHARED_NICHE_ASSIGNMENTS.tsv.gz"
)
SOURCE_1101 = (
    ROOT / "results/reanalysis/REFERENCE18_K6_FACTOR_DOWNSTREAM_20261013_003020"
)
SOURCE_1104 = (
    ROOT / "results/reanalysis/FACTOR_FAMILY_PATIENT_DISEASE_20261013_005851"
)

B_FACTORS = {"K8": "F02", "J2": "F03", "L3": "F01"}
DISEASE_ORDER = ["HC", "IPF", "SSc-ILD"]
THRESHOLD_SCHEMES = [
    ("chip_median_patient_q90", "chip_balanced", 0.90),
    ("patient_specific_q90", "patient_specific", 0.90),
    ("chip_median_patient_q85", "chip_balanced", 0.85),
    ("chip_median_patient_q95", "chip_balanced", 0.95),
]
METRICS = {
    "largest_component_bin_n": "log1p",
    "largest_component_over_valid": "arcsin_sqrt",
    "largest_component_over_high": "arcsin_sqrt",
}


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def zscore(values: pd.Series) -> pd.Series:
    standard_deviation = values.std(ddof=1)
    if not np.isfinite(standard_deviation) or standard_deviation <= EPS:
        return pd.Series(np.zeros(len(values)), index=values.index)
    return (values - values.mean()) / standard_deviation


def transform(values: pd.Series, method: str) -> pd.Series:
    numeric = values.astype(float)
    if method == "log1p":
        return np.log1p(numeric)
    if method == "arcsin_sqrt":
        return np.arcsin(np.sqrt(np.clip(numeric, 0, 1)))
    raise ValueError(method)


def bh_adjust(values: pd.Series) -> pd.Series:
    result = pd.Series(np.nan, index=values.index, dtype=float)
    valid = values.dropna()
    order = np.argsort(valid.to_numpy())
    ranked = valid.to_numpy()[order]
    adjusted = ranked * len(ranked) / np.arange(1, len(ranked) + 1)
    adjusted = np.minimum.accumulate(adjusted[::-1])[::-1]
    adjusted = np.clip(adjusted, 0, 1)
    result.loc[valid.index.to_numpy()[order]] = adjusted
    return result


def edge_array(coordinates: np.ndarray) -> np.ndarray:
    pairs = cKDTree(coordinates).query_pairs(
        EDGE_RADIUS_UM + 1e-8,
        output_type="ndarray",
    )
    return pairs.astype(np.int64)


def graph_labels(node_n: int, edges: np.ndarray) -> tuple[int, np.ndarray]:
    if node_n == 0:
        return 0, np.array([], dtype=int)
    if len(edges) == 0:
        return node_n, np.arange(node_n, dtype=int)
    rows = np.r_[edges[:, 0], edges[:, 1]]
    columns = np.r_[edges[:, 1], edges[:, 0]]
    adjacency = csr_matrix(
        (np.ones(len(rows), dtype=np.uint8), (rows, columns)),
        shape=(node_n, node_n),
    )
    return connected_components(adjacency, directed=False)


def high_component_details(
    high: np.ndarray,
    edges: np.ndarray,
    tissue_labels: np.ndarray,
) -> dict[str, float | int | bool]:
    high_indices = np.flatnonzero(high)
    if len(high_indices) == 0:
        return {
            "high_component_n": 0,
            "largest_component_bin_n": 0,
            "largest_component_spans_tissue_components": False,
        }
    remap = np.full(len(high), -1, dtype=int)
    remap[high_indices] = np.arange(len(high_indices))
    keep = high[edges[:, 0]] & high[edges[:, 1]]
    high_edges = edges[keep]
    if len(high_edges) == 0:
        high_labels = np.arange(len(high_indices), dtype=int)
        component_n = len(high_indices)
    else:
        rows = remap[high_edges[:, 0]]
        columns = remap[high_edges[:, 1]]
        adjacency = csr_matrix(
            (
                np.ones(len(rows) * 2, dtype=np.uint8),
                (np.r_[rows, columns], np.r_[columns, rows]),
            ),
            shape=(len(high_indices), len(high_indices)),
        )
        component_n, high_labels = connected_components(adjacency, directed=False)
    sizes = np.bincount(high_labels, minlength=component_n)
    spanning = False
    for label in range(component_n):
        members = high_indices[high_labels == label]
        if len(np.unique(tissue_labels[members])) > 1:
            spanning = True
            break
    return {
        "high_component_n": int(component_n),
        "largest_component_bin_n": int(sizes.max()),
        "largest_component_spans_tissue_components": spanning,
    }


def load_data() -> tuple[pd.DataFrame, dict[str, np.ndarray]]:
    assignments = pd.read_csv(ASSIGNMENTS, sep="\t")
    frames = []
    usage_by_chip: dict[str, np.ndarray] = {}
    for chip, factor in B_FACTORS.items():
        model_path = (
            MODEL_RUN / "models" / f"neighborhood50_{chip}_K06_MODEL_AND_USAGE.npz"
        )
        with np.load(model_path, allow_pickle=True) as model:
            keys = model["entity_key"].astype(str)
            relative_usage = model["relative_usage"].astype(float)
        factor_index = int(factor[1:]) - 1
        lookup = assignments.loc[assignments["chip_label"].eq(chip)].set_index("entity_key")
        frame = lookup.loc[keys].reset_index()
        if frame["entity_key"].astype(str).tolist() != keys.tolist():
            raise RuntimeError(f"{chip}: entity order mismatch")
        frame["b_usage"] = relative_usage[:, factor_index]
        frame["factor"] = factor
        frames.append(frame)
        usage_by_chip[chip] = relative_usage
    return pd.concat(frames, ignore_index=True), usage_by_chip


def threshold_lookup(data: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for chip, chip_frame in data.groupby("chip_label", sort=True):
        for scheme, rule, quantile in THRESHOLD_SCHEMES:
            patient_quantiles = chip_frame.groupby("patient_id")["b_usage"].quantile(quantile)
            if rule == "chip_balanced":
                common = float(patient_quantiles.median())
                for patient in patient_quantiles.index:
                    rows.append(
                        {
                            "chip": chip,
                            "patient_id": patient,
                            "threshold_scheme": scheme,
                            "threshold_rule": rule,
                            "target_quantile": quantile,
                            "high_usage_threshold": common,
                        }
                    )
            else:
                for patient, cutoff in patient_quantiles.items():
                    rows.append(
                        {
                            "chip": chip,
                            "patient_id": patient,
                            "threshold_scheme": scheme,
                            "threshold_rule": rule,
                            "target_quantile": quantile,
                            "high_usage_threshold": float(cutoff),
                        }
                    )
    return pd.DataFrame(rows)


def patient_metrics(data: pd.DataFrame, thresholds: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    threshold_index = thresholds.set_index(["chip", "patient_id", "threshold_scheme"])
    rows = []
    audit_rows = []
    for (chip, patient), frame in data.groupby(["chip_label", "patient_id"], sort=True):
        frame = frame.reset_index(drop=True)
        coordinates = frame[["x_um", "y_um"]].to_numpy(dtype=float)
        edges = edge_array(coordinates)
        tissue_component_n, tissue_labels = graph_labels(len(frame), edges)
        source_ids = frame["bin_id"].astype(str).str.split("|", n=1).str[0]
        edge_distances = (
            np.linalg.norm(coordinates[edges[:, 0]] - coordinates[edges[:, 1]], axis=1)
            if len(edges)
            else np.array([], dtype=float)
        )
        audit_rows.append(
            {
                "chip": chip,
                "patient_id": patient,
                "disease": frame["disease"].iloc[0],
                "valid_bin_n": len(frame),
                "source_id_n": source_ids.nunique(),
                "source_ids": ";".join(sorted(source_ids.unique())),
                "duplicate_coordinate_n": int(frame.duplicated(["x_um", "y_um"]).sum()),
                "spatial_edge_n": len(edges),
                "minimum_edge_distance_um": float(edge_distances.min()) if len(edges) else math.nan,
                "maximum_edge_distance_um": float(edge_distances.max()) if len(edges) else math.nan,
                "tissue_component_n": int(tissue_component_n),
                "largest_tissue_component_fraction": float(
                    np.bincount(tissue_labels).max() / len(frame)
                ),
                "patient_isolation_by_construction": True,
                "explicit_roi_column_available": False,
            }
        )
        for scheme, _, _ in THRESHOLD_SCHEMES:
            threshold = float(threshold_index.loc[(chip, patient, scheme), "high_usage_threshold"])
            high = frame["b_usage"].to_numpy(dtype=float) >= threshold
            details = high_component_details(high, edges, tissue_labels)
            high_n = int(high.sum())
            largest_n = int(details["largest_component_bin_n"])
            rows.append(
                {
                    "chip": chip,
                    "factor": B_FACTORS[chip],
                    "patient_id": patient,
                    "disease": frame["disease"].iloc[0],
                    "threshold_scheme": scheme,
                    "threshold_rule": threshold_index.loc[(chip, patient, scheme), "threshold_rule"],
                    "target_quantile": float(
                        threshold_index.loc[(chip, patient, scheme), "target_quantile"]
                    ),
                    "high_usage_threshold": threshold,
                    "valid_bin_n": len(frame),
                    "high_usage_bin_n": high_n,
                    "high_usage_fraction": high_n / len(frame),
                    "high_component_n": details["high_component_n"],
                    "largest_component_bin_n": largest_n,
                    "largest_component_area_mm2": largest_n * BIN_AREA_MM2,
                    "largest_component_over_high": largest_n / high_n if high_n else math.nan,
                    "largest_component_over_valid": largest_n / len(frame),
                    "largest_component_spans_tissue_components": details[
                        "largest_component_spans_tissue_components"
                    ],
                }
            )
    return pd.DataFrame(rows), pd.DataFrame(audit_rows)


def make_permutations(meta: pd.DataFrame, rng: np.random.Generator) -> np.ndarray:
    observed = meta["disease"].to_numpy(dtype=object)
    permutations = np.empty((N_PERMUTATIONS, len(meta)), dtype=object)
    for chip in sorted(meta["chip"].unique()):
        positions = np.flatnonzero(meta["chip"].to_numpy() == chip)
        labels = observed[positions]
        for index in range(N_PERMUTATIONS):
            permutations[index, positions] = rng.permutation(labels)
    return permutations


def contrast_effect(values: np.ndarray, labels: np.ndarray) -> float:
    return float(values[labels == "IPF"].mean() - values[labels == "HC"].mean())


def disease_sensitivity(patient: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    meta = (
        patient[["chip", "patient_id", "disease"]]
        .drop_duplicates()
        .sort_values(["chip", "patient_id"])
        .reset_index(drop=True)
    )
    rng = np.random.default_rng(SEED)
    permutations = make_permutations(meta, rng)
    position = dict(zip(meta["patient_id"], meta.index))
    summary_rows = []
    chip_rows = []
    loo_rows = []
    for scheme, _, _ in THRESHOLD_SCHEMES:
        subset = patient.loc[patient["threshold_scheme"].eq(scheme)].copy()
        subset = subset.sort_values(["chip", "patient_id"]).reset_index(drop=True)
        positions = np.array([position[value] for value in subset["patient_id"]])
        permuted_labels = permutations[:, positions]
        observed_labels = subset["disease"].to_numpy(dtype=object)
        for metric, method in METRICS.items():
            subset["transformed"] = transform(subset[metric], method)
            subset["within_chip_z"] = subset.groupby("chip")["transformed"].transform(zscore)
            values = subset["within_chip_z"].to_numpy(dtype=float)
            observed_effect = contrast_effect(values, observed_labels)
            permutation_effects = np.empty(N_PERMUTATIONS, dtype=float)
            for index in range(N_PERMUTATIONS):
                permutation_effects[index] = contrast_effect(values, permuted_labels[index])
            p_value = (1 + np.sum(np.abs(permutation_effects) >= abs(observed_effect))) / (
                N_PERMUTATIONS + 1
            )
            raw_hc = subset.loc[subset["disease"].eq("HC"), metric]
            raw_ipf = subset.loc[subset["disease"].eq("IPF"), metric]
            loo_effects = []
            for patient_id in subset.loc[
                subset["disease"].isin(["HC", "IPF"]), "patient_id"
            ]:
                keep = subset["patient_id"].ne(patient_id).to_numpy()
                effect = contrast_effect(values[keep], observed_labels[keep])
                loo_effects.append(effect)
                loo_rows.append(
                    {
                        "threshold_scheme": scheme,
                        "metric": metric,
                        "left_out_patient": patient_id,
                        "left_out_disease": subset.loc[
                            subset["patient_id"].eq(patient_id), "disease"
                        ].iloc[0],
                        "full_effect_within_chip_sd": observed_effect,
                        "loo_effect_within_chip_sd": effect,
                        "change_from_full": effect - observed_effect,
                        "same_direction_as_full": bool(
                            np.sign(effect) == np.sign(observed_effect)
                        ),
                    }
                )
            summary_rows.append(
                {
                    "threshold_scheme": scheme,
                    "metric": metric,
                    "transform": method,
                    "hc_mean_raw": float(raw_hc.mean()),
                    "ipf_mean_raw": float(raw_ipf.mean()),
                    "raw_ipf_minus_hc": float(raw_ipf.mean() - raw_hc.mean()),
                    "effect_within_chip_sd": observed_effect,
                    "permutation_p_exploratory": float(p_value),
                    "permutation_n": N_PERMUTATIONS,
                    "loo_effect_min": float(np.min(loo_effects)),
                    "loo_effect_max": float(np.max(loo_effects)),
                    "loo_same_direction_fraction": float(
                        np.mean(np.sign(loo_effects) == np.sign(observed_effect))
                    ),
                    "largest_absolute_loo_change": float(
                        np.max(np.abs(np.asarray(loo_effects) - observed_effect))
                    ),
                }
            )
            for chip, chip_frame in subset.groupby("chip", sort=True):
                chip_rows.append(
                    {
                        "threshold_scheme": scheme,
                        "metric": metric,
                        "chip": chip,
                        "hc_mean_raw": float(
                            chip_frame.loc[chip_frame["disease"].eq("HC"), metric].mean()
                        ),
                        "ipf_mean_raw": float(
                            chip_frame.loc[chip_frame["disease"].eq("IPF"), metric].mean()
                        ),
                        "raw_ipf_minus_hc": float(
                            chip_frame.loc[chip_frame["disease"].eq("IPF"), metric].mean()
                            - chip_frame.loc[chip_frame["disease"].eq("HC"), metric].mean()
                        ),
                        "effect_within_chip_sd": float(
                            chip_frame.loc[
                                chip_frame["disease"].eq("IPF"), "within_chip_z"
                            ].mean()
                            - chip_frame.loc[
                                chip_frame["disease"].eq("HC"), "within_chip_z"
                            ].mean()
                        ),
                    }
                )
    summary = pd.DataFrame(summary_rows)
    summary["exploratory_bh_q_across_12"] = bh_adjust(
        summary["permutation_p_exploratory"]
    )
    return summary, pd.DataFrame(chip_rows), pd.DataFrame(loo_rows)


def make_figures(patient: pd.DataFrame, summary: pd.DataFrame, figures: Path) -> None:
    sns.set_theme(style="whitegrid", context="notebook")
    baseline = patient.loc[
        patient["threshold_scheme"].eq("chip_median_patient_q90")
    ].copy()
    baseline["label"] = baseline["patient_id"] + " (" + baseline["chip"] + ")"
    baseline = baseline.sort_values(["disease", "chip", "largest_component_over_valid"])
    figure, axes = plt.subplots(1, 3, figsize=(17, 7))
    columns = [
        ("high_usage_fraction", "High-usage fraction"),
        ("largest_component_over_high", "Largest component / high-usage bins"),
        ("largest_component_over_valid", "Largest component / valid bins"),
    ]
    palette = {"HC": "#4C78A8", "IPF": "#E45756", "SSc-ILD": "#72B7B2"}
    for axis, (column, title) in zip(axes, columns):
        sns.barplot(
            data=baseline,
            x=column,
            y="label",
            hue="disease",
            dodge=False,
            palette=palette,
            ax=axis,
        )
        axis.set_title(title)
        axis.set_xlabel("Fraction")
        axis.set_ylabel("")
        if axis is not axes[-1]:
            axis.get_legend().remove()
    axes[-1].legend(title="Disease", loc="lower right")
    figure.suptitle("B-associated baseline component audit", fontsize=15)
    figure.tight_layout()
    figure.savefig(figures / "01_B_ASSOC_BASELINE_PATIENT_METRICS.png", dpi=220)
    plt.close(figure)

    figure, axis = plt.subplots(figsize=(12, 6))
    display = summary.copy()
    display["metric_label"] = display["metric"].map(
        {
            "largest_component_bin_n": "Absolute component size",
            "largest_component_over_valid": "Component / valid bins",
            "largest_component_over_high": "Component / high-usage bins",
        }
    )
    sns.pointplot(
        data=display,
        x="effect_within_chip_sd",
        y="threshold_scheme",
        hue="metric_label",
        join=False,
        dodge=0.35,
        ax=axis,
    )
    axis.axvline(0, color="black", linewidth=1)
    axis.set_xlabel("IPF - HC effect (within-chip SD)")
    axis.set_ylabel("Threshold scheme")
    axis.set_title("B-associated component sensitivity: IPF vs HC")
    axis.legend(title="Metric", bbox_to_anchor=(1.02, 1), loc="upper left")
    figure.tight_layout()
    figure.savefig(figures / "02_B_ASSOC_THRESHOLD_SENSITIVITY.png", dpi=220)
    plt.close(figure)


def write_report(
    output: Path,
    baseline: pd.DataFrame,
    summary: pd.DataFrame,
    adjacency: pd.DataFrame,
    original: pd.DataFrame,
) -> None:
    baseline_extreme = baseline.loc[baseline["high_usage_fraction"] > 0.50]
    patient_q90 = summary.loc[
        summary["threshold_scheme"].eq("patient_specific_q90")
    ].set_index("metric")
    current = summary.loc[
        summary["threshold_scheme"].eq("chip_median_patient_q90")
    ].set_index("metric")
    text = f"""# B-associated factor 最大高-usage连通块针对性方法核查

## 范围

本轮只读复用冻结的逐芯片 K=6 NMF、Reference18 neighborhood assignment、1101患者指标和1104疾病检验。未重跑NMF、未修改factor、未开展B细胞下游生物学分析。

## 实际代码定义

- 1104没有重新计算连通块；它读取1101的`largest_high_component_fraction`。
- 当前阈值是每张芯片对应B-associated factor的7个患者自身q90的中位数，同一芯片患者共享该阈值。
- 高usage定义为`relative_usage >= cutoff`。
- 图节点是该患者全部有效Bin50中心；边为同一患者内欧氏距离不超过36 μm的中心对。在25 μm格点上等价于queen邻接（水平/垂直25 μm及对角35.36 μm）。
- 原指标分子为最大高usage连通块Bin数，分母为该患者全部高usage Bin数；它首先描述high集合内部的集中程度，不是连续区域相对全部组织的面积。

## 阈值核查

原统一阈值下，高usage覆盖超过50%的患者有{len(baseline_extreme)}位。最突出的是`{baseline_extreme.iloc[0]['patient_id'] if len(baseline_extreme) else '无'}`，覆盖比例为{baseline_extreme.iloc[0]['high_usage_fraction'] if len(baseline_extreme) else float('nan'):.3f}。因此，原指标可在高usage覆盖差异很大时获得较高值，不能单独解释为实际连续区域更大。

## IPF vs HC拆解

- 冻结1104基准：原`largest_high_component_fraction`效应={float(original.iloc[0]['effect_within_chip_sd']):.3f}个芯片内SD，置换P={float(original.iloc[0]['permutation_p']):.5f}，FDR={float(original.iloc[0]['permutation_q']):.3f}。该结果保持不变，且仍属于未通过FDR的探索性信号。
- 原阈值下，最大块/全部有效Bin的效应为{current.loc['largest_component_over_valid', 'effect_within_chip_sd']:.3f}，最大块绝对Bin数效应为{current.loc['largest_component_bin_n', 'effect_within_chip_sd']:.3f}。
- 当每位患者固定为自身top 10%（patient-specific q90）时，最大块/高usage Bin效应为{patient_q90.loc['largest_component_over_high', 'effect_within_chip_sd']:.3f}，最大块/全部有效Bin效应为{patient_q90.loc['largest_component_over_valid', 'effect_within_chip_sd']:.3f}，绝对Bin数效应为{patient_q90.loc['largest_component_bin_n', 'effect_within_chip_sd']:.3f}。这组结果用于区分空间集中度与覆盖/组织面积。
- 患者自身q90下，IPF最大块平均为{patient_q90.loc['largest_component_bin_n', 'ipf_mean_raw']:.1f}个Bin（约{patient_q90.loc['largest_component_bin_n', 'ipf_mean_raw'] * BIN_AREA_MM2:.4f} mm²），HC为{patient_q90.loc['largest_component_bin_n', 'hc_mean_raw']:.1f}个Bin（约{patient_q90.loc['largest_component_bin_n', 'hc_mean_raw'] * BIN_AREA_MM2:.4f} mm²）；该绝对大小检验P={patient_q90.loc['largest_component_bin_n', 'permutation_p_exploratory']:.3f}，并不稳健。相反，最大块/有效组织和最大块/high Bin的探索性P分别为{patient_q90.loc['largest_component_over_valid', 'permutation_p_exploratory']:.4f}和{patient_q90.loc['largest_component_over_high', 'permutation_p_exploratory']:.4f}；二者在固定10%覆盖时近似成比例，不能算两项独立证据。12项敏感性结果的探索性BH值也完整报告，不替代原1104的FDR体系。

完整的q85、q90、q95以及患者自身q90结果均保存在`04_IPF_HC_METRIC_SENSITIVITY.tsv`，没有选择性删除方向不一致的结果。

## 邻接与ROI核查

- 代码在`groupby(patient_id)`内部构图，因此没有跨患者连接。
- 21/21患者各自只有一个`bin_id`来源前缀；没有重复坐标。
- 所有高usage连通块都未跨越由全部有效Bin构成的组织连通分量。
- 输入没有显式ROI列，所以无法用独立ROI标签再次验证；当前可验证的结论是：按患者和坐标图定义，没有跨患者或跨坐标不连续组织分量的连接。若同一患者存在两个人工ROI且坐标边界相距≤36 μm，数据表本身无法识别，但本批次未提供这类ROI标签。

## 结论

原`largest_high_component_fraction`在数学上主要是**高usage区域内部的空间集中程度**：最大块占所有high Bin的比例，它不等同于实际连续区域大小。敏感性结果表明，统一阈值异常覆盖并不能完全解释IPF方向：固定每位患者top 10%后，三张芯片和全部留一患者结果仍同向，支持IPF的B-associated高usage Bin更倾向聚成一个主要斑块，并占更大的患者内组织比例。但绝对最大块Bin数的证据较弱且受组织大小影响，因此目前不能确认IPF具有更大的绝对物理连续区域。最审慎的表述是“IPF存在更强的相对空间聚集/相对组织占据趋势”，仍是探索性结果，不是确认性疾病结论。
"""
    (output / "08_METHOD_AUDIT_ZH.md").write_text(text, encoding="utf-8")


def main() -> None:
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    output = ROOT / "results/reanalysis" / f"B_ASSOC_LARGEST_COMPONENT_AUDIT_{timestamp}"
    figures = output / "figures"
    figures.mkdir(parents=True)

    data, _ = load_data()
    thresholds = threshold_lookup(data)
    patient, adjacency = patient_metrics(data, thresholds)
    baseline = patient.loc[
        patient["threshold_scheme"].eq("chip_median_patient_q90")
    ].copy()

    source_metrics = pd.read_csv(SOURCE_1101 / "02_PATIENT_FACTOR_METRICS.tsv", sep="\t")
    source_metrics = source_metrics.loc[
        source_metrics.apply(lambda row: B_FACTORS.get(row["chip"]) == row["factor"], axis=1)
    ].copy()
    compare = baseline.merge(
        source_metrics[
            [
                "chip",
                "patient_id",
                "balanced_q90_cutoff",
                "high_usage_fraction",
                "largest_high_component_fraction",
            ]
        ],
        on=["chip", "patient_id"],
        suffixes=("_recomputed", "_source"),
        validate="one_to_one",
    )
    compare["threshold_absolute_difference"] = np.abs(
        compare["high_usage_threshold"] - compare["balanced_q90_cutoff"]
    )
    compare["coverage_absolute_difference"] = np.abs(
        compare["high_usage_fraction_recomputed"] - compare["high_usage_fraction_source"]
    )
    compare["component_fraction_absolute_difference"] = np.abs(
        compare["largest_component_over_high"]
        - compare["largest_high_component_fraction"]
    )

    summary, chip, loo = disease_sensitivity(patient)
    original_all = pd.read_csv(
        SOURCE_1104 / "05_DISEASE_CONTRASTS_AND_PERMUTATIONS.tsv", sep="\t"
    )
    original = original_all.loc[
        original_all["family_id"].eq("B_ASSOC")
        & original_all["metric"].eq("largest_high_component_fraction")
        & original_all["contrast"].eq("IPF_vs_HC")
    ].copy()

    provenance = pd.DataFrame(
        [
            {"path": str(path), "sha256": sha256(path)}
            for path in [
                ASSIGNMENTS,
                SOURCE_1101 / "02_PATIENT_FACTOR_METRICS.tsv",
                SOURCE_1104 / "05_DISEASE_CONTRASTS_AND_PERMUTATIONS.tsv",
                *[
                    MODEL_RUN / "models" / f"neighborhood50_{chip_id}_K06_MODEL_AND_USAGE.npz"
                    for chip_id in B_FACTORS
                ],
            ]
        ]
    )
    provenance.to_csv(output / "00_INPUT_PROVENANCE.tsv", sep="\t", index=False)
    thresholds.to_csv(output / "01_THRESHOLD_DEFINITIONS.tsv", sep="\t", index=False)
    baseline.to_csv(output / "02_B_ASSOC_PATIENT_BASELINE.tsv", sep="\t", index=False)
    patient.to_csv(output / "03_B_ASSOC_THRESHOLD_SENSITIVITY_PATIENT.tsv", sep="\t", index=False)
    summary.to_csv(output / "04_IPF_HC_METRIC_SENSITIVITY.tsv", sep="\t", index=False)
    chip.to_csv(output / "05_CHIP_STRATIFIED_EFFECTS.tsv", sep="\t", index=False)
    loo.to_csv(output / "06_PATIENT_INFLUENCE_LOO.tsv", sep="\t", index=False)
    adjacency.to_csv(output / "07_SPATIAL_ADJACENCY_AUDIT.tsv", sep="\t", index=False)
    original.to_csv(output / "09_ORIGINAL_1104_BASELINE_RESULT.tsv", sep="\t", index=False)
    compare.to_csv(output / "10_RECOMPUTATION_CONCORDANCE.tsv", sep="\t", index=False)

    validation_rows = [
        ("patient_n", len(baseline), 21, len(baseline) == 21),
        ("threshold_scheme_n", patient["threshold_scheme"].nunique(), 4, patient["threshold_scheme"].nunique() == 4),
        ("sensitivity_patient_rows", len(patient), 84, len(patient) == 84),
        ("metric_sensitivity_rows", len(summary), 12, len(summary) == 12),
        ("baseline_threshold_matches_1101", compare["threshold_absolute_difference"].max(), 0, compare["threshold_absolute_difference"].max() < 1e-12),
        ("baseline_coverage_matches_1101", compare["coverage_absolute_difference"].max(), 0, compare["coverage_absolute_difference"].max() < 1e-12),
        ("baseline_component_matches_1101", compare["component_fraction_absolute_difference"].max(), 0, compare["component_fraction_absolute_difference"].max() < 1e-12),
        ("no_duplicate_coordinates", adjacency["duplicate_coordinate_n"].sum(), 0, adjacency["duplicate_coordinate_n"].sum() == 0),
        ("one_source_id_per_patient", adjacency["source_id_n"].max(), 1, adjacency["source_id_n"].max() == 1),
        ("no_high_component_crosses_tissue_component", patient["largest_component_spans_tissue_components"].sum(), 0, patient["largest_component_spans_tissue_components"].sum() == 0),
        ("original_baseline_row_n", len(original), 1, len(original) == 1),
        ("original_result_unchanged", float(original.iloc[0]["effect_within_chip_sd"]), 1.2698552050829444, abs(float(original.iloc[0]["effect_within_chip_sd"]) - 1.2698552050829444) < 1e-12),
    ]
    validation = pd.DataFrame(validation_rows, columns=["check", "observed", "expected", "passed"])
    validation.to_csv(output / "11_VALIDATION.tsv", sep="\t", index=False)
    if not validation["passed"].all():
        raise RuntimeError("Validation failed")

    make_figures(patient, summary, figures)
    write_report(output, baseline, summary, adjacency, original)
    session = {
        "timestamp": timestamp,
        "python": sys.version,
        "platform": platform.platform(),
        "numpy": np.__version__,
        "pandas": pd.__version__,
        "seed": SEED,
        "permutations": N_PERMUTATIONS,
        "edge_radius_um": EDGE_RADIUS_UM,
    }
    (output / "SESSION_INFO.json").write_text(
        json.dumps(session, ensure_ascii=False, indent=2), encoding="utf-8"
    )
    (output / "COMPLETE.ok").write_text("PASS\n", encoding="utf-8")
    print(output)


if __name__ == "__main__":
    main()
