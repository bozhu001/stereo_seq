#!/usr/bin/env python3
"""Uniform component sensitivity audit for four frozen primary NMF families.

This is a non-selective extension of the approved B-associated audit. It
reuses frozen K=6 relative usages and applies identical thresholds, component
definitions, disease contrasts, within-chip permutations, and influence
checks to B-, Ciliated-, Fibroblast-, and Lymphatic-associated factor families.
No NMF/RCTD model is refit and no factor identity is changed.
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
SOURCE_1101 = ROOT / "results/reanalysis/REFERENCE18_K6_FACTOR_DOWNSTREAM_20261013_003020"
SOURCE_1104 = ROOT / "results/reanalysis/FACTOR_FAMILY_PATIENT_DISEASE_20261013_005851"

FAMILY_FACTORS = {
    "B-associated": {"K8": "F02", "J2": "F03", "L3": "F01"},
    "Ciliated-associated": {"K8": "F03", "J2": "F04", "L3": "F02"},
    "Fibroblast-associated": {"K8": "F04", "J2": "F06", "L3": "F05"},
    "Lymphatic-associated": {"K8": "F01", "J2": "F02", "L3": "F04"},
}
FAMILY_IDS = {
    "B-associated": "B_ASSOC",
    "Ciliated-associated": "CILIATED",
    "Fibroblast-associated": "FIBROBLAST",
    "Lymphatic-associated": "LYMPHATIC",
}
DISEASE_CODES = {"HC": 0, "IPF": 1, "SSc-ILD": 2}
CONTRASTS = [
    ("IPF_vs_HC", "IPF", "HC"),
    ("SSc-ILD_vs_HC", "SSc-ILD", "HC"),
    ("SSc-ILD_vs_IPF", "SSc-ILD", "IPF"),
]
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
    return cKDTree(coordinates).query_pairs(
        EDGE_RADIUS_UM + 1e-8,
        output_type="ndarray",
    ).astype(np.int64)


def graph_labels(node_n: int, edges: np.ndarray) -> tuple[int, np.ndarray]:
    if len(edges) == 0:
        return node_n, np.arange(node_n, dtype=int)
    rows = np.r_[edges[:, 0], edges[:, 1]]
    columns = np.r_[edges[:, 1], edges[:, 0]]
    adjacency = csr_matrix(
        (np.ones(len(rows), dtype=np.uint8), (rows, columns)),
        shape=(node_n, node_n),
    )
    return connected_components(adjacency, directed=False)


def component_details(
    high: np.ndarray,
    edges: np.ndarray,
    tissue_labels: np.ndarray,
) -> tuple[int, int, bool]:
    high_indices = np.flatnonzero(high)
    if len(high_indices) == 0:
        return 0, 0, False
    remap = np.full(len(high), -1, dtype=int)
    remap[high_indices] = np.arange(len(high_indices))
    keep = high[edges[:, 0]] & high[edges[:, 1]]
    high_edges = edges[keep]
    if len(high_edges) == 0:
        component_n = len(high_indices)
        labels = np.arange(component_n, dtype=int)
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
        component_n, labels = connected_components(adjacency, directed=False)
    sizes = np.bincount(labels, minlength=component_n)
    spanning = any(
        len(np.unique(tissue_labels[high_indices[labels == label]])) > 1
        for label in range(component_n)
    )
    return int(component_n), int(sizes.max()), bool(spanning)


def load_long_usage() -> pd.DataFrame:
    assignments = pd.read_csv(ASSIGNMENTS, sep="\t")
    pieces = []
    for chip in ("K8", "J2", "L3"):
        model_path = MODEL_RUN / "models" / f"neighborhood50_{chip}_K06_MODEL_AND_USAGE.npz"
        with np.load(model_path, allow_pickle=True) as model:
            keys = model["entity_key"].astype(str)
            relative_usage = model["relative_usage"].astype(float)
        lookup = assignments.loc[assignments["chip_label"].eq(chip)].set_index("entity_key")
        metadata = lookup.loc[keys].reset_index()
        if metadata["entity_key"].astype(str).tolist() != keys.tolist():
            raise RuntimeError(f"{chip}: entity order mismatch")
        for family, mapping in FAMILY_FACTORS.items():
            factor = mapping[chip]
            frame = metadata.copy()
            frame["family_label"] = family
            frame["family_id"] = FAMILY_IDS[family]
            frame["factor"] = factor
            frame["usage"] = relative_usage[:, int(factor[1:]) - 1]
            pieces.append(frame)
    return pd.concat(pieces, ignore_index=True)


def build_thresholds(data: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for (family, chip), frame in data.groupby(["family_label", "chip_label"], sort=True):
        for scheme, rule, quantile in THRESHOLD_SCHEMES:
            patient_quantiles = frame.groupby("patient_id")["usage"].quantile(quantile)
            common = float(patient_quantiles.median()) if rule == "chip_balanced" else math.nan
            for patient, patient_cutoff in patient_quantiles.items():
                rows.append(
                    {
                        "family_id": FAMILY_IDS[family],
                        "family_label": family,
                        "chip": chip,
                        "factor": FAMILY_FACTORS[family][chip],
                        "patient_id": patient,
                        "threshold_scheme": scheme,
                        "threshold_rule": rule,
                        "target_quantile": quantile,
                        "high_usage_threshold": (
                            common if rule == "chip_balanced" else float(patient_cutoff)
                        ),
                    }
                )
    return pd.DataFrame(rows)


def build_patient_metrics(
    data: pd.DataFrame,
    thresholds: pd.DataFrame,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    threshold_index = thresholds.set_index(
        ["family_label", "chip", "patient_id", "threshold_scheme"]
    )
    metric_rows = []
    adjacency_rows = []
    graph_cache: dict[tuple[str, str], tuple[np.ndarray, np.ndarray, pd.DataFrame]] = {}
    for (chip, patient), base in data.loc[
        data["family_label"].eq("B-associated")
    ].groupby(["chip_label", "patient_id"], sort=True):
        base = base.reset_index(drop=True)
        coordinates = base[["x_um", "y_um"]].to_numpy(dtype=float)
        edges = edge_array(coordinates)
        tissue_n, tissue_labels = graph_labels(len(base), edges)
        source_ids = base["bin_id"].astype(str).str.split("|", n=1).str[0]
        distances = (
            np.linalg.norm(coordinates[edges[:, 0]] - coordinates[edges[:, 1]], axis=1)
            if len(edges)
            else np.array([], dtype=float)
        )
        graph_cache[(chip, patient)] = (edges, tissue_labels, base)
        adjacency_rows.append(
            {
                "chip": chip,
                "patient_id": patient,
                "disease": base["disease"].iloc[0],
                "valid_bin_n": len(base),
                "source_id_n": source_ids.nunique(),
                "duplicate_coordinate_n": int(base.duplicated(["x_um", "y_um"]).sum()),
                "spatial_edge_n": len(edges),
                "minimum_edge_distance_um": float(distances.min()) if len(distances) else math.nan,
                "maximum_edge_distance_um": float(distances.max()) if len(distances) else math.nan,
                "tissue_component_n": tissue_n,
                "largest_tissue_component_fraction": float(
                    np.bincount(tissue_labels).max() / len(base)
                ),
                "patient_isolation_by_construction": True,
                "explicit_roi_column_available": False,
            }
        )

    for (family, chip, patient), frame in data.groupby(
        ["family_label", "chip_label", "patient_id"], sort=True
    ):
        frame = frame.reset_index(drop=True)
        edges, tissue_labels, base = graph_cache[(chip, patient)]
        if frame["entity_key"].tolist() != base["entity_key"].tolist():
            raise RuntimeError("Cached spatial graph entity order mismatch")
        values = frame["usage"].to_numpy(dtype=float)
        for scheme, _, _ in THRESHOLD_SCHEMES:
            key = (family, chip, patient, scheme)
            threshold_row = threshold_index.loc[key]
            cutoff = float(threshold_row["high_usage_threshold"])
            high = values >= cutoff
            component_n, largest_n, spanning = component_details(
                high, edges, tissue_labels
            )
            high_n = int(high.sum())
            metric_rows.append(
                {
                    "family_id": FAMILY_IDS[family],
                    "family_label": family,
                    "chip": chip,
                    "factor": FAMILY_FACTORS[family][chip],
                    "patient_id": patient,
                    "disease": frame["disease"].iloc[0],
                    "threshold_scheme": scheme,
                    "threshold_rule": threshold_row["threshold_rule"],
                    "target_quantile": float(threshold_row["target_quantile"]),
                    "high_usage_threshold": cutoff,
                    "valid_bin_n": len(frame),
                    "high_usage_bin_n": high_n,
                    "high_usage_fraction": high_n / len(frame),
                    "high_component_n": component_n,
                    "largest_component_bin_n": largest_n,
                    "largest_component_area_mm2": largest_n * BIN_AREA_MM2,
                    "largest_component_over_high": largest_n / high_n if high_n else math.nan,
                    "largest_component_over_valid": largest_n / len(frame),
                    "largest_component_spans_tissue_components": spanning,
                }
            )
    return pd.DataFrame(metric_rows), pd.DataFrame(adjacency_rows)


def make_permutations(meta: pd.DataFrame) -> np.ndarray:
    rng = np.random.default_rng(SEED)
    observed = meta["disease"].map(DISEASE_CODES).to_numpy(dtype=np.uint8)
    permutations = np.empty((N_PERMUTATIONS, len(meta)), dtype=np.uint8)
    for chip in sorted(meta["chip"].unique()):
        positions = np.flatnonzero(meta["chip"].to_numpy() == chip)
        labels = observed[positions]
        for index in range(N_PERMUTATIONS):
            permutations[index, positions] = rng.permutation(labels)
    return permutations


def disease_tests(
    patient: pd.DataFrame,
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    meta = (
        patient[["chip", "patient_id", "disease"]]
        .drop_duplicates()
        .sort_values(["chip", "patient_id"])
        .reset_index(drop=True)
    )
    permutations = make_permutations(meta)
    position = dict(zip(meta["patient_id"], meta.index))
    summary_rows = []
    chip_rows = []
    loo_rows = []
    for family in FAMILY_FACTORS:
        for scheme, _, _ in THRESHOLD_SCHEMES:
            subset = patient.loc[
                patient["family_label"].eq(family)
                & patient["threshold_scheme"].eq(scheme)
            ].copy()
            subset = subset.sort_values(["chip", "patient_id"]).reset_index(drop=True)
            positions = np.array([position[value] for value in subset["patient_id"]])
            permuted = permutations[:, positions]
            observed_codes = subset["disease"].map(DISEASE_CODES).to_numpy(dtype=np.uint8)
            for metric, method in METRICS.items():
                subset["transformed"] = transform(subset[metric], method)
                subset["within_chip_z"] = subset.groupby("chip")["transformed"].transform(zscore)
                values = subset["within_chip_z"].to_numpy(dtype=float)
                for contrast, numerator, denominator in CONTRASTS:
                    numerator_code = DISEASE_CODES[numerator]
                    denominator_code = DISEASE_CODES[denominator]
                    observed_effect = float(
                        values[observed_codes == numerator_code].mean()
                        - values[observed_codes == denominator_code].mean()
                    )
                    numerator_mask = permuted == numerator_code
                    denominator_mask = permuted == denominator_code
                    permutation_effects = (
                        (numerator_mask @ values) / numerator_mask.sum(axis=1)
                        - (denominator_mask @ values) / denominator_mask.sum(axis=1)
                    )
                    p_value = (
                        1 + np.sum(np.abs(permutation_effects) >= abs(observed_effect))
                    ) / (N_PERMUTATIONS + 1)
                    loo_effects = []
                    contrast_patients = subset.loc[
                        subset["disease"].isin([numerator, denominator]), "patient_id"
                    ]
                    for patient_id in contrast_patients:
                        keep = subset["patient_id"].ne(patient_id).to_numpy()
                        kept_values = values[keep]
                        kept_codes = observed_codes[keep]
                        effect = float(
                            kept_values[kept_codes == numerator_code].mean()
                            - kept_values[kept_codes == denominator_code].mean()
                        )
                        loo_effects.append(effect)
                        loo_rows.append(
                            {
                                "family_id": FAMILY_IDS[family],
                                "family_label": family,
                                "threshold_scheme": scheme,
                                "metric": metric,
                                "contrast": contrast,
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
                    chip_effects = []
                    for chip, chip_frame in subset.groupby("chip", sort=True):
                        effect = float(
                            chip_frame.loc[
                                chip_frame["disease"].eq(numerator), "within_chip_z"
                            ].mean()
                            - chip_frame.loc[
                                chip_frame["disease"].eq(denominator), "within_chip_z"
                            ].mean()
                        )
                        chip_effects.append(effect)
                        chip_rows.append(
                            {
                                "family_id": FAMILY_IDS[family],
                                "family_label": family,
                                "threshold_scheme": scheme,
                                "metric": metric,
                                "contrast": contrast,
                                "chip": chip,
                                "effect_within_chip_sd": effect,
                            }
                        )
                    summary_rows.append(
                        {
                            "family_id": FAMILY_IDS[family],
                            "family_label": family,
                            "threshold_scheme": scheme,
                            "metric": metric,
                            "transform": method,
                            "contrast": contrast,
                            "effect_within_chip_sd": observed_effect,
                            "permutation_p_exploratory": float(p_value),
                            "permutation_n": N_PERMUTATIONS,
                            "chip_effect_min": float(np.min(chip_effects)),
                            "chip_effect_max": float(np.max(chip_effects)),
                            "chip_same_direction_fraction": float(
                                np.mean(np.sign(chip_effects) == np.sign(observed_effect))
                            ),
                            "all_three_chips_same_direction": bool(
                                len(set(np.sign(chip_effects))) == 1
                            ),
                            "loo_effect_min": float(np.min(loo_effects)),
                            "loo_effect_max": float(np.max(loo_effects)),
                            "loo_same_direction_fraction": float(
                                np.mean(np.sign(loo_effects) == np.sign(observed_effect))
                            ),
                        }
                    )
    summary = pd.DataFrame(summary_rows)
    summary["exploratory_bh_q_global_144"] = bh_adjust(
        summary["permutation_p_exploratory"]
    )
    summary["exploratory_bh_q_within_contrast_48"] = summary.groupby(
        "contrast", group_keys=False
    )["permutation_p_exploratory"].apply(bh_adjust)
    return summary, pd.DataFrame(chip_rows), pd.DataFrame(loo_rows)


def summarize_robustness(
    patient: pd.DataFrame,
    summary: pd.DataFrame,
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    group_rows = []
    for keys, frame in patient.groupby(
        ["family_id", "family_label", "threshold_scheme", "chip", "disease"],
        sort=True,
    ):
        family_id, family, scheme, chip, disease = keys
        for metric in METRICS:
            values = frame[metric].astype(float)
            group_rows.append(
                {
                    "family_id": family_id,
                    "family_label": family,
                    "threshold_scheme": scheme,
                    "chip": chip,
                    "disease": disease,
                    "metric": metric,
                    "patient_n": len(values),
                    "mean": float(values.mean()),
                    "median": float(values.median()),
                    "minimum": float(values.min()),
                    "maximum": float(values.max()),
                }
            )

    robustness_rows = []
    for keys, frame in summary.groupby(
        ["family_id", "family_label", "metric", "contrast"], sort=True
    ):
        family_id, family, metric, contrast = keys
        effects = frame["effect_within_chip_sd"].to_numpy(dtype=float)
        signs = np.sign(effects)
        robustness_rows.append(
            {
                "family_id": family_id,
                "family_label": family,
                "metric": metric,
                "contrast": contrast,
                "threshold_scheme_n": len(frame),
                "effect_min": float(effects.min()),
                "effect_max": float(effects.max()),
                "positive_threshold_scheme_n": int(np.sum(effects > 0)),
                "negative_threshold_scheme_n": int(np.sum(effects < 0)),
                "same_direction_all_thresholds": bool(len(set(signs)) == 1),
                "all_three_chips_same_direction_scheme_n": int(
                    frame["all_three_chips_same_direction"].sum()
                ),
                "loo_full_direction_scheme_n": int(
                    np.isclose(frame["loo_same_direction_fraction"], 1.0).sum()
                ),
                "nominal_p_lt_0_05_scheme_n": int(
                    (frame["permutation_p_exploratory"] < 0.05).sum()
                ),
                "minimum_permutation_p": float(
                    frame["permutation_p_exploratory"].min()
                ),
                "minimum_global_q": float(
                    frame["exploratory_bh_q_global_144"].min()
                ),
            }
        )
    robustness = pd.DataFrame(robustness_rows)

    primary = summary.loc[
        summary["threshold_scheme"].eq("patient_specific_q90")
        & summary["metric"].eq("largest_component_over_high")
    ].copy()
    recommendation_rows = []
    for row in primary.itertuples():
        if row.all_three_chips_same_direction and row.loo_same_direction_fraction == 1.0:
            if row.permutation_p_exploratory < 0.05:
                classification = "Nominal cross-chip exploratory candidate"
            elif abs(row.effect_within_chip_sd) >= 0.8:
                classification = "Cross-chip directional candidate"
            else:
                classification = "Secondary descriptive candidate"
        else:
            classification = "Unresolved across chips"
        recommendation_rows.append(
            {
                "family_id": row.family_id,
                "family_label": row.family_label,
                "contrast": row.contrast,
                "recommended_variable": (
                    "patient_specific_q90_largest_component_over_high"
                ),
                "classification": classification,
                "effect_within_chip_sd": row.effect_within_chip_sd,
                "permutation_p_exploratory": row.permutation_p_exploratory,
                "global_q_144": row.exploratory_bh_q_global_144,
                "all_three_chips_same_direction": row.all_three_chips_same_direction,
                "loo_same_direction_fraction": row.loo_same_direction_fraction,
                "confirmed_disease_association": False,
                "interpretation_limit": (
                    "Relative top-decile concentration only; not absolute physical area."
                ),
            }
        )
    return pd.DataFrame(group_rows), robustness, pd.DataFrame(recommendation_rows)


def make_figures(
    patient: pd.DataFrame,
    summary: pd.DataFrame,
    figures: Path,
) -> None:
    sns.set_theme(style="whitegrid", context="notebook")
    ipf_hc = summary.loc[summary["contrast"].eq("IPF_vs_HC")].copy()
    metric_labels = {
        "largest_component_bin_n": "Absolute component size",
        "largest_component_over_valid": "Component / valid bins",
        "largest_component_over_high": "Component / high-usage bins",
    }
    figure, axes = plt.subplots(2, 2, figsize=(17, 11), sharex=True)
    for axis, family in zip(axes.ravel(), FAMILY_FACTORS):
        frame = ipf_hc.loc[ipf_hc["family_label"].eq(family)].copy()
        frame["metric_label"] = frame["metric"].map(metric_labels)
        sns.pointplot(
            data=frame,
            x="effect_within_chip_sd",
            y="threshold_scheme",
            hue="metric_label",
            join=False,
            dodge=0.35,
            ax=axis,
        )
        axis.axvline(0, color="black", linewidth=1)
        axis.set_title(family)
        axis.set_xlabel("IPF - HC effect (within-chip SD)")
        axis.set_ylabel("Threshold scheme")
        if axis is not axes[0, 1]:
            axis.get_legend().remove()
    axes[0, 1].legend(title="Metric", bbox_to_anchor=(1.02, 1), loc="upper left")
    figure.suptitle("Primary factor-family component sensitivity", fontsize=16)
    figure.tight_layout()
    figure.savefig(figures / "01_PRIMARY_FAMILY_THRESHOLD_SENSITIVITY.png", dpi=220)
    plt.close(figure)

    patient_q90 = patient.loc[patient["threshold_scheme"].eq("patient_specific_q90")].copy()
    figure, axes = plt.subplots(2, 2, figsize=(15, 10))
    for axis, family in zip(axes.ravel(), FAMILY_FACTORS):
        frame = patient_q90.loc[patient_q90["family_label"].eq(family)]
        sns.boxplot(
            data=frame,
            x="disease",
            y="largest_component_over_high",
            order=["HC", "IPF", "SSc-ILD"],
            color="#d9e6f2",
            ax=axis,
        )
        sns.stripplot(
            data=frame,
            x="disease",
            y="largest_component_over_high",
            order=["HC", "IPF", "SSc-ILD"],
            hue="chip",
            dodge=False,
            size=7,
            ax=axis,
        )
        axis.set_title(family)
        axis.set_xlabel("")
        axis.set_ylabel("Largest component / top-decile bins")
        if axis is not axes[0, 1]:
            axis.get_legend().remove()
    axes[0, 1].legend(title="Chip", bbox_to_anchor=(1.02, 1), loc="upper left")
    figure.suptitle("Patient-specific q90 component concentration", fontsize=16)
    figure.tight_layout()
    figure.savefig(figures / "02_PATIENT_Q90_COMPONENT_CONCENTRATION.png", dpi=220)
    plt.close(figure)


def write_report(
    output: Path,
    summary: pd.DataFrame,
    original: pd.DataFrame,
    recommendation: pd.DataFrame,
    extreme_coverage: pd.DataFrame,
) -> None:
    ipf_hc = summary.loc[
        summary["contrast"].eq("IPF_vs_HC")
        & summary["threshold_scheme"].eq("patient_specific_q90")
        & summary["metric"].eq("largest_component_over_high")
    ].sort_values("family_label")
    lines = [
        "# 四个primary factor-family连通块统一敏感性核查",
        "",
        "## 范围",
        "",
        "本轮将已批准的B-associated核查规则不加选择地应用于B、Ciliated、Fibroblast和Lymphatic四个三芯片primary family。未重跑NMF/RCTD，未修改factor-family映射，未开展factor–factor分析。",
        "",
        "## 固定方法",
        "",
        "- 阈值：芯片内患者q85/q90/q95中位数，以及患者自身q90。",
        "- 指标：最大块绝对Bin数、最大块/全部有效Bin、最大块/全部high Bin。",
        "- 邻接：患者内36 μm queen-like邻接；不跨患者或坐标图组织分量。",
        "- 推断：患者为单位；指标在family×chip内转换并标准化；疾病标签仅在芯片内置换100,000次。",
        "- 完整报告IPF-HC、SSc-ILD-HC和SSc-ILD-IPF，探索性BH同时按全部144项和每个contrast内48项计算。",
        "",
        "## patient-specific q90下IPF vs HC集中度",
        "",
        "| Family | Effect (within-chip SD) | P | Global q | 3-chip same direction | LOO direction fraction |",
        "|---|---:|---:|---:|---|---:|",
    ]
    for row in ipf_hc.itertuples():
        lines.append(
            f"| {row.family_label} | {row.effect_within_chip_sd:.3f} | "
            f"{row.permutation_p_exploratory:.4f} | "
            f"{row.exploratory_bh_q_global_144:.3f} | "
            f"{row.all_three_chips_same_direction} | {row.loo_same_direction_fraction:.3f} |"
        )
    lines.extend(
        [
            "",
            "## 统一多重校正与保留原则",
            "",
            f"全部{len(summary)}项探索性敏感性检验中，通过global BH q<0.05的结果为{int((summary['exploratory_bh_q_global_144'] < 0.05).sum())}项。因此没有任何疾病关联被升级为确认性结论。patient-specific q90的相对集中度可作为后续患者级空间变量，但其疾病方向按`13_FINAL_VARIABLE_RECOMMENDATION.tsv`分为nominal cross-chip exploratory、cross-chip directional、secondary descriptive或unresolved。",
            "",
            "## 统一阈值覆盖异常",
            "",
            f"原chip-median patient q90规则下，高usage覆盖超过50%的family-patient组合共有{len(extreme_coverage)}个，完整列于`14_EXTREME_BASELINE_COVERAGE.tsv`。这类异常不限于B-associated，因此芯片统一阈值不宜作为跨患者面积比较的唯一规则；患者自身q90应保留为主要拓扑敏感性。",
            "",
            "## 解释原则",
            "",
            "`largest_component_over_high`衡量top集合内部的集中程度；`largest_component_over_valid`衡量患者内相对组织占据；`largest_component_bin_n`才对应绝对物理块大小。患者自身q90固定约10% high，因此前两个比例近似成比例，不能视为两项独立证据。",
            "",
            "原1104结果原样保存在`09_ORIGINAL_1104_PRIMARY_COMPONENT_RESULTS.tsv`，本轮敏感性结果不替代其P值和FDR。只有阈值敏感、三芯片方向、LOO和三种分母共同支持时，才建议冻结为后续空间变量；其余保留为描述性或unresolved。",
            "",
            f"原1104基准行数：{len(original)}。完整144项结果见`04_ALL_DISEASE_COMPONENT_SENSITIVITY.tsv`。",
            "",
            f"最终推荐行数：{len(recommendation)}（4 families × 3 contrasts）；所有`confirmed_disease_association`均为False。",
        ]
    )
    (output / "08_METHOD_AND_RESULTS_ZH.md").write_text(
        "\n".join(lines) + "\n", encoding="utf-8"
    )


def main() -> None:
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    output = ROOT / "results/reanalysis" / f"PRIMARY_FACTOR_COMPONENT_SENSITIVITY_{timestamp}"
    figures = output / "figures"
    figures.mkdir(parents=True)

    data = load_long_usage()
    thresholds = build_thresholds(data)
    patient, adjacency = build_patient_metrics(data, thresholds)
    summary, chip, loo = disease_tests(patient)
    group_summary, robustness, recommendation = summarize_robustness(patient, summary)

    source = pd.read_csv(SOURCE_1101 / "02_PATIENT_FACTOR_METRICS.tsv", sep="\t")
    map_rows = [
        (FAMILY_IDS[family], family, chip_id, factor)
        for family, mapping in FAMILY_FACTORS.items()
        for chip_id, factor in mapping.items()
    ]
    mapping = pd.DataFrame(
        map_rows, columns=["family_id", "family_label", "chip", "factor"]
    )
    source = source.merge(mapping, on=["chip", "factor"], how="inner")
    baseline = patient.loc[
        patient["threshold_scheme"].eq("chip_median_patient_q90")
    ].copy()
    concordance = baseline.merge(
        source[
            [
                "family_label",
                "chip",
                "patient_id",
                "balanced_q90_cutoff",
                "high_usage_fraction",
                "largest_high_component_fraction",
            ]
        ],
        on=["family_label", "chip", "patient_id"],
        suffixes=("_recomputed", "_source"),
        validate="one_to_one",
    )
    concordance["threshold_absolute_difference"] = np.abs(
        concordance["high_usage_threshold"] - concordance["balanced_q90_cutoff"]
    )
    concordance["coverage_absolute_difference"] = np.abs(
        concordance["high_usage_fraction_recomputed"]
        - concordance["high_usage_fraction_source"]
    )
    concordance["component_fraction_absolute_difference"] = np.abs(
        concordance["largest_component_over_high"]
        - concordance["largest_high_component_fraction"]
    )

    original_all = pd.read_csv(
        SOURCE_1104 / "05_DISEASE_CONTRASTS_AND_PERMUTATIONS.tsv", sep="\t"
    )
    original = original_all.loc[
        original_all["family_id"].isin(FAMILY_IDS.values())
        & original_all["metric"].eq("largest_high_component_fraction")
    ].copy()
    extreme_coverage = baseline.loc[baseline["high_usage_fraction"] > 0.50].copy()

    provenance_paths = [
        ASSIGNMENTS,
        SOURCE_1101 / "02_PATIENT_FACTOR_METRICS.tsv",
        SOURCE_1104 / "05_DISEASE_CONTRASTS_AND_PERMUTATIONS.tsv",
        *[
            MODEL_RUN / "models" / f"neighborhood50_{chip_id}_K06_MODEL_AND_USAGE.npz"
            for chip_id in ("K8", "J2", "L3")
        ],
    ]
    pd.DataFrame(
        [{"path": str(path), "sha256": sha256(path)} for path in provenance_paths]
    ).to_csv(output / "00_INPUT_PROVENANCE.tsv", sep="\t", index=False)
    mapping.to_csv(output / "01_FROZEN_PRIMARY_FAMILY_MAP.tsv", sep="\t", index=False)
    thresholds.to_csv(output / "02_THRESHOLD_DEFINITIONS.tsv", sep="\t", index=False)
    patient.to_csv(output / "03_PATIENT_COMPONENT_SENSITIVITY.tsv", sep="\t", index=False)
    summary.to_csv(output / "04_ALL_DISEASE_COMPONENT_SENSITIVITY.tsv", sep="\t", index=False)
    group_summary.to_csv(
        output / "04A_CHIP_DISEASE_RAW_SUMMARIES.tsv", sep="\t", index=False
    )
    chip.to_csv(output / "05_CHIP_STRATIFIED_EFFECTS.tsv", sep="\t", index=False)
    loo.to_csv(output / "06_PATIENT_INFLUENCE_LOO.tsv", sep="\t", index=False)
    adjacency.to_csv(output / "07_SPATIAL_ADJACENCY_AUDIT.tsv", sep="\t", index=False)
    original.to_csv(
        output / "09_ORIGINAL_1104_PRIMARY_COMPONENT_RESULTS.tsv", sep="\t", index=False
    )
    concordance.to_csv(output / "10_BASELINE_RECOMPUTATION_CONCORDANCE.tsv", sep="\t", index=False)
    robustness.to_csv(
        output / "12_CROSS_THRESHOLD_ROBUSTNESS.tsv", sep="\t", index=False
    )
    recommendation.to_csv(
        output / "13_FINAL_VARIABLE_RECOMMENDATION.tsv", sep="\t", index=False
    )
    extreme_coverage.to_csv(
        output / "14_EXTREME_BASELINE_COVERAGE.tsv", sep="\t", index=False
    )

    validations = [
        ("family_n", patient["family_label"].nunique(), 4),
        ("patient_n", patient["patient_id"].nunique(), 21),
        ("patient_metric_rows", len(patient), 336),
        ("disease_test_rows", len(summary), 144),
        ("chip_effect_rows", len(chip), 432),
        ("baseline_concordance_rows", len(concordance), 84),
        ("threshold_max_abs_diff", concordance["threshold_absolute_difference"].max(), 0),
        ("coverage_max_abs_diff", concordance["coverage_absolute_difference"].max(), 0),
        ("component_max_abs_diff", concordance["component_fraction_absolute_difference"].max(), 0),
        ("duplicate_coordinate_n", adjacency["duplicate_coordinate_n"].sum(), 0),
        ("maximum_source_id_n", adjacency["source_id_n"].max(), 1),
        ("cross_tissue_component_count", patient["largest_component_spans_tissue_components"].sum(), 0),
        ("original_result_rows", len(original), 12),
        ("raw_group_summary_rows", len(group_summary), 432),
        ("cross_threshold_robustness_rows", len(robustness), 36),
        ("final_recommendation_rows", len(recommendation), 12),
        (
            "confirmed_disease_association_count",
            recommendation["confirmed_disease_association"].sum(),
            0,
        ),
        ("extreme_baseline_coverage_rows", len(extreme_coverage), 6),
    ]
    validation_rows = []
    for check, observed, expected in validations:
        tolerance = 1e-12 if "diff" in check else 0
        passed = abs(float(observed) - float(expected)) <= tolerance
        validation_rows.append(
            {"check": check, "observed": observed, "expected": expected, "passed": passed}
        )
    validation = pd.DataFrame(validation_rows)
    validation.to_csv(output / "11_VALIDATION.tsv", sep="\t", index=False)
    if not validation["passed"].all():
        raise RuntimeError("Validation failed")

    make_figures(patient, summary, figures)
    write_report(output, summary, original, recommendation, extreme_coverage)
    session = {
        "timestamp": timestamp,
        "python": sys.version,
        "platform": platform.platform(),
        "numpy": np.__version__,
        "pandas": pd.__version__,
        "seed": SEED,
        "permutations": N_PERMUTATIONS,
        "edge_radius_um": EDGE_RADIUS_UM,
        "global_exploratory_test_n": len(summary),
    }
    (output / "SESSION_INFO.json").write_text(
        json.dumps(session, ensure_ascii=False, indent=2), encoding="utf-8"
    )
    (output / "COMPLETE.ok").write_text("PASS\n", encoding="utf-8")
    print(output)


if __name__ == "__main__":
    main()
