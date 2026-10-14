"""Close out stage-1 cross-chip edgeR evidence without refitting edgeR.

This script reuses frozen family and per-chip edgeR tables.  It derives an
approximate standard error from the one-degree-of-freedom edgeR QL-F statistic
and performs fixed-effect and DerSimonian-Laird random-effect meta-analysis
across mapped chip-specific factors.  The existing joint family edgeR model
remains the primary model; meta-analysis is a replication sensitivity layer.
"""

from __future__ import annotations

import argparse
import json
import math
import platform
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import chi2, norm


ENDPOINT = "paired_high_minus_control"
CHIPS = ("K8", "J2", "L3")
FDR_CUTOFF = 0.05
MIN_ABS_LOG2FC = 0.5


def bh_adjust(values: pd.Series) -> pd.Series:
    """Benjamini-Hochberg adjustment preserving missing values and index."""
    result = pd.Series(np.nan, index=values.index, dtype=float)
    valid = values.notna() & np.isfinite(values)
    if not valid.any():
        return result
    p = values.loc[valid].astype(float)
    order = np.argsort(p.to_numpy())
    ranked = p.to_numpy()[order]
    adjusted = ranked * len(ranked) / np.arange(1, len(ranked) + 1)
    adjusted = np.minimum.accumulate(adjusted[::-1])[::-1]
    adjusted = np.clip(adjusted, 0, 1)
    restored = np.empty_like(adjusted)
    restored[order] = adjusted
    result.loc[valid] = restored
    return result


def standard_error(log2_fc: float, ql_f: float) -> float:
    """Recover the effect-scale SE implied by a 1-df QL-F statistic."""
    if not np.isfinite(log2_fc) or not np.isfinite(ql_f) or ql_f <= 0:
        return np.nan
    return abs(log2_fc) / math.sqrt(ql_f)


def meta_row(frame: pd.DataFrame) -> dict[str, float | int]:
    usable = frame.loc[
        frame["log2_fc"].notna()
        & frame["derived_se"].notna()
        & np.isfinite(frame["derived_se"])
        & frame["derived_se"].gt(0)
    ].copy()
    k = len(usable)
    if k < 2:
        return {"meta_chip_n": k}

    effect = usable["log2_fc"].to_numpy(dtype=float)
    variance = np.square(usable["derived_se"].to_numpy(dtype=float))
    weight = 1.0 / variance
    fixed = float(np.sum(weight * effect) / np.sum(weight))
    fixed_se = float(math.sqrt(1.0 / np.sum(weight)))
    fixed_z = fixed / fixed_se
    fixed_p = float(2 * norm.sf(abs(fixed_z)))
    q = float(np.sum(weight * np.square(effect - fixed)))
    q_df = k - 1
    q_p = float(chi2.sf(q, q_df))
    i2 = float(max(0.0, (q - q_df) / q) * 100) if q > 0 else 0.0
    denominator = float(np.sum(weight) - np.sum(np.square(weight)) / np.sum(weight))
    tau2 = float(max(0.0, (q - q_df) / denominator)) if denominator > 0 else 0.0
    random_weight = 1.0 / (variance + tau2)
    random = float(np.sum(random_weight * effect) / np.sum(random_weight))
    random_se = float(math.sqrt(1.0 / np.sum(random_weight)))
    random_z = random / random_se
    random_p = float(2 * norm.sf(abs(random_z)))
    signs = np.sign(effect)

    return {
        "meta_chip_n": k,
        "meta_concordant_n": int(max(np.sum(signs > 0), np.sum(signs < 0))),
        "fixed_log2_fc": fixed,
        "fixed_se": fixed_se,
        "fixed_ci_low": fixed - 1.96 * fixed_se,
        "fixed_ci_high": fixed + 1.96 * fixed_se,
        "fixed_p_value": fixed_p,
        "heterogeneity_q": q,
        "heterogeneity_df": q_df,
        "heterogeneity_p_value": q_p,
        "i2_percent": i2,
        "tau2_dl": tau2,
        "random_log2_fc": random,
        "random_se": random_se,
        "random_ci_low": random - 1.96 * random_se,
        "random_ci_high": random + 1.96 * random_se,
        "random_p_value": random_p,
    }


def classify(row: pd.Series) -> str:
    def safe_int(name: str) -> int:
        value = row.get(name, np.nan)
        return int(value) if pd.notna(value) else 0

    observed = safe_int("meta_chip_n")
    all_same = observed >= 2 and safe_int("joint_direction_chip_n") == observed
    joint_global = bool(row.get("joint_global_coding_fdr", 1) < FDR_CUTOFF)
    meta_global = bool(row.get("fixed_global_meta_fdr", 1) < FDR_CUTOFF)
    meta_within = bool(row.get("fixed_within_family_contrast_fdr", 1) < FDR_CUTOFF)
    chip_global = safe_int("chip_global_fdr_n")
    material_opposite = safe_int("material_opposite_chip_n")
    loo_available = bool(row.get("loo_available", False))
    loo_flips = row.get("loo_sign_flip_n", np.nan)
    stable_loo = loo_available and pd.notna(loo_flips) and float(loo_flips) == 0
    i2 = row.get("i2_percent", np.nan)
    low_heterogeneity = pd.notna(i2) and float(i2) < 50
    high_heterogeneity = pd.notna(i2) and float(i2) >= 75

    if all_same and joint_global and meta_global and stable_loo and low_heterogeneity:
        return "Cross-chip robust"
    if all_same and (joint_global or meta_global) and (stable_loo or not loo_available) and not high_heterogeneity:
        return "Cross-chip supported"
    if material_opposite > 0 or high_heterogeneity or (loo_available and float(loo_flips) > 0):
        return "Heterogeneous or patient-sensitive"
    if all_same and (joint_global or meta_global or meta_within or chip_global >= 1):
        return "Directionally consistent"
    if chip_global >= 1:
        return "Chip-contained support"
    return "Unresolved / within-model only"


def load_inputs(source: Path) -> tuple[pd.DataFrame, ...]:
    family_map = pd.read_csv(source / "01_FACTOR_FAMILY_MAP.tsv", sep="\t")
    family = pd.read_csv(source / "22_FAMILY_CODING_RNA_DISEASE_TESTS.tsv.gz", sep="\t")
    chip = pd.read_csv(source / "23_CHIP_FACTOR_CODING_RNA_DISEASE_TESTS.tsv.gz", sep="\t")
    candidates = pd.read_csv(source / "24_FAMILY_CODING_RNA_CANDIDATES.tsv", sep="\t")
    loo = pd.read_csv(source / "19_GLOBAL_SIGNAL_CROSS_CHIP_LOO_AUDIT.tsv", sep="\t")
    family = family.loc[family["endpoint"].eq(ENDPOINT) & family["coding_rna"].astype(bool)].copy()
    chip = chip.loc[chip["endpoint"].eq(ENDPOINT) & chip["coding_rna"].astype(bool)].copy()
    candidates = candidates.loc[candidates["endpoint"].eq(ENDPOINT)].copy()
    return family_map, family, chip, candidates, loo


def build_meta_universe(
    family_map: pd.DataFrame, chip: pd.DataFrame
) -> tuple[pd.DataFrame, pd.DataFrame]:
    chip = chip.copy()
    chip["derived_se"] = [standard_error(a, b) for a, b in zip(chip["log2_fc"], chip["ql_f"])]
    long_parts: list[pd.DataFrame] = []
    for family_id, mappings in family_map.groupby("family_id", sort=False):
        mapping_lookup = dict(zip(mappings["chip"], mappings["factor"]))
        pieces = []
        for chip_name, factor in mapping_lookup.items():
            unit_id = f"{chip_name}::{factor}"
            part = chip.loc[chip["unit_id"].eq(unit_id)].copy()
            part["family_id"] = family_id
            part["family_label"] = mappings["family_label"].iloc[0]
            part["family_scope"] = mappings["family_scope"].iloc[0]
            part["chip"] = chip_name
            part["factor"] = factor
            pieces.append(part)
        family_long = pd.concat(pieces, ignore_index=True)
        long_parts.append(family_long)
    long = pd.concat(long_parts, ignore_index=True)
    long = long.loc[
        long["derived_se"].notna()
        & np.isfinite(long["derived_se"])
        & long["derived_se"].gt(0)
    ].copy()
    long["variance"] = np.square(long["derived_se"])
    long["weight"] = 1.0 / long["variance"]
    long["weight_effect"] = long["weight"] * long["log2_fc"]
    long["weight_effect2"] = long["weight"] * np.square(long["log2_fc"])
    long["weight2"] = np.square(long["weight"])
    long["positive"] = long["log2_fc"].gt(0).astype(int)
    long["negative"] = long["log2_fc"].lt(0).astype(int)
    group_cols = [
        "family_id",
        "family_label",
        "family_scope",
        "comparison",
        "ensembl_id",
        "gene_symbol",
    ]
    meta = long.groupby(group_cols, as_index=False).agg(
        meta_chip_n=("chip", "nunique"),
        sum_weight=("weight", "sum"),
        sum_weight_effect=("weight_effect", "sum"),
        sum_weight_effect2=("weight_effect2", "sum"),
        sum_weight2=("weight2", "sum"),
        positive_n=("positive", "sum"),
        negative_n=("negative", "sum"),
    )
    meta = meta.loc[meta["meta_chip_n"].ge(2)].copy()
    meta["mapped_chip_n"] = meta["family_id"].map(family_map.groupby("family_id").size())
    meta["meta_concordant_n"] = meta[["positive_n", "negative_n"]].max(axis=1)
    meta["fixed_log2_fc"] = meta["sum_weight_effect"] / meta["sum_weight"]
    meta["fixed_se"] = np.sqrt(1.0 / meta["sum_weight"])
    meta["fixed_ci_low"] = meta["fixed_log2_fc"] - 1.96 * meta["fixed_se"]
    meta["fixed_ci_high"] = meta["fixed_log2_fc"] + 1.96 * meta["fixed_se"]
    meta["fixed_p_value"] = 2 * norm.sf(np.abs(meta["fixed_log2_fc"] / meta["fixed_se"]))
    meta["heterogeneity_q"] = (
        meta["sum_weight_effect2"]
        - np.square(meta["sum_weight_effect"]) / meta["sum_weight"]
    ).clip(lower=0)
    meta["heterogeneity_df"] = meta["meta_chip_n"] - 1
    meta["heterogeneity_p_value"] = chi2.sf(
        meta["heterogeneity_q"], meta["heterogeneity_df"]
    )
    meta["i2_percent"] = np.where(
        meta["heterogeneity_q"].gt(0),
        np.maximum(
            0,
            (meta["heterogeneity_q"] - meta["heterogeneity_df"])
            / meta["heterogeneity_q"],
        )
        * 100,
        0,
    )
    dl_denominator = meta["sum_weight"] - meta["sum_weight2"] / meta["sum_weight"]
    meta["tau2_dl"] = np.where(
        dl_denominator.gt(0),
        np.maximum(0, (meta["heterogeneity_q"] - meta["heterogeneity_df"]) / dl_denominator),
        0,
    )

    long = long.merge(
        meta[group_cols + ["tau2_dl"]], on=group_cols, how="inner"
    )
    long["random_weight"] = 1.0 / (long["variance"] + long["tau2_dl"])
    long["random_weight_effect"] = long["random_weight"] * long["log2_fc"]
    random_summary = long.groupby(group_cols, as_index=False).agg(
        sum_random_weight=("random_weight", "sum"),
        sum_random_weight_effect=("random_weight_effect", "sum"),
    )
    meta = meta.merge(random_summary, on=group_cols, how="left")
    meta["random_log2_fc"] = meta["sum_random_weight_effect"] / meta["sum_random_weight"]
    meta["random_se"] = np.sqrt(1.0 / meta["sum_random_weight"])
    meta["random_ci_low"] = meta["random_log2_fc"] - 1.96 * meta["random_se"]
    meta["random_ci_high"] = meta["random_log2_fc"] + 1.96 * meta["random_se"]
    meta["random_p_value"] = 2 * norm.sf(np.abs(meta["random_log2_fc"] / meta["random_se"]))
    meta["fixed_within_family_contrast_fdr"] = meta.groupby(
        ["family_id", "comparison"], group_keys=False
    )["fixed_p_value"].transform(bh_adjust)
    meta["random_within_family_contrast_fdr"] = meta.groupby(
        ["family_id", "comparison"], group_keys=False
    )["random_p_value"].transform(bh_adjust)
    meta["fixed_global_meta_fdr"] = bh_adjust(meta["fixed_p_value"])
    meta["random_global_meta_fdr"] = bh_adjust(meta["random_p_value"])
    meta = meta.drop(
        columns=[
            "sum_weight",
            "sum_weight_effect",
            "sum_weight_effect2",
            "sum_weight2",
            "positive_n",
            "negative_n",
            "sum_random_weight",
            "sum_random_weight_effect",
        ]
    )
    return meta, long


def make_candidate_closeout(
    family_map: pd.DataFrame,
    family: pd.DataFrame,
    chip_long: pd.DataFrame,
    candidates: pd.DataFrame,
    loo: pd.DataFrame,
    meta: pd.DataFrame,
) -> pd.DataFrame:
    keys = ["family_id", "comparison", "ensembl_id", "gene_symbol"]
    candidates = candidates.rename(columns={"unit_id": "family_id"}).copy()
    candidates = candidates.drop_duplicates(keys)
    joint = family.rename(
        columns={
            "unit_id": "family_id",
            "log2_fc": "joint_log2_fc",
            "p_value": "joint_p_value",
            "within_model_coding_fdr": "joint_within_model_coding_fdr",
            "global_coding_fdr_within_level_endpoint": "joint_global_coding_fdr",
        }
    )
    joint_cols = keys + [
        "joint_log2_fc",
        "joint_p_value",
        "joint_within_model_coding_fdr",
        "joint_global_coding_fdr",
    ]
    base = candidates[keys].merge(joint[joint_cols], on=keys, how="left")
    base = base.merge(
        family_map[["family_id", "family_label", "family_scope"]].drop_duplicates(),
        on="family_id",
        how="left",
    )
    base = base.merge(meta, on=keys + ["family_label", "family_scope"], how="left", suffixes=("", "_meta"))

    joint_sign = np.sign(base.set_index(keys)["joint_log2_fc"])
    chip_rows = []
    for row in base.itertuples(index=False):
        mapping = family_map.loc[family_map["family_id"].eq(row.family_id)]
        observed = 0
        supportive = 0
        opposite = 0
        neutral = 0
        chip_global = 0
        chip_within = 0
        rec = {"family_id": row.family_id, "comparison": row.comparison, "ensembl_id": row.ensembl_id, "gene_symbol": row.gene_symbol}
        for chip_name in CHIPS:
            m = mapping.loc[mapping["chip"].eq(chip_name)]
            if m.empty:
                rec[f"{chip_name}_factor"] = ""
                rec[f"{chip_name}_log2_fc"] = np.nan
                rec[f"{chip_name}_se"] = np.nan
                rec[f"{chip_name}_ci_low"] = np.nan
                rec[f"{chip_name}_ci_high"] = np.nan
                rec[f"{chip_name}_global_coding_fdr"] = np.nan
                rec[f"{chip_name}_status"] = "N/A"
                continue
            factor = m["factor"].iloc[0]
            hit = chip_long.loc[
                chip_long["unit_id"].eq(f"{chip_name}::{factor}")
                & chip_long["comparison"].eq(row.comparison)
                & chip_long["ensembl_id"].eq(row.ensembl_id)
            ]
            rec[f"{chip_name}_factor"] = factor
            if hit.empty:
                rec[f"{chip_name}_log2_fc"] = np.nan
                rec[f"{chip_name}_se"] = np.nan
                rec[f"{chip_name}_ci_low"] = np.nan
                rec[f"{chip_name}_ci_high"] = np.nan
                rec[f"{chip_name}_global_coding_fdr"] = np.nan
                rec[f"{chip_name}_status"] = "Filtered"
                continue
            value = hit.iloc[0]
            effect = float(value["log2_fc"])
            se = float(value["derived_se"])
            observed += 1
            if abs(effect) < MIN_ABS_LOG2FC:
                neutral += 1
            elif np.sign(effect) == np.sign(row.joint_log2_fc):
                supportive += 1
            else:
                opposite += 1
            global_fdr = float(value["global_coding_fdr_within_level_endpoint"])
            within_fdr = float(value["within_model_coding_fdr"])
            chip_global += int(global_fdr < FDR_CUTOFF and abs(effect) >= MIN_ABS_LOG2FC)
            chip_within += int(within_fdr < FDR_CUTOFF and abs(effect) >= MIN_ABS_LOG2FC)
            status = "G" if global_fdr < FDR_CUTOFF and abs(effect) >= MIN_ABS_LOG2FC else (
                "W" if within_fdr < FDR_CUTOFF and abs(effect) >= MIN_ABS_LOG2FC else "ns"
            )
            rec[f"{chip_name}_log2_fc"] = effect
            rec[f"{chip_name}_se"] = se
            rec[f"{chip_name}_ci_low"] = effect - 1.96 * se
            rec[f"{chip_name}_ci_high"] = effect + 1.96 * se
            rec[f"{chip_name}_global_coding_fdr"] = global_fdr
            rec[f"{chip_name}_status"] = status
        rec.update(
            {
                "observed_chip_n": observed,
                "joint_direction_chip_n": supportive + neutral,
                "material_supportive_chip_n": supportive,
                "material_opposite_chip_n": opposite,
                "material_neutral_chip_n": neutral,
                "chip_global_fdr_n": chip_global,
                "chip_within_model_fdr_n": chip_within,
            }
        )
        chip_rows.append(rec)
    base = base.merge(pd.DataFrame(chip_rows), on=keys, how="left")

    loo_cols = [
        "family_id",
        "comparison",
        "gene_symbol",
        "loo_min_effect",
        "loo_max_effect",
        "loo_sign_flip_n",
        "most_influential_patient",
        "effect_without_most_influential",
        "max_abs_loo_shift",
    ]
    base = base.merge(loo[loo_cols], on=["family_id", "comparison", "gene_symbol"], how="left")
    base["loo_available"] = base["loo_sign_flip_n"].notna()
    base["all_observed_directionally_compatible"] = (
        base["observed_chip_n"].ge(2)
        & base["joint_direction_chip_n"].eq(base["observed_chip_n"])
    )
    base["joint_primary_supported"] = base["joint_global_coding_fdr"].lt(FDR_CUTOFF)
    base["fixed_meta_sensitivity_supported"] = base["fixed_global_meta_fdr"].lt(FDR_CUTOFF)
    base["high_heterogeneity_flag"] = base["i2_percent"].ge(75)
    base["patient_loo_instability_flag"] = (
        base["loo_available"] & base["loo_sign_flip_n"].gt(0)
    )
    base["evidence_tier"] = base.apply(classify, axis=1)
    order = {
        "Cross-chip robust": 0,
        "Cross-chip supported": 1,
        "Directionally consistent": 2,
        "Chip-contained support": 3,
        "Heterogeneous or patient-sensitive": 4,
        "Unresolved / within-model only": 5,
    }
    base["tier_order"] = base["evidence_tier"].map(order)
    base = base.sort_values(["tier_order", "family_scope", "family_label", "comparison", "joint_global_coding_fdr", "gene_symbol"])
    return base.drop(columns="tier_order")


def plot_joint_vs_meta(closeout: pd.DataFrame, path: Path) -> None:
    fig, ax = plt.subplots(figsize=(8.5, 7.2))
    scatter = ax.scatter(
        closeout["joint_log2_fc"],
        closeout["fixed_log2_fc"],
        c=closeout["i2_percent"],
        cmap="viridis",
        vmin=0,
        vmax=100,
        s=70,
        edgecolor="white",
        linewidth=0.8,
    )
    lim = max(1.0, np.nanmax(np.abs(closeout[["joint_log2_fc", "fixed_log2_fc"]].to_numpy()))) * 1.08
    ax.plot([-lim, lim], [-lim, lim], linestyle="--", color="#718096", linewidth=1)
    ax.axhline(0, color="#a0aec0", linewidth=0.8)
    ax.axvline(0, color="#a0aec0", linewidth=0.8)
    for row in closeout.itertuples(index=False):
        if row.joint_global_coding_fdr < 0.05 or row.fixed_global_meta_fdr < 0.05:
            ax.annotate(row.gene_symbol, (row.joint_log2_fc, row.fixed_log2_fc), xytext=(3, 3), textcoords="offset points", fontsize=8)
    ax.set_xlim(-lim, lim)
    ax.set_ylim(-lim, lim)
    ax.set_xlabel("Joint family edgeR log2FC")
    ax.set_ylabel("Fixed-effect chip meta-analysis log2FC")
    ax.set_title("Joint model versus chip-wise meta-analysis")
    cbar = fig.colorbar(scatter, ax=ax)
    cbar.set_label("I² heterogeneity (%)")
    fig.tight_layout()
    fig.savefig(path, dpi=220)
    plt.close(fig)


def plot_candidate_forest(closeout: pd.DataFrame, path: Path) -> None:
    ordered = closeout.sort_values(["family_label", "comparison", "gene_symbol"]).reset_index(drop=True)
    labels = [f"{r.family_label} | {r.comparison.replace('_vs_', ' vs ')} | {r.gene_symbol}" for r in ordered.itertuples()]
    fig_height = max(8, 0.48 * len(ordered) + 2)
    fig, ax = plt.subplots(figsize=(11.5, fig_height))
    y = np.arange(len(ordered))
    offsets = {"K8": -0.18, "J2": 0.0, "L3": 0.18}
    colors = {"K8": "#2b6cb0", "J2": "#38a169", "L3": "#dd6b20"}
    for chip_name in CHIPS:
        effect = ordered[f"{chip_name}_log2_fc"].to_numpy(float)
        low = ordered[f"{chip_name}_ci_low"].to_numpy(float)
        high = ordered[f"{chip_name}_ci_high"].to_numpy(float)
        valid = np.isfinite(effect) & np.isfinite(low) & np.isfinite(high)
        ax.errorbar(
            effect[valid],
            y[valid] + offsets[chip_name],
            xerr=np.vstack((effect[valid] - low[valid], high[valid] - effect[valid])),
            fmt="o",
            color=colors[chip_name],
            markersize=4.5,
            linewidth=1,
            capsize=2,
            label=chip_name,
        )
    ax.errorbar(
        ordered["fixed_log2_fc"],
        y + 0.34,
        xerr=np.vstack((ordered["fixed_log2_fc"] - ordered["fixed_ci_low"], ordered["fixed_ci_high"] - ordered["fixed_log2_fc"])),
        fmt="D",
        color="#1a202c",
        markersize=4.5,
        linewidth=1.2,
        capsize=2,
        label="Fixed meta",
    )
    ax.axvline(0, color="#718096", linewidth=1)
    ax.set_yticks(y)
    ax.set_yticklabels(labels, fontsize=8)
    ax.invert_yaxis()
    ax.set_xlabel("log2 fold change (positive favors first disease)")
    ax.set_title("Chip-specific effects and fixed-effect meta estimates")
    ax.legend(ncol=4, loc="lower right")
    fig.tight_layout()
    fig.savefig(path, dpi=220)
    plt.close(fig)


def write_report(closeout: pd.DataFrame, meta_sig: pd.DataFrame, output: Path) -> None:
    tier = closeout["evidence_tier"].value_counts().to_dict()
    joint_sig = int((closeout["joint_global_coding_fdr"] < 0.05).sum())
    meta_sig_n = int((closeout["fixed_global_meta_fdr"] < 0.05).sum())
    all_direction = int((closeout["joint_direction_chip_n"] == closeout["observed_chip_n"]).sum())
    lines = [
        "# 第一阶段 edgeR 跨芯片证据收口（分析草案）",
        "",
        "## 范围与主终点",
        "",
        "- 仅复用冻结的 K=6 factor-family 映射、family edgeR、单芯片 edgeR 和既有患者 LOO；未重跑 RCTD、NMF、pseudobulk 或 edgeR。",
        "- 主终点固定为 `paired_high_minus_control`，即同一患者 factor-high 区域减 matched control 区域后的疾病差异。",
        "- 联合 family edgeR（控制芯片）仍为主分析；芯片级 meta-analysis 仅用于复现与异质性敏感性评估。",
        "",
        "## 新增方法",
        "",
        "- 由单芯片 edgeR 的 1-df QL-F 统计量近似恢复 effect-scale SE：`abs(log2FC) / sqrt(QL-F)`。",
        "- 计算 inverse-variance fixed-effect meta、DerSimonian-Laird random-effect meta、Cochran Q、I² 与 tau²。",
        "- 因每个 family 只有 2–3 张芯片，random-effect 和 I² 均不作为单独淘汰标准。",
        "- Meta FDR 在全部可估计 protein-coding family × comparison × gene 检验中统一校正；同时保留 family/contrast 内 FDR。",
        "",
        "## 数据概览",
        "",
        f"- 既有候选信号：{len(closeout)}。",
        f"- Joint global coding FDR < 0.05：{joint_sig}。",
        f"- Fixed meta global FDR < 0.05：{meta_sig_n}。",
        f"- 所有可估计芯片方向与 Joint 一致或效应中性：{all_direction}/{len(closeout)}。",
        f"- 全 meta 检验中 fixed global FDR < 0.05：{len(meta_sig)}。",
        "",
        "## 描述性证据等级",
        "",
    ]
    for name, count in tier.items():
        lines.append(f"- {name}: {count}")
    lines.extend(
        [
            "",
            "证据等级不是新的筛除门槛。联合模型、单芯片结果、meta、异质性和患者 LOO 应并列阅读；B-associated 等患者敏感信号继续保留，但不得表述为三芯片独立显著。",
            "",
            "## 关键限制",
            "",
            "- 芯片级 SE 来自已保存 QL-F 的近似反推，并非重新拟合 edgeR 后直接导出的协方差矩阵。",
            "- 仅 2–3 张芯片时，I²、tau² 和 random-effect 估计很不稳定。",
            "- Meta 与 joint family 模型使用相同患者数据，是互补敏感性分析，不是独立验证队列。",
            "- 候选表来自既有 edgeR 候选规则；完整 meta 结果另存，避免只报告被 Joint 预筛选的基因。",
        ]
    )
    (output / "06_FIRST_STAGE_CLOSEOUT_DRAFT_ZH.md").write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--project-root", type=Path, required=True)
    parser.add_argument("--source", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()

    root = args.project_root.resolve()
    source = (root / args.source).resolve() if not args.source.is_absolute() else args.source.resolve()
    output = (root / args.output).resolve() if not args.output.is_absolute() else args.output.resolve()
    figures = output / "figures"
    figures.mkdir(parents=True, exist_ok=False)

    family_map, family, chip, candidates, loo = load_inputs(source)
    meta, chip_long = build_meta_universe(family_map, chip)
    closeout = make_candidate_closeout(family_map, family, chip_long, candidates, loo, meta)
    meta_sig = meta.loc[meta["fixed_global_meta_fdr"].lt(FDR_CUTOFF)].copy()

    manifest = pd.DataFrame(
        {
            "input": [
                "01_FACTOR_FAMILY_MAP.tsv",
                "22_FAMILY_CODING_RNA_DISEASE_TESTS.tsv.gz",
                "23_CHIP_FACTOR_CODING_RNA_DISEASE_TESTS.tsv.gz",
                "24_FAMILY_CODING_RNA_CANDIDATES.tsv",
                "19_GLOBAL_SIGNAL_CROSS_CHIP_LOO_AUDIT.tsv",
            ],
            "source_directory": str(source),
        }
    )
    manifest.to_csv(output / "00_INPUT_MANIFEST.tsv", sep="\t", index=False)
    meta.to_csv(output / "01_META_ALL_CODING_TESTS.tsv.gz", sep="\t", index=False, compression="gzip")
    closeout.to_csv(output / "02_CANDIDATE_EVIDENCE_CLOSEOUT.tsv", sep="\t", index=False)
    meta_sig.to_csv(output / "03_META_SIGNIFICANT_SENSITIVITY.tsv", sep="\t", index=False)
    closeout.groupby(["family_scope", "evidence_tier"], dropna=False).size().reset_index(name="candidate_n").to_csv(
        output / "04_EVIDENCE_TIER_SUMMARY.tsv", sep="\t", index=False
    )

    qa = pd.DataFrame(
        [
            ("endpoint_is_paired", True, ENDPOINT),
            ("candidate_count_matches_existing", len(closeout) == len(candidates.drop_duplicates(["unit_id", "comparison", "ensembl_id", "gene_symbol"])), len(closeout)),
            ("all_meta_have_at_least_two_chips", bool(meta["meta_chip_n"].ge(2).all()), int(meta["meta_chip_n"].min())),
            ("all_meta_se_positive", bool((meta["fixed_se"] > 0).all()), int(meta["fixed_se"].isna().sum())),
            ("all_q_nonnegative", bool(meta["heterogeneity_q"].ge(0).all()), float(meta["heterogeneity_q"].min())),
            ("all_i2_in_range", bool(meta["i2_percent"].between(0, 100).all()), f"{meta['i2_percent'].min():.3f} to {meta['i2_percent'].max():.3f}"),
            ("no_upstream_refit", True, "read-only frozen edgeR tables"),
        ],
        columns=["check", "passed", "detail"],
    )
    qa.to_csv(output / "05_QA.tsv", sep="\t", index=False)
    plot_joint_vs_meta(closeout, figures / "01_JOINT_VS_FIXED_META.png")
    plot_candidate_forest(closeout, figures / "02_CANDIDATE_CHIP_META_FOREST.png")
    write_report(closeout, meta_sig, output)

    session = {
        "python": platform.python_version(),
        "pandas": pd.__version__,
        "numpy": np.__version__,
        "scipy": __import__("scipy").__version__,
        "matplotlib": __import__("matplotlib").__version__,
        "source": str(source),
        "output": str(output),
        "candidate_n": int(len(closeout)),
        "meta_test_n": int(len(meta)),
        "meta_global_fdr_n": int(len(meta_sig)),
    }
    (output / "SESSION_INFO.json").write_text(json.dumps(session, indent=2), encoding="utf-8")
    (output / "ANALYSIS_COMPLETE.ok").write_text("complete\n", encoding="utf-8")
    print(output)


if __name__ == "__main__":
    main()
