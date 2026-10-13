"""Build a disease-blind external Reference18 atlas around frozen NMF factors.

The analysis reuses frozen per-chip K=6 neighborhood NMF usages and frozen
per-Bin50 Reference18 proportions.  It never refits RCTD or NMF.  The primary
external compartment is the 100--200 um annulus around each eligible Bin50.
Inference is performed from patient-level correlations after 200 um spatial
block aggregation; disease labels are carried only as provenance.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import platform
import sys
import time
from dataclasses import dataclass
from datetime import datetime
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scipy
from scipy.spatial import cKDTree
from scipy.stats import rankdata, spearmanr, wilcoxon


SEED = 20261013
CHIPS = ("K8", "J2", "L3")
RINGS = ((0.0, 50.0), (50.0, 100.0), (100.0, 200.0), (200.0, 300.0))
BLOCK_SIZES = (150.0, 200.0, 250.0)
PRIMARY_RING = "100-200"
PRIMARY_BLOCK = 200.0
MIN_RING_BINS = {"0-50": 3, "50-100": 8, "100-200": 20, "200-300": 30}
MIN_SECTORS = {"0-50": 2, "50-100": 3, "100-200": 4, "200-300": 4}
MIN_CENTERS_PER_BLOCK = 5
MIN_BLOCKS_PER_PATIENT = 8

FACTOR_MAP = {
    "Ciliated": {"K8": 3, "J2": 4, "L3": 2},
    "Fibroblast": {"K8": 4, "J2": 6, "L3": 5},
    "Lymphatic": {"K8": 1, "J2": 2, "L3": 4},
}


@dataclass
class ChipData:
    chip: str
    frame: pd.DataFrame
    state_names: list[str]
    usage: np.ndarray


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def bh_adjust(values: pd.Series) -> np.ndarray:
    p = values.to_numpy(dtype=float)
    valid = np.isfinite(p)
    answer = np.full(len(p), np.nan)
    if not valid.any():
        return answer
    pv = p[valid]
    order = np.argsort(pv)
    ranked = pv[order]
    adjusted = ranked * len(ranked) / np.arange(1, len(ranked) + 1)
    adjusted = np.minimum.accumulate(adjusted[::-1])[::-1]
    restored = np.empty_like(adjusted)
    restored[order] = np.minimum(adjusted, 1.0)
    answer[np.flatnonzero(valid)] = restored
    return answer


def safe_spearman(x: np.ndarray, y: np.ndarray) -> float:
    mask = np.isfinite(x) & np.isfinite(y)
    if mask.sum() < 5 or np.nanstd(x[mask]) == 0 or np.nanstd(y[mask]) == 0:
        return np.nan
    return float(spearmanr(x[mask], y[mask]).statistic)


def residualize_rank(values: np.ndarray, covariates: np.ndarray) -> np.ndarray:
    y = rankdata(values, method="average").astype(float)
    x = np.asarray(covariates, dtype=float)
    x = np.column_stack([np.ones(len(x)), x])
    scale = np.nanstd(x[:, 1:], axis=0)
    scale[scale == 0] = 1.0
    x[:, 1:] = (x[:, 1:] - np.nanmean(x[:, 1:], axis=0)) / scale
    coef = np.linalg.lstsq(x, y, rcond=None)[0]
    return y - x @ coef


def load_chip_data(project: Path, chip: str) -> tuple[ChipData, list[Path]]:
    frozen_root = project / "results" / "reanalysis" / (
        "REFERENCE18_UMI50_VOLTRON_PER_CHIP_NICHE_20261005_153036"
    )
    model_root = project / "results" / "reanalysis" / (
        "REFERENCE18_PER_CHIP_DUAL_NMF_20261012_204500_attempt4_alpha0"
    ) / "models"
    frame_path = frozen_root / chip / "01_FROZEN_REFERENCE18_INPUT.tsv.gz"
    model_path = model_root / f"neighborhood50_{chip}_K06_MODEL_AND_USAGE.npz"
    frame = pd.read_csv(frame_path, sep="\t")
    weight_columns = [column for column in frame if column.startswith("weight__")]
    state_names = [column.removeprefix("weight__") for column in weight_columns]
    with np.load(model_path, allow_pickle=True) as model:
        keys = model["entity_key"].astype(str)
        usage = model["relative_usage"].astype(np.float64)
    index = pd.Series(np.arange(len(frame)), index=frame["entity_key"].astype(str))
    if index.index.has_duplicates:
        raise RuntimeError(f"{chip}: duplicate entity_key in frozen input")
    locations = index.reindex(keys)
    if locations.isna().any():
        raise RuntimeError(f"{chip}: model keys missing from frozen input")
    eligible = frame.iloc[locations.astype(int).to_numpy()].copy().reset_index(drop=True)
    if eligible["entity_key"].astype(str).tolist() != keys.tolist():
        raise RuntimeError(f"{chip}: entity alignment failed")
    if usage.shape != (len(eligible), 6):
        raise RuntimeError(f"{chip}: unexpected usage shape {usage.shape}")
    if not np.allclose(frame[weight_columns].sum(axis=1), 1.0, atol=1e-5):
        raise RuntimeError(f"{chip}: Reference18 rows do not sum to one")
    eligible["model_row"] = np.arange(len(eligible))
    return ChipData(chip, frame, state_names, usage), [frame_path, model_path]


def annular_composition(
    all_frame: pd.DataFrame,
    centers: pd.DataFrame,
    weight_columns: list[str],
) -> tuple[dict[str, np.ndarray], list[dict[str, object]]]:
    n = len(centers)
    outputs = {
        f"{int(inner)}-{int(outer)}": np.full((n, len(weight_columns)), np.nan, dtype=np.float32)
        for inner, outer in RINGS
    }
    counts = {name: np.zeros(n, dtype=np.int32) for name in outputs}
    sectors = {name: np.zeros(n, dtype=np.int8) for name in outputs}
    qc_rows: list[dict[str, object]] = []
    group_columns = ["patient_id", "roi_id", "tissue_component_id"]
    all_groups = {key: group for key, group in all_frame.groupby(group_columns, sort=False)}
    for key, center_group in centers.groupby(group_columns, sort=False):
        tissue = all_groups.get(key)
        if tissue is None:
            raise RuntimeError(f"Missing tissue group {key}")
        tissue_xy = tissue[["x_um", "y_um"]].to_numpy(dtype=float)
        tissue_w = tissue[weight_columns].to_numpy(dtype=np.float64)
        center_xy = center_group[["x_um", "y_um"]].to_numpy(dtype=float)
        center_rows = center_group["model_row"].to_numpy(dtype=int)
        tree = cKDTree(tissue_xy)
        for start in range(0, len(center_group), 512):
            stop = min(start + 512, len(center_group))
            neighbors = tree.query_ball_point(center_xy[start:stop], 300.0 + 1e-8)
            for local_index, neighbor_index in enumerate(neighbors):
                row = center_rows[start + local_index]
                idx = np.asarray(neighbor_index, dtype=int)
                delta = tissue_xy[idx] - center_xy[start + local_index]
                distance = np.sqrt(np.square(delta).sum(axis=1))
                angle = np.mod(np.arctan2(delta[:, 1], delta[:, 0]), 2 * np.pi)
                for inner, outer in RINGS:
                    name = f"{int(inner)}-{int(outer)}"
                    mask = (distance > inner + 1e-8) & (distance <= outer + 1e-8)
                    if not mask.any():
                        continue
                    selected = idx[mask]
                    outputs[name][row] = tissue_w[selected].mean(axis=0).astype(np.float32)
                    counts[name][row] = int(mask.sum())
                    sectors[name][row] = int(
                        np.unique(np.floor(angle[mask] / (np.pi / 4)).astype(int)).size
                    )
    for name in outputs:
        centers[f"ring_n__{name}"] = counts[name]
        centers[f"ring_sectors__{name}"] = sectors[name]
        centers[f"ring_valid__{name}"] = (
            (counts[name] >= MIN_RING_BINS[name])
            & (sectors[name] >= MIN_SECTORS[name])
        )
        for (patient, roi, component), group in centers.groupby(group_columns, sort=False):
            idx = group.index.to_numpy()
            qc_rows.append(
                {
                    "chip": str(group["chip_label"].iloc[0]),
                    "patient_id": patient,
                    "roi_id": roi,
                    "tissue_component_id": component,
                    "ring_um": name,
                    "eligible_centers": len(group),
                    "valid_centers": int(centers.loc[idx, f"ring_valid__{name}"].sum()),
                    "valid_fraction": float(centers.loc[idx, f"ring_valid__{name}"].mean()),
                    "median_ring_bins": float(np.median(counts[name][idx])),
                    "median_angular_sectors": float(np.median(sectors[name][idx])),
                }
            )
    return outputs, qc_rows


def block_associations(
    chip_data: ChipData,
    centers: pd.DataFrame,
    annuli: dict[str, np.ndarray],
) -> tuple[list[dict[str, object]], list[dict[str, object]]]:
    association_rows: list[dict[str, object]] = []
    contrast_rows: list[dict[str, object]] = []
    for family, mapping in FACTOR_MAP.items():
        factor = mapping[chip_data.chip]
        family_centers = centers.copy()
        family_centers["usage"] = chip_data.usage[:, factor - 1]
        for ring_name, composition in annuli.items():
            valid_column = f"ring_valid__{ring_name}"
            for patient, patient_frame in family_centers.groupby("patient_id", sort=True):
                rows = patient_frame.index.to_numpy()
                valid = patient_frame[valid_column].to_numpy(dtype=bool)
                rows = rows[valid]
                if len(rows) < 50:
                    continue
                patient_valid = family_centers.loc[rows].copy()
                patient_comp = composition[rows]
                q90 = float(patient_valid["usage"].quantile(0.90))
                q40 = float(patient_valid["usage"].quantile(0.40))
                q60 = float(patient_valid["usage"].quantile(0.60))
                high = patient_valid["usage"].to_numpy() >= q90
                middle = patient_valid["usage"].between(q40, q60, inclusive="both").to_numpy()
                for state_index, state in enumerate(chip_data.state_names):
                    contrast_rows.append(
                        {
                            "chip": chip_data.chip,
                            "patient_id": patient,
                            "disease_provenance_only": patient_valid["disease"].iloc[0],
                            "factor_family": family,
                            "factor": f"F{factor:02d}",
                            "ring_um": ring_name,
                            "state": state,
                            "high_q90_n": int(high.sum()),
                            "middle_q40_q60_n": int(middle.sum()),
                            "high_minus_middle": float(
                                patient_comp[high, state_index].mean()
                                - patient_comp[middle, state_index].mean()
                            ),
                        }
                    )
                for block_size in BLOCK_SIZES:
                    work = patient_valid.copy()
                    work["block_x"] = np.floor(work["x_um"] / block_size).astype(int)
                    work["block_y"] = np.floor(work["y_um"] / block_size).astype(int)
                    work["row_position"] = np.arange(len(work))
                    block_keys = ["roi_id", "tissue_component_id", "block_x", "block_y"]
                    block_meta = work.groupby(block_keys, sort=False).agg(
                        center_n=("entity_key", "size"),
                        usage=("usage", "mean"),
                        nUMI=("nUMI", "mean"),
                        nGene=("nGene", "mean"),
                        neighbor_n=("ring_n__0-50", "mean"),
                        ring_n=(f"ring_n__{ring_name}", "mean"),
                    ).reset_index()
                    block_ids = work.groupby(block_keys, sort=False).ngroup().to_numpy()
                    comp_block = pd.DataFrame(patient_comp, columns=chip_data.state_names)
                    comp_block["block_id"] = block_ids
                    comp_means = comp_block.groupby("block_id", sort=False).mean().reset_index(drop=True)
                    block_meta = block_meta.reset_index(drop=True)
                    keep = block_meta["center_n"].to_numpy() >= MIN_CENTERS_PER_BLOCK
                    block_meta = block_meta.loc[keep].reset_index(drop=True)
                    comp_means = comp_means.loc[keep].reset_index(drop=True)
                    if len(block_meta) < MIN_BLOCKS_PER_PATIENT:
                        continue
                    covariates = np.column_stack(
                        [
                            np.log1p(block_meta["nUMI"].to_numpy()),
                            np.log1p(block_meta["nGene"].to_numpy()),
                            block_meta["neighbor_n"].to_numpy(),
                            np.log1p(block_meta["ring_n"].to_numpy()),
                        ]
                    )
                    usage_values = block_meta["usage"].to_numpy(dtype=float)
                    usage_residual = residualize_rank(usage_values, covariates)
                    for state in chip_data.state_names:
                        values = comp_means[state].to_numpy(dtype=float)
                        raw_rho = safe_spearman(usage_values, values)
                        state_residual = residualize_rank(values, covariates)
                        adjusted_rho = safe_spearman(usage_residual, state_residual)
                        association_rows.append(
                            {
                                "chip": chip_data.chip,
                                "patient_id": patient,
                                "disease_provenance_only": patient_valid["disease"].iloc[0],
                                "factor_family": family,
                                "factor": f"F{factor:02d}",
                                "ring_um": ring_name,
                                "block_size_um": int(block_size),
                                "state": state,
                                "valid_centers": len(patient_valid),
                                "spatial_blocks": len(block_meta),
                                "raw_spearman_rho": raw_rho,
                                "technical_adjusted_spearman_rho": adjusted_rho,
                            }
                        )
    return association_rows, contrast_rows


def summarize_primary(associations: pd.DataFrame) -> pd.DataFrame:
    primary = associations.loc[
        associations["ring_um"].eq(PRIMARY_RING)
        & associations["block_size_um"].eq(int(PRIMARY_BLOCK))
    ].copy()
    rows: list[dict[str, object]] = []
    for (family, state), group in primary.groupby(["factor_family", "state"], sort=True):
        values = group["technical_adjusted_spearman_rho"].dropna().to_numpy(dtype=float)
        raw = group["raw_spearman_rho"].dropna().to_numpy(dtype=float)
        if len(values) == 0:
            continue
        direction = 1.0 if np.median(values) >= 0 else -1.0
        try:
            p_value = float(wilcoxon(values, alternative="two-sided").pvalue)
        except ValueError:
            p_value = 1.0
        chip_medians = group.groupby("chip")["technical_adjusted_spearman_rho"].median()
        rows.append(
            {
                "factor_family": family,
                "state": state,
                "patients": len(values),
                "median_adjusted_rho": float(np.median(values)),
                "median_raw_rho": float(np.median(raw)) if len(raw) else np.nan,
                "patient_direction_fraction": float(np.mean(np.sign(values) == direction)),
                "wilcoxon_p": p_value,
                "chip_median_K8": chip_medians.get("K8", np.nan),
                "chip_median_J2": chip_medians.get("J2", np.nan),
                "chip_median_L3": chip_medians.get("L3", np.nan),
                "chips_same_direction": int(np.sum(np.sign(chip_medians) == direction)),
                "raw_adjusted_same_direction": bool(
                    np.sign(np.median(raw)) == np.sign(np.median(values))
                ) if len(raw) else False,
            }
        )
    summary = pd.DataFrame(rows)
    summary["bh_q_54_primary_tests"] = bh_adjust(summary["wilcoxon_p"])
    sensitivity = associations.loc[associations["ring_um"].eq(PRIMARY_RING)].groupby(
        ["factor_family", "state", "block_size_um"]
    )["technical_adjusted_spearman_rho"].median().unstack("block_size_um")
    sensitivity.columns = [f"median_rho_block_{int(column)}um" for column in sensitivity.columns]
    summary = summary.merge(sensitivity.reset_index(), on=["factor_family", "state"], how="left")
    block_columns = [f"median_rho_block_{int(size)}um" for size in BLOCK_SIZES]
    summary["block_scale_same_direction"] = summary.apply(
        lambda row: bool(
            all(np.sign(row[column]) == np.sign(row["median_adjusted_rho"]) for column in block_columns)
        ),
        axis=1,
    )
    summary["direction_interpretation"] = np.where(
        summary["median_adjusted_rho"] > 0,
        "external enrichment/continuity",
        "relative compositional depletion; closure-sensitive",
    )
    summary["self_state_continuity"] = (
        ((summary["factor_family"] == "Ciliated") & (summary["state"] == "Ciliated epithelial"))
        | ((summary["factor_family"] == "Fibroblast") & (summary["state"] == "Fibroblast"))
        | ((summary["factor_family"] == "Lymphatic") & (summary["state"] == "Lymphatic endothelial"))
    )
    def classify(row: pd.Series) -> str:
        if (
            row["bh_q_54_primary_tests"] < 0.05
            and row["chips_same_direction"] == 3
            and row["patient_direction_fraction"] >= 2 / 3
            and row["raw_adjusted_same_direction"]
            and row["block_scale_same_direction"]
        ):
            return "Cross-chip reproducible external association"
        if (
            row["wilcoxon_p"] < 0.05
            and row["chips_same_direction"] >= 2
            and row["patient_direction_fraction"] >= 0.60
            and row["block_scale_same_direction"]
        ):
            return "Partially reproducible external association"
        return "Unresolved external association"
    summary["evidence_class"] = summary.apply(classify, axis=1)
    return summary.sort_values(["factor_family", "bh_q_54_primary_tests", "state"])


def create_figures(summary: pd.DataFrame, contrasts: pd.DataFrame, output: Path) -> None:
    figure_dir = output / "figures"
    figure_dir.mkdir(parents=True, exist_ok=True)
    state_order = summary.groupby("state")["median_adjusted_rho"].apply(lambda x: np.nanmax(np.abs(x))).sort_values(ascending=False).index
    matrix = summary.pivot(index="factor_family", columns="state", values="median_adjusted_rho").reindex(columns=state_order)
    fig, ax = plt.subplots(figsize=(15, 4.3))
    image = ax.imshow(matrix.to_numpy(), aspect="auto", cmap="RdBu_r", vmin=-0.6, vmax=0.6)
    ax.set_xticks(np.arange(len(matrix.columns)), matrix.columns, rotation=55, ha="right", fontsize=8)
    ax.set_yticks(np.arange(len(matrix.index)), matrix.index)
    ax.set_title("External Reference18 associations (100–200 µm; patient-level median adjusted rho)")
    fig.colorbar(image, ax=ax, label="Median adjusted Spearman rho", shrink=0.8)
    fig.tight_layout()
    fig.savefig(figure_dir / "01_EXTERNAL_ASSOCIATION_HEATMAP.png", dpi=220)
    plt.close(fig)

    top = summary.assign(abs_rho=summary["median_adjusted_rho"].abs()).sort_values(
        ["factor_family", "abs_rho"], ascending=[True, False]
    ).groupby("factor_family").head(5)
    selected = set(zip(top["factor_family"], top["state"]))
    radial = contrasts.groupby(["factor_family", "state", "ring_um"])["high_minus_middle"].median().reset_index()
    radial = radial.loc[[pair in selected for pair in zip(radial["factor_family"], radial["state"])]]
    ring_order = [f"{int(a)}-{int(b)}" for a, b in RINGS]
    fig, axes = plt.subplots(1, 3, figsize=(16, 4.8), sharey=True)
    for ax, family in zip(axes, FACTOR_MAP):
        family_data = radial.loc[radial["factor_family"].eq(family)]
        for state, group in family_data.groupby("state"):
            values = group.set_index("ring_um")["high_minus_middle"].reindex(ring_order)
            ax.plot(ring_order, values, marker="o", label=state)
        ax.axhline(0, color="black", linewidth=0.7)
        ax.axvspan(1.75, 2.25, color="#eeeeee", zorder=-2)
        ax.set_title(family)
        ax.set_xlabel("Annulus (µm)")
        ax.tick_params(axis="x", rotation=35)
        ax.legend(fontsize=7, frameon=False)
    axes[0].set_ylabel("Median q90-minus-middle Reference18 proportion")
    fig.suptitle("Radial external composition around high-usage factor centers")
    fig.tight_layout()
    fig.savefig(figure_dir / "02_RADIAL_EXTERNAL_COMPOSITION.png", dpi=220)
    plt.close(fig)

    candidate = summary.loc[summary["evidence_class"].ne("Unresolved external association")].copy()
    if candidate.empty:
        candidate = summary.assign(abs_rho=summary["median_adjusted_rho"].abs()).nlargest(12, "abs_rho")
    candidate["label"] = candidate["factor_family"] + " | " + candidate["state"]
    candidate = candidate.sort_values("median_adjusted_rho")
    fig, ax = plt.subplots(figsize=(9, max(4, 0.35 * len(candidate))))
    y = np.arange(len(candidate))
    ax.barh(y, candidate["median_adjusted_rho"], color=np.where(candidate["median_adjusted_rho"] >= 0, "#b2182b", "#2166ac"))
    ax.set_yticks(y, candidate["label"], fontsize=8)
    ax.axvline(0, color="black", linewidth=0.8)
    ax.set_xlabel("Median adjusted Spearman rho")
    ax.set_title("Candidate reproducible external microenvironment relationships")
    fig.tight_layout()
    fig.savefig(figure_dir / "03_CANDIDATE_EXTERNAL_RELATIONSHIPS.png", dpi=220)
    plt.close(fig)


def write_report(
    output: Path,
    summary: pd.DataFrame,
    qc: pd.DataFrame,
    associations: pd.DataFrame,
) -> None:
    confirmed = summary.loc[summary["evidence_class"].str.startswith("Cross-chip")]
    partial = summary.loc[summary["evidence_class"].str.startswith("Partially")]
    primary_assoc = associations.loc[
        associations["ring_um"].eq(PRIMARY_RING)
        & associations["block_size_um"].eq(int(PRIMARY_BLOCK))
    ]
    lines = [
        "# 三类非B NMF factor外部微环境图谱（疾病盲法）",
        "",
        "## 范围",
        "",
        "本分析只复用冻结的逐芯片K=6 Neighborhood NMF连续usage与逐Bin50 Reference18比例。",
        "未重新运行RCTD或NMF，未开展疾病比较，未分析B-associated factor，也未进行factor-factor接触分析。",
        "",
        "## 主方法",
        "",
        "- 主外环为100–200 μm；0–50 μm仅作内部正对照，50–100及200–300 μm为径向敏感性。",
        "- 外环只在同一patient、ROI及tissue_component内汇总，避免跨组织间隙错误连接。",
        "- 每个外环至少20个Bin50并覆盖8个角扇区中的至少4个，才进入100–200 μm主分析。",
        "- 先按150/200/250 μm空间块聚合，再在患者内计算factor usage与外环细胞状态比例的相关。",
        "- 主结果为200 μm块的患者级技术调整Spearman rho；调整nUMI、nGene、原50 μm有效邻居数和外环Bin数。",
        "- 21位患者是复现单位；疾病标签仅保留作溯源，未进入发现、筛选或检验。",
        "- 正向关系解释为外环富集/连续性；负向关系仅表示组成比例相对降低，受总和为1的闭合约束影响，不能直接解释为细胞排斥。",
        "",
        "## 数据完整性",
        "",
        f"- 100–200 μm主环有效中心中位覆盖率：{qc.loc[qc['ring_um'].eq(PRIMARY_RING), 'valid_fraction'].median():.3f}。",
        f"- 主分析患者×factor×state相关记录：{len(primary_assoc):,}。",
        f"- BH校正后跨芯片可复现关系：{len(confirmed)}；部分复现关系：{len(partial)}。",
        "",
        "## 证据分级",
        "",
    ]
    for label, data in (("跨芯片可复现", confirmed), ("部分复现", partial)):
        lines.append(f"### {label}")
        lines.append("")
        if data.empty:
            lines.append("无。")
        else:
            for row in data.itertuples():
                lines.append(
                    f"- {row.factor_family} → {row.state}: median adjusted rho={row.median_adjusted_rho:.3f}, "
                    f"患者同向比例={row.patient_direction_fraction:.2f}, 3个芯片同向数={row.chips_same_direction}, "
                    f"P={row.wilcoxon_p:.4g}, BH q={row.bh_q_54_primary_tests:.4g}."
                )
        lines.append("")
    lines.extend(
        [
            "## 解释边界",
            "",
            "这些结果描述的是factor连续usage较高处外围Reference18组成的重复关系，不等同于细胞直接接触、配体受体通讯或疾病效应。",
            "空间分块降低了重叠邻域造成的伪重复，但不是严格的组织形状保持空间零模型；因此未通过跨患者/跨芯片标准的关系保留为unresolved。",
        ]
    )
    (output / "08_METHOD_AND_RESULTS_ZH.md").write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--project-root", type=Path, default=Path(__file__).resolve().parents[1])
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    project = args.project_root.resolve()
    output = args.output.resolve()
    output.mkdir(parents=True, exist_ok=False)
    started = time.time()
    np.random.seed(SEED)

    input_paths: list[Path] = []
    all_qc: list[dict[str, object]] = []
    all_associations: list[dict[str, object]] = []
    all_contrasts: list[dict[str, object]] = []
    factor_rows: list[dict[str, object]] = []
    state_names_reference: list[str] | None = None
    for chip in CHIPS:
        chip_data, paths = load_chip_data(project, chip)
        input_paths.extend(paths)
        if state_names_reference is None:
            state_names_reference = chip_data.state_names
        elif chip_data.state_names != state_names_reference:
            raise RuntimeError("Reference18 state order differs by chip")
        model_path = paths[1]
        with np.load(model_path, allow_pickle=True) as model:
            keys = model["entity_key"].astype(str)
        frame_index = chip_data.frame.set_index(chip_data.frame["entity_key"].astype(str), drop=False)
        centers = frame_index.loc[keys].copy().reset_index(drop=True)
        centers["model_row"] = np.arange(len(centers))
        weight_columns = [f"weight__{state}" for state in chip_data.state_names]
        annuli, qc_rows = annular_composition(chip_data.frame, centers, weight_columns)
        association_rows, contrast_rows = block_associations(chip_data, centers, annuli)
        all_qc.extend(qc_rows)
        all_associations.extend(association_rows)
        all_contrasts.extend(contrast_rows)
        for family, mapping in FACTOR_MAP.items():
            factor_rows.append({"factor_family": family, "chip": chip, "factor": f"F{mapping[chip]:02d}"})

    provenance = pd.DataFrame({"path": [str(path) for path in input_paths], "sha256": [sha256(path) for path in input_paths]})
    provenance.to_csv(output / "00_INPUT_PROVENANCE.tsv", sep="\t", index=False)
    pd.DataFrame(factor_rows).to_csv(output / "01_FROZEN_NON_B_FACTOR_MAP.tsv", sep="\t", index=False)
    qc = pd.DataFrame(all_qc)
    associations = pd.DataFrame(all_associations)
    contrasts = pd.DataFrame(all_contrasts)
    qc.to_csv(output / "02_ANNULUS_COVERAGE_QC.tsv", sep="\t", index=False)
    patient_qc = qc.groupby(["chip", "patient_id", "ring_um"], as_index=False).agg(
        eligible_centers=("eligible_centers", "sum"),
        valid_centers=("valid_centers", "sum"),
        tissue_components=("tissue_component_id", "nunique"),
    )
    patient_qc["valid_fraction"] = patient_qc["valid_centers"] / patient_qc["eligible_centers"]
    patient_qc.to_csv(output / "02A_PATIENT_ANNULUS_COVERAGE_QC.tsv", sep="\t", index=False)
    associations.to_csv(output / "03_PATIENT_BLOCK_ASSOCIATIONS.tsv.gz", sep="\t", index=False, compression="gzip")
    contrasts.to_csv(output / "04_PATIENT_Q90_RADIAL_CONTRASTS.tsv.gz", sep="\t", index=False, compression="gzip")
    summary = summarize_primary(associations)
    summary.to_csv(output / "05_PRIMARY_EXTERNAL_ASSOCIATION_SUMMARY.tsv", sep="\t", index=False)
    candidates = summary.loc[summary["evidence_class"].ne("Unresolved external association")].copy()
    candidates.to_csv(output / "06_FROZEN_RELATIONSHIP_CANDIDATES.tsv", sep="\t", index=False)
    unresolved = summary.loc[summary["evidence_class"].eq("Unresolved external association")].copy()
    unresolved.to_csv(output / "07_UNRESOLVED_RELATIONSHIPS.tsv", sep="\t", index=False)
    create_figures(summary, contrasts, output)
    write_report(output, summary, qc, associations)

    validation = pd.DataFrame(
        [
            ("three_non_b_families_only", set(associations["factor_family"]) == set(FACTOR_MAP), "error"),
            ("b_family_absent", not associations["factor_family"].str.contains("^B$").any(), "error"),
            ("all_18_states", associations["state"].nunique() == 18, "error"),
            ("all_21_patients_primary", associations.loc[(associations["ring_um"] == PRIMARY_RING) & (associations["block_size_um"] == 200), "patient_id"].nunique() == 21, "error"),
            ("all_three_chips_primary", associations.loc[(associations["ring_um"] == PRIMARY_RING) & (associations["block_size_um"] == 200), "chip"].nunique() == 3, "error"),
            ("disease_not_used_in_summary_groups", "disease_provenance_only" not in summary.columns, "error"),
            ("primary_test_count_54", len(summary) == 54, "error"),
            ("finite_primary_rho", np.isfinite(summary["median_adjusted_rho"]).all(), "error"),
            ("candidate_and_unresolved_partition", len(candidates) + len(unresolved) == len(summary), "error"),
        ],
        columns=["check", "passed", "severity"],
    )
    validation.to_csv(output / "09_VALIDATION.tsv", sep="\t", index=False)
    if not validation.loc[validation["severity"].eq("error"), "passed"].all():
        raise RuntimeError("Validation failed; inspect 09_VALIDATION.tsv")
    session = {
        "created_at": datetime.now().isoformat(),
        "elapsed_seconds": time.time() - started,
        "seed": SEED,
        "python": sys.version,
        "platform": platform.platform(),
        "numpy": np.__version__,
        "pandas": pd.__version__,
        "scipy": scipy.__version__,
        "primary_ring_um": PRIMARY_RING,
        "primary_block_um": PRIMARY_BLOCK,
    }
    (output / "SESSION_INFO.json").write_text(json.dumps(session, indent=2), encoding="utf-8")
    (output / "COMPLETE.ok").write_text("complete\n", encoding="utf-8")
    print(json.dumps({"output": str(output), "elapsed_seconds": session["elapsed_seconds"], "candidates": len(candidates)}, indent=2))


if __name__ == "__main__":
    main()
