#!/usr/bin/env python3
"""Render the frozen K=7 spatial niche atlas and Chinese findings reports."""

from __future__ import annotations

import argparse
import math
from pathlib import Path

import matplotlib
import numpy as np
import pandas as pd

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib.backends.backend_pdf import PdfPages


ROOT = Path(__file__).resolve().parents[1]
FAMILIES = (
    "LEC-associated", "Blood-endothelial", "AT1-associated",
    "Perivascular-stromal", "B-enriched/mixed-immune",
    "Ciliated-associated", "Fibroblast-associated",
)
SHORT = {
    "LEC-associated": "LEC", "Blood-endothelial": "Blood endo",
    "AT1-associated": "AT1", "Perivascular-stromal": "Perivascular",
    "B-enriched/mixed-immune": "B / mixed immune",
    "Ciliated-associated": "Ciliated", "Fibroblast-associated": "Fibroblast",
}
COLORS = dict(zip(FAMILIES, sns.color_palette("tab10", 7)))


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output-dir", required=True, type=Path)
    return parser.parse_args()


def scatter(ax, x, y, values=None, *, cmap="viridis", vmax=None, color=None, size=2.0):
    if values is None:
        ax.scatter(x, y, s=size, color=color or "#666666", linewidths=0, rasterized=True)
    else:
        ax.scatter(x, y, s=size, c=values, cmap=cmap, vmin=0, vmax=vmax,
                   linewidths=0, rasterized=True)
    ax.set_aspect("equal")
    if not ax.yaxis_inverted():
        ax.invert_yaxis()
    ax.set_xticks([])
    ax.set_yticks([])


def save_page(fig, pdf: PdfPages, png: Path | None = None, rect=(0.0, 0.0, 1.0, 0.95)):
    fig.tight_layout(rect=rect, h_pad=1.5)
    pdf.savefig(fig, bbox_inches="tight")
    if png is not None:
        fig.savefig(png, dpi=240, bbox_inches="tight", facecolor="white")
    plt.close(fig)


def cosine(left: np.ndarray, right: np.ndarray) -> float:
    denominator = np.linalg.norm(left) * np.linalg.norm(right)
    return float(np.dot(left, right) / denominator) if denominator > 0 else np.nan


def update_reproducibility(output: Path, composition: pd.DataFrame) -> pd.DataFrame:
    reproducibility = pd.read_csv(output / "07_PATIENT_REPRODUCIBILITY.tsv", sep="\t")
    reproducibility = reproducibility.drop(
        columns=["median_loading_hotspot_cosine", "minimum_loading_hotspot_cosine"],
        errors="ignore",
    )
    loading = composition.loc[composition.record_type.eq("factor_loading")]
    observed = composition.loc[
        composition.record_type.eq("observed_composition") & composition.region.eq("q95_hotspot")
    ]
    states = sorted(observed.state.dropna().unique())
    rows = []
    for key, group in observed.groupby(["family_label", "chip", "patient_id"]):
        family, chip, patient = key
        obs = group.set_index("state").reindex(states).value.to_numpy(float)
        load = loading.loc[
            loading.family_label.eq(family) & loading.chip.eq(chip)
        ].set_index("state").reindex(states).value.to_numpy(float)
        rows.append((family, patient, cosine(obs, load)))
    cosines = pd.DataFrame(rows, columns=["family_label", "patient_id", "loading_hotspot_composition_cosine"])
    summary = cosines.groupby("family_label").loading_hotspot_composition_cosine.agg(["median", "min"]).reset_index()
    summary.columns = ["family_label", "median_loading_hotspot_cosine", "minimum_loading_hotspot_cosine"]
    reproducibility = reproducibility.merge(summary, on="family_label", how="left")
    reproducibility.to_csv(output / "07_PATIENT_REPRODUCIBILITY.tsv", sep="\t", index=False)
    cosines.to_csv(output / "work/PATIENT_FACTOR_LOADING_COMPOSITION_COSINES.tsv", sep="\t", index=False)
    return reproducibility


def write_reports(output: Path, reproducibility: pd.DataFrame, classification: pd.DataFrame,
                  metrics: pd.DataFrame, marker: pd.DataFrame, relationships: pd.DataFrame):
    q95 = metrics.loc[(metrics.record_type == "threshold_summary") & (metrics.threshold == "q95")]
    expected_flag = marker.expected_for_family.astype(str).str.lower().eq("true")
    expected = marker.loc[(marker.record_type == "marker_gene") & expected_flag]
    marker_summary = expected.groupby("family_label").agg(
        marker_delta=("hotspot_minus_background_log2cpm", "median"),
        marker_positive_fraction=("hotspot_minus_background_log2cpm", lambda x: float(np.mean(x > 0))),
    )
    close = relationships.loc[relationships.distance_band_um == "0-50"]
    close_summary = close.groupby(["family_a", "family_b"]).log2_enrichment.median().sort_values(ascending=False)
    top_pairs = close_summary.head(5)

    methods = [
        "# K=7 空间组织生态位表征：方法与发现",
        "",
        "## 分析边界",
        "",
        "本轮完全复用冻结的 K=7 NMF factor、relative usage、跨芯片 family map、Reference18 neighborhood composition、原始 Bin50 counts 以及既有 edgeR/GSEA。未重新运行 RCTD、NMF、edgeR 或 GSEA，也未进行新的疾病差异检验。",
        "",
        "方法框架参考 Mayr et al. (2024) 对去卷积细胞组成、NMF niche 与空间邻接的分层解释，但本项目保持自身的 Bin50、Reference18 和逐芯片 factor 定义；composition-derived factor 不被自动当作经过组织学验证的 niche。论文：https://pmc.ncbi.nlm.nih.gov/articles/PMC11313858/",
        "",
        "## 空间定义",
        "",
        "- 连续图使用每位患者内部的直接 K=7 relative usage；跨芯片色标不共享。",
        "- q90/q95/q99 使用患者内精确 top 10%/5%/1% 排名，避免大量并列值扩大热点。",
        "- 邻接为同患者 25 μm 网格的四/八邻域；不跨患者连边。组织不相连片段由全组织图自动区分。",
        "- 连通组件面积按每个 Bin50 为 625 μm²；紧凑度为 4πA/P²；边缘比例以全组织八邻域不足 8 个定义。",
        "- 0-50、50-100、100-200 μm 关系为描述性 pair enrichment；同时单列同位置 hotspot overlap，避免把混合 usage 当作物理邻近。",
        "",
        "## 核心结果",
        "",
    ]
    for row in reproducibility.itertuples(index=False):
        category = classification.loc[classification.family_label.eq(row.family_label), "classification"].iloc[0]
        methods.append(
            f"- **{row.family_label}**：q95 形成预设连通结构的患者 {int(row.structured_patient_n)}/21；"
            f"最大连通块占比中位数 {row.median_largest_component_fraction:.3f}；"
            f"热点组成与对应 loading 的 cosine 中位数 {row.median_loading_hotspot_cosine:.3f}；"
            f"预期 marker 的患者阳性比例 {row.expected_marker_positive_patient_fraction:.2f}。分类：{category}。"
        )
    methods.extend([
        "",
        "## 解释限制",
        "",
        "Reference18 composition 与 NMF loading 的一致性属于模型内部解释，不是独立验证。marker 支持来自原始计数的描述性定位，但仍可能反映细胞组成和测序深度。50 μm center-excluded neighborhood 输入本身会产生空间平滑，因此只有重复出现、具有连通形态且有原始表达支持的结构才进入候选 niche；本轮结果不支持直接命名 TLS、fibrotic niche 或疾病特异性 niche。",
    ])
    (output / "09_METHODS_AND_FINDINGS_ZH.md").write_text("\n".join(methods) + "\n", encoding="utf-8")

    next_steps = [
        "# 下一步建议",
        "",
        "1. 优先对 C 类结构逐患者叠加 ssDNA/组织学图像复核边界；当前最值得看的是 B-mixed、Ciliated、Fibroblast 与 Perivascular 在少数患者中的大连通块。",
        "2. 对候选区域进行独立的原位 marker 验证；组成因子与 Reference18 的一致性不能替代验证。",
        "3. 在组织结构冻结之前，不新增疾病显著性检验，也不根据疾病结果调整 q90/q95/q99。",
        "4. q99 仅用于极高 usage 核心；若热点少于 20 个 Bin50 或高度孤立，不将其视为完整 niche。",
        "5. 后续疾病比较必须以患者为重复并分层审计芯片异质性；这将在本任务完成后恢复上一项 joint vs chip-stratified edgeR/GSEA 工作。",
        "",
        "## 需要优先人工复核的空间关系（q95，0-50 μm，中位 log2 enrichment）",
        "",
    ]
    for (left, right), value in top_pairs.items():
        next_steps.append(f"- {left} ↔ {right}: {value:.3f}")
    next_steps.extend([
        "",
        "这些关系仅表示共定位/邻近，不表示细胞通讯或因果作用。",
        "",
        "## 与 SSc-ILD 研究的关系",
        "",
        "免疫-血管、成纤维-周血管以及上皮结构可能具有后续 SSc-ILD 研究价值，但本轮是疾病盲法表征，尚未证明任何结构具有 SSc-ILD 特异性。",
    ])
    (output / "10_NEXT_STEPS_ZH.md").write_text("\n".join(next_steps) + "\n", encoding="utf-8")


def main() -> None:
    args = parse_args()
    output = args.output_dir.resolve()
    figures = output / "figures"
    figures.mkdir(exist_ok=True)
    data = np.load(output / "work/K7_SPATIAL_FRAME.npz", allow_pickle=True)
    frame = pd.DataFrame({
        "patient_id": data["patient_id"].astype(str), "chip": data["chip"].astype(str),
        "disease": data["disease"].astype(str), "x": data["x_um"], "y": data["y_um"],
        "nUMI": data["nUMI"], "nGene": data["nGene"],
    })
    usage = data["usage"].astype(float)
    composition = pd.read_csv(output / "04_FACTOR_REFERENCE18_COMPOSITION.tsv", sep="\t")
    metrics = pd.read_csv(output / "03_FACTOR_HOTSPOT_AND_COMPONENT_METRICS.tsv", sep="\t", low_memory=False)
    marker = pd.read_csv(output / "05_FACTOR_MARKER_AND_PROGRAM_SUPPORT.tsv", sep="\t", low_memory=False)
    relationships = pd.read_csv(output / "06_FACTOR_SPATIAL_RELATIONSHIPS.tsv", sep="\t")
    classification = pd.read_csv(output / "08_CANDIDATE_NICHE_CLASSIFICATION.tsv", sep="\t")
    reproducibility = update_reproducibility(output, composition)
    write_reports(output, reproducibility, classification, metrics, marker, relationships)

    patient_meta = frame[["patient_id", "chip", "disease"]].drop_duplicates().sort_values(["chip", "disease", "patient_id"])
    chip_scales = {}
    for chip in frame.chip.unique():
        take = frame.chip.eq(chip).to_numpy()
        for index, family in enumerate(FAMILIES):
            chip_scales[(chip, family)] = float(np.quantile(usage[take, index], 0.99))

    pdf_path = output / "02_PATIENT_FACTOR_SPATIAL_ATLAS.pdf"
    with PdfPages(pdf_path) as pdf:
        fig = plt.figure(figsize=(16, 10))
        fig.text(0.07, 0.86, "K=7 Spatial Tissue Niche Characterization", fontsize=28, weight="bold")
        fig.text(0.07, 0.78, "21 patients | 7 factor families | patient-contained physical graphs", fontsize=16)
        fig.text(0.07, 0.66, "Continuous usage, q90/q95/q99 hotspot sensitivity, composition, markers, and spatial architecture", fontsize=14)
        fig.text(0.07, 0.52, "Important: composition factor is not automatically a validated biological niche", fontsize=18, color="#b22222", weight="bold")
        fig.text(0.07, 0.40, "All scales are chip/family-specific. No disease testing was performed.", fontsize=13)
        fig.text(0.07, 0.12, "Generated from frozen K=7 inputs; coordinates never cross patient boundaries.", fontsize=11)
        plt.axis("off")
        save_page(fig, pdf)

        for patient_row in patient_meta.itertuples(index=False):
            patient = patient_row.patient_id
            take = frame.patient_id.eq(patient).to_numpy()
            current = frame.loc[take]
            current_usage = usage[take]
            fig, axes = plt.subplots(2, 5, figsize=(22, 10))
            scatter(axes[0, 0], current.x, current.y, color="#555555", size=3)
            axes[0, 0].set_title("Tissue bins")
            scatter(axes[0, 1], current.x, current.y, np.log1p(current.nUMI), cmap="magma",
                    vmax=float(np.quantile(np.log1p(current.nUMI), 0.99)), size=3)
            axes[0, 1].set_title("log1p nUMI")
            scatter(axes[0, 2], current.x, current.y, np.log1p(current.nGene), cmap="magma",
                    vmax=float(np.quantile(np.log1p(current.nGene), 0.99)), size=3)
            axes[0, 2].set_title("log1p nGene")
            for family_index, family in enumerate(FAMILIES):
                ax = axes.flat[family_index + 3]
                scatter(ax, current.x, current.y, current_usage[:, family_index], cmap="viridis",
                        vmax=chip_scales[(patient_row.chip, family)], size=3)
                ax.set_title(f"{SHORT[family]}\n{patient_row.chip}-specific scale")
            fig.suptitle(f"{patient} | {patient_row.disease} | {patient_row.chip}", fontsize=18)
            safe = patient.replace("/", "_")
            save_page(fig, pdf, figures / f"PATIENT_{safe}_CONTINUOUS_FACTOR_ATLAS.png")

            fig, axes = plt.subplots(7, 3, figsize=(15, 27))
            for family_index, family in enumerate(FAMILIES):
                values = current_usage[:, family_index]
                for threshold_index, quantile in enumerate((0.90, 0.95, 0.99)):
                    count = max(1, int(math.ceil((1 - quantile) * len(values))))
                    order = np.lexsort((np.arange(len(values)), -values))
                    high = np.zeros(len(values), bool)
                    high[order[:count]] = True
                    ax = axes[family_index, threshold_index]
                    scatter(ax, current.x, current.y, color="#dddddd", size=2)
                    scatter(ax, current.x.to_numpy()[high], current.y.to_numpy()[high],
                            color=COLORS[family], size=5)
                    ax.set_title(f"{SHORT[family]} q{int(quantile*100)} | n={count}")
            fig.suptitle(f"Patient-specific hotspot sensitivity | {patient} | {patient_row.disease} | {patient_row.chip}", fontsize=18)
            save_page(fig, pdf, figures / f"PATIENT_{safe}_HOTSPOT_Q90_Q95_Q99.png")

        loading = composition.loc[composition.record_type.eq("factor_loading")]
        states = list(data["state_names"].astype(str))
        load_matrix = loading.pivot_table(index="family_label", columns="state", values="value", aggfunc="mean").reindex(FAMILIES)[states]
        observed = composition.loc[(composition.record_type == "observed_composition") & (composition.region == "q95_hotspot")]
        obs_matrix = observed.pivot_table(index="family_label", columns="state", values="value", aggfunc="mean").reindex(FAMILIES)[states]
        fig, axes = plt.subplots(2, 1, figsize=(20, 10))
        sns.heatmap(load_matrix, cmap="mako", ax=axes[0], cbar_kws={"label": "Mean normalized loading"})
        axes[0].set_title("Frozen K=7 factor loadings (mean across chip-specific factors)")
        sns.heatmap(obs_matrix, cmap="mako", ax=axes[1], cbar_kws={"label": "Mean Reference18 weight"})
        axes[1].set_title("Observed q95 hotspot neighborhood composition (mean across patients)")
        for ax in axes:
            ax.set_xlabel("")
            ax.set_ylabel("")
        save_page(fig, pdf, figures / "SUMMARY_FACTOR_REFERENCE18_COMPOSITION.png")

        expected_flag = marker.expected_for_family.astype(str).str.lower().eq("true")
        expected = marker.loc[(marker.record_type == "marker_gene") & expected_flag]
        marker_matrix = expected.pivot_table(index="family_label", columns="gene_symbol",
                                             values="hotspot_minus_background_log2cpm", aggfunc="median").reindex(FAMILIES)
        fig, ax = plt.subplots(figsize=(18, 7))
        sns.heatmap(marker_matrix, cmap="vlag", center=0, ax=ax,
                    cbar_kws={"label": "Median hotspot - background log2(CPM+1)"})
        ax.set_title("Raw marker expression support across patients")
        ax.set_xlabel("Marker gene")
        ax.set_ylabel("Factor family")
        save_page(fig, pdf, figures / "SUMMARY_FACTOR_MARKER_SUPPORT.png")

        close = relationships.loc[relationships.distance_band_um.eq("0-50")]
        rel = close.groupby(["family_a", "family_b"]).log2_enrichment.median()
        matrix = pd.DataFrame(0.0, index=FAMILIES, columns=FAMILIES)
        for (left, right), value in rel.items():
            matrix.loc[left, right] = value
            matrix.loc[right, left] = value
        for index in range(len(matrix)):
            matrix.iat[index, index] = np.nan
        fig, ax = plt.subplots(figsize=(11, 9))
        sns.heatmap(matrix, cmap="vlag", center=0, annot=True, fmt=".2f", ax=ax,
                    cbar_kws={"label": "Median log2 enrichment at 0-50 um"})
        ax.set_title("q95 hotspot physical proximity (descriptive)")
        save_page(fig, pdf, figures / "SUMMARY_FACTOR_SPATIAL_RELATIONSHIPS.png")

        fig, axes = plt.subplots(1, 2, figsize=(16, 6))
        axes[0].barh([SHORT[x] for x in reproducibility.family_label], reproducibility.structured_patient_n,
                     color=[COLORS[x] for x in reproducibility.family_label])
        axes[0].set_xlim(0, 21)
        axes[0].set_xlabel("Patients with a connected q95 structure")
        axes[0].set_title("Spatial recurrence")
        axes[1].barh([SHORT[x] for x in reproducibility.family_label],
                     reproducibility.median_loading_hotspot_cosine,
                     color=[COLORS[x] for x in reproducibility.family_label])
        axes[1].set_xlim(0, 1)
        axes[1].set_xlabel("Median loading-hotspot composition cosine")
        axes[1].set_title("Model-internal composition agreement")
        save_page(fig, pdf, figures / "SUMMARY_FACTOR_REPRODUCIBILITY.png")

        q95 = metrics.loc[(metrics.record_type == "threshold_summary") & (metrics.threshold == "q95")]

        patient_scores = (
            q95.groupby("patient_id", as_index=False)
            .agg(mean_largest_component_fraction=("largest_component_fraction", "mean"))
            .sort_values("mean_largest_component_fraction", ascending=False)
            .head(6)
        )
        fig, axes = plt.subplots(2, 3, figsize=(18, 11))
        for ax, score_row in zip(axes.flat, patient_scores.itertuples(index=False)):
            take = frame.patient_id.eq(score_row.patient_id).to_numpy()
            current = frame.loc[take]
            current_usage = usage[take]
            dominant = np.argmax(current_usage, axis=1)
            point_colors = np.asarray([COLORS[FAMILIES[index]] for index in dominant])
            ax.scatter(
                current.x,
                current.y,
                s=3,
                c=point_colors,
                linewidths=0,
                rasterized=True,
            )
            ax.set_aspect("equal")
            ax.invert_yaxis()
            ax.set_xticks([])
            ax.set_yticks([])
            meta = patient_meta.loc[patient_meta.patient_id.eq(score_row.patient_id)].iloc[0]
            ax.set_title(
                f"{score_row.patient_id} | {meta.disease} | {meta.chip}\n"
                f"mean max-component fraction={score_row.mean_largest_component_fraction:.2f}"
            )
        handles = [
            plt.Line2D([0], [0], marker="o", linestyle="", color=COLORS[family],
                       label=SHORT[family], markersize=7)
            for family in FAMILIES
        ]
        fig.legend(handles=handles, loc="lower center", ncol=4, frameon=False)
        fig.suptitle("Representative patient tissue atlas: dominant K=7 factor per bin", fontsize=17)
        save_page(
            fig,
            pdf,
            figures / "SUMMARY_REPRESENTATIVE_PATIENT_TISSUE_ATLAS.png",
            rect=(0.0, 0.08, 1.0, 0.94),
        )

        fig, axes = plt.subplots(2, 4, figsize=(20, 10))
        for family_index, family in enumerate(FAMILIES):
            best = q95.loc[q95.family_label.eq(family)].sort_values(
                ["largest_component_bin_n", "largest_component_fraction"], ascending=False
            ).iloc[0]
            take = frame.patient_id.eq(best.patient_id).to_numpy()
            current = frame.loc[take]
            values = usage[take, family_index]
            count = max(1, int(math.ceil(0.05 * len(values))))
            order = np.lexsort((np.arange(len(values)), -values))
            high = np.zeros(len(values), bool)
            high[order[:count]] = True
            ax = axes.flat[family_index]
            scatter(ax, current.x, current.y, color="#dddddd", size=2)
            scatter(ax, current.x.to_numpy()[high], current.y.to_numpy()[high], color=COLORS[family], size=7)
            ax.set_title(f"{SHORT[family]} | {best.patient_id}\nmax component={int(best.largest_component_bin_n)} bins")
        axes.flat[-1].axis("off")
        fig.suptitle("Best-connected q95 example per factor family (not disease-selected)", fontsize=17)
        save_page(fig, pdf, figures / "SUMMARY_REPRESENTATIVE_CONNECTED_HOTSPOTS.png")

    validation = pd.read_csv(output / "99_VALIDATION.tsv", sep="\t")
    extra = pd.DataFrame([
        ("atlas_pdf_exists", True, pdf_path.exists()),
        ("atlas_pdf_nonempty", True, pdf_path.stat().st_size > 100_000),
        ("patient_continuous_png_n", 21, len(list(figures.glob("PATIENT_*_CONTINUOUS_FACTOR_ATLAS.png")))),
        ("patient_hotspot_png_n", 21, len(list(figures.glob("PATIENT_*_HOTSPOT_Q90_Q95_Q99.png")))),
        ("summary_png_n", 6, len(list(figures.glob("SUMMARY_*.png")))),
        ("methods_report_exists", True, (output / "09_METHODS_AND_FINDINGS_ZH.md").exists()),
        ("next_steps_exists", True, (output / "10_NEXT_STEPS_ZH.md").exists()),
    ], columns=["check", "expected", "observed"])
    extra["passed"] = extra.expected.astype(str).eq(extra.observed.astype(str))
    validation = validation.loc[~validation.check.isin(extra.check)]
    validation = pd.concat([validation, extra], ignore_index=True)
    validation.to_csv(output / "99_VALIDATION.tsv", sep="\t", index=False)
    if not validation.passed.all():
        raise RuntimeError("Render validation failed")
    (output / "COMPLETE.ok").write_text("COMPLETE\n", encoding="utf-8")
    print(pdf_path)


if __name__ == "__main__":
    main()
