#!/usr/bin/env python3
"""Finalize and audit step-7 frozen-K3 patient transcriptomics.

Consumes only the completed pseudobulk/edgeR tables. It does not refit models.
The purpose is to distinguish primary discoveries from findings that remain
supported after excluding patient/niche libraries containing fewer than 20 bins.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from datetime import datetime
from pathlib import Path

import matplotlib
import numpy as np
import pandas as pd

matplotlib.use("Agg")
import matplotlib.pyplot as plt


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_OUT = (
    ROOT
    / "results/reanalysis/ALL21_K3_NICHE_PSEUDOBULK_DEG_20261011_201008"
)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, default=DEFAULT_OUT)
    args = parser.parse_args()
    output = args.output.resolve()

    deg = pd.read_csv(output / "10_NICHE_SPECIFIC_EDGER_DEG.tsv", sep="\t")
    programs = pd.read_csv(output / "09_PROGRAM_DISEASE_TESTS.tsv", sep="\t")
    primary = deg.loc[deg["analysis"].eq("primary_all21")].copy()
    sensitivity = deg.loc[
        deg["analysis"].eq("sensitivity_exclude_lt20_bins")
    ].copy()
    keys = ["niche", "comparison", "ensembl_id"]
    merged = primary.merge(
        sensitivity,
        on=keys,
        suffixes=("_primary", "_sensitivity"),
        validate="one_to_one",
    )
    merged["primary_supported"] = (
        (merged["within_contrast_fdr_primary"] < 0.05)
        & (merged["log2_fc_primary"].abs() >= 0.5)
    )
    merged["same_direction"] = (
        merged["log2_fc_primary"] * merged["log2_fc_sensitivity"] > 0
    )
    merged["sensitivity_supported"] = (
        (merged["within_contrast_fdr_sensitivity"] < 0.05)
        & (merged["log2_fc_sensitivity"].abs() >= 0.5)
    )
    merged["robust_to_lt20_bin_exclusion"] = (
        merged["primary_supported"]
        & merged["same_direction"]
        & merged["sensitivity_supported"]
    )
    audit = merged.loc[merged["primary_supported"]].copy()
    audit = audit.sort_values(
        ["robust_to_lt20_bin_exclusion", "within_contrast_fdr_primary"],
        ascending=[False, True],
    )
    audit.to_csv(
        output / "16_PRIMARY_DEG_SENSITIVITY_AUDIT.tsv", sep="\t", index=False
    )

    summary = (
        audit.groupby(["niche", "comparison"], observed=True)
        .agg(
            primary_supported_genes=("ensembl_id", "size"),
            same_direction_genes=("same_direction", "sum"),
            sensitivity_supported_genes=("sensitivity_supported", "sum"),
            robust_genes=("robust_to_lt20_bin_exclusion", "sum"),
        )
        .reset_index()
    )
    all_combinations = pd.MultiIndex.from_product(
        [
            ["CN1", "CN2", "CN3"],
            ["IPF_vs_HC", "SSc-ILD_vs_HC", "SSc-ILD_vs_IPF"],
        ],
        names=["niche", "comparison"],
    ).to_frame(index=False)
    summary = all_combinations.merge(summary, how="left").fillna(0)
    for column in summary.columns[2:]:
        summary[column] = summary[column].astype(int)
    summary.to_csv(output / "17_DEG_ROBUSTNESS_SUMMARY.tsv", sep="\t", index=False)

    fig, axes = plt.subplots(3, 3, figsize=(12, 11), constrained_layout=True)
    niches = ["CN1", "CN2", "CN3"]
    comparisons = ["IPF_vs_HC", "SSc-ILD_vs_HC", "SSc-ILD_vs_IPF"]
    for row, niche in enumerate(niches):
        for col, comparison in enumerate(comparisons):
            ax = axes[row, col]
            sub = audit.loc[
                audit["niche"].eq(niche) & audit["comparison"].eq(comparison)
            ]
            if sub.empty:
                ax.text(0.5, 0.5, "No primary FDR finding", ha="center", va="center")
                ax.set_xlim(-2, 2)
                ax.set_ylim(-2, 2)
            else:
                colors = np.where(
                    sub["robust_to_lt20_bin_exclusion"], "#B2182B", "#9E9E9E"
                )
                ax.scatter(
                    sub["log2_fc_primary"],
                    sub["log2_fc_sensitivity"],
                    c=colors,
                    s=25,
                    alpha=0.8,
                    edgecolors="none",
                )
                low = min(sub["log2_fc_primary"].min(), sub["log2_fc_sensitivity"].min())
                high = max(sub["log2_fc_primary"].max(), sub["log2_fc_sensitivity"].max())
                pad = max(0.15, (high - low) * 0.1)
                ax.plot([low - pad, high + pad], [low - pad, high + pad], "--", color="black", lw=0.7)
                ax.set_xlim(low - pad, high + pad)
                ax.set_ylim(low - pad, high + pad)
                ax.text(
                    0.03,
                    0.97,
                    f"robust {int(sub['robust_to_lt20_bin_exclusion'].sum())}/{len(sub)}",
                    transform=ax.transAxes,
                    ha="left",
                    va="top",
                    fontsize=9,
                )
            ax.axhline(0, color="grey", lw=0.5)
            ax.axvline(0, color="grey", lw=0.5)
            ax.set_title(f"{niche} | {comparison}", fontsize=10)
            if row == 2:
                ax.set_xlabel("Primary log2FC")
            if col == 0:
                ax.set_ylabel("Sensitivity log2FC")
    fig.suptitle(
        "Primary DEG effect concordance after excluding <20-bin libraries",
        fontsize=14,
    )
    fig.savefig(
        output / "figures/FIG03_PRIMARY_VS_LOW_BIN_EXCLUSION_DEG.png", dpi=180
    )
    plt.close(fig)

    primary_program = programs.loc[programs["analysis"].eq("primary_all21")]
    program_order = [
        "ECM organization",
        "Collagen formation",
        "Epithelial differentiation",
        "Angiogenesis",
        "Inflammatory response",
        "Interferon-gamma response",
        "T-cell activation",
        "B-cell activation",
    ]
    effect = primary_program.assign(
        row_label=lambda x: x["niche"] + " | " + x["program_label"]
    ).pivot(index="row_label", columns="comparison", values="effect_sd_units")
    ordered_rows = [f"{n} | {p}" for n in niches for p in program_order]
    effect = effect.reindex(index=ordered_rows, columns=comparisons)
    limit = max(0.6, float(np.nanmax(np.abs(effect.to_numpy()))))
    fig, ax = plt.subplots(figsize=(8.5, 10), constrained_layout=True)
    image = ax.imshow(effect.to_numpy(), cmap="RdBu_r", vmin=-limit, vmax=limit, aspect="auto")
    ax.set_xticks(range(len(comparisons)), comparisons, rotation=20, ha="right")
    ax.set_yticks(range(len(ordered_rows)), ordered_rows, fontsize=8)
    ax.set_title("K3 niche-specific fixed-program disease effects\npatient-level; chip-adjusted; SD units")
    colorbar = fig.colorbar(image, ax=ax, shrink=0.65)
    colorbar.set_label("Disease effect (SD units)")
    for boundary in (7.5, 15.5):
        ax.axhline(boundary, color="black", linewidth=1.2)
    fig.savefig(output / "figures/FIG01_PROGRAM_DISEASE_EFFECT_HEATMAP.png", dpi=180)
    plt.close(fig)

    fig, axes = plt.subplots(3, 3, figsize=(13, 11), constrained_layout=True)
    for row, niche in enumerate(niches):
        for col, comparison in enumerate(comparisons):
            ax = axes[row, col]
            sub = primary.loc[
                primary["niche"].eq(niche) & primary["comparison"].eq(comparison)
            ].copy()
            sub["minus_log10_fdr"] = -np.log10(
                np.maximum(sub["within_contrast_fdr"], np.finfo(float).tiny)
            )
            significant = (
                sub["within_contrast_fdr"].lt(0.05)
                & sub["log2_fc"].abs().ge(0.5)
            )
            ax.scatter(
                sub.loc[~significant, "log2_fc"],
                sub.loc[~significant, "minus_log10_fdr"],
                s=3,
                color="#808080",
                alpha=0.28,
                rasterized=True,
            )
            ax.scatter(
                sub.loc[significant, "log2_fc"],
                sub.loc[significant, "minus_log10_fdr"],
                s=10,
                color="#B2182B",
                alpha=0.9,
            )
            for _, hit in sub.loc[significant].nsmallest(3, "within_contrast_fdr").iterrows():
                label = hit["gene_symbol"] if pd.notna(hit["gene_symbol"]) else hit["ensembl_id"]
                ax.annotate(
                    str(label),
                    (hit["log2_fc"], hit["minus_log10_fdr"]),
                    xytext=(3, 3),
                    textcoords="offset points",
                    fontsize=7,
                )
            ax.axhline(-np.log10(0.05), color="grey", linestyle="--", linewidth=0.7)
            ax.axvline(-0.5, color="grey", linestyle="--", linewidth=0.7)
            ax.axvline(0.5, color="grey", linestyle="--", linewidth=0.7)
            ax.set_ylim(bottom=0, top=max(1.5, float(sub["minus_log10_fdr"].max()) * 1.1))
            ax.set_title(f"{niche} | {comparison}", fontsize=10)
            if row == 2:
                ax.set_xlabel("log2 fold change")
            if col == 0:
                ax.set_ylabel("-log10 within-contrast BH-FDR")
            ax.text(
                0.03,
                0.97,
                f"FDR+effect: {int(significant.sum())}",
                transform=ax.transAxes,
                va="top",
                fontsize=8,
            )
    fig.suptitle("Frozen K3 niche-specific disease DEG", fontsize=14)
    fig.savefig(output / "figures/FIG02_NICHE_DISEASE_DEG_VOLCANO.png", dpi=180)
    plt.close(fig)

    robust = audit.loc[audit["robust_to_lt20_bin_exclusion"]].copy()
    robust_global_both = int(
        (
            robust["global_fdr_primary"].lt(0.05)
            & robust["global_fdr_sensitivity"].lt(0.05)
        ).sum()
    )
    summary_text = f"""# Step 7 — frozen K3 niche-specific patient transcriptomics

## Scope

- Frozen disease-blind CN1/CN2/CN3 assignments; no RCTD, neighborhood, or K rerun.
- 126,199 evaluable Bin50 raw integer-count profiles were summed into 63 patient-by-niche pseudobulks.
- Patient is the biological replicate. Each niche was modeled separately with `~ chip + disease`.
- Primary includes all 21 patients. Sensitivity excludes the two patient/niche libraries with fewer than 20 bins.
- Eight fixed programs come from versioned local Reactome 2022, GO BP 2025, and Hallmark 2020 GMT files.

## Main results

- Fixed programs: **{int((primary_program['global_fdr'] < 0.05).sum())} / {len(primary_program)}** niche-by-contrast tests passed global BH-FDR < 0.05.
- Gene-level Primary: **{int(audit.shape[0])}** genes passed within-contrast BH-FDR < 0.05 and |log2FC| >= 0.5.
- Low-bin sensitivity: **{int(robust.shape[0])} / {int(audit.shape[0])}** Primary genes retained the same direction, BH-FDR < 0.05, and |log2FC| >= 0.5 after excluding the two <20-bin libraries.
- Of those robust genes, **{robust_global_both}** passed the more conservative pooled global FDR < 0.05 in both analyses.
- Robust findings occurred in CN2 SSc-ILD vs HC ({int(((robust['niche'] == 'CN2') & (robust['comparison'] == 'SSc-ILD_vs_HC')).sum())}), CN3 IPF vs HC ({int(((robust['niche'] == 'CN3') & (robust['comparison'] == 'IPF_vs_HC')).sum())}), and CN3 SSc-ILD vs HC ({int(((robust['niche'] == 'CN3') & (robust['comparison'] == 'SSc-ILD_vs_HC')).sum())}).
- Most robust CN3 genes are cilia/airway-associated and are lower in disease than HC. Because CN3 itself is ciliated/immune/basal enriched, this is most safely interpreted as a reproducible epithelial-composition/state signal inside a mixed niche, not as a new pure cell type or a fibrosis pathway activation claim.

## Files to read first

- `09_PROGRAM_DISEASE_TESTS.tsv`: fixed-program effects and multiplicity correction.
- `16_PRIMARY_DEG_SENSITIVITY_AUDIT.tsv`: all Primary FDR findings and low-bin robustness.
- `17_DEG_ROBUSTNESS_SUMMARY.tsv`: counts by niche and disease contrast.
- `figures/FIG01_PROGRAM_DISEASE_EFFECT_HEATMAP.png`: program-level effect overview.
- `figures/FIG03_PRIMARY_VS_LOW_BIN_EXCLUSION_DEG.png`: robustness of gene-level effects.

## Interpretation boundary

A niche-level transcriptomic difference is an association in mixed Bin50 tissue. It is not evidence of cell-cell interaction, pathway activation, causality, or a pure cell state. The fixed program analysis is the confirmatory layer; individual DEG findings are a discovery layer and require biological annotation and, ideally, orthogonal validation.
"""
    (output / "15_SUMMARY.md").write_text(summary_text, encoding="utf-8")

    manifest_rows = []
    for path in sorted(output.glob("*.tsv")) + sorted((output / "figures").glob("*.png")):
        manifest_rows.append(
            {
                "path": str(path.relative_to(output)),
                "bytes": path.stat().st_size,
                "sha256": sha256(path),
            }
        )
    pd.DataFrame(manifest_rows).to_csv(
        output / "18_OUTPUT_MANIFEST.tsv", sep="\t", index=False
    )
    run = {
        "completed_at": datetime.now().astimezone().isoformat(timespec="seconds"),
        "primary_supported_genes": int(audit.shape[0]),
        "robust_to_low_bin_exclusion": int(robust.shape[0]),
        "primary_global_program_findings": int(
            (primary_program["global_fdr"] < 0.05).sum()
        ),
        "status": "ANALYSIS_COMPLETE_PENDING_USER_APPROVAL",
    }
    (output / "FINAL_STATUS.json").write_text(
        json.dumps(run, indent=2), encoding="utf-8"
    )
    print(output)


if __name__ == "__main__":
    main()
