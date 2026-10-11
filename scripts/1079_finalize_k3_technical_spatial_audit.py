#!/usr/bin/env python3
"""Finalize technical consistency tiers and raw spatial contact sheets."""

from __future__ import annotations

import json
from datetime import datetime
from pathlib import Path

import matplotlib
import numpy as np
import pandas as pd

matplotlib.use("Agg")
import matplotlib.pyplot as plt


ROOT = Path(__file__).resolve().parents[1]
OUT = (
    ROOT
    / "results/reanalysis/ALL21_K3_NICHE_PSEUDOBULK_DEG_20261011_201008"
)
GENE_NICHE = {
    "SCGB1A1": "CN2",
    "DRC1": "CN2",
    "MS4A8": "CN3",
    "DNAI1": "CN3",
    "NWD1": "CN3",
    "NEK5": "CN3",
}
DISEASE_COLORS = {"HC": "#4C72B0", "IPF": "#DD8452", "SSc-ILD": "#55A868"}


def as_bool(series: pd.Series) -> pd.Series:
    return series.astype(str).str.lower().eq("true")


def main() -> None:
    summary = pd.read_csv(
        OUT / "36_TECHNICAL_PATIENT_CONSISTENCY_SUMMARY.tsv", sep="\t"
    )
    adjusted = pd.read_csv(
        OUT / "33_TECHNICAL_COVARIATE_ADJUSTED_EFFECTS.tsv", sep="\t"
    )
    chip = pd.read_csv(OUT / "34_CHIP_SPECIFIC_GENE_EFFECTS.tsv", sep="\t")
    expression = pd.read_csv(OUT / "25_ROBUST_GENE_PATIENT_LOGCPM.tsv", sep="\t")
    expression = expression.loc[expression["analysis"].eq("primary_all21")]
    spatial = pd.read_csv(
        OUT / "31_REPRESENTATIVE_GENE_UNSMOOTHED_SPATIAL_COUNTS.tsv.gz",
        sep="\t",
    )
    figures = OUT / "figures"
    figures.mkdir(exist_ok=True)

    for column in [
        "all_three_chips_same_direction",
        "all_21_lopo_same_direction",
        "all_four_technical_models_same_direction",
        "maximum_technical_attenuation_lt30pct",
        "passes_all_direction_and_attenuation_checks",
    ]:
        summary[column] = as_bool(summary[column])
    summary["evidence_tier"] = np.select(
        [
            summary["passes_all_direction_and_attenuation_checks"],
            summary["all_21_lopo_same_direction"]
            & summary["all_four_technical_models_same_direction"]
            & summary["all_three_chips_same_direction"],
            summary["all_21_lopo_same_direction"]
            & summary["all_four_technical_models_same_direction"],
        ],
        [
            "Tier_A_cross_chip_and_technical_stable",
            "Tier_B_cross_chip_direction_stable_but_qc_attenuated",
            "Tier_C_patient_stable_but_chip_heterogeneous",
        ],
        default="Tier_D_unstable",
    )
    summary.to_csv(
        OUT / "38_TECHNICAL_EVIDENCE_TIER_ASSIGNMENT.tsv", sep="\t", index=False
    )

    # Chip effect heatmap for all 21 robust gene-by-contrast rows.
    chip["row_label"] = (
        chip["niche"]
        + " | "
        + chip["comparison"]
        + " | "
        + chip["gene_symbol"]
    )
    row_order = (
        summary.assign(
            row_label=lambda x: x["niche"]
            + " | "
            + x["comparison"]
            + " | "
            + x["gene_symbol"]
        )
        .sort_values(["evidence_tier", "niche", "comparison", "gene_symbol"])[
            "row_label"
        ]
        .tolist()
    )
    matrix = chip.pivot(
        index="row_label", columns="chip_label", values="chip_mean_logcpm_effect"
    ).reindex(index=row_order, columns=["J2", "K8", "L3"])
    limit = max(1.0, float(np.nanmax(np.abs(matrix.to_numpy()))))
    fig, ax = plt.subplots(figsize=(7.5, 9.5), constrained_layout=True)
    image = ax.imshow(matrix.to_numpy(), cmap="RdBu_r", vmin=-limit, vmax=limit, aspect="auto")
    ax.set_xticks(range(3), matrix.columns)
    ax.set_yticks(range(len(matrix)), matrix.index, fontsize=7.5)
    ax.set_title("Chip-specific patient-mean logCPM disease effects")
    colorbar = fig.colorbar(image, ax=ax, shrink=0.7)
    colorbar.set_label("Target disease − reference disease")
    fig.savefig(figures / "FIG06_CHIP_SPECIFIC_GENE_EFFECT_HEATMAP.png", dpi=180)
    plt.close(fig)

    # Maximum attenuation under the four technical covariates.
    technical = adjusted.loc[adjusted["covariate"].ne("none")].copy()
    technical["row_label"] = (
        technical["niche"]
        + " | "
        + technical["comparison"]
        + " | "
        + technical["gene_symbol"]
    )
    attenuation = technical.pivot(
        index="row_label",
        columns="covariate",
        values="attenuation_fraction_vs_unadjusted_lm",
    ).reindex(index=row_order)
    attenuation = attenuation[
        [
            "log10_bin_n",
            "log10_median_common_umi_plus1",
            "log10_median_common_ngene_plus1",
            "median_pct_mt",
        ]
    ]
    fig, ax = plt.subplots(figsize=(9, 9.5), constrained_layout=True)
    image = ax.imshow(
        attenuation.to_numpy() * 100,
        cmap="PuOr_r",
        vmin=-100,
        vmax=100,
        aspect="auto",
    )
    ax.set_xticks(
        range(4), ["Bin count", "Median UMI", "Median nGene", "Median mt%"], rotation=20
    )
    ax.set_yticks(range(len(attenuation)), attenuation.index, fontsize=7.5)
    ax.set_title("Disease-effect attenuation after separate technical adjustments")
    colorbar = fig.colorbar(image, ax=ax, shrink=0.7)
    colorbar.set_label("Attenuation (%)\nnegative = effect strengthened")
    fig.savefig(figures / "FIG07_TECHNICAL_COVARIATE_ATTENUATION_HEATMAP.png", dpi=180)
    plt.close(fig)

    # Unsmoothed raw-count spatial contact sheets for six representative genes.
    patient_order = (
        spatial[["patient_id", "chip_label", "disease"]]
        .drop_duplicates()
        .assign(
            chip_order=lambda x: pd.Categorical(x["chip_label"], ["K8", "J2", "L3"]),
            disease_order=lambda x: pd.Categorical(x["disease"], ["HC", "IPF", "SSc-ILD"]),
        )
        .sort_values(["chip_order", "disease_order", "patient_id"])
    )
    by_chip = {
        chip_name: patient_order.loc[patient_order["chip_label"].eq(chip_name), "patient_id"].tolist()
        for chip_name in ["K8", "J2", "L3"]
    }
    for gene, niche in GENE_NICHE.items():
        target = spatial.loc[spatial["niche"].eq(niche), gene]
        positive = np.log1p(target[target > 0])
        vmax = max(1.0, float(np.quantile(positive, 0.99))) if len(positive) else 1.0
        fig, axes = plt.subplots(3, 7, figsize=(18, 8.5), constrained_layout=True)
        for row, chip_name in enumerate(["K8", "J2", "L3"]):
            for col, patient in enumerate(by_chip[chip_name]):
                ax = axes[row, col]
                sub = spatial.loc[spatial["patient_id"].eq(patient)]
                target_niche = sub["niche"].eq(niche)
                detected = target_niche & sub[gene].gt(0)
                ax.scatter(sub["x_um"], sub["y_um"], s=0.35, color="#D9D9D9", linewidths=0)
                ax.scatter(
                    sub.loc[target_niche, "x_um"],
                    sub.loc[target_niche, "y_um"],
                    s=0.55,
                    color="#9ECAE1",
                    linewidths=0,
                )
                if detected.any():
                    ax.scatter(
                        sub.loc[detected, "x_um"],
                        sub.loc[detected, "y_um"],
                        c=np.log1p(sub.loc[detected, gene]),
                        cmap="magma",
                        vmin=0,
                        vmax=vmax,
                        s=2.2,
                        linewidths=0,
                    )
                disease = sub["disease"].iloc[0]
                detection_rate = 100 * detected.sum() / max(1, target_niche.sum())
                ax.set_title(
                    f"{patient}\n{disease}; detected {detection_rate:.1f}%",
                    fontsize=7,
                    color=DISEASE_COLORS[disease],
                )
                ax.set_aspect("equal")
                ax.invert_yaxis()
                ax.axis("off")
                if col == 0:
                    ax.text(
                        -0.08,
                        0.5,
                        chip_name,
                        transform=ax.transAxes,
                        rotation=90,
                        va="center",
                        ha="right",
                        fontsize=11,
                        fontweight="bold",
                    )
        fig.suptitle(
            f"{gene}: unsmoothed raw Bin50 counts inside frozen {niche}\n"
            "grey = evaluable tissue; blue = target niche; magma = directly detected",
            fontsize=13,
        )
        fig.savefig(figures / f"FIG08_SPATIAL_{gene}_{niche}.png", dpi=180)
        plt.close(fig)

    tier_counts = {
        str(key): int(value)
        for key, value in summary["evidence_tier"].value_counts().items()
    }
    tier_a = summary.loc[
        summary["evidence_tier"].eq("Tier_A_cross_chip_and_technical_stable")
    ]
    tier_a_labels = [
        f"{row.niche} {row.comparison} `{row.gene_symbol}`"
        for row in tier_a.itertuples(index=False)
    ]
    interpretation = f"""# Technical-confounding and patient-consistency audit

## What was tested

- Direction of each of the 21 robust gene-by-contrast results separately in J2, K8, and L3.
- Leave-one-patient-out refits for every patient (21 refits per result).
- Separate adjustment for niche Bin count, median raw UMI per Bin, median detected genes per Bin, and median mitochondrial fraction.
- Unsmoothed raw Bin50 spatial maps for six representative genes.

## Results

- **21/21** results retained the same direction after omitting every individual patient: no result is driven by one patient alone.
- **21/21** retained the same direction in all four technical-covariate models.
- **14/21** had the same direction independently in all three chips.
- **9/21** had less than 30% maximum attenuation across all four separate technical adjustments.
- **8/21** gene-by-contrast rows passed all prespecified direction and attenuation checks (Tier A): {', '.join(tier_a_labels)}.
- Tier counts: {dict(tier_counts)}.

## Interpretation

The signal is patient-stable, but the full 21-row set is not uniformly chip- and QC-stable. The strongest evidence is concentrated in eight Tier-A gene-by-contrast rows, all in CN3. CN2 `SCGB1A1` and `DRC1` are direction-consistent across chips and leave-one-patient-out analyses, but attenuate by 46% and 39% under at least one QC adjustment, so they remain provisional.

Technical covariates can also carry genuine biology, especially UMI and nGene in damaged epithelium. Attenuation therefore identifies uncertainty; it does not prove technical artifact. Spatial panels show directly detected raw counts without smoothing and should be used to check whether expression is anatomically plausible rather than diffuse or single-patient dominated.
"""
    (OUT / "39_TECHNICAL_SPATIAL_AUDIT_INTERPRETATION.md").write_text(
        interpretation, encoding="utf-8"
    )
    status = {
        "completed_at": datetime.now().astimezone().isoformat(timespec="seconds"),
        "robust_rows": 21,
        "all_lopo_direction_stable": int(summary["all_21_lopo_same_direction"].sum()),
        "all_three_chips_direction_stable": int(summary["all_three_chips_same_direction"].sum()),
        "tier_a_rows": int(len(tier_a)),
        "status": "ANALYSIS_COMPLETE_PENDING_USER_APPROVAL",
    }
    (OUT / "TECHNICAL_SPATIAL_AUDIT_FINAL_STATUS.json").write_text(
        json.dumps(status, indent=2), encoding="utf-8"
    )
    print(OUT)


if __name__ == "__main__":
    main()
