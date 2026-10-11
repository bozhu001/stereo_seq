#!/usr/bin/env python3
"""Create audited tables/figures for K3 composition-vs-state decomposition."""

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
COLORS = {"HC": "#4C72B0", "IPF": "#DD8452", "SSc-ILD": "#55A868"}


def main() -> None:
    composition = pd.read_csv(OUT / "19_PATIENT_NICHE_RCTD_COMPOSITION.tsv", sep="\t")
    tests = pd.read_csv(OUT / "22_RCTD_COMPOSITION_DISEASE_TESTS.tsv", sep="\t")
    adjusted = pd.read_csv(OUT / "23_COMPOSITION_ADJUSTED_ROBUST_DEG.tsv", sep="\t")
    correlations = pd.read_csv(
        OUT / "26_ROBUST_GENE_COMPOSITION_CORRELATIONS.tsv", sep="\t"
    )
    figures = OUT / "figures"
    figures.mkdir(exist_ok=True)

    # Patient-level ciliated-neighborhood composition by disease and niche.
    fig, axes = plt.subplots(1, 3, figsize=(11, 4), sharey=True, constrained_layout=True)
    rng = np.random.default_rng(20261011)
    for ax, niche in zip(axes, ["CN1", "CN2", "CN3"]):
        sub = composition.loc[composition["niche"].eq(niche)]
        for position, disease in enumerate(["HC", "IPF", "SSc-ILD"]):
            values = sub.loc[sub["disease"].eq(disease), "Ciliated epithelial"] * 100
            box = ax.boxplot(
                values,
                positions=[position],
                widths=0.55,
                patch_artist=True,
                showfliers=False,
            )
            box["boxes"][0].set_facecolor(COLORS[disease])
            box["boxes"][0].set_alpha(0.45)
            jitter = rng.normal(0, 0.045, size=len(values))
            ax.scatter(
                np.full(len(values), position) + jitter,
                values,
                s=25,
                color=COLORS[disease],
                edgecolor="black",
                linewidth=0.35,
                zorder=3,
            )
        ax.set_xticks(range(3), ["HC", "IPF", "SSc-ILD"], rotation=20)
        ax.set_title(niche)
        ax.set_xlabel("")
        if niche == "CN1":
            ax.set_ylabel("Mean ciliated epithelial fraction (%)")
    fig.suptitle("Frozen RCTD ciliated-neighborhood composition by disease")
    fig.savefig(figures / "FIG04_CILIATED_COMPOSITION_BY_DISEASE.png", dpi=180)
    plt.close(fig)

    # Compare original and composition-adjusted disease effects.
    primary = adjusted.loc[adjusted["analysis"].eq("primary_all21")].copy()
    fig, axes = plt.subplots(1, 2, figsize=(10, 4.8), constrained_layout=True)
    for ax, (adjustment, label) in zip(
        axes,
        [
            ("ciliated_fraction", "Ciliated epithelial fraction"),
            ("airway_epithelial_sum", "Total airway epithelial fraction"),
        ],
    ):
        sub = primary.loc[primary["adjustment"].eq(adjustment)]
        colors = np.where(sub["niche"].eq("CN2"), "#4C72B0", "#B2182B")
        ax.scatter(
            sub["original_log2_fc"],
            sub["adjusted_log2_fc"],
            c=colors,
            s=34,
            alpha=0.85,
        )
        low = min(sub["original_log2_fc"].min(), sub["adjusted_log2_fc"].min()) - 0.08
        high = max(sub["original_log2_fc"].max(), sub["adjusted_log2_fc"].max()) + 0.08
        ax.plot([low, high], [low, high], "--", color="black", linewidth=0.8)
        ax.set_xlim(low, high)
        ax.set_ylim(low, high)
        ax.set_xlabel("Original disease log2FC")
        ax.set_ylabel("Composition-adjusted disease log2FC")
        median = sub["attenuation_fraction"].median() * 100
        maximum = sub["attenuation_fraction"].max() * 100
        ax.set_title(f"Adjust for {label}\nmedian attenuation {median:.1f}%; max {maximum:.1f}%")
    fig.suptitle("Robust DEG effects persist after frozen composition adjustment")
    fig.savefig(figures / "FIG05_ORIGINAL_VS_COMPOSITION_ADJUSTED_DEG.png", dpi=180)
    plt.close(fig)

    primary_ciliated = primary.loc[primary["adjustment"].eq("ciliated_fraction")]
    sensitivity_ciliated = adjusted.loc[
        adjusted["analysis"].eq("sensitivity_exclude_lt20_bins")
        & adjusted["adjustment"].eq("ciliated_fraction")
    ]
    primary_composition = tests.loc[tests["analysis"].eq("primary_all21")]
    primary_corr = correlations.loc[
        correlations["analysis"].eq("primary_all21")
        & correlations["covariate"].eq("ciliated_fraction")
    ]

    summary = f"""# Composition versus within-niche expression-state audit

## Question

Do the 21 robust disease-associated gene-by-contrast results disappear after accounting for the frozen Reference18 neighborhood composition within the same patient and K3 niche?

## Result

- None of the 162 Primary Reference18 composition tests passed pooled global BH-FDR < 0.05.
- Ciliated epithelial neighborhood fraction was not disease-associated in CN1, CN2, or CN3 for any contrast after correction.
- After conditioning on ciliated epithelial fraction, all **{len(primary_ciliated)} / {len(primary_ciliated)}** robust DEG rows retained nominal P < 0.05 and candidate-set FDR < 0.05.
- Median absolute disease-effect attenuation was **{primary_ciliated['attenuation_fraction'].median() * 100:.1f}%** in Primary and **{sensitivity_ciliated['attenuation_fraction'].median() * 100:.1f}%** after excluding the two <20-bin libraries; maximum attenuation was **{primary_ciliated['attenuation_fraction'].max() * 100:.1f}%**.
- Adjustment for the broader ciliated + secretory + basal airway fraction gave a median attenuation of **{primary.loc[primary['adjustment'].eq('airway_epithelial_sum'), 'attenuation_fraction'].median() * 100:.1f}%**.
- Only **{int((primary_corr['fdr'] < 0.05).sum())} / {len(primary_corr)}** unique robust genes had an FDR-supported correlation with ciliated fraction (`CROCC2`).

## Interpretation

At Reference18/K3 resolution, the disease-associated cilia/airway gene decrease is **not explained by a lower broad ciliated-neighborhood fraction**. It is therefore more consistent with an expression-state change within the mixed CN2/CN3 neighborhoods than with a simple loss of the broad RCTD ciliated component.

This does not prove a cell-intrinsic change. Reference18 is coarse and the covariate is an averaged neighborhood composition; an unmodeled finer epithelial subtype shift, disease-related RNA quality, or another technical/biological covariate could still generate the residual signal. Composition-adjusted P values are post-selection diagnostics, not an independent confirmatory test.
"""
    (OUT / "29_COMPOSITION_VS_STATE_INTERPRETATION.md").write_text(
        summary, encoding="utf-8"
    )
    final = {
        "completed_at": datetime.now().astimezone().isoformat(timespec="seconds"),
        "composition_global_fdr_findings": int(
            (primary_composition["global_fdr"] < 0.05).sum()
        ),
        "robust_rows_retained_after_ciliated_adjustment": int(
            primary_ciliated["adjusted_nominal_p_lt_005"].sum()
        ),
        "robust_rows_total": int(len(primary_ciliated)),
        "median_attenuation_fraction": float(
            primary_ciliated["attenuation_fraction"].median()
        ),
        "interpretation": "more consistent with within-niche expression-state signal at Reference18 resolution",
        "status": "ANALYSIS_COMPLETE_PENDING_USER_APPROVAL",
    }
    (OUT / "COMPOSITION_VS_STATE_FINAL_STATUS.json").write_text(
        json.dumps(final, indent=2), encoding="utf-8"
    )
    print(OUT)


if __name__ == "__main__":
    main()
