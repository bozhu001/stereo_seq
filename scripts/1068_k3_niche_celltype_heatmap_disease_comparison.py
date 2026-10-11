#!/usr/bin/env python3
"""Create boss-facing niche composition and disease-comparison deliverables."""

from __future__ import annotations

import argparse
import hashlib
import json
import platform
import sys
from pathlib import Path

import h5py
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--atlas-dir", required=True, type=Path)
    parser.add_argument("--output-dir", required=True, type=Path)
    parser.add_argument("--k", required=True, type=int)
    return parser.parse_args()


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def main() -> None:
    args = parse_args()
    atlas_dir = args.atlas_dir.resolve()
    output = args.output_dir.resolve()
    output.mkdir(parents=True, exist_ok=False)
    figures = output / "figures"
    figures.mkdir()

    assignment_path = atlas_dir / "02_ALL21_SHARED_NICHE_ASSIGNMENTS.tsv.gz"
    h5_path = atlas_dir / "12_SHARED_NICHE_ASSAY.h5"
    patient_path = atlas_dir / "06_PATIENT_NICHE_ABUNDANCE.tsv"
    disease_path = atlas_dir / "09_NICHE_DISEASE_TESTS.tsv"
    input_paths = [assignment_path, h5_path, patient_path, disease_path]
    if not all(path.exists() for path in input_paths):
        raise FileNotFoundError("Required atlas outputs are missing")

    assignments = pd.read_csv(assignment_path, sep="\t")
    patients = pd.read_csv(patient_path, sep="\t")
    disease_tests = pd.read_csv(disease_path, sep="\t")
    with h5py.File(h5_path, "r") as handle:
        composition = handle["composition"][:]
        state_names = [value.decode() if isinstance(value, bytes) else str(value) for value in handle["state_names"][:]]
        h5_keys = [value.decode() if isinstance(value, bytes) else str(value) for value in handle["entity_key"][:]]

    if h5_keys != assignments["entity_key"].astype(str).tolist():
        raise RuntimeError("HDF5 and assignment entity order differ")
    evaluable = assignments["niche"].ne("NOT_EVALUABLE").to_numpy()
    niches = sorted(assignments.loc[evaluable, "niche"].unique())
    expected_niches = [f"CN{index}" for index in range(1, args.k + 1)]
    if niches != expected_niches:
        raise RuntimeError(f"Expected exploratory K={args.k} labels, observed {niches}")

    mean_rows = []
    for niche in niches:
        mask = assignments["niche"].eq(niche).to_numpy()
        values = composition[mask]
        for index, state in enumerate(state_names):
            mean_rows.append(
                {
                    "niche": niche,
                    "cell_type": state,
                    "mean_rctd_fraction": float(np.nanmean(values[:, index])),
                    "median_rctd_fraction": float(np.nanmedian(values[:, index])),
                    "bin_n": int(mask.sum()),
                }
            )
    composition_long = pd.DataFrame(mean_rows)
    raw_matrix = composition_long.pivot(index="niche", columns="cell_type", values="mean_rctd_fraction").reindex(index=niches, columns=state_names)
    enrichment = raw_matrix.copy()
    for column in enrichment.columns:
        sd = enrichment[column].std(ddof=0)
        enrichment[column] = 0.0 if sd == 0 else (enrichment[column] - enrichment[column].mean()) / sd
    enrichment_long = enrichment.stack().rename("celltype_z_across_niches").reset_index()
    composition_long = composition_long.merge(enrichment_long, on=["niche", "cell_type"], how="left")

    fraction_columns = [f"fraction_{niche}" for niche in niches]
    patient_long = patients.melt(
        id_vars=[
            "patient_id",
            "chip_label",
            "disease",
            "eligible_bin_n",
            "evaluable_bin_n",
            "median_nUMI",
        ],
        value_vars=fraction_columns,
        var_name="niche",
        value_name="patient_niche_fraction",
    )
    patient_long["niche"] = patient_long["niche"].str.replace("fraction_", "", regex=False)
    disease_summary = (
        patient_long.groupby(["disease", "niche"], observed=True)["patient_niche_fraction"]
        .agg(patient_n="size", mean="mean", median="median", sd="std", minimum="min", maximum="max")
        .reset_index()
    )

    composition_long.to_csv(output / "01_K3_NICHE_CELLTYPE_COMPOSITION.tsv", sep="\t", index=False)
    patient_long.to_csv(output / "02_PATIENT_NICHE_FRACTIONS_LONG.tsv", sep="\t", index=False)
    disease_summary.to_csv(output / "03_DISEASE_NICHE_SUMMARY.tsv", sep="\t", index=False)
    disease_tests.to_csv(output / "04_PATIENT_LEVEL_DISEASE_EFFECTS.tsv", sep="\t", index=False)
    pd.DataFrame(
        {
            "input_file": [str(path) for path in input_paths],
            "sha256": [sha256(path) for path in input_paths],
        }
    ).to_csv(output / "00_INPUT_MANIFEST.tsv", sep="\t", index=False)

    sns.set_theme(style="white")
    fig, axes = plt.subplots(2, 1, figsize=(16, 7.5), constrained_layout=True)
    sns.heatmap(
        raw_matrix * 100,
        cmap="YlGnBu",
        annot=True,
        fmt=".1f",
        linewidths=0.5,
        cbar_kws={"label": "Mean RCTD fraction (%)"},
        ax=axes[0],
    )
    axes[0].set_title(f"A. Exploratory K={args.k} niche × cell-type composition (raw RCTD fractions)")
    axes[0].set_xlabel("")
    axes[0].set_ylabel("Niche")
    sns.heatmap(
        enrichment,
        cmap="vlag",
        center=0,
        vmin=-1.5,
        vmax=1.5,
        annot=True,
        fmt=".1f",
        linewidths=0.5,
        cbar_kws={"label": "Cell-type-wise z score across niches"},
        ax=axes[1],
    )
    axes[1].set_title("B. Relative cell-type enrichment used for niche interpretation")
    axes[1].set_xlabel("Reference18 cell type")
    axes[1].set_ylabel("Niche")
    axes[1].tick_params(axis="x", rotation=55)
    for label in axes[1].get_xticklabels():
        label.set_horizontalalignment("right")
    fig.savefig(figures / f"FIG01_K{args.k}_NICHE_X_CELLTYPE_HEATMAP.png", dpi=300, bbox_inches="tight")
    fig.savefig(figures / f"FIG01_K{args.k}_NICHE_X_CELLTYPE_HEATMAP.pdf", bbox_inches="tight")
    plt.close(fig)

    disease_order = ["HC", "IPF", "SSc-ILD"]
    fig, axes = plt.subplots(1, len(niches), figsize=(12, 4), sharey=True, constrained_layout=True)
    palette = {"HC": "#4C72B0", "IPF": "#DD8452", "SSc-ILD": "#55A868"}
    for axis, niche in zip(axes, niches):
        frame = patient_long[patient_long["niche"].eq(niche)]
        sns.boxplot(data=frame, x="disease", y="patient_niche_fraction", order=disease_order, palette=palette, showfliers=False, ax=axis)
        sns.stripplot(data=frame, x="disease", y="patient_niche_fraction", order=disease_order, color="black", size=5, jitter=0.12, ax=axis)
        axis.set_title(niche)
        axis.set_xlabel("")
        axis.set_ylabel("Patient niche fraction" if niche == niches[0] else "")
    fig.suptitle(f"Exploratory K={args.k} niche abundance by disease (one point per patient)")
    fig.savefig(figures / "FIG02_PATIENT_NICHE_FRACTION_BY_DISEASE.png", dpi=300, bbox_inches="tight")
    fig.savefig(figures / "FIG02_PATIENT_NICHE_FRACTION_BY_DISEASE.pdf", bbox_inches="tight")
    plt.close(fig)

    contrast_order = ["IPF_vs_HC", "SSc-ILD_vs_HC", "SSc-ILD_vs_IPF"]
    effect = disease_tests.pivot(index="contrast", columns="niche", values="log_odds_effect").reindex(index=contrast_order, columns=niches)
    q_values = disease_tests.pivot(index="contrast", columns="niche", values="permutation_fdr_global").reindex(index=contrast_order, columns=niches)
    annotations = effect.copy().astype(object)
    for row in contrast_order:
        for niche in niches:
            annotations.loc[row, niche] = f"{effect.loc[row, niche]:+.2f}\nq={q_values.loc[row, niche]:.2f}"
    fig, ax = plt.subplots(figsize=(6.5, 4.2))
    bound = max(1.0, float(np.nanmax(np.abs(effect.to_numpy()))))
    sns.heatmap(
        effect,
        cmap="vlag",
        center=0,
        vmin=-bound,
        vmax=bound,
        annot=annotations,
        fmt="",
        linewidths=0.8,
        cbar_kws={"label": "Corrected patient-level log odds effect"},
        ax=ax,
    )
    ax.set_title("Disease differences in niche abundance\nchip-adjusted; within-chip permutation FDR")
    ax.set_xlabel("Niche")
    ax.set_ylabel("Contrast")
    fig.tight_layout()
    fig.savefig(figures / "FIG03_DISEASE_EFFECT_HEATMAP.png", dpi=300, bbox_inches="tight")
    fig.savefig(figures / "FIG03_DISEASE_EFFECT_HEATMAP.pdf", bbox_inches="tight")
    plt.close(fig)

    validation = pd.DataFrame(
        [
            ("assignment_h5_order", True, len(h5_keys)),
            (f"k{args.k}_niches", niches == expected_niches, ";".join(niches)),
            ("patients_21", patients["patient_id"].nunique() == 21, patients["patient_id"].nunique()),
            ("composition_rows", len(composition_long) == 18 * args.k, len(composition_long)),
            ("patient_fraction_rows", len(patient_long) == 21 * args.k, len(patient_long)),
            ("finite_effects", len(disease_tests) == 3 * args.k and np.isfinite(disease_tests["log_odds_effect"]).all(), len(disease_tests)),
        ],
        columns=["check", "pass", "observed"],
    )
    validation.to_csv(output / "05_FINAL_VALIDATION.tsv", sep="\t", index=False)
    if not validation["pass"].all():
        raise RuntimeError("Final validation failed")

    summary = [
        f"# Exploratory K={args.k} niche composition and disease comparison",
        "",
        "- The niche x cell-type heatmap uses uncorrected raw RCTD fractions.",
        f"- K={args.k} assignments were learned disease-blind in chip-centered CLR space.",
        "- Disease comparisons use one observation per patient, include chip as a covariate, and use within-chip permutation FDR.",
        f"- K={args.k} is a user-requested exploratory resolution and did not pass the predeclared held-chip transfer threshold; all biological conclusions remain exploratory.",
        f"- Significant niche-disease contrasts at permutation FDR < 0.05: {(disease_tests['permutation_fdr_global'] < 0.05).sum()} / {len(disease_tests)}.",
    ]
    (output / "06_SUMMARY.md").write_text("\n".join(summary) + "\n", encoding="utf-8")
    session = {
        "python": sys.version,
        "platform": platform.platform(),
        "numpy": np.__version__,
        "pandas": pd.__version__,
    }
    (output / "SESSION_INFO.json").write_text(json.dumps(session, indent=2), encoding="utf-8")
    command = f'"{sys.executable}" "{Path(__file__).resolve()}" --atlas-dir "{atlas_dir}" --output-dir "{output}" --k "{args.k}"'
    (output / "RUN_LOG.txt").write_text(command + "\n", encoding="utf-8")
    (output / "COMPLETE.ok").write_text("PASS\n", encoding="utf-8")
    print(output)


if __name__ == "__main__":
    main()
