#!/usr/bin/env python3
"""Create the final English summary figure for harmonized K=7 robustness.

Inputs are persisted disease-blind K2-K12 rank audits, primary K=7 loadings,
state-removal sensitivity models, the independent factor-matching review, and
the K8 20,000-iteration convergence extension. No NMF model is refitted.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import platform
import sys
from pathlib import Path

import matplotlib
import numpy as np
import pandas as pd

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import seaborn as sns


ROOT = Path(__file__).resolve().parents[1]
RANK_RUN = ROOT / (
    "results/reanalysis/"
    "REFERENCE18_NEIGHBORHOOD_NMF_K2_K12_RANK_SELECTION_20261016_110753"
)
PRIMARY_RUN = ROOT / (
    "results/reanalysis/"
    "REFERENCE18_NEIGHBORHOOD_NMF_K7_K12_EXTENSION_20261015_attempt5"
)
SENSITIVITY_RUN = ROOT / (
    "results/reanalysis/"
    "REFERENCE18_NEIGHBORHOOD_NMF_K7_STATE_SENSITIVITY_20261016_210302_attempt5"
)
REVIEW_RUN = ROOT / (
    "results/reanalysis/"
    "REFERENCE18_NEIGHBORHOOD_NMF_K7_SENSITIVITY_EVIDENCE_REVIEW_20261016_212300"
)
CONVERGENCE_RUN = ROOT / (
    "results/reanalysis/"
    "REFERENCE18_NEIGHBORHOOD_NMF_K8_K7_CONVERGENCE_EXTENSION_20261016_212100"
)
CHIPS = ("K8", "J2", "L3")
FAMILY_ORDER = (
    "LEC-associated",
    "Blood-endothelial",
    "AT1-associated",
    "Perivascular-stromal",
    "B-enriched/mixed-immune",
    "Ciliated-associated",
    "Fibroblast-associated",
)
STATE_ORDER = (
    "Alveolar/resident macrophage",
    "Inflammatory monocyte/macrophage",
    "Dendritic cell",
    "T",
    "NK",
    "Mast/Basophil",
    "B",
    "Plasma",
    "Ciliated epithelial",
    "AT2/alveolar transitional",
    "Secretory/other airway epithelial",
    "AT1",
    "Basal-like epithelial",
    "Blood endothelial",
    "Lymphatic endothelial",
    "Perivascular/mesothelial stromal",
    "Fibroblast",
    "Myofibroblast",
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output-dir", type=Path, required=True)
    return parser.parse_args()


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def add_panel_letter(axis, letter: str) -> None:
    axis.text(
        -0.10,
        1.08,
        letter,
        transform=axis.transAxes,
        fontsize=16,
        fontweight="bold",
        va="top",
        ha="left",
    )


def load_inputs():
    paths = {
        "rank": RANK_RUN / "01_K2_K12_RANK_AUDIT.tsv",
        "cross": SENSITIVITY_RUN / "07_CROSS_CHIP_FAMILY_MAP.tsv",
        "independent": REVIEW_RUN / "01_INDEPENDENT_BEST_MATCHES.tsv",
        "retained": REVIEW_RUN / "02_RETAINED_FAMILY_SUMMARY.tsv",
        "b": REVIEW_RUN / "03_B_SUCCESSOR_SUMMARY.tsv",
        "decisions": REVIEW_RUN / "04_DECISIONS.tsv",
        "convergence": CONVERGENCE_RUN / "01_CONVERGENCE_EXTENSION.tsv",
    }
    rank = pd.read_csv(paths["rank"], sep="\t")
    cross = pd.read_csv(paths["cross"], sep="\t")
    independent = pd.read_csv(paths["independent"], sep="\t")
    retained = pd.read_csv(paths["retained"], sep="\t")
    b_summary = pd.read_csv(paths["b"], sep="\t")
    decisions = pd.read_csv(paths["decisions"], sep="\t")
    convergence = pd.read_csv(paths["convergence"], sep="\t")

    loadings = {}
    for chip in CHIPS:
        path = PRIMARY_RUN / "models" / f"neighborhood50_{chip}_K07_LOADINGS.tsv"
        paths[f"loading_{chip}"] = path
        loadings[chip] = pd.read_csv(path, sep="\t", index_col=0)
    return (
        paths,
        rank,
        cross,
        independent,
        retained,
        b_summary,
        decisions,
        convergence,
        loadings,
    )


def family_centroids(cross: pd.DataFrame, loadings: dict[str, pd.DataFrame]) -> pd.DataFrame:
    main = cross.loc[cross["scenario"].eq("main18")].copy()
    rows = []
    for family in FAMILY_ORDER:
        row = main.loc[main["family_label"].eq(family)].iloc[0]
        vectors = []
        for chip in CHIPS:
            factor = row[f"{chip}_factor"]
            vectors.append(loadings[chip].loc[factor, list(STATE_ORDER)].to_numpy(float))
        centroid = np.mean(vectors, axis=0)
        centroid /= centroid.sum()
        rows.append(centroid)
    return pd.DataFrame(rows, index=FAMILY_ORDER, columns=STATE_ORDER)


def main() -> None:
    args = parse_args()
    output = args.output_dir.resolve()
    output.mkdir(parents=True, exist_ok=False)
    (
        paths,
        rank,
        cross,
        independent,
        retained,
        b_summary,
        decisions,
        convergence,
        loadings,
    ) = load_inputs()
    if not decisions["pass"].all():
        raise RuntimeError("Evidence review contains a failed decision")
    centroids = family_centroids(cross, loadings)

    sns.set_theme(style="whitegrid", context="notebook")
    colors = {"K8": "#4C78A8", "J2": "#F58518", "L3": "#54A24B"}
    fig = plt.figure(figsize=(20, 17), constrained_layout=True)
    grid = fig.add_gridspec(3, 2, height_ratios=(1.0, 1.22, 1.0))

    # A: rank stability with K=7 fixed for harmonized comparison.
    axis = fig.add_subplot(grid[0, 0])
    for chip in CHIPS:
        frame = rank.loc[rank["chip"].eq(chip)].sort_values("rank")
        axis.plot(
            frame["rank"],
            frame["restart_fraction_all_factors_ge_0_80"],
            marker="o",
            color=colors[chip],
            label=f"{chip}: restart",
        )
        axis.plot(
            frame["rank"],
            frame["lopo_fraction_all_factors_ge_0_80"],
            marker="s",
            linestyle="--",
            color=colors[chip],
            alpha=0.78,
            label=f"{chip}: LOPO",
        )
    axis.axvline(7, color="#B22222", linewidth=2, linestyle=":")
    axis.axhline(0.80, color="grey", linewidth=1, linestyle=":")
    axis.set(xlabel="NMF rank K", ylabel="Fraction passing cosine ≥ 0.80", ylim=(-0.03, 1.04))
    axis.set_xticks(range(2, 13))
    axis.set_title("Harmonized K=7 is supported in all three chips", fontweight="bold")
    axis.legend(ncol=2, fontsize=8, loc="lower left")
    add_panel_letter(axis, "A")

    # B: the seven cross-chip centroid composition signatures.
    axis = fig.add_subplot(grid[0, 1])
    sns.heatmap(
        centroids,
        cmap="mako",
        vmin=0,
        vmax=float(np.ceil(centroids.to_numpy().max() * 10) / 10),
        ax=axis,
        cbar_kws={"label": "Mean normalized state loading"},
    )
    axis.set_title("Seven shared K=7 composition signatures", fontweight="bold")
    axis.set_xlabel("Reference18 state")
    axis.set_ylabel("")
    axis.tick_params(axis="x", rotation=65, labelsize=7)
    axis.tick_params(axis="y", rotation=0, labelsize=8)
    add_panel_letter(axis, "B")

    # C: explicit cross-chip reproducibility for each primary family.
    axis = fig.add_subplot(grid[1, 0])
    main_cross = cross.loc[cross["scenario"].eq("main18")].set_index("family_label")
    pairwise = main_cross.loc[
        FAMILY_ORDER, ["K8_J2_cosine", "K8_L3_cosine", "J2_L3_cosine"]
    ].copy()
    pairwise.columns = ["K8–J2", "K8–L3", "J2–L3"]
    sns.heatmap(
        pairwise,
        cmap="YlGnBu",
        vmin=0.80,
        vmax=1.0,
        annot=True,
        fmt=".2f",
        ax=axis,
        cbar_kws={"label": "Loading cosine"},
    )
    axis.set_title("Every primary factor family is shared across chips", fontweight="bold")
    axis.set_xlabel("")
    axis.set_ylabel("")
    axis.tick_params(axis="x", rotation=0)
    axis.tick_params(axis="y", rotation=0, labelsize=8)
    add_panel_letter(axis, "C")

    # D: retained-factor composition and spatial persistence under state removal.
    axis = fig.add_subplot(grid[1, 1])
    display = retained.copy()
    display["model"] = display["scenario"].map(
        {"no_basal_lec16": "Remove Basal + LEC", "no_b17": "Remove B"}
    )
    display["row"] = display["model"] + " | " + display["chip"]
    matrix = display.set_index("row")[[
        "minimum_best_loading_cosine",
        "minimum_patient_usage_spearman",
        "minimum_q90_hotspot_jaccard",
    ]]
    matrix.columns = [
        "Min loading\ncosine",
        "Min median-patient\nscore Spearman",
        "Min median-patient\nQ90 Jaccard",
    ]
    desired_rows = [
        f"{model} | {chip}"
        for model in ("Remove Basal + LEC", "Remove B")
        for chip in CHIPS
    ]
    matrix = matrix.reindex(desired_rows)
    sns.heatmap(
        matrix,
        cmap="viridis",
        vmin=0,
        vmax=1,
        annot=True,
        fmt=".2f",
        ax=axis,
        cbar_kws={"label": "Persistence metric"},
    )
    axis.set_title("The six retained families remain compositionally and spatially stable", fontweight="bold")
    axis.set_xlabel("")
    axis.set_ylabel("")
    axis.tick_params(axis="x", rotation=0, labelsize=8)
    axis.tick_params(axis="y", rotation=0, labelsize=8)
    add_panel_letter(axis, "D")

    # E: the deleted B-enriched factor has no coherent single-factor successor.
    axis = fig.add_subplot(grid[2, 0])
    positions = np.arange(len(CHIPS))
    b_ordered = b_summary.set_index("chip").loc[list(CHIPS)]
    width = 0.34
    axis.bar(
        positions - width / 2,
        b_ordered["maximum_candidate_patient_usage_spearman"],
        width,
        label="Best median-patient score Spearman",
        color="#E45756",
    )
    axis.bar(
        positions + width / 2,
        b_ordered["maximum_candidate_q90_hotspot_jaccard"],
        width,
        label="Best median-patient Q90 Jaccard",
        color="#72B7B2",
    )
    axis.axhline(0.30, color="#E45756", linestyle=":", linewidth=1)
    axis.axhline(0.10, color="#72B7B2", linestyle=":", linewidth=1)
    axis.set_xticks(positions, CHIPS)
    axis.set_ylim(0, 0.36)
    axis.set_ylabel("Maximum across all seven B-removed factors")
    axis.set_title("The B-enriched spatial axis is not relabeled as one T/NK factor", fontweight="bold")
    axis.legend(fontsize=8, loc="upper right")
    add_panel_letter(axis, "E")

    # F: higher iteration allowance confirms the K8 sensitivity solutions.
    axis = fig.add_subplot(grid[2, 1])
    conv = convergence.copy()
    conv["label"] = conv["scenario"].map(
        {"no_basal_lec16": "Remove Basal + LEC", "no_b17": "Remove B"}
    )
    x = np.arange(len(conv))
    axis.bar(
        x - width / 2,
        conv["convergence_fraction"],
        width,
        color="#4C78A8",
        label="20,000-step convergence fraction",
    )
    axis.bar(
        x + width / 2,
        conv["matched_5k_to_20k_min_cosine"],
        width,
        color="#F2CF5B",
        label="Minimum 5k-to-20k factor cosine",
    )
    axis.axhline(0.95, color="grey", linestyle=":", linewidth=1)
    axis.set_xticks(x, conv["label"])
    axis.set_ylim(0, 1.04)
    axis.set_ylabel("Fraction or cosine")
    axis.set_title("Extended K8 optimization confirms the same solutions", fontweight="bold")
    axis.legend(fontsize=8, loc="lower right")
    add_panel_letter(axis, "F")

    fig.suptitle(
        "Harmonized K=7 Neighborhood NMF is cross-chip reproducible and state-removal robust",
        fontsize=20,
        fontweight="bold",
    )
    footer = (
        "Primary model: frozen 18-state, centre-excluded 50-µm Neighborhood NMF; "
        "50 restarts; 495 balanced observations per patient; disease labels unused. "
        "B is interpreted as B-enriched/mixed-immune; LEC remains sensitivity-qualified."
    )
    fig.text(0.5, -0.005, footer, ha="center", va="top", fontsize=9)
    png = output / "REFERENCE18_K7_ROBUSTNESS_SUMMARY.png"
    pdf = output / "REFERENCE18_K7_ROBUSTNESS_SUMMARY.pdf"
    svg = output / "REFERENCE18_K7_ROBUSTNESS_SUMMARY.svg"
    fig.savefig(png, dpi=300, bbox_inches="tight", facecolor="white")
    fig.savefig(pdf, bbox_inches="tight", facecolor="white")
    fig.savefig(svg, bbox_inches="tight", facecolor="white")
    plt.close(fig)

    frozen = pd.DataFrame(
        [
            ("primary_cross_chip_rank", "K=7 for K8, J2, and L3"),
            ("higher_rank_models", "deferred; revisit only if common K=7 fails a downstream requirement"),
            ("primary_families", ";".join(FAMILY_ORDER)),
            ("B_label", "B-enriched/mixed-immune; retained as a primary immune-focused neighborhood factor, not pure B-cell abundance"),
            ("LEC_label", "LEC-associated; upstream RCTD evidence is weak and interpretation is sensitivity-qualified"),
            ("T_NK_interpretation", "signals are present but do not form separate primary K=7 factors; absence of a factor is not absence of biology"),
            ("gene_expression_use", "independent corroboration for biological naming; not required to define the RCTD-composition NMF factors"),
            ("disease_use_in_discovery", "none"),
        ],
        columns=["decision", "frozen_value"],
    )
    frozen.to_csv(output / "01_FROZEN_K7_INTERPRETATION.tsv", sep="\t", index=False)
    pd.DataFrame(
        [{"path": str(path.resolve()), "sha256": sha256(path)} for path in paths.values()]
    ).to_csv(output / "00_INPUT_MANIFEST.tsv", sep="\t", index=False)
    checks = pd.DataFrame(
        [
            ("all_evidence_decisions_pass", bool(decisions["pass"].all()), str(decisions["pass"].tolist())),
            ("seven_primary_families", centroids.shape == (7, 18), str(centroids.shape)),
            ("all_cross_chip_cosines_ge_0_80", bool((pairwise.to_numpy() >= 0.80).all()), f"minimum={pairwise.to_numpy().min():.6f}"),
            ("all_three_formats_nonempty", all(path.exists() and path.stat().st_size > 0 for path in (png, pdf, svg)), str([path.stat().st_size for path in (png, pdf, svg)])),
            ("disease_blind", True, "No disease field used in plotting decisions"),
        ],
        columns=["check", "pass", "detail"],
    )
    checks.to_csv(output / "02_VALIDATION.tsv", sep="\t", index=False)
    if not checks["pass"].all():
        raise RuntimeError("Figure validation failed")
    status = {
        "status": "FIGURE_COMPLETE_PENDING_VISUAL_QA",
        "png": str(png),
        "pdf": str(pdf),
        "svg": str(svg),
    }
    (output / "FINAL_STATUS.json").write_text(json.dumps(status, indent=2), encoding="utf-8")
    (output / "RUN_LOG.txt").write_text(
        f"Script: scripts/{Path(__file__).name}\n"
        f"Command: {sys.executable} scripts/{Path(__file__).name} --output-dir {output}\n"
        f"Python: {sys.version.replace(chr(10), ' ')}\n"
        f"Platform: {platform.platform()}\n"
        "NMF refit: no\nDisease used: no\nStatus: PASS\n",
        encoding="utf-8",
    )
    print(json.dumps(status, indent=2))


if __name__ == "__main__":
    main()
