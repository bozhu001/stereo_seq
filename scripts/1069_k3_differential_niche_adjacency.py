#!/usr/bin/env python3
"""Patient-level differential adjacency analysis for the frozen exploratory K3 atlas."""

from __future__ import annotations

import argparse
import hashlib
import json
import platform
import sys
from itertools import combinations
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from scipy.spatial import cKDTree
from scipy.stats import t as student_t


SEED = 20261011
ADJACENCY_RADIUS_UM = 36.0
PERMUTATIONS = 10_000
NICHES = ("CN1", "CN2", "CN3")
DISEASES = ("HC", "IPF", "SSc-ILD")
CONTRASTS = ("IPF_vs_HC", "SSc-ILD_vs_HC", "SSc-ILD_vs_IPF")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--assignment-file", required=True, type=Path)
    parser.add_argument("--output-dir", required=True, type=Path)
    return parser.parse_args()


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def bh_adjust(values: np.ndarray) -> np.ndarray:
    values = np.asarray(values, dtype=float)
    order = np.argsort(values)
    ranked = values[order]
    adjusted = np.minimum.accumulate(
        (ranked * len(ranked) / np.arange(1, len(ranked) + 1))[::-1]
    )[::-1]
    output = np.empty_like(adjusted)
    output[order] = np.minimum(adjusted, 1)
    return output


def design_matrix(chip: np.ndarray, disease: np.ndarray) -> np.ndarray:
    return np.column_stack(
        [
            np.ones(len(chip)),
            (chip == "J2").astype(float),
            (chip == "L3").astype(float),
            (disease == "IPF").astype(float),
            (disease == "SSc-ILD").astype(float),
        ]
    )


def fit_hc3(y: np.ndarray, x: np.ndarray) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    xtx_inv = np.linalg.pinv(x.T @ x)
    beta = xtx_inv @ x.T @ y
    residual = y - x @ beta
    leverage = np.sum((x @ xtx_inv) * x, axis=1)
    adjusted = residual / np.maximum(1 - leverage, 1e-8)
    meat = x.T @ ((adjusted**2)[:, None] * x)
    covariance = xtx_inv @ meat @ xtx_inv
    se = np.sqrt(np.maximum(np.diag(covariance), 0))
    statistic = np.divide(beta, se, out=np.zeros_like(beta), where=se > 0)
    p = 2 * student_t.sf(np.abs(statistic), df=max(1, len(y) - x.shape[1]))
    return beta, se, p


def contrasts(beta: np.ndarray) -> dict[str, float]:
    return {
        "IPF_vs_HC": float(beta[3]),
        "SSc-ILD_vs_HC": float(beta[4]),
        "SSc-ILD_vs_IPF": float(beta[4] - beta[3]),
    }


def patient_adjacency(frame: pd.DataFrame) -> tuple[pd.DataFrame, int]:
    frame = frame[frame["niche"].isin(NICHES)].reset_index(drop=True)
    coordinates = frame[["x_um", "y_um"]].to_numpy(float)
    edges = cKDTree(coordinates).query_pairs(
        r=ADJACENCY_RADIUS_UM + 1e-8,
        output_type="ndarray",
    )
    labels = frame["niche"].to_numpy()
    node_n = len(frame)
    edge_n = len(edges)
    rows = []
    for left, right in combinations(NICHES, 2):
        observed = int(
            np.sum(
                ((labels[edges[:, 0]] == left) & (labels[edges[:, 1]] == right))
                | ((labels[edges[:, 0]] == right) & (labels[edges[:, 1]] == left))
            )
        )
        left_n = int(np.sum(labels == left))
        right_n = int(np.sum(labels == right))
        expected = edge_n * 2 * left_n * right_n / (node_n * (node_n - 1))
        enrichment = np.log2((observed + 0.5) / (expected + 0.5))
        rows.append(
            {
                "niche_pair": f"{left}--{right}",
                "left_niche": left,
                "right_niche": right,
                "observed_edge_n": observed,
                "expected_edge_n": expected,
                "log2_adjacency_enrichment": enrichment,
                "edge_n": edge_n,
                "evaluable_bin_n": node_n,
            }
        )
    return pd.DataFrame(rows), edge_n


def main() -> None:
    args = parse_args()
    assignment_path = args.assignment_file.resolve()
    output = args.output_dir.resolve()
    output.mkdir(parents=True, exist_ok=False)
    figures = output / "figures"
    figures.mkdir()

    data = pd.read_csv(assignment_path, sep="\t")
    observed_niches = sorted(data.loc[data["niche"] != "NOT_EVALUABLE", "niche"].unique())
    if observed_niches != list(NICHES):
        raise RuntimeError(f"Expected frozen K3 niches, observed {observed_niches}")

    patient_rows = []
    for patient_id, frame in data.groupby("patient_id", sort=True):
        adjacency, _ = patient_adjacency(frame)
        adjacency.insert(0, "patient_id", patient_id)
        adjacency.insert(1, "chip_label", frame["chip_label"].iloc[0])
        adjacency.insert(2, "disease", frame["disease"].iloc[0])
        patient_rows.append(adjacency)
    patient_table = pd.concat(patient_rows, ignore_index=True)

    result_rows = []
    rng = np.random.default_rng(SEED)
    patient_metadata = data.drop_duplicates("patient_id")[["patient_id", "chip_label", "disease"]]
    for niche_pair, frame in patient_table.groupby("niche_pair", sort=True):
        frame = patient_metadata.merge(frame, on=["patient_id", "chip_label", "disease"], how="left")
        y = frame["log2_adjacency_enrichment"].to_numpy(float)
        chip = frame["chip_label"].to_numpy()
        disease = frame["disease"].to_numpy()
        x = design_matrix(chip, disease)
        beta, se, p = fit_hc3(y, x)
        observed_effects = contrasts(beta)
        permuted_effects = {name: np.empty(PERMUTATIONS) for name in CONTRASTS}
        chip_indices = [np.flatnonzero(chip == value) for value in sorted(np.unique(chip))]
        for permutation in range(PERMUTATIONS):
            permuted_disease = disease.copy()
            for indices in chip_indices:
                permuted_disease[indices] = rng.permutation(permuted_disease[indices])
            permuted_beta, _, _ = fit_hc3(y, design_matrix(chip, permuted_disease))
            for name, value in contrasts(permuted_beta).items():
                permuted_effects[name][permutation] = value
        for name in CONTRASTS:
            if name == "IPF_vs_HC":
                coefficient_index = 3
                contrast_se = se[3]
                contrast_p = p[3]
            elif name == "SSc-ILD_vs_HC":
                coefficient_index = 4
                contrast_se = se[4]
                contrast_p = p[4]
            else:
                coefficient_index = -1
                contrast_se = np.nan
                contrast_p = np.nan
            effect = observed_effects[name]
            permutation_p = (1 + np.sum(np.abs(permuted_effects[name]) >= abs(effect))) / (PERMUTATIONS + 1)
            result_rows.append(
                {
                    "niche_pair": niche_pair,
                    "contrast": name,
                    "difference_in_log2_adjacency_enrichment": effect,
                    "hc3_se": contrast_se,
                    "hc3_p": contrast_p,
                    "within_chip_permutation_p": permutation_p,
                    "coefficient_index": coefficient_index,
                }
            )
    results = pd.DataFrame(result_rows)
    results["permutation_fdr"] = bh_adjust(results["within_chip_permutation_p"].to_numpy())
    results["supported_finding"] = results["permutation_fdr"] < 0.05

    disease_summary = (
        patient_table.groupby(["disease", "niche_pair"], observed=True)["log2_adjacency_enrichment"]
        .agg(patient_n="size", mean="mean", median="median", sd="std")
        .reset_index()
    )
    patient_table.to_csv(output / "01_PATIENT_K3_NICHE_ADJACENCY.tsv", sep="\t", index=False)
    disease_summary.to_csv(output / "02_DISEASE_ADJACENCY_SUMMARY.tsv", sep="\t", index=False)
    results.to_csv(output / "03_DIFFERENTIAL_NICHE_ADJACENCY.tsv", sep="\t", index=False)
    pd.DataFrame(
        {"input_file": [str(assignment_path)], "sha256": [sha256(assignment_path)]}
    ).to_csv(output / "00_INPUT_MANIFEST.tsv", sep="\t", index=False)

    sns.set_theme(style="whitegrid")
    pair_order = ["CN1--CN2", "CN1--CN3", "CN2--CN3"]
    disease_order = ["HC", "IPF", "SSc-ILD"]
    palette = {"HC": "#4C72B0", "IPF": "#DD8452", "SSc-ILD": "#55A868"}
    fig, axes = plt.subplots(1, 3, figsize=(13, 4.5), sharey=True, constrained_layout=True)
    for axis, niche_pair in zip(axes, pair_order):
        frame = patient_table[patient_table["niche_pair"].eq(niche_pair)]
        sns.boxplot(
            data=frame,
            x="disease",
            y="log2_adjacency_enrichment",
            order=disease_order,
            hue="disease",
            palette=palette,
            legend=False,
            showfliers=False,
            ax=axis,
        )
        sns.stripplot(
            data=frame,
            x="disease",
            y="log2_adjacency_enrichment",
            order=disease_order,
            color="black",
            jitter=0.12,
            size=5,
            ax=axis,
        )
        axis.axhline(0, color="grey", linestyle="--", linewidth=1)
        axis.set_title(niche_pair)
        axis.set_xlabel("")
        axis.set_ylabel("Log2 adjacency enrichment" if niche_pair == pair_order[0] else "")
    fig.suptitle("K=3 niche–niche spatial adjacency by disease (one point per patient)")
    fig.savefig(figures / "FIG01_PATIENT_K3_NICHE_ADJACENCY_BY_DISEASE.png", dpi=300, bbox_inches="tight")
    fig.savefig(figures / "FIG01_PATIENT_K3_NICHE_ADJACENCY_BY_DISEASE.pdf", bbox_inches="tight")
    plt.close(fig)

    contrast_order = list(CONTRASTS)
    effect = results.pivot(index="contrast", columns="niche_pair", values="difference_in_log2_adjacency_enrichment").reindex(index=contrast_order, columns=pair_order)
    fdr = results.pivot(index="contrast", columns="niche_pair", values="permutation_fdr").reindex(index=contrast_order, columns=pair_order)
    annotation = effect.copy().astype(object)
    for contrast in contrast_order:
        for niche_pair in pair_order:
            annotation.loc[contrast, niche_pair] = f"{effect.loc[contrast, niche_pair]:+.2f}\nq={fdr.loc[contrast, niche_pair]:.2f}"
    limit = max(0.5, float(np.nanmax(np.abs(effect.to_numpy()))))
    fig, ax = plt.subplots(figsize=(7, 4.5))
    sns.heatmap(
        effect,
        cmap="vlag",
        center=0,
        vmin=-limit,
        vmax=limit,
        annot=annotation,
        fmt="",
        linewidths=0.8,
        cbar_kws={"label": "Difference in log2 adjacency enrichment"},
        ax=ax,
    )
    ax.set_title("Disease-associated changes in K=3 niche adjacency\nchip-adjusted; within-chip permutation FDR")
    ax.set_xlabel("Niche pair")
    ax.set_ylabel("Contrast")
    fig.tight_layout()
    fig.savefig(figures / "FIG02_DIFFERENTIAL_K3_NICHE_ADJACENCY_HEATMAP.png", dpi=300, bbox_inches="tight")
    fig.savefig(figures / "FIG02_DIFFERENTIAL_K3_NICHE_ADJACENCY_HEATMAP.pdf", bbox_inches="tight")
    plt.close(fig)

    validation = pd.DataFrame(
        [
            ("patients_21", patient_table["patient_id"].nunique() == 21, patient_table["patient_id"].nunique()),
            ("three_inter_niche_pairs", patient_table["niche_pair"].nunique() == 3, patient_table["niche_pair"].nunique()),
            ("patient_pair_rows_63", len(patient_table) == 63, len(patient_table)),
            ("disease_tests_9", len(results) == 9, len(results)),
            ("finite_enrichment", np.isfinite(patient_table["log2_adjacency_enrichment"]).all(), len(patient_table)),
            ("finite_effects", np.isfinite(results["difference_in_log2_adjacency_enrichment"]).all(), len(results)),
        ],
        columns=["check", "pass", "observed"],
    )
    validation.to_csv(output / "04_FINAL_VALIDATION.tsv", sep="\t", index=False)
    if not validation["pass"].all():
        raise RuntimeError("Validation failed")

    summary = [
        "# K3 differential niche adjacency",
        "",
        "- Frozen input: common exploratory K=3 assignments; K=4 is excluded.",
        "- Spatial edges: exact Euclidean distance <=36 um, within patient only.",
        "- Per-patient adjacency is log2(observed/expected), where expected preserves patient-specific niche counts.",
        "- Disease effects are patient-level, chip-adjusted, with 10,000 within-chip disease-label permutations.",
        f"- Supported differential adjacency findings at permutation FDR <0.05: {int(results['supported_finding'].sum())}/{len(results)}.",
    ]
    (output / "05_SUMMARY.md").write_text("\n".join(summary) + "\n", encoding="utf-8")
    (output / "SESSION_INFO.json").write_text(
        json.dumps(
            {
                "python": sys.version,
                "platform": platform.platform(),
                "numpy": np.__version__,
                "pandas": pd.__version__,
            },
            indent=2,
        ),
        encoding="utf-8",
    )
    command = f'"{sys.executable}" "{Path(__file__).resolve()}" --assignment-file "{assignment_path}" --output-dir "{output}"'
    (output / "RUN_LOG.txt").write_text(command + "\n", encoding="utf-8")
    (output / "COMPLETE.ok").write_text("PASS\n", encoding="utf-8")
    print(output)


if __name__ == "__main__":
    main()
