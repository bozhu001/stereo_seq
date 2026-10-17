#!/usr/bin/env python3
"""Stress-test the harmonized per-chip Neighborhood Reference18 K=7 model.

The persisted 18-state K=7 model is the primary model. Two disease-blind
sensitivity models are fitted independently on the same frozen centre-excluded
50-um neighborhood assay: (1) remove Basal-like epithelial and Lymphatic
endothelial; (2) remove B. Every refit keeps K=7, the original patient-balanced
reference, solver, 50-restart medoid rule, and leave-one-patient-out audit.
Factor persistence is assessed from common-state loading cosine, same-Bin50
usage correlation, and patient-wise top-decile spatial overlap.
"""

from __future__ import annotations

import argparse
import hashlib
import importlib.util
import json
import platform
import sys
import time
from pathlib import Path

import matplotlib
import numpy as np
import pandas as pd
from joblib import Parallel, delayed
from scipy.optimize import linear_sum_assignment
from scipy.stats import spearmanr

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import seaborn as sns


ROOT = Path(__file__).resolve().parents[1]
BASE_SCRIPT = ROOT / "scripts/1093_reference18_per_chip_dual_nmf_rank_audit.py"
PRIMARY_RUN = ROOT / (
    "results/reanalysis/"
    "REFERENCE18_NEIGHBORHOOD_NMF_K7_K12_EXTENSION_20261015_attempt5"
)
CHIPS = ("K8", "J2", "L3")
ROUTE = "neighborhood50"
RANK = 7
SCENARIOS = {
    "main18": (),
    "no_basal_lec16": ("Basal-like epithelial", "Lymphatic endothelial"),
    "no_b17": ("B",),
}
FAMILY_ORDER = (
    "LEC-associated",
    "Blood-endothelial",
    "AT1-associated",
    "Perivascular-stromal",
    "B-enriched/mixed-immune",
    "Ciliated-associated",
    "Fibroblast-associated",
)


def load_base():
    spec = importlib.util.spec_from_file_location("dual_nmf_base", BASE_SCRIPT)
    if spec is None or spec.loader is None:
        raise RuntimeError("Could not load script 1093")
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--n-jobs", type=int, default=4)
    parser.add_argument("--main-restarts", type=int, default=50)
    parser.add_argument("--lopo-restarts", type=int, default=3)
    parser.add_argument("--max-iter", type=int, default=5_000)
    return parser.parse_args()


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def atomic_table(frame: pd.DataFrame, path: Path, **kwargs) -> None:
    temporary = path.with_name(f".{path.name}.tmp")
    frame.to_csv(temporary, **kwargs)
    temporary.replace(path)


def atomic_npz(path: Path, arrays: dict[str, np.ndarray]) -> None:
    temporary = path.with_name(f".{path.name}.tmp")
    with temporary.open("wb") as handle:
        np.savez_compressed(handle, **arrays)
    temporary.replace(path)


def subset_matrix(
    matrix: np.ndarray, states: list[str], removed: tuple[str, ...]
) -> tuple[np.ndarray, list[str]]:
    keep = [index for index, state in enumerate(states) if state not in removed]
    selected = matrix[:, keep].copy()
    totals = selected.sum(axis=1, keepdims=True)
    if np.any(totals <= 0):
        raise RuntimeError(f"Zero row after removing states: {removed}")
    selected /= totals
    return selected, [states[index] for index in keep]


def cosine_matrix(first: np.ndarray, second: np.ndarray) -> np.ndarray:
    first = first / np.maximum(np.linalg.norm(first, axis=1, keepdims=True), 1e-12)
    second = second / np.maximum(np.linalg.norm(second, axis=1, keepdims=True), 1e-12)
    return first @ second.T


def aligned_loadings(
    loadings: np.ndarray,
    source_states: list[str],
    common_states: list[str],
) -> np.ndarray:
    index = {state: column for column, state in enumerate(source_states)}
    aligned = loadings[:, [index[state] for state in common_states]].copy()
    aligned /= np.maximum(aligned.sum(axis=1, keepdims=True), 1e-12)
    return aligned


def patient_spatial_metrics(
    first: np.ndarray,
    second: np.ndarray,
    patients: np.ndarray,
) -> tuple[float, float, float]:
    correlations = []
    jaccards = []
    for patient in sorted(np.unique(patients)):
        mask = patients == patient
        rho = spearmanr(first[mask], second[mask]).statistic
        if np.isfinite(rho):
            correlations.append(float(rho))
        first_cut = np.quantile(first[mask], 0.90)
        second_cut = np.quantile(second[mask], 0.90)
        first_hot = mask & (first >= first_cut)
        second_hot = mask & (second >= second_cut)
        union = np.count_nonzero(first_hot | second_hot)
        if union:
            jaccards.append(
                np.count_nonzero(first_hot & second_hot) / union
            )
    global_rho = float(spearmanr(first, second).statistic)
    return (
        global_rho,
        float(np.median(correlations)),
        float(np.median(jaccards)),
    )


def factor_label(vector: np.ndarray, states: list[str]) -> str:
    dominant = states[int(np.argmax(vector))]
    labels = {
        "Lymphatic endothelial": "LEC-associated",
        "Blood endothelial": "Blood-endothelial",
        "AT1": "AT1-associated",
        "Perivascular/mesothelial stromal": "Perivascular-stromal",
        "B": "B-enriched/mixed-immune",
        "Ciliated epithelial": "Ciliated-associated",
        "Fibroblast": "Fibroblast-associated",
    }
    return labels.get(dominant, f"{dominant}-associated")


def load_primary(dataset) -> dict[str, object]:
    model_path = PRIMARY_RUN / "models" / (
        f"neighborhood50_{dataset.chip}_K07_MODEL_AND_USAGE.npz"
    )
    loading_path = PRIMARY_RUN / "models" / (
        f"neighborhood50_{dataset.chip}_K07_LOADINGS.tsv"
    )
    with np.load(model_path, allow_pickle=True) as archive:
        usage = archive["relative_usage"].astype(np.float64)
        entity = archive["entity_key"].astype(str)
        dominant = archive["dominant_factor"].astype(np.int64)
        components = archive["components"].astype(np.float64)
    loading_frame = pd.read_csv(loading_path, sep="\t", index_col=0)
    if loading_frame.columns.tolist() != dataset.state_names:
        raise RuntimeError(f"{dataset.chip}: primary loading state order differs")
    if not np.array_equal(
        entity, dataset.metadata["entity_key"].astype(str).to_numpy()
    ):
        raise RuntimeError(f"{dataset.chip}: primary entity order differs")
    return {
        "states": dataset.state_names,
        "loadings": loading_frame.to_numpy(float),
        "usage": usage,
        "dominant": dominant,
        "components": components,
        "model_path": model_path,
        "loading_path": loading_path,
    }


def fit_scenario(base, dataset, scenario, removed, args, scenario_index):
    matrix, states = subset_matrix(dataset.matrix, dataset.state_names, removed)
    patients = dataset.metadata["patient_id"].astype(str).to_numpy()
    reference = base.balanced_indices(
        patients, base.ROUTES.index(ROUTE), base.CHIPS.index(dataset.chip)
    )
    seed_base = (
        base.SEED
        + 80_000_000
        + scenario_index * 10_000_000
        + base.CHIPS.index(dataset.chip) * 1_000_000
        + RANK * 10_000
    )
    fits = Parallel(n_jobs=args.n_jobs, prefer="threads")(
        delayed(base.fit_nmf)(
            matrix[reference], RANK, seed_base + restart, args.max_iter
        )
        for restart in range(args.main_restarts)
    )
    medoid_index, restart = base.select_medoid(fits)
    components = fits[medoid_index].components
    usage, loadings, projection_iter = base.project_fixed(
        matrix, components, args.max_iter
    )
    lopo = base.lopo_audit(
        matrix,
        patients,
        reference,
        RANK,
        components,
        args.lopo_restarts,
        args.max_iter,
        seed_base + 500_000,
        args.n_jobs,
    )
    edges = base.spatial_edges(dataset.metadata)
    qc = base.factor_qc(usage, loadings, dataset.metadata, edges, states)
    error = float(np.linalg.norm(matrix - usage @ loadings) / np.linalg.norm(matrix))
    audit = {
        "scenario": scenario,
        "chip": dataset.chip,
        "rank": RANK,
        "states_n": len(states),
        "removed_states": ";".join(removed),
        "input_n": len(matrix),
        "balanced_reference_n": len(reference),
        "medoid_seed": fits[medoid_index].seed,
        "projection_n_iter": projection_iter,
        "normalized_reconstruction_error": error,
        "restart_median_factor_cosine": float(
            restart["mean_factor_cosine_to_medoid"].median()
        ),
        "restart_min_factor_cosine": float(
            restart["min_factor_cosine_to_medoid"].min()
        ),
        "restart_fraction_all_factors_ge_0_80": float(
            (restart["min_factor_cosine_to_medoid"] >= 0.80).mean()
        ),
        "main_convergence_fraction": float(
            (restart["n_iter"] < args.max_iter).mean()
        ),
        "lopo_median_factor_cosine": float(
            lopo["mean_matched_factor_cosine"].median()
        ),
        "lopo_min_factor_cosine": float(
            lopo["min_matched_factor_cosine"].min()
        ),
        "lopo_fraction_all_factors_ge_0_80": float(
            (lopo["min_matched_factor_cosine"] >= 0.80).mean()
        ),
        "minimum_high_usage_patient_n": int(qc["high_usage_patient_n"].min()),
        "maximum_patient_contribution": float(qc["max_patient_contribution"].max()),
        "maximum_abs_usage_log_umi_spearman": float(
            qc["usage_log_umi_spearman"].abs().max()
        ),
        "median_spatial_edge_usage_spearman": float(
            qc["spatial_edge_usage_spearman"].median()
        ),
        "dominant_factor_edge_agreement": float(
            np.mean(usage.argmax(1)[edges[:, 0]] == usage.argmax(1)[edges[:, 1]])
        ),
    }
    restart.insert(0, "chip", dataset.chip)
    restart.insert(0, "scenario", scenario)
    lopo.insert(0, "chip", dataset.chip)
    lopo.insert(0, "scenario", scenario)
    qc.insert(0, "chip", dataset.chip)
    qc.insert(0, "scenario", scenario)
    return {
        "states": states,
        "loadings": loadings,
        "usage": usage,
        "dominant": usage.argmax(1) + 1,
        "components": components,
        "restart": restart,
        "lopo": lopo,
        "qc": qc,
        "audit": audit,
    }


def primary_audit_row(chip: str) -> dict[str, object]:
    audit = pd.read_csv(PRIMARY_RUN / "03_RANK_AUDIT.tsv", sep="\t")
    row = audit.loc[audit["chip"].eq(chip) & audit["rank"].eq(RANK)].iloc[0]
    restart = pd.read_csv(PRIMARY_RUN / "04_RESTART_STABILITY.tsv.gz", sep="\t")
    restart = restart.loc[restart["chip"].eq(chip) & restart["rank"].eq(RANK)]
    lopo = pd.read_csv(PRIMARY_RUN / "05_LOPO_PATIENT.tsv.gz", sep="\t")
    lopo = lopo.loc[lopo["chip"].eq(chip) & lopo["rank"].eq(RANK)]
    return {
        "scenario": "main18",
        "chip": chip,
        "rank": RANK,
        "states_n": 18,
        "removed_states": "",
        "input_n": int(row["input_n"]),
        "balanced_reference_n": int(row["balanced_reference_n"]),
        "medoid_seed": int(row["medoid_seed"]),
        "projection_n_iter": int(row["projection_n_iter"]),
        "normalized_reconstruction_error": float(
            row["normalized_reconstruction_error"]
        ),
        "restart_median_factor_cosine": float(
            row["restart_median_factor_cosine"]
        ),
        "restart_min_factor_cosine": float(row["restart_min_factor_cosine"]),
        "restart_fraction_all_factors_ge_0_80": float(
            (restart["min_factor_cosine_to_medoid"] >= 0.80).mean()
        ),
        "main_convergence_fraction": float(row["main_convergence_fraction"]),
        "lopo_median_factor_cosine": float(row["lopo_median_factor_cosine"]),
        "lopo_min_factor_cosine": float(row["lopo_min_factor_cosine"]),
        "lopo_fraction_all_factors_ge_0_80": float(
            (lopo["min_matched_factor_cosine"] >= 0.80).mean()
        ),
        "minimum_high_usage_patient_n": int(row["minimum_high_usage_patient_n"]),
        "maximum_patient_contribution": float(row["maximum_patient_contribution"]),
        "maximum_abs_usage_log_umi_spearman": float(
            row["maximum_abs_usage_log_umi_spearman"]
        ),
        "median_spatial_edge_usage_spearman": float(
            row["median_spatial_edge_usage_spearman"]
        ),
        "dominant_factor_edge_agreement": float(
            row["dominant_factor_edge_agreement"]
        ),
    }


def match_primary_to_sensitivity(
    chip: str,
    primary: dict[str, object],
    sensitivity: dict[str, object],
    scenario: str,
    patients: np.ndarray,
) -> pd.DataFrame:
    common = list(sensitivity["states"])
    primary_aligned = aligned_loadings(
        primary["loadings"], primary["states"], common
    )
    sensitivity_aligned = aligned_loadings(
        sensitivity["loadings"], sensitivity["states"], common
    )
    similarity = cosine_matrix(primary_aligned, sensitivity_aligned)
    rows, columns = linear_sum_assignment(-similarity)
    output = []
    for left, right in zip(rows, columns):
        global_rho, median_patient_rho, median_jaccard = patient_spatial_metrics(
            primary["usage"][:, left], sensitivity["usage"][:, right], patients
        )
        order = np.argsort(sensitivity_aligned[right])[::-1][:5]
        output.append(
            {
                "chip": chip,
                "scenario": scenario,
                "primary_factor": f"F{left + 1:02d}",
                "primary_family": factor_label(
                    primary["loadings"][left], primary["states"]
                ),
                "sensitivity_factor": f"F{right + 1:02d}",
                "common_state_loading_cosine": float(similarity[left, right]),
                "same_bin_global_usage_spearman": global_rho,
                "median_patient_usage_spearman": median_patient_rho,
                "median_patient_q90_hotspot_jaccard": median_jaccard,
                "sensitivity_top_states": ";".join(common[index] for index in order),
                "sensitivity_top_loadings": ";".join(
                    f"{sensitivity_aligned[right, index]:.6f}" for index in order
                ),
            }
        )
    return pd.DataFrame(output)


def cross_chip_family_map(
    models: dict[tuple[str, str], dict[str, object]], scenario: str
) -> pd.DataFrame:
    states = list(models[(scenario, "K8")]["states"])
    anchor = models[(scenario, "K8")]["loadings"]
    matches = {"K8": np.arange(RANK, dtype=int)}
    cosines = {"K8": np.ones(RANK)}
    for chip in ("J2", "L3"):
        current = models[(scenario, chip)]["loadings"]
        similarity = cosine_matrix(anchor, current)
        rows, columns = linear_sum_assignment(-similarity)
        mapping = np.full(RANK, -1, dtype=int)
        scores = np.full(RANK, np.nan)
        mapping[rows] = columns
        scores[rows] = similarity[rows, columns]
        matches[chip] = mapping
        cosines[chip] = scores
    output = []
    for anchor_index in range(RANK):
        vectors = [
            models[(scenario, chip)]["loadings"][matches[chip][anchor_index]]
            for chip in CHIPS
        ]
        centroid = np.mean(vectors, axis=0)
        order = np.argsort(centroid)[::-1][:5]
        output.append(
            {
                "scenario": scenario,
                "family_id": f"CF{anchor_index + 1:02d}",
                "family_label": factor_label(anchor[anchor_index], states),
                "K8_factor": f"F{anchor_index + 1:02d}",
                "J2_factor": f"F{matches['J2'][anchor_index] + 1:02d}",
                "L3_factor": f"F{matches['L3'][anchor_index] + 1:02d}",
                "K8_J2_cosine": float(cosines["J2"][anchor_index]),
                "K8_L3_cosine": float(cosines["L3"][anchor_index]),
                "J2_L3_cosine": float(
                    cosine_matrix(vectors[1][None, :], vectors[2][None, :])[0, 0]
                ),
                "centroid_top_states": ";".join(states[index] for index in order),
                "centroid_top_loadings": ";".join(
                    f"{centroid[index]:.6f}" for index in order
                ),
            }
        )
    return pd.DataFrame(output)


def b_factor_trajectory(
    matches: pd.DataFrame,
    models: dict[tuple[str, str], dict[str, object]],
) -> pd.DataFrame:
    rows = []
    for chip in CHIPS:
        primary = models[("main18", chip)]
        b_index = primary["states"].index("B")
        factor_index = int(np.argmax(primary["loadings"][:, b_index]))
        factor = f"F{factor_index + 1:02d}"
        row = matches.loc[
            matches["chip"].eq(chip)
            & matches["scenario"].eq("no_b17")
            & matches["primary_factor"].eq(factor)
        ].iloc[0]
        sensitivity = models[("no_b17", chip)]
        sensitivity_index = int(str(row["sensitivity_factor"])[1:]) - 1
        vector = sensitivity["loadings"][sensitivity_index]
        state_index = {state: index for index, state in enumerate(sensitivity["states"])}
        immune_states = [
            "Alveolar/resident macrophage",
            "Inflammatory monocyte/macrophage",
            "Dendritic cell",
            "T",
            "NK",
            "Mast/Basophil",
            "Plasma",
        ]
        rows.append(
            {
                **row.to_dict(),
                "primary_B_loading": float(
                    primary["loadings"][factor_index, b_index]
                ),
                "matched_factor_T_loading": float(vector[state_index["T"]]),
                "matched_factor_NK_loading": float(vector[state_index["NK"]]),
                "matched_factor_DC_loading": float(
                    vector[state_index["Dendritic cell"]]
                ),
                "matched_factor_Mast_loading": float(
                    vector[state_index["Mast/Basophil"]]
                ),
                "matched_factor_Plasma_loading": float(
                    vector[state_index["Plasma"]]
                ),
                "matched_factor_total_nonB_immune_loading": float(
                    sum(vector[state_index[state]] for state in immune_states)
                ),
            }
        )
    return pd.DataFrame(rows)


def b_factor_all_candidates(
    models: dict[tuple[str, str], dict[str, object]],
    metadata: dict[str, pd.DataFrame],
) -> pd.DataFrame:
    """Score every B-removed factor against the primary B-enriched factor."""
    rows = []
    for chip in CHIPS:
        primary = models[("main18", chip)]
        sensitivity = models[("no_b17", chip)]
        b_index = primary["states"].index("B")
        primary_index = int(np.argmax(primary["loadings"][:, b_index]))
        common = list(sensitivity["states"])
        primary_vector = aligned_loadings(
            primary["loadings"][[primary_index]], primary["states"], common
        )[0]
        sensitivity_loadings = aligned_loadings(
            sensitivity["loadings"], sensitivity["states"], common
        )
        loading_cosines = cosine_matrix(
            primary_vector[None, :], sensitivity_loadings
        )[0]
        patients = metadata[chip]["patient_id"].astype(str).to_numpy()
        state_index = {state: index for index, state in enumerate(common)}
        immune_states = [
            "Alveolar/resident macrophage",
            "Inflammatory monocyte/macrophage",
            "Dendritic cell",
            "T",
            "NK",
            "Mast/Basophil",
            "Plasma",
        ]
        for candidate_index in range(RANK):
            global_rho, patient_rho, jaccard = patient_spatial_metrics(
                primary["usage"][:, primary_index],
                sensitivity["usage"][:, candidate_index],
                patients,
            )
            vector = sensitivity_loadings[candidate_index]
            order = np.argsort(vector)[::-1][:5]
            rows.append(
                {
                    "chip": chip,
                    "primary_B_factor": f"F{primary_index + 1:02d}",
                    "candidate_no_B_factor": f"F{candidate_index + 1:02d}",
                    "common_state_loading_cosine": float(
                        loading_cosines[candidate_index]
                    ),
                    "same_bin_global_usage_spearman": global_rho,
                    "median_patient_usage_spearman": patient_rho,
                    "median_patient_q90_hotspot_jaccard": jaccard,
                    "T_loading": float(vector[state_index["T"]]),
                    "NK_loading": float(vector[state_index["NK"]]),
                    "DC_loading": float(vector[state_index["Dendritic cell"]]),
                    "Mast_loading": float(vector[state_index["Mast/Basophil"]]),
                    "Plasma_loading": float(vector[state_index["Plasma"]]),
                    "total_nonB_immune_loading": float(
                        sum(vector[state_index[state]] for state in immune_states)
                    ),
                    "top_states": ";".join(common[index] for index in order),
                }
            )
    candidates = pd.DataFrame(rows)
    candidates["best_by_loading_cosine"] = candidates.groupby("chip")[
        "common_state_loading_cosine"
    ].transform("max").eq(candidates["common_state_loading_cosine"])
    candidates["best_by_patient_usage_spearman"] = candidates.groupby("chip")[
        "median_patient_usage_spearman"
    ].transform("max").eq(candidates["median_patient_usage_spearman"])
    return candidates


def make_figures(
    audit: pd.DataFrame,
    matches: pd.DataFrame,
    cross_chip: pd.DataFrame,
    b_trajectory: pd.DataFrame,
    output: Path,
) -> None:
    metrics = [
        "restart_fraction_all_factors_ge_0_80",
        "lopo_fraction_all_factors_ge_0_80",
        "main_convergence_fraction",
    ]
    long = audit.melt(
        id_vars=["scenario", "chip"], value_vars=metrics,
        var_name="metric", value_name="value",
    )
    fig, axes = plt.subplots(1, 3, figsize=(15, 4.6), constrained_layout=True)
    for axis, chip in zip(axes, CHIPS):
        sns.barplot(
            data=long.loc[long["chip"].eq(chip)],
            x="scenario", y="value", hue="metric", ax=axis,
        )
        axis.axhline(0.80, color="black", linestyle=":", linewidth=1)
        axis.set_ylim(0, 1.05)
        axis.set_title(chip)
        axis.tick_params(axis="x", rotation=20)
        if chip != "L3":
            axis.legend_.remove()
    axes[-1].legend(fontsize=7, loc="lower right")
    fig.suptitle("K=7 numerical robustness under Reference18 state removal")
    fig.savefig(output / "figures/FIG01_NUMERICAL_ROBUSTNESS.png", dpi=220)
    plt.close(fig)

    display = matches.copy()
    display["row"] = display["chip"] + " | " + display["scenario"]
    display["column"] = display["primary_family"]
    for metric, filename, title, vmin in (
        ("common_state_loading_cosine", "FIG02_LOADING_PERSISTENCE.png", "Matched loading cosine", 0),
        ("median_patient_usage_spearman", "FIG03_SPATIAL_SCORE_PERSISTENCE.png", "Median within-patient Bin50 usage Spearman", -1),
        ("median_patient_q90_hotspot_jaccard", "FIG04_HOTSPOT_PERSISTENCE.png", "Median patient Q90 hotspot Jaccard", 0),
    ):
        matrix = display.pivot(index="row", columns="column", values=metric)
        ordered_columns = [column for column in FAMILY_ORDER if column in matrix]
        matrix = matrix.reindex(columns=ordered_columns)
        fig, axis = plt.subplots(figsize=(13, 5.5), constrained_layout=True)
        sns.heatmap(
            matrix, vmin=vmin, vmax=1, cmap="viridis", annot=True,
            fmt=".2f", ax=axis,
        )
        axis.set_title(title)
        axis.set_xlabel("")
        axis.set_ylabel("")
        axis.tick_params(axis="x", rotation=35)
        fig.savefig(output / f"figures/{filename}", dpi=220)
        plt.close(fig)

    cross = cross_chip.copy()
    cross["minimum_pairwise_cosine"] = cross[
        ["K8_J2_cosine", "K8_L3_cosine", "J2_L3_cosine"]
    ].min(axis=1)
    matrix = cross.pivot(
        index="scenario", columns="family_id", values="minimum_pairwise_cosine"
    )
    matrix = matrix.reindex(columns=[f"CF{index + 1:02d}" for index in range(RANK)])
    fig, axis = plt.subplots(figsize=(13, 3.8), constrained_layout=True)
    sns.heatmap(matrix, vmin=0, vmax=1, cmap="viridis", annot=True, fmt=".2f", ax=axis)
    axis.set_title("Minimum pairwise cross-chip loading cosine")
    axis.set_xlabel("")
    axis.set_ylabel("")
    axis.tick_params(axis="x", rotation=35)
    fig.savefig(output / "figures/FIG05_CROSS_CHIP_FAMILY_ROBUSTNESS.png", dpi=220)
    plt.close(fig)

    fig, axes = plt.subplots(1, 3, figsize=(14, 4.3), constrained_layout=True)
    for axis, (_, row) in zip(axes, b_trajectory.iterrows()):
        values = [
            row["matched_factor_T_loading"],
            row["matched_factor_NK_loading"],
            row["matched_factor_DC_loading"],
            row["matched_factor_Mast_loading"],
            row["matched_factor_Plasma_loading"],
        ]
        axis.bar(["T", "NK", "DC", "Mast", "Plasma"], values)
        axis.set_ylim(0, max(0.5, max(values) * 1.15))
        axis.set_title(
            f"{row['chip']}: {row['primary_factor']} -> {row['sensitivity_factor']}\n"
            f"cos={row['common_state_loading_cosine']:.2f}, "
            f"rho={row['median_patient_usage_spearman']:.2f}"
        )
        axis.tick_params(axis="x", rotation=30)
    fig.suptitle("Where the primary B-enriched factor maps after removing B")
    fig.savefig(output / "figures/FIG06_B_FACTOR_TRAJECTORY.png", dpi=220)
    plt.close(fig)


def main() -> None:
    args = parse_args()
    output = args.output_dir.resolve()
    output.mkdir(parents=True, exist_ok=False)
    (output / "models").mkdir()
    (output / "figures").mkdir()
    started = time.time()
    base = load_base()

    datasets, input_audit, input_paths = base.load_datasets()
    datasets = [dataset for dataset in datasets if dataset.route == ROUTE]
    input_audit = input_audit.loc[input_audit["route"].eq(ROUTE)].copy()
    models: dict[tuple[str, str], dict[str, object]] = {}
    audit_rows = []
    restart_tables = []
    lopo_tables = []
    qc_tables = []

    for dataset in datasets:
        primary = load_primary(dataset)
        models[("main18", dataset.chip)] = primary
        audit_rows.append(primary_audit_row(dataset.chip))
        for scenario_index, scenario in enumerate(
            ("no_basal_lec16", "no_b17"), start=1
        ):
            print(
                f"[{time.strftime('%H:%M:%S')}] {dataset.chip}/{scenario}: fitting",
                flush=True,
            )
            fitted = fit_scenario(
                base,
                dataset,
                scenario,
                SCENARIOS[scenario],
                args,
                scenario_index,
            )
            models[(scenario, dataset.chip)] = fitted
            audit_rows.append(fitted["audit"])
            restart_tables.append(fitted["restart"])
            lopo_tables.append(fitted["lopo"])
            qc_tables.append(fitted["qc"])
            stem = f"{scenario}_{dataset.chip}_K07"
            atomic_npz(
                output / "models" / f"{stem}_MODEL_AND_USAGE.npz",
                {
                    "components": fitted["components"].astype(np.float32),
                    "normalized_loadings": fitted["loadings"].astype(np.float32),
                    "relative_usage": fitted["usage"].astype(np.float32),
                    "entity_key": dataset.metadata["entity_key"].astype(str).to_numpy(),
                    "dominant_factor": fitted["dominant"].astype(np.int8),
                    "state_names": np.asarray(fitted["states"], dtype=str),
                },
            )
            atomic_table(
                pd.DataFrame(
                    fitted["loadings"],
                    index=[f"F{index + 1:02d}" for index in range(RANK)],
                    columns=fitted["states"],
                ),
                output / "models" / f"{stem}_LOADINGS.tsv",
                sep="\t",
                index_label="factor",
            )
            print(
                f"[{time.strftime('%H:%M:%S')}] {dataset.chip}/{scenario}: complete",
                flush=True,
            )

    comparisons = []
    metadata = {dataset.chip: dataset.metadata for dataset in datasets}
    for chip in CHIPS:
        patients = metadata[chip]["patient_id"].astype(str).to_numpy()
        for scenario in ("no_basal_lec16", "no_b17"):
            comparisons.append(
                match_primary_to_sensitivity(
                    chip,
                    models[("main18", chip)],
                    models[(scenario, chip)],
                    scenario,
                    patients,
                )
            )
    matches = pd.concat(comparisons, ignore_index=True)
    cross_chip = pd.concat(
        [cross_chip_family_map(models, scenario) for scenario in SCENARIOS],
        ignore_index=True,
    )
    b_trajectory = b_factor_trajectory(matches, models)
    b_candidates = b_factor_all_candidates(models, metadata)
    audit = pd.DataFrame(audit_rows)

    method = pd.DataFrame(
        [
            ("primary model", "persisted per-chip 18-state Neighborhood50 NMF K=7"),
            ("harmonized rank", "K=7 for K8, J2, and L3"),
            ("sensitivity A", "remove Basal-like epithelial and Lymphatic endothelial; renormalize rows; refit K=7"),
            ("sensitivity B", "remove B; renormalize rows; refit K=7"),
            ("fit reference", f"{base.PER_PATIENT} rows per patient; 7 patients per chip"),
            ("main restarts", str(args.main_restarts)),
            ("LOPO restarts", str(args.lopo_restarts)),
            ("NMF objective", "unpenalized Frobenius; random init; coordinate descent"),
            ("factor matching", "Hungarian assignment on common-state loading cosine"),
            ("spatial persistence", "same-Bin50 usage Spearman and within-patient Q90 hotspot Jaccard"),
            ("disease labels used", "no"),
        ],
        columns=["parameter", "locked_value"],
    )
    manifest_paths = list(input_paths)
    for chip in CHIPS:
        manifest_paths.extend(
            [
                models[("main18", chip)]["model_path"],
                models[("main18", chip)]["loading_path"],
            ]
        )
    atomic_table(method, output / "00_METHOD_LOCK.tsv", sep="\t", index=False)
    atomic_table(
        pd.DataFrame(
            [
                {"path": str(path.resolve()), "sha256": sha256(path)}
                for path in manifest_paths
            ]
        ),
        output / "00_INPUT_MANIFEST.tsv",
        sep="\t",
        index=False,
    )
    atomic_table(input_audit, output / "01_INPUT_AUDIT.tsv", sep="\t", index=False)
    atomic_table(audit, output / "02_MODEL_ROBUSTNESS.tsv", sep="\t", index=False)
    atomic_table(
        pd.concat(restart_tables, ignore_index=True),
        output / "03_SENSITIVITY_RESTARTS.tsv.gz",
        sep="\t", index=False, compression="gzip",
    )
    atomic_table(
        pd.concat(lopo_tables, ignore_index=True),
        output / "04_SENSITIVITY_LOPO.tsv.gz",
        sep="\t", index=False, compression="gzip",
    )
    atomic_table(
        pd.concat(qc_tables, ignore_index=True),
        output / "05_SENSITIVITY_FACTOR_QC.tsv",
        sep="\t", index=False,
    )
    atomic_table(matches, output / "06_PRIMARY_TO_SENSITIVITY_MATCHES.tsv", sep="\t", index=False)
    atomic_table(cross_chip, output / "07_CROSS_CHIP_FAMILY_MAP.tsv", sep="\t", index=False)
    atomic_table(b_trajectory, output / "08_B_FACTOR_TRAJECTORY.tsv", sep="\t", index=False)
    atomic_table(b_candidates, output / "09_B_FACTOR_ALL_CANDIDATES.tsv", sep="\t", index=False)
    make_figures(audit, matches, cross_chip, b_trajectory, output)

    checks = pd.DataFrame(
        [
            ("three_chips_loaded", len(datasets) == 3, str(len(datasets))),
            ("nine_models_available", len(models) == 9, str(len(models))),
            ("six_sensitivity_models_fitted", len(audit) == 9, str(len(audit) - 3)),
            ("all_model_rows_sum_to_one", all(np.allclose(item["loadings"].sum(1), 1, atol=1e-6) for item in models.values()), "checked"),
            ("all_usage_rows_sum_to_one", all(np.allclose(item["usage"].sum(1), 1, atol=1e-5) for item in models.values()), "checked"),
            ("all_primary_factors_matched", len(matches) == 42, str(len(matches))),
            ("three_B_trajectories", len(b_trajectory) == 3, str(len(b_trajectory))),
            ("twenty_one_B_candidate_matches", len(b_candidates) == 21, str(len(b_candidates))),
            ("all_cross_chip_families_present", len(cross_chip) == 21, str(len(cross_chip))),
            ("no_nonfinite_match_metric", bool(np.isfinite(matches[["common_state_loading_cosine", "same_bin_global_usage_spearman", "median_patient_usage_spearman", "median_patient_q90_hotspot_jaccard"]].to_numpy()).all()), "checked"),
            ("disease_blind", True, "disease labels not supplied to NMF or matching"),
        ],
        columns=["check", "pass", "detail"],
    )
    atomic_table(checks, output / "10_VALIDATION.tsv", sep="\t", index=False)
    if not checks["pass"].all():
        raise RuntimeError(
            "Validation failed: "
            + ";".join(checks.loc[~checks["pass"], "check"])
        )

    status = {
        "status": "ANALYSIS_COMPLETE_PENDING_EVIDENCE_REVIEW",
        "primary": "18 states; common K=7",
        "sensitivity_models": list(SCENARIOS)[1:],
        "chips": list(CHIPS),
        "models_available": len(models),
        "elapsed_seconds": time.time() - started,
    }
    (output / "FINAL_STATUS.json").write_text(
        json.dumps(status, indent=2), encoding="utf-8"
    )
    (output / "COMPLETE.ok").write_text("PASS\n", encoding="utf-8")
    (output / "RUN_LOG.txt").write_text(
        f"Script: scripts/{Path(__file__).name}\n"
        f"Command: {sys.executable} scripts/{Path(__file__).name} "
        f"--output-dir {output} --n-jobs {args.n_jobs} "
        f"--main-restarts {args.main_restarts} "
        f"--lopo-restarts {args.lopo_restarts} --max-iter {args.max_iter}\n"
        f"Python: {sys.version.replace(chr(10), ' ')}\n"
        f"Platform: {platform.platform()}\n"
        "Disease used: no\nStatus: PASS\n",
        encoding="utf-8",
    )
    print(json.dumps(status, indent=2), flush=True)


if __name__ == "__main__":
    main()
