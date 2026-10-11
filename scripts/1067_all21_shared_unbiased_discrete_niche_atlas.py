#!/usr/bin/env python3
"""Build and test a shared All21 unbiased discrete niche atlas.

The workflow consumes the three frozen Reference18 UMI>=50 chip inputs, builds
one common center-excluded 50-um neighborhood-composition space, selects among
K=3/4/5 without using disease labels, freezes a shared niche atlas, and only
then performs patient-level disease comparisons. Historical chip-specific
labels and aggregate/B-lineage labels are never read.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import platform
import sys
from itertools import combinations
from pathlib import Path

import h5py
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from scipy.spatial import cKDTree
from scipy.stats import t as student_t
from sklearn.cluster import KMeans
from sklearn.metrics import adjusted_rand_score, silhouette_score
from sklearn.metrics.cluster import contingency_matrix
from sklearn.metrics.pairwise import pairwise_distances_argmin_min


SEED = 20261011
K_VALUES = (3, 4, 5)
REPEATS = 20
MAX_REFERENCE_PER_PATIENT = 500
RADIUS_UM = 50.0
ADJACENCY_RADIUS_UM = 36.0
PERMUTATIONS = 10_000
DISEASES = ("HC", "IPF", "SSc-ILD")
CONTRASTS = ("IPF_vs_HC", "SSc-ILD_vs_HC", "SSc-ILD_vs_IPF")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input-root", required=True, type=Path)
    parser.add_argument("--output-dir", required=True, type=Path)
    parser.add_argument(
        "--feature-transform",
        choices=("raw_clr", "chip_centered_clr"),
        default="raw_clr",
    )
    parser.add_argument("--forced-k", type=int, choices=K_VALUES)
    return parser.parse_args()


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def bh_adjust(p_values: np.ndarray) -> np.ndarray:
    values = np.asarray(p_values, dtype=float)
    output = np.full(values.shape, np.nan)
    valid = np.isfinite(values)
    p = values[valid]
    if not len(p):
        return output
    order = np.argsort(p)
    ranked = p[order]
    adjusted = np.minimum.accumulate(
        (ranked * len(ranked) / np.arange(1, len(ranked) + 1))[::-1]
    )[::-1]
    restored = np.empty_like(adjusted)
    restored[order] = np.minimum(adjusted, 1.0)
    output[valid] = restored
    return output


def voltron_clr(composition: np.ndarray) -> np.ndarray:
    x = np.asarray(composition, dtype=np.float64)
    x = x / x.sum(axis=1, keepdims=True)
    positive_log_sum = np.where(x > 0, np.log1p(x), 0.0).sum(axis=1)
    denominator = np.exp(positive_log_sum / x.shape[1])
    output = np.log1p(x / denominator[:, None])
    if not np.isfinite(output).all():
        raise RuntimeError("Non-finite CLR values")
    return output.astype(np.float32)


def read_inputs(root: Path) -> tuple[pd.DataFrame, list[str], list[Path]]:
    paths = [root / chip / "01_FROZEN_REFERENCE18_INPUT.tsv.gz" for chip in ("K8", "J2", "L3")]
    frames = [pd.read_csv(path, sep="\t") for path in paths]
    columns = [column for column in frames[0].columns if column.startswith("weight__")]
    if len(columns) != 18:
        raise RuntimeError("Expected 18 Reference18 weight columns")
    for frame in frames[1:]:
        if [column for column in frame.columns if column.startswith("weight__")] != columns:
            raise RuntimeError("Reference18 state order differs across chips")
    data = pd.concat(frames, ignore_index=True)
    return data, [column.removeprefix("weight__") for column in columns], paths


def build_neighborhood_composition(
    data: pd.DataFrame, weight_columns: list[str]
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    weights = data[weight_columns].to_numpy(dtype=np.float64)
    composition = np.full(weights.shape, np.nan, dtype=np.float32)
    neighbor_n = np.zeros(len(data), dtype=np.int32)
    spatial_edges = []
    for patient_id, frame in data.groupby("patient_id", sort=True):
        global_indices = frame.index.to_numpy()
        coordinates = frame[["x_um", "y_um"]].to_numpy(dtype=float)
        tree = cKDTree(coordinates)
        neighborhoods = tree.query_ball_point(coordinates, r=RADIUS_UM + 1e-8, workers=-1)
        for local_index, neighbors in enumerate(neighborhoods):
            neighbor_array = np.asarray(neighbors, dtype=int)
            neighbor_array = neighbor_array[neighbor_array != local_index]
            if len(neighbor_array):
                global_index = global_indices[local_index]
                total = weights[global_indices[neighbor_array]].sum(axis=0)
                composition[global_index] = total / total.sum()
                neighbor_n[global_index] = len(neighbor_array)
        pairs = tree.query_pairs(r=ADJACENCY_RADIUS_UM + 1e-8, output_type="ndarray")
        if len(pairs):
            spatial_edges.append(global_indices[pairs])
    edges = np.vstack(spatial_edges).astype(np.int64)
    return composition, neighbor_n, edges


def balanced_indices(
    patient: np.ndarray,
    seed: int,
    eligible: np.ndarray,
    reference_per_patient: int,
) -> np.ndarray:
    rng = np.random.default_rng(seed)
    selected = []
    for value in sorted(np.unique(patient)):
        choices = np.flatnonzero((patient == value) & eligible)
        if len(choices) < reference_per_patient:
            raise RuntimeError(
                f"Patient {value} has fewer than {reference_per_patient} evaluable bins"
            )
        selected.append(rng.choice(choices, reference_per_patient, replace=False))
    return np.sort(np.concatenate(selected))


def cramer_v(labels: np.ndarray, groups: np.ndarray) -> float:
    table = pd.crosstab(pd.Series(labels), pd.Series(groups)).to_numpy(dtype=float)
    expected = table.sum(axis=1, keepdims=True) @ table.sum(axis=0, keepdims=True) / table.sum()
    chi2 = np.divide((table - expected) ** 2, expected, out=np.zeros_like(expected), where=expected > 0).sum()
    return float(np.sqrt((chi2 / table.sum()) / max(1, min(table.shape) - 1)))


def consensus_for_k(
    features: np.ndarray,
    patient: np.ndarray,
    chip: np.ndarray,
    eligible: np.ndarray,
    edges: np.ndarray,
    k: int,
    reference_per_patient: int,
) -> tuple[dict, np.ndarray, np.ndarray, pd.DataFrame]:
    eligible_indices = np.flatnonzero(eligible)
    features_eligible = features[eligible]
    labels_all = []
    centers_all = []
    fit_rows = []
    for repeat in range(REPEATS):
        ref = balanced_indices(
            patient,
            SEED + 1000 * k + repeat,
            eligible,
            reference_per_patient,
        )
        model = KMeans(
            n_clusters=k,
            random_state=SEED + 10_000 * k + repeat,
            n_init=25,
            max_iter=500,
            algorithm="lloyd",
        ).fit(features[ref])
        labels = model.predict(features_eligible)
        labels_all.append(labels)
        centers_all.append(model.cluster_centers_)
        fit_rows.append(
            {
                "k": k,
                "repeat": repeat + 1,
                "sample_seed": SEED + 1000 * k + repeat,
                "kmeans_seed": SEED + 10_000 * k + repeat,
                "inertia_balanced_fit": model.inertia_,
                "converged_iteration_n": model.n_iter_,
            }
        )
    pairwise = np.eye(REPEATS)
    for left, right in combinations(range(REPEATS), 2):
        value = adjusted_rand_score(labels_all[left], labels_all[right])
        pairwise[left, right] = value
        pairwise[right, left] = value
    mean_ari = (pairwise.sum(axis=1) - 1) / (REPEATS - 1)
    medoid = int(np.argmax(mean_ari))
    labels = labels_all[medoid]
    centers = centers_all[medoid]
    sizes = np.bincount(labels, minlength=k)
    patient_table = pd.crosstab(pd.Series(patient[eligible_indices]), pd.Series(labels))
    patient_fractions = patient_table.div(patient_table.sum(axis=1), axis=0)
    coverage = ((patient_table >= 20) & (patient_fractions >= 0.01)).sum(axis=0)
    chip_table = pd.crosstab(pd.Series(chip[eligible_indices]), pd.Series(labels))
    chip_coverage = (chip_table > 0).sum(axis=0)
    contribution = patient_table.div(patient_table.sum(axis=0), axis=1)
    edge_mask = eligible[edges[:, 0]] & eligible[edges[:, 1]]
    mapped = np.full(len(eligible), -1, dtype=int)
    mapped[eligible_indices] = labels
    spatial_agreement = np.mean(mapped[edges[edge_mask, 0]] == mapped[edges[edge_mask, 1]])
    silhouette_rng = np.random.default_rng(SEED + k)
    silhouette_idx = silhouette_rng.choice(
        len(features_eligible), min(5000, len(features_eligible)), replace=False
    )
    silhouette = silhouette_score(features_eligible[silhouette_idx], labels[silhouette_idx])
    total_ss = np.sum((features_eligible - features_eligible.mean(axis=0)) ** 2)
    _, distances = pairwise_distances_argmin_min(features_eligible, centers)
    within_ss = np.sum(distances**2)

    cross_chip_values = []
    for held_chip in sorted(np.unique(chip)):
        held = eligible & (chip == held_chip)
        train = eligible & (chip != held_chip)
        for repeat in range(5):
            rng = np.random.default_rng(SEED + 100_000 * k + 100 * repeat + len(held_chip))
            selected = []
            for value in sorted(np.unique(patient[train])):
                choices = np.flatnonzero(train & (patient == value))
                selected.append(rng.choice(choices, reference_per_patient, replace=False))
            ref = np.concatenate(selected)
            model = KMeans(
                n_clusters=k,
                random_state=SEED + 200_000 * k + repeat,
                n_init=25,
                max_iter=500,
            ).fit(features[ref])
            predicted = model.predict(features[held])
            cross_chip_values.append(
                {
                    "k": k,
                    "held_out_chip": held_chip,
                    "repeat": repeat + 1,
                    "ari_vs_global_labels": adjusted_rand_score(mapped[held], predicted),
                }
            )
    cross_chip = pd.DataFrame(cross_chip_values)
    metrics = {
        "k": k,
        "medoid_repeat": medoid + 1,
        "repeat_ari_mean": np.mean(pairwise[np.triu_indices(REPEATS, 1)]),
        "repeat_ari_median": np.median(pairwise[np.triu_indices(REPEATS, 1)]),
        "repeat_ari_min": np.min(pairwise[np.triu_indices(REPEATS, 1)]),
        "cross_chip_ari_mean": cross_chip["ari_vs_global_labels"].mean(),
        "cross_chip_ari_median": cross_chip["ari_vs_global_labels"].median(),
        "cross_chip_ari_min": cross_chip["ari_vs_global_labels"].min(),
        "silhouette_mean": silhouette,
        "variance_explained_fraction": 1 - within_ss / total_ss,
        "spatial_edge_agreement": spatial_agreement,
        "min_niche_n": sizes.min(),
        "min_niche_fraction": sizes.min() / sizes.sum(),
        "min_patient_coverage_n": coverage.min(),
        "min_chip_coverage_n": chip_coverage.min(),
        "max_single_patient_contribution": contribution.max().max(),
        "chip_cramer_v": cramer_v(labels, chip[eligible_indices]),
    }
    fit_audit = pd.DataFrame(fit_rows)
    fit_audit["is_medoid_repeat"] = fit_audit["repeat"].eq(medoid + 1)
    return metrics, mapped, centers, pd.concat([fit_audit, cross_chip], ignore_index=True, sort=False)


def zscore(series: pd.Series) -> pd.Series:
    if series.std(ddof=0) == 0:
        return pd.Series(np.zeros(len(series)), index=series.index)
    return (series - series.mean()) / series.std(ddof=0)


def relabel_by_composition(
    labels: np.ndarray, composition: np.ndarray, state_names: list[str], eligible: np.ndarray
) -> tuple[np.ndarray, pd.DataFrame]:
    old_values = sorted(np.unique(labels[eligible]))
    rows = []
    for old in old_values:
        mean = composition[(labels == old) & eligible].mean(axis=0)
        top = np.argsort(mean)[::-1]
        rows.append(
            {
                "old_label": old,
                "dominant_state_index": int(top[0]),
                "dominant_weight": float(mean[top[0]]),
                "second_state_index": int(top[1]),
                "top3_states": ";".join(np.asarray(state_names)[top[:3]]),
                **{f"mean_{state}": float(mean[i]) for i, state in enumerate(state_names)},
            }
        )
    summary = pd.DataFrame(rows).sort_values(
        ["dominant_state_index", "second_state_index", "dominant_weight"],
        ascending=[True, True, False],
    )
    mapping = {old: new for new, old in enumerate(summary["old_label"], start=1)}
    relabeled = np.full(len(labels), -1, dtype=int)
    for old, new in mapping.items():
        relabeled[labels == old] = new
    summary["niche_number"] = summary["old_label"].map(mapping)
    summary["niche"] = summary["niche_number"].map(lambda x: f"CN{x}")
    summary = summary.sort_values("niche_number").drop(
        columns=["old_label", "dominant_state_index", "second_state_index"]
    )
    return relabeled, summary


def design_matrix(chip: np.ndarray, disease: np.ndarray, depth: np.ndarray | None = None) -> np.ndarray:
    columns = [
        np.ones(len(chip)),
        (chip == "J2").astype(float),
        (chip == "L3").astype(float),
        (disease == "IPF").astype(float),
        (disease == "SSc-ILD").astype(float),
    ]
    if depth is not None:
        columns.append((depth - depth.mean()) / depth.std(ddof=1))
    return np.column_stack(columns)


def contrast_vectors(p: int) -> dict[str, np.ndarray]:
    output = {}
    for name, ipf, ssc in (
        ("IPF_vs_HC", 1.0, 0.0),
        ("SSc-ILD_vs_HC", 0.0, 1.0),
        ("SSc-ILD_vs_IPF", -1.0, 1.0),
    ):
        vector = np.zeros(p)
        vector[3] = ipf
        vector[4] = ssc
        output[name] = vector
    return output


def fit_hc3(y: np.ndarray, x: np.ndarray) -> dict[str, dict[str, np.ndarray]]:
    xtx_inv = np.linalg.inv(x.T @ x)
    beta = xtx_inv @ x.T @ y
    residual = y - x @ beta
    leverage = np.sum((x @ xtx_inv) * x, axis=1)
    residual_sq = (residual / (1 - leverage[:, None])) ** 2
    df = x.shape[0] - x.shape[1]
    output = {}
    for name, contrast in contrast_vectors(x.shape[1]).items():
        effect = contrast @ beta
        influence = x @ xtx_inv @ contrast
        variance = (influence[:, None] ** 2 * residual_sq).sum(axis=0)
        se = np.sqrt(np.maximum(variance, 0))
        statistic = np.divide(effect, se, out=np.zeros_like(effect), where=se > 0)
        output[name] = {
            "effect": effect,
            "se": se,
            "p": 2 * student_t.sf(np.abs(statistic), df=df),
        }
    return output


def permuted_contrast_p(
    y: np.ndarray, chip: np.ndarray, disease: np.ndarray, observed: dict[str, dict[str, np.ndarray]]
) -> dict[str, np.ndarray]:
    rng = np.random.default_rng(SEED + 999)
    exceed = {name: np.zeros(y.shape[1], dtype=int) for name in CONTRASTS}
    chip_indices = [np.flatnonzero(chip == value) for value in np.unique(chip)]
    for _ in range(PERMUTATIONS):
        permuted = disease.copy()
        for indices in chip_indices:
            permuted[indices] = rng.permutation(permuted[indices])
        x = design_matrix(chip, permuted)
        beta = np.linalg.inv(x.T @ x) @ x.T @ y
        for name, contrast in contrast_vectors(x.shape[1]).items():
            effect = contrast @ beta
            exceed[name] += np.abs(effect) >= np.abs(observed[name]["effect"]) - 1e-12
    return {name: (value + 1) / (PERMUTATIONS + 1) for name, value in exceed.items()}


def multivariate_f(y: np.ndarray, x_reduced: np.ndarray, x_full: np.ndarray) -> float:
    beta_reduced = np.linalg.lstsq(x_reduced, y, rcond=None)[0]
    beta_full = np.linalg.lstsq(x_full, y, rcond=None)[0]
    sse_reduced = np.sum((y - x_reduced @ beta_reduced) ** 2)
    sse_full = np.sum((y - x_full @ beta_full) ** 2)
    df1 = x_full.shape[1] - x_reduced.shape[1]
    df2 = len(y) - x_full.shape[1]
    return float(((sse_reduced - sse_full) / df1) / (sse_full / df2))


def global_composition_tests(patient_table: pd.DataFrame, niche_columns: list[str]) -> pd.DataFrame:
    rows = []
    comparisons = [("overall", list(DISEASES))] + [
        ("IPF_vs_HC", ["HC", "IPF"]),
        ("SSc-ILD_vs_HC", ["HC", "SSc-ILD"]),
        ("SSc-ILD_vs_IPF", ["IPF", "SSc-ILD"]),
    ]
    rng = np.random.default_rng(SEED + 555)
    for name, groups in comparisons:
        frame = patient_table[patient_table["disease"].isin(groups)].copy()
        fractions = frame[niche_columns].to_numpy(dtype=float)
        clr = np.log(fractions + 1e-5)
        clr -= clr.mean(axis=1, keepdims=True)
        chip = frame["chip_label"].to_numpy()
        disease = frame["disease"].to_numpy()
        x_reduced = np.column_stack(
            [np.ones(len(frame)), (chip == "J2").astype(float), (chip == "L3").astype(float)]
        )
        if name == "overall":
            x_full = design_matrix(chip, disease)
        else:
            target = groups[1]
            x_full = np.column_stack([x_reduced, (disease == target).astype(float)])
        observed = multivariate_f(clr, x_reduced, x_full)
        exceed = 0
        chip_indices = [np.flatnonzero(chip == value) for value in np.unique(chip)]
        for _ in range(PERMUTATIONS):
            permuted = disease.copy()
            for indices in chip_indices:
                permuted[indices] = rng.permutation(permuted[indices])
            if name == "overall":
                perm_full = design_matrix(chip, permuted)
            else:
                perm_full = np.column_stack([x_reduced, (permuted == groups[1]).astype(float)])
            exceed += multivariate_f(clr, x_reduced, perm_full) >= observed - 1e-12
        rows.append(
            {
                "comparison": name,
                "patient_n": len(frame),
                "pseudo_f": observed,
                "within_chip_permutation_p": (exceed + 1) / (PERMUTATIONS + 1),
            }
        )
    result = pd.DataFrame(rows)
    result["fdr"] = bh_adjust(result["within_chip_permutation_p"].to_numpy())
    return result


def main() -> None:
    args = parse_args()
    input_root = args.input_root.resolve()
    output = args.output_dir.resolve()
    output.mkdir(parents=True, exist_ok=False)
    figures = output / "figures"
    figures.mkdir()

    data, state_names, input_paths = read_inputs(input_root)
    weight_columns = [f"weight__{state}" for state in state_names]
    if len(data) != 126_449 or data["entity_key"].nunique() != len(data):
        raise RuntimeError("Unexpected All21 input rows or duplicate entity keys")
    patient_design = data.drop_duplicates("patient_id")[["patient_id", "chip_label", "disease"]]
    expected = pd.DataFrame(
        [[2, 2, 3]] * 3, index=["J2", "K8", "L3"], columns=list(DISEASES)
    )
    observed = pd.crosstab(patient_design["chip_label"], patient_design["disease"]).reindex(
        index=expected.index, columns=expected.columns, fill_value=0
    )
    if not observed.equals(expected):
        raise RuntimeError("Chip/disease design is not the expected balanced 2/2/3 structure")

    composition, neighbor_n, edges = build_neighborhood_composition(data, weight_columns)
    eligible = neighbor_n > 0
    evaluable_by_patient = (
        pd.Series(eligible, index=data.index)
        .groupby(data["patient_id"])
        .sum()
        .astype(int)
    )
    reference_per_patient = min(MAX_REFERENCE_PER_PATIENT, int(evaluable_by_patient.min()))
    features = np.zeros_like(composition, dtype=np.float32)
    features[eligible] = voltron_clr(composition[eligible])
    chip_order = sorted(data["chip_label"].unique())
    chip_offsets = np.zeros((len(chip_order), features.shape[1]), dtype=float)
    if args.feature_transform == "chip_centered_clr":
        patient_centroids = []
        for patient_id, frame in data[eligible].groupby("patient_id", sort=True):
            patient_centroids.append(
                {
                    "patient_id": patient_id,
                    "chip_label": frame["chip_label"].iloc[0],
                    "centroid": features[frame.index].mean(axis=0),
                }
            )
        grand_centroid = np.mean([row["centroid"] for row in patient_centroids], axis=0)
        for chip_index, chip_label in enumerate(chip_order):
            chip_centroid = np.mean(
                [
                    row["centroid"]
                    for row in patient_centroids
                    if row["chip_label"] == chip_label
                ],
                axis=0,
            )
            chip_offsets[chip_index] = chip_centroid - grand_centroid
            chip_mask = eligible & data["chip_label"].eq(chip_label).to_numpy()
            features[chip_mask] -= chip_offsets[chip_index]

    all_metrics = []
    labels_by_k = {}
    centers_by_k = {}
    audits = []
    for k in K_VALUES:
        metrics, labels, centers, audit = consensus_for_k(
            features,
            data["patient_id"].to_numpy(),
            data["chip_label"].to_numpy(),
            eligible,
            edges,
            k,
            reference_per_patient,
        )
        all_metrics.append(metrics)
        labels_by_k[k] = labels
        centers_by_k[k] = centers
        audits.append(audit)
    k_table = pd.DataFrame(all_metrics)
    k_table["selection_eligible"] = (
        (k_table["repeat_ari_median"] >= 0.80)
        & (k_table["cross_chip_ari_median"] >= 0.50)
        & (k_table["min_niche_fraction"] >= 0.05)
        & (k_table["min_patient_coverage_n"] >= 18)
        & (k_table["min_chip_coverage_n"] == 3)
        & (k_table["max_single_patient_contribution"] <= 0.35)
    )
    k_table["selection_score"] = (
        zscore(k_table["repeat_ari_median"])
        + zscore(k_table["cross_chip_ari_median"])
        + zscore(k_table["silhouette_mean"])
        + 0.5 * zscore(k_table["spatial_edge_agreement"])
        + 0.5 * zscore(k_table["variance_explained_fraction"])
        - 0.25 * (k_table["k"] - min(K_VALUES))
    )
    candidates = k_table[k_table["selection_eligible"]]
    if args.forced_k is not None:
        selected_k = args.forced_k
        selection_status = "user_requested_exploratory_fixed_K"
    elif len(candidates):
        selected_k = int(candidates.sort_values(["selection_score", "k"], ascending=[False, True]).iloc[0]["k"])
        selection_status = "selected_among_predeclared_eligible_K"
    else:
        selected_k = int(k_table.sort_values(["selection_score", "k"], ascending=[False, True]).iloc[0]["k"])
        selection_status = "descriptive_fallback_no_K_passed_all_predeclared_criteria"
    k_table["selected"] = k_table["k"].eq(selected_k)
    k_table["selection_status"] = selection_status

    labels, niche_composition = relabel_by_composition(
        labels_by_k[selected_k], composition, state_names, eligible
    )
    data["official_full_radius_neighbor_n"] = neighbor_n
    data["niche_number"] = np.where(eligible, labels, np.nan)
    data["niche"] = np.where(eligible, np.char.add("CN", labels.astype(str)), "NOT_EVALUABLE")
    niche_names = [f"CN{i}" for i in range(1, selected_k + 1)]

    patient_rows = []
    for patient_id, frame in data.groupby("patient_id", sort=True):
        evaluable = frame[frame["niche"] != "NOT_EVALUABLE"]
        counts = evaluable["niche"].value_counts()
        row = {
            "patient_id": patient_id,
            "chip_label": frame["chip_label"].iloc[0],
            "disease": frame["disease"].iloc[0],
            "eligible_bin_n": len(frame),
            "evaluable_bin_n": len(evaluable),
            "median_nUMI": float(frame["nUMI"].median()),
        }
        for niche in niche_names:
            row[f"count_{niche}"] = int(counts.get(niche, 0))
            row[f"fraction_{niche}"] = counts.get(niche, 0) / max(1, len(evaluable))
        patient_rows.append(row)
    patient_table = pd.DataFrame(patient_rows)

    fraction_columns = [f"fraction_{niche}" for niche in niche_names]
    count_columns = [f"count_{niche}" for niche in niche_names]
    totals = patient_table["evaluable_bin_n"].to_numpy(dtype=float)
    counts = patient_table[count_columns].to_numpy(dtype=float)
    fractions_corrected = (counts + 0.5) / (totals[:, None] + 1.0)
    y = np.log(fractions_corrected / (1 - fractions_corrected))
    chip = patient_table["chip_label"].to_numpy()
    disease = patient_table["disease"].to_numpy()
    primary_fit = fit_hc3(y, design_matrix(chip, disease))
    permutation_p = permuted_contrast_p(y, chip, disease, primary_fit)
    depth = np.log1p(patient_table["median_nUMI"].to_numpy(dtype=float))
    depth_fit = fit_hc3(y, design_matrix(chip, disease, depth))
    result_rows = []
    for contrast in CONTRASTS:
        for index, niche in enumerate(niche_names):
            result_rows.append(
                {
                    "niche": niche,
                    "contrast": contrast,
                    "log_odds_effect": primary_fit[contrast]["effect"][index],
                    "odds_ratio": np.exp(primary_fit[contrast]["effect"][index]),
                    "hc3_se": primary_fit[contrast]["se"][index],
                    "hc3_p": primary_fit[contrast]["p"][index],
                    "permutation_p": permutation_p[contrast][index],
                    "depth_adjusted_log_odds_effect": depth_fit[contrast]["effect"][index],
                    "depth_adjusted_hc3_p": depth_fit[contrast]["p"][index],
                    "depth_adjusted_same_direction": np.sign(depth_fit[contrast]["effect"][index])
                    == np.sign(primary_fit[contrast]["effect"][index]),
                }
            )
    disease_results = pd.DataFrame(result_rows)
    disease_results["hc3_fdr_global"] = bh_adjust(disease_results["hc3_p"].to_numpy())
    disease_results["permutation_fdr_global"] = bh_adjust(
        disease_results["permutation_p"].to_numpy()
    )
    disease_results["supported_finding"] = (
        (disease_results["hc3_fdr_global"] < 0.05)
        & (disease_results["permutation_fdr_global"] < 0.05)
        & disease_results["depth_adjusted_same_direction"]
    )
    global_tests = global_composition_tests(patient_table, fraction_columns)

    chip_effect_rows = []
    for chip_value in sorted(patient_table["chip_label"].unique()):
        frame = patient_table[patient_table["chip_label"] == chip_value]
        for contrast, group_a, group_b in (
            ("IPF_vs_HC", "HC", "IPF"),
            ("SSc-ILD_vs_HC", "HC", "SSc-ILD"),
            ("SSc-ILD_vs_IPF", "IPF", "SSc-ILD"),
        ):
            for niche in niche_names:
                a = frame.loc[frame["disease"] == group_a, f"fraction_{niche}"]
                b = frame.loc[frame["disease"] == group_b, f"fraction_{niche}"]
                chip_effect_rows.append(
                    {
                        "chip_label": chip_value,
                        "contrast": contrast,
                        "niche": niche,
                        "group_a_mean_fraction": a.mean(),
                        "group_b_mean_fraction": b.mean(),
                        "fraction_difference": b.mean() - a.mean(),
                    }
                )
    chip_effects = pd.DataFrame(chip_effect_rows)

    component_rows = []
    for patient_id, frame in data[data["niche"] != "NOT_EVALUABLE"].groupby("patient_id", sort=True):
        frame_indices = frame.index.to_numpy()
        lookup = {value: index for index, value in enumerate(frame_indices)}
        edge_mask = np.isin(edges[:, 0], frame_indices) & np.isin(edges[:, 1], frame_indices)
        patient_edges = edges[edge_mask]
        for niche in niche_names:
            niche_global = frame.index[frame["niche"] == niche].to_numpy()
            niche_set = set(niche_global.tolist())
            adjacency = {value: [] for value in niche_global}
            for left, right in patient_edges:
                if left in niche_set and right in niche_set:
                    adjacency[left].append(right)
                    adjacency[right].append(left)
            visited = set()
            sizes = []
            for start in niche_global:
                if start in visited:
                    continue
                stack = [start]
                visited.add(start)
                size = 0
                while stack:
                    current = stack.pop()
                    size += 1
                    for neighbor in adjacency[current]:
                        if neighbor not in visited:
                            visited.add(neighbor)
                            stack.append(neighbor)
                sizes.append(size)
            component_rows.append(
                {
                    "patient_id": patient_id,
                    "chip_label": frame["chip_label"].iloc[0],
                    "disease": frame["disease"].iloc[0],
                    "niche": niche,
                    "niche_bin_n": len(niche_global),
                    "component_n": len(sizes),
                    "component_ge5_n": int(np.sum(np.asarray(sizes) >= 5)),
                    "largest_component_n": max(sizes) if sizes else 0,
                    "largest_component_fraction": max(sizes) / len(niche_global) if sizes else 0,
                }
            )
    components = pd.DataFrame(component_rows)

    method = pd.DataFrame(
        {
            "parameter": [
                "input", "reference", "umi_min", "patients", "radius_um", "center_bin",
                "radius_implementation", "transform", "fit_balance", "k_candidates", "repeats",
                "selection_rule", "disease_used_for_discovery", "primary_disease_model",
                "permutations", "seed",
            ],
            "value": [
                str(input_root), "six-dataset Reference18", 50, 21, RADIUS_UM, "excluded",
                "exact Euclidean radius within patient/ROI; no RANN k-cap; no cross-patient edges",
                (
                    "VoltRon-compatible zero-safe CLR of normalized neighbor-weight sums"
                    if args.feature_transform == "raw_clr"
                    else "zero-safe CLR followed by disease-blind equal-patient chip mean centering"
                ),
                f"{reference_per_patient} evaluable Bin50 per patient per repeat (minimum available, capped at {MAX_REFERENCE_PER_PATIENT})",
                ";".join(map(str, K_VALUES)), REPEATS,
                "predeclared stability, held-chip transfer, size, patient/chip coverage, single-patient dominance, spatial agreement, parsimony",
                "FALSE", "patient corrected-logit niche fraction ~ chip + disease; HC3",
                f"{PERMUTATIONS} disease-label permutations within chip", SEED,
            ],
        }
    )

    method.to_csv(output / "00_METHOD_AND_INPUT_LOCK.tsv", sep="\t", index=False)
    pd.DataFrame(
        {
            "input_file": [str(path.resolve()) for path in input_paths],
            "sha256": [sha256(path) for path in input_paths],
        }
    ).to_csv(output / "01_INPUT_MANIFEST.tsv", sep="\t", index=False)
    data[[
        "entity_key", "bin_id", "patient_id", "chip_label", "disease", "x_um", "y_um",
        "nUMI", "nGene", "official_full_radius_neighbor_n", "niche_number", "niche",
    ]].to_csv(output / "02_ALL21_SHARED_NICHE_ASSIGNMENTS.tsv.gz", sep="\t", index=False, compression="gzip")
    k_table.to_csv(output / "03_K_SELECTION_AUDIT.tsv", sep="\t", index=False)
    pd.concat(audits, ignore_index=True).to_csv(
        output / "04_K_REPEAT_AND_HELD_CHIP_AUDIT.tsv.gz", sep="\t", index=False, compression="gzip"
    )
    niche_composition.to_csv(output / "05_SELECTED_K_NICHE_COMPOSITION.tsv", sep="\t", index=False)
    patient_table.to_csv(output / "06_PATIENT_NICHE_ABUNDANCE.tsv", sep="\t", index=False)
    components.to_csv(output / "07_PATIENT_NICHE_SPATIAL_COMPONENTS.tsv", sep="\t", index=False)
    global_tests.to_csv(output / "08_GLOBAL_COMPOSITION_DISEASE_TESTS.tsv", sep="\t", index=False)
    disease_results.to_csv(output / "09_NICHE_DISEASE_TESTS.tsv", sep="\t", index=False)
    chip_effects.to_csv(output / "10_CHIP_SPECIFIC_EFFECTS.tsv", sep="\t", index=False)
    pd.DataFrame(chip_offsets, index=chip_order, columns=state_names).rename_axis(
        "chip_label"
    ).reset_index().to_csv(output / "11_CHIP_CLR_OFFSETS.tsv", sep="\t", index=False)

    string_dtype = h5py.string_dtype(encoding="utf-8")
    with h5py.File(output / "12_SHARED_NICHE_ASSAY.h5", "w") as handle:
        handle.create_dataset(
            "entity_key",
            data=data["entity_key"].astype(str).to_numpy(dtype=object),
            dtype=string_dtype,
        )
        handle.create_dataset(
            "state_names",
            data=np.asarray(state_names, dtype=object),
            dtype=string_dtype,
        )
        handle.create_dataset(
            "composition",
            data=composition.astype(np.float32),
            compression="gzip",
        )
        handle.create_dataset(
            "clr",
            data=features.astype(np.float32),
            compression="gzip",
        )
        handle.create_dataset(
            "neighbor_n",
            data=data["official_full_radius_neighbor_n"].to_numpy(dtype=np.int32),
            compression="gzip",
        )
        handle.create_dataset(
            "eligible",
            data=eligible.astype(np.uint8),
            compression="gzip",
        )
        handle.create_dataset(
            "selected_niche_number",
            data=labels.astype(np.int16),
            compression="gzip",
        )
        for k, candidate_labels in labels_by_k.items():
            handle.create_dataset(
                f"candidate_k{k}_label",
                data=candidate_labels.astype(np.int16),
                compression="gzip",
            )

    sns.set_theme(style="whitegrid")
    metric_long = k_table.melt(
        id_vars=["k", "selected"],
        value_vars=[
            "repeat_ari_median", "cross_chip_ari_median", "silhouette_mean",
            "spatial_edge_agreement", "variance_explained_fraction",
        ],
        var_name="metric",
        value_name="value",
    )
    fig, ax = plt.subplots(figsize=(10, 5))
    sns.lineplot(data=metric_long, x="k", y="value", hue="metric", marker="o", ax=ax)
    ax.axvline(selected_k, color="black", linestyle="--", label=f"selected K={selected_k}")
    ax.set_title("Disease-blind shared-niche K audit")
    ax.set_xticks(K_VALUES)
    fig.tight_layout()
    fig.savefig(figures / "FIG01_K_SELECTION_AUDIT.png", dpi=220)
    fig.savefig(figures / "FIG01_K_SELECTION_AUDIT.pdf")
    plt.close(fig)

    heat = niche_composition.set_index("niche")[[f"mean_{state}" for state in state_names]]
    heat.columns = state_names
    fig, ax = plt.subplots(figsize=(14, 2.5 + selected_k * 0.45))
    sns.heatmap(heat, cmap="viridis", ax=ax)
    ax.set_title(f"Shared All21 niche composition (K={selected_k})")
    fig.tight_layout()
    fig.savefig(figures / "FIG02_NICHE_REFERENCE18_COMPOSITION.png", dpi=220)
    fig.savefig(figures / "FIG02_NICHE_REFERENCE18_COMPOSITION.pdf")
    plt.close(fig)

    palette = dict(zip(niche_names, sns.color_palette("tab10", selected_k)))
    fig, axes = plt.subplots(7, 3, figsize=(15, 29))
    for ax, (patient_id, frame) in zip(axes.flat, data.groupby("patient_id", sort=True)):
        for niche in niche_names:
            subset = frame[frame["niche"] == niche]
            ax.scatter(subset["x_um"], subset["y_um"], s=2, color=palette[niche], rasterized=True)
        ax.set_title(f"{patient_id} | {frame['disease'].iloc[0]} | {frame['chip_label'].iloc[0]}", fontsize=9)
        ax.set_aspect("equal")
        ax.invert_yaxis()
        ax.axis("off")
    fig.suptitle(f"Shared All21 spatial niche atlas (K={selected_k})", y=0.995)
    fig.tight_layout()
    fig.savefig(figures / "FIG03_ALL21_SPATIAL_NICHE_ATLAS.png", dpi=220)
    fig.savefig(figures / "FIG03_ALL21_SPATIAL_NICHE_ATLAS.pdf")
    plt.close(fig)

    proportion_long = patient_table.melt(
        id_vars=["patient_id", "chip_label", "disease"],
        value_vars=fraction_columns,
        var_name="niche",
        value_name="fraction",
    )
    proportion_long["niche"] = proportion_long["niche"].str.removeprefix("fraction_")
    patient_order = patient_table.sort_values(["chip_label", "disease", "patient_id"])["patient_id"]
    pivot = proportion_long.pivot(index="patient_id", columns="niche", values="fraction").loc[patient_order]
    ax = pivot.plot(kind="bar", stacked=True, figsize=(15, 5), color=[palette[c] for c in pivot.columns])
    ax.set_ylabel("Niche fraction")
    ax.set_title("Shared niche composition by patient")
    ax.legend(title="Niche", bbox_to_anchor=(1.01, 1), loc="upper left")
    plt.tight_layout()
    plt.savefig(figures / "FIG04_PATIENT_NICHE_PROPORTIONS.png", dpi=220)
    plt.savefig(figures / "FIG04_PATIENT_NICHE_PROPORTIONS.pdf")
    plt.close()

    fig, axes = plt.subplots(1, selected_k, figsize=(4 * selected_k, 4), sharey=False)
    if selected_k == 1:
        axes = [axes]
    for ax, niche in zip(axes, niche_names):
        sns.boxplot(data=proportion_long[proportion_long["niche"] == niche], x="disease", y="fraction", ax=ax, color="white")
        sns.stripplot(data=proportion_long[proportion_long["niche"] == niche], x="disease", y="fraction", hue="chip_label", ax=ax, dodge=False)
        ax.set_title(niche)
        if ax.get_legend() is not None:
            ax.get_legend().remove()
    fig.suptitle("Patient-level shared niche fractions by disease")
    fig.tight_layout()
    fig.savefig(figures / "FIG05_DISEASE_PATIENT_DISTRIBUTIONS.png", dpi=220)
    fig.savefig(figures / "FIG05_DISEASE_PATIENT_DISTRIBUTIONS.pdf")
    plt.close(fig)

    validation_rows = [
        ("input_rows_126449", len(data) == 126449, len(data)),
        ("unique_entities", data["entity_key"].nunique() == len(data), data["entity_key"].nunique()),
        ("patients_21", data["patient_id"].nunique() == 21, data["patient_id"].nunique()),
        ("balanced_chip_disease", observed.equals(expected), observed.to_json()),
        ("states_18", len(state_names) == 18, len(state_names)),
        ("weights_sum_one", np.allclose(data[weight_columns].sum(axis=1), 1), data[weight_columns].sum(axis=1).min()),
        ("no_cross_patient_edges", np.all(data.loc[edges[:, 0], "patient_id"].to_numpy() == data.loc[edges[:, 1], "patient_id"].to_numpy()), len(edges)),
        ("nonempty_neighborhoods", eligible.sum() > 0, int(eligible.sum())),
        ("all_k_completed", set(k_table["k"]) == set(K_VALUES), k_table["k"].tolist()),
        ("selected_one_k", int(k_table["selected"].sum()) == 1, selected_k),
        ("disease_blind_selection", True, "disease not passed to clustering or K metrics"),
        ("all_evaluable_assigned", data.loc[eligible, "niche"].ne("NOT_EVALUABLE").all(), int(eligible.sum())),
        ("patient_fractions_sum_one", np.allclose(patient_table[fraction_columns].sum(axis=1), 1), patient_table[fraction_columns].sum(axis=1).min()),
        ("finite_disease_results", np.isfinite(disease_results.select_dtypes(include=np.number)).all().all(), len(disease_results)),
    ]
    validation = pd.DataFrame(validation_rows, columns=["check", "pass", "observed"])
    validation.to_csv(output / "15_FINAL_VALIDATION.tsv", sep="\t", index=False)
    if not validation["pass"].all():
        raise RuntimeError(f"Validation failed: {validation.loc[~validation['pass'], 'check'].tolist()}")

    selected_row = k_table[k_table["selected"]].iloc[0]
    summary = [
        "# All21 shared unbiased discrete niche atlas",
        "",
        "## Scope",
        "",
        "This is the primary boss-requested workflow: one disease-blind shared discrete niche atlas, followed by patient-level disease comparison.",
        "Historical chip-specific labels, B-lineage masks and aggregate labels were not used.",
        "",
        "## Frozen definition",
        "",
        f"- Input: {len(data):,} Reference18 UMI>=50 Bin50 from 21 patients.",
        f"- Center-excluded exact 50-um neighborhoods; evaluable bins: {eligible.sum():,}; non-evaluable isolated bins: {(~eligible).sum():,}.",
        f"- Feature transform: {args.feature_transform}.",
        f"- Disease-blind candidates: K={','.join(map(str, K_VALUES))}; {REPEATS} patient-balanced repeats per K.",
        f"- Selected K: {selected_k} ({selection_status}).",
        f"- Selected-K repeat median ARI: {selected_row['repeat_ari_median']:.3f}; held-chip median ARI: {selected_row['cross_chip_ari_median']:.3f}.",
        "",
        "## Disease comparison",
        "",
        f"- Global composition tests with FDR<0.05: {(global_tests['fdr'] < 0.05).sum()} / {len(global_tests)}.",
        f"- Niche-contrast tests with global permutation FDR<0.05: {(disease_results['permutation_fdr_global'] < 0.05).sum()} / {len(disease_results)}.",
        f"- Fully supported niche findings: {disease_results['supported_finding'].sum()}.",
        "",
        "Disease labels were introduced only after K and shared niche assignments were frozen.",
    ]
    (output / "16_SUMMARY.md").write_text("\n".join(summary) + "\n", encoding="utf-8")
    (output / "SESSION_INFO.json").write_text(
        json.dumps({"python": sys.version, "platform": platform.platform(), "numpy": np.__version__, "pandas": pd.__version__}, indent=2),
        encoding="utf-8",
    )
    command = (
        f'"{sys.executable}" "{Path(__file__).resolve()}" '
        f'--input-root "{input_root}" --output-dir "{output}" '
        f'--feature-transform "{args.feature_transform}"'
        + (f' --forced-k "{args.forced_k}"' if args.forced_k is not None else "")
    )
    (output / "RUN_LOG.txt").write_text(command + "\n", encoding="utf-8")
    (output / "COMPLETE.ok").write_text("PASS\n", encoding="utf-8")
    print(output)


if __name__ == "__main__":
    main()
