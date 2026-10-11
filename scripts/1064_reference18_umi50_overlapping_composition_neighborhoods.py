#!/usr/bin/env python3
"""Patient-level disease testing of overlapping composition neighborhoods.

This analysis reuses the frozen Reference18 UMI>=50, 50-um neighborhood
composition assay. It does not rerun RCTD, rebuild spatial neighborhoods, or
select a discrete K. Disease-blind, patient-balanced reference bins define
overlapping k-nearest-neighbor regions in the common 18-state composition
space. Abundance is then tested with patients as independent replicates.
"""

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
from scipy.stats import t as student_t
from sklearn.neighbors import NearestNeighbors


SEED = 20261011
PERMUTATIONS = 10_000
REFERENCE_PER_PATIENT = 500
REFERENCE_K = 100
ANCHOR_PROPORTION = 0.10
PRIMARY_K = 500
SENSITIVITY_K = (300, 800)
MIN_BINS_PER_SUPPORTED_PATIENT = 5
MIN_SUPPORTED_PATIENTS = 6
MIN_CHIPS_WITH_TWO_PATIENTS = 2
MAX_PATIENT_FRACTION = 0.80
CONTRASTS = ("IPF_vs_HC", "SSc-ILD_vs_HC", "SSc-ILD_vs_IPF")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input-h5", required=True, type=Path)
    parser.add_argument("--output-dir", required=True, type=Path)
    return parser.parse_args()


def decode(values: np.ndarray) -> np.ndarray:
    return np.asarray(
        [x.decode("utf-8") if isinstance(x, bytes) else str(x) for x in values],
        dtype=object,
    )


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
    adjusted = np.minimum(adjusted, 1.0)
    restored = np.empty_like(adjusted)
    restored[order] = adjusted
    output[valid] = restored
    return output


def voltron_clr(composition: np.ndarray) -> np.ndarray:
    """Reproduce the frozen workflow's VoltRon zero-safe CLR transform."""
    x = np.asarray(composition, dtype=np.float64)
    totals = x.sum(axis=1, keepdims=True)
    if np.any(~np.isfinite(totals)) or np.any(totals <= 0):
        raise RuntimeError("Invalid composition row total")
    x = x / totals
    positive_log_sum = np.where(x > 0, np.log1p(x), 0.0).sum(axis=1)
    denominator = np.exp(positive_log_sum / x.shape[1])
    transformed = np.log1p(x / denominator[:, None])
    if not np.isfinite(transformed).all():
        raise RuntimeError("Non-finite CLR output")
    return transformed.astype(np.float32)


def design_matrix(chip: np.ndarray, disease: np.ndarray, depth: np.ndarray | None = None) -> np.ndarray:
    columns = [
        np.ones(len(chip)),
        (chip == "Y40105J2").astype(float),
        (chip == "Y40105L3").astype(float),
        (disease == "IPF").astype(float),
        (disease == "SSc-ILD").astype(float),
    ]
    if depth is not None:
        columns.append((depth - depth.mean()) / depth.std(ddof=1))
    return np.column_stack(columns)


def contrast_vectors(p: int) -> dict[str, np.ndarray]:
    vectors = {}
    for name, ipf, ssc in (
        ("IPF_vs_HC", 1.0, 0.0),
        ("SSc-ILD_vs_HC", 0.0, 1.0),
        ("SSc-ILD_vs_IPF", -1.0, 1.0),
    ):
        vector = np.zeros(p)
        vector[3] = ipf
        vector[4] = ssc
        vectors[name] = vector
    return vectors


def fit_hc3(y: np.ndarray, x: np.ndarray) -> dict[str, dict[str, np.ndarray]]:
    if np.linalg.matrix_rank(x) != x.shape[1]:
        raise RuntimeError("Disease/chip design matrix is rank deficient")
    xtx_inv = np.linalg.inv(x.T @ x)
    beta = xtx_inv @ x.T @ y
    residual = y - x @ beta
    leverage = np.sum((x @ xtx_inv) * x, axis=1)
    scaled_residual_sq = (residual / (1.0 - leverage[:, None])) ** 2
    df = x.shape[0] - x.shape[1]
    result = {}
    for name, contrast in contrast_vectors(x.shape[1]).items():
        effect = contrast @ beta
        influence = x @ xtx_inv @ contrast
        variance = (influence[:, None] ** 2 * scaled_residual_sq).sum(axis=0)
        se = np.sqrt(np.maximum(variance, 0.0))
        statistic = np.divide(effect, se, out=np.zeros_like(effect), where=se > 0)
        p_value = 2.0 * student_t.sf(np.abs(statistic), df=df)
        result[name] = {"effect": effect, "se": se, "p": p_value}
    return result


def permuted_p_values(
    y: np.ndarray,
    chip: np.ndarray,
    disease: np.ndarray,
    observed: dict[str, dict[str, np.ndarray]],
    seed: int,
) -> dict[str, np.ndarray]:
    rng = np.random.default_rng(seed)
    exceed = {name: np.zeros(y.shape[1], dtype=np.int32) for name in CONTRASTS}
    chip_indices = [np.flatnonzero(chip == value) for value in np.unique(chip)]
    for _ in range(PERMUTATIONS):
        permuted = disease.copy()
        for indices in chip_indices:
            permuted[indices] = rng.permutation(permuted[indices])
        x = design_matrix(chip, permuted)
        beta = np.linalg.inv(x.T @ x) @ x.T @ y
        for name, contrast in contrast_vectors(x.shape[1]).items():
            effect = contrast @ beta
            exceed[name] += np.abs(effect) >= (np.abs(observed[name]["effect"]) - 1e-12)
    return {name: (counts + 1.0) / (PERMUTATIONS + 1.0) for name, counts in exceed.items()}


def patient_balanced_reference(patient: np.ndarray) -> np.ndarray:
    rng = np.random.default_rng(SEED)
    selected = []
    for value in sorted(np.unique(patient)):
        indices = np.flatnonzero(patient == value)
        if len(indices) < REFERENCE_PER_PATIENT:
            raise RuntimeError(f"Patient {value} has fewer than {REFERENCE_PER_PATIENT} bins")
        selected.append(rng.choice(indices, REFERENCE_PER_PATIENT, replace=False))
    return np.sort(np.concatenate(selected))


def refined_anchor_indices(reference_features: np.ndarray) -> tuple[np.ndarray, int]:
    rng = np.random.default_rng(SEED + 1)
    graph = NearestNeighbors(n_neighbors=REFERENCE_K, metric="euclidean", n_jobs=4)
    graph.fit(reference_features)
    n_initial = int(np.ceil(len(reference_features) * ANCHOR_PROPORTION))
    initial = np.sort(rng.choice(len(reference_features), n_initial, replace=False))
    neighbors = graph.kneighbors(reference_features[initial], return_distance=False)
    medians = np.median(reference_features[neighbors], axis=1)
    refined = graph.kneighbors(medians, n_neighbors=1, return_distance=False).ravel()
    return np.unique(refined), n_initial


def membership_counts(
    full_features: np.ndarray,
    centers: np.ndarray,
    patient_codes: np.ndarray,
    k: int,
) -> tuple[np.ndarray, np.ndarray]:
    model = NearestNeighbors(n_neighbors=k, metric="euclidean", n_jobs=4)
    model.fit(full_features)
    indices = model.kneighbors(centers, return_distance=False)
    counts = np.vstack(
        [np.bincount(patient_codes[row], minlength=patient_codes.max() + 1) for row in indices]
    )
    return indices, counts


def support_filter(counts: np.ndarray, patient_chip: np.ndarray) -> pd.DataFrame:
    supported = counts >= MIN_BINS_PER_SUPPORTED_PATIENT
    patient_n = supported.sum(axis=1)
    chip_n = np.zeros(len(counts), dtype=int)
    for chip in np.unique(patient_chip):
        chip_n += (supported[:, patient_chip == chip].sum(axis=1) >= 2).astype(int)
    max_fraction = counts.max(axis=1) / counts.sum(axis=1)
    passed = (
        (patient_n >= MIN_SUPPORTED_PATIENTS)
        & (chip_n >= MIN_CHIPS_WITH_TWO_PATIENTS)
        & (max_fraction <= MAX_PATIENT_FRACTION)
    )
    return pd.DataFrame(
        {
            "supported_patient_n": patient_n,
            "chips_with_at_least_two_supported_patients": chip_n,
            "max_patient_fraction": max_fraction,
            "support_filter_pass": passed,
        }
    )


def loo_direction_fraction(
    y: np.ndarray,
    chip: np.ndarray,
    disease: np.ndarray,
    observed: dict[str, dict[str, np.ndarray]],
) -> dict[str, np.ndarray]:
    same = {name: np.zeros(y.shape[1], dtype=int) for name in CONTRASTS}
    valid = {name: np.zeros(y.shape[1], dtype=int) for name in CONTRASTS}
    for omitted in range(len(chip)):
        keep = np.arange(len(chip)) != omitted
        x = design_matrix(chip[keep], disease[keep])
        if np.linalg.matrix_rank(x) != x.shape[1]:
            continue
        beta = np.linalg.inv(x.T @ x) @ x.T @ y[keep]
        for name, contrast in contrast_vectors(x.shape[1]).items():
            effect = contrast @ beta
            obs = observed[name]["effect"]
            nonzero = (effect != 0) & (obs != 0)
            same[name] += ((np.sign(effect) == np.sign(obs)) & nonzero).astype(int)
            valid[name] += nonzero.astype(int)
    return {
        name: np.divide(same[name], valid[name], out=np.zeros(y.shape[1]), where=valid[name] > 0)
        for name in CONTRASTS
    }


def analyze_scale(
    counts: np.ndarray,
    patients: pd.DataFrame,
    neighborhood_ids: np.ndarray,
    k: int,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    support = support_filter(counts, patients["chip_id"].to_numpy())
    support.insert(0, "neighborhood_id", neighborhood_ids)
    support.insert(1, "membership_k", k)
    passed = support["support_filter_pass"].to_numpy()
    totals = patients["valid_bin_n"].to_numpy(dtype=float)
    proportions = (counts.T + 0.5) / (totals[:, None] + 1.0)
    y_all = np.log(proportions / (1.0 - proportions))
    y = y_all[:, passed]
    chip = patients["chip_id"].to_numpy()
    disease = patients["disease"].to_numpy()
    x = design_matrix(chip, disease)
    primary = fit_hc3(y, x)
    permutation = permuted_p_values(y, chip, disease, primary, SEED + k)
    depth = np.log1p(patients["median_nUMI"].to_numpy(dtype=float))
    depth_fit = fit_hc3(y, design_matrix(chip, disease, depth))
    loo = loo_direction_fraction(y, chip, disease, primary)
    kept_ids = neighborhood_ids[passed]
    rows = []
    for contrast in CONTRASTS:
        for j, neighborhood_id in enumerate(kept_ids):
            effect = primary[contrast]["effect"][j]
            depth_effect = depth_fit[contrast]["effect"][j]
            rows.append(
                {
                    "neighborhood_id": neighborhood_id,
                    "membership_k": k,
                    "contrast": contrast,
                    "log_odds_effect": effect,
                    "odds_ratio": np.exp(effect),
                    "hc3_se": primary[contrast]["se"][j],
                    "hc3_p": primary[contrast]["p"][j],
                    "permutation_p": permutation[contrast][j],
                    "depth_adjusted_log_odds_effect": depth_effect,
                    "depth_adjusted_hc3_p": depth_fit[contrast]["p"][j],
                    "depth_adjusted_same_direction": np.sign(depth_effect) == np.sign(effect),
                    "loo_direction_fraction": loo[contrast][j],
                }
            )
    result = pd.DataFrame(rows)
    result["hc3_fdr_within_contrast"] = result.groupby("contrast")["hc3_p"].transform(
        lambda x: bh_adjust(x.to_numpy())
    )
    result["permutation_fdr_within_contrast"] = result.groupby("contrast")[
        "permutation_p"
    ].transform(lambda x: bh_adjust(x.to_numpy()))
    result["hc3_fdr_global"] = bh_adjust(result["hc3_p"].to_numpy())
    result["permutation_fdr_global"] = bh_adjust(result["permutation_p"].to_numpy())
    result["supported_finding"] = (
        (result["hc3_fdr_global"] < 0.05)
        & (result["permutation_fdr_global"] < 0.05)
        & result["depth_adjusted_same_direction"]
        & (result["loo_direction_fraction"] >= 0.80)
    )
    return result, support


def make_figures(
    output: Path,
    patients: pd.DataFrame,
    primary_counts: np.ndarray,
    results: pd.DataFrame,
    compositions: pd.DataFrame,
) -> None:
    sns.set_theme(style="whitegrid")
    primary = results[results["membership_k"] == PRIMARY_K].copy()
    fig, axes = plt.subplots(1, 2, figsize=(12, 4.5))
    sns.barplot(data=patients, x="chip_short", y="valid_bin_n", hue="disease", ax=axes[0])
    axes[0].set_title("Patient input bin counts")
    axes[0].set_xlabel("Chip")
    axes[0].set_ylabel("Valid Bin50")
    support_n = (primary_counts >= MIN_BINS_PER_SUPPORTED_PATIENT).sum(axis=1)
    axes[1].hist(support_n, bins=np.arange(0.5, 22.5, 1), color="#4472C4")
    axes[1].axvline(MIN_SUPPORTED_PATIENTS, color="black", linestyle="--")
    axes[1].set_title(f"Primary k={PRIMARY_K}: patient support")
    axes[1].set_xlabel("Patients with >=5 member bins")
    axes[1].set_ylabel("Composition neighborhoods")
    fig.tight_layout()
    fig.savefig(output / "FIG01_DESIGN_AND_SUPPORT_QC.png", dpi=220)
    fig.savefig(output / "FIG01_DESIGN_AND_SUPPORT_QC.pdf")
    plt.close(fig)

    fig, axes = plt.subplots(1, 3, figsize=(15, 4.5), sharey=True)
    for ax, contrast in zip(axes, CONTRASTS):
        data = primary[primary["contrast"] == contrast]
        ax.scatter(data["log_odds_effect"], -np.log10(data["permutation_p"]), s=12, alpha=0.65)
        ax.axhline(-np.log10(0.05), color="grey", linestyle="--")
        ax.axvline(0, color="grey", linewidth=0.8)
        ax.set_title(contrast)
        ax.set_xlabel("Patient-level log odds effect")
    axes[0].set_ylabel("-log10 within-chip permutation p")
    fig.tight_layout()
    fig.savefig(output / "FIG02_PRIMARY_DISEASE_EFFECTS.png", dpi=220)
    fig.savefig(output / "FIG02_PRIMARY_DISEASE_EFFECTS.pdf")
    plt.close(fig)

    top_ids = (
        primary.sort_values("permutation_p")
        .drop_duplicates("neighborhood_id")
        .head(20)["neighborhood_id"]
        .to_numpy()
    )
    comp = compositions.set_index("neighborhood_id").loc[top_ids]
    comp = comp[[c for c in comp.columns if c.startswith("mean_")]]
    comp.columns = [c.removeprefix("mean_") for c in comp.columns]
    fig, ax = plt.subplots(figsize=(13, max(6, len(comp) * 0.32)))
    sns.heatmap(comp, cmap="viridis", ax=ax)
    ax.set_title("Composition of top primary exploratory neighborhoods")
    ax.set_xlabel("Reference18 state")
    ax.set_ylabel("Composition neighborhood")
    fig.tight_layout()
    fig.savefig(output / "FIG03_TOP_NEIGHBORHOOD_COMPOSITION.png", dpi=220)
    fig.savefig(output / "FIG03_TOP_NEIGHBORHOOD_COMPOSITION.pdf")
    plt.close(fig)


def main() -> None:
    args = parse_args()
    input_h5 = args.input_h5.resolve()
    output = args.output_dir.resolve()
    output.mkdir(parents=True, exist_ok=False)

    with h5py.File(input_h5, "r") as handle:
        composition = handle["radius_50um/niche_composition"][:]
        state_names = decode(handle["state_names"][:])
        patient = decode(handle["metadata/patient_id"][:])
        sample = decode(handle["metadata/sample_id"][:])
        chip = decode(handle["metadata/chip_id"][:])
        disease = decode(handle["metadata/disease"][:])
        n_umi = handle["metadata/nUMI"][:]

    if composition.shape != (126_449, 18):
        raise RuntimeError(f"Unexpected composition shape: {composition.shape}")
    if not np.allclose(composition.sum(axis=1), 1.0, atol=2e-5):
        raise RuntimeError("Composition rows do not sum to one")
    if not np.array_equal(patient, sample):
        raise RuntimeError("patient_id and sample_id differ")

    patient_order = np.asarray(sorted(np.unique(patient)), dtype=object)
    patient_code_map = {value: i for i, value in enumerate(patient_order)}
    patient_codes = np.asarray([patient_code_map[value] for value in patient], dtype=int)
    patient_rows = []
    for value in patient_order:
        idx = np.flatnonzero(patient == value)
        patient_rows.append(
            {
                "patient_id": value,
                "sample_id": sample[idx[0]],
                "chip_id": chip[idx[0]],
                "chip_short": chip[idx[0]][-2:],
                "disease": disease[idx[0]],
                "valid_bin_n": len(idx),
                "median_nUMI": float(np.median(n_umi[idx])),
            }
        )
    patients = pd.DataFrame(patient_rows)
    expected = pd.DataFrame(
        [[2, 2, 3]] * 3,
        index=["Y40102K8", "Y40105J2", "Y40105L3"],
        columns=["HC", "IPF", "SSc-ILD"],
    )
    observed = pd.crosstab(patients["chip_id"], patients["disease"]).reindex(
        index=expected.index, columns=expected.columns, fill_value=0
    )
    if not observed.equals(expected):
        raise RuntimeError(f"Unexpected chip/disease design:\n{observed}")

    features = voltron_clr(composition)
    reference_indices = patient_balanced_reference(patient)
    anchor_reference_indices, n_initial = refined_anchor_indices(features[reference_indices])
    anchor_full_indices = reference_indices[anchor_reference_indices]
    centers = features[anchor_full_indices]
    neighborhood_ids = np.asarray(
        [f"CN{i:04d}" for i in range(1, len(centers) + 1)], dtype=object
    )

    all_results = []
    all_support = []
    scale_counts = {}
    scale_indices = {}
    for k in (PRIMARY_K,) + SENSITIVITY_K:
        indices, counts = membership_counts(features, centers, patient_codes, k)
        scale_indices[k] = indices
        scale_counts[k] = counts
        results, support = analyze_scale(counts, patients, neighborhood_ids, k)
        all_results.append(results)
        all_support.append(support)
    results = pd.concat(all_results, ignore_index=True)
    support = pd.concat(all_support, ignore_index=True)

    primary_indices = scale_indices[PRIMARY_K]
    global_mean = composition.mean(axis=0)
    composition_rows = []
    for i, neighborhood_id in enumerate(neighborhood_ids):
        mean_composition = composition[primary_indices[i]].mean(axis=0)
        enrichment = mean_composition / np.maximum(global_mean, 1e-12)
        top_mean = np.argsort(mean_composition)[::-1][:3]
        top_enrichment = np.argsort(enrichment)[::-1][:3]
        row = {
            "neighborhood_id": neighborhood_id,
            "anchor_bin_index": int(anchor_full_indices[i]),
            "top3_mean_states": ";".join(state_names[top_mean]),
            "top3_enriched_states": ";".join(state_names[top_enrichment]),
        }
        row.update({f"mean_{state}": mean_composition[j] for j, state in enumerate(state_names)})
        row.update({f"enrichment_{state}": enrichment[j] for j, state in enumerate(state_names)})
        composition_rows.append(row)
    compositions = pd.DataFrame(composition_rows)
    results = results.merge(
        compositions[["neighborhood_id", "top3_mean_states", "top3_enriched_states"]],
        on="neighborhood_id",
        how="left",
    )

    count_rows = []
    for k, counts in scale_counts.items():
        for i, neighborhood_id in enumerate(neighborhood_ids):
            for j, row in patients.iterrows():
                count_rows.append(
                    {
                        "membership_k": k,
                        "neighborhood_id": neighborhood_id,
                        "patient_id": row["patient_id"],
                        "chip_id": row["chip_id"],
                        "disease": row["disease"],
                        "member_bin_n": int(counts[i, j]),
                        "valid_bin_n": int(row["valid_bin_n"]),
                        "member_fraction": counts[i, j] / row["valid_bin_n"],
                    }
                )
    patient_counts = pd.DataFrame(count_rows)

    primary = results[results["membership_k"] == PRIMARY_K]
    sensitivity_summary = []
    for k in SENSITIVITY_K:
        other = results[results["membership_k"] == k]
        for contrast in CONTRASTS:
            left = primary[primary["contrast"] == contrast][
                ["neighborhood_id", "log_odds_effect"]
            ]
            right = other[other["contrast"] == contrast][
                ["neighborhood_id", "log_odds_effect"]
            ]
            joined = left.merge(right, on="neighborhood_id", suffixes=("_primary", "_sensitivity"))
            sensitivity_summary.append(
                {
                    "sensitivity_k": k,
                    "contrast": contrast,
                    "shared_tested_neighborhood_n": len(joined),
                    "spearman_effect_correlation": joined[
                        ["log_odds_effect_primary", "log_odds_effect_sensitivity"]
                    ].corr(method="spearman").iloc[0, 1],
                    "direction_agreement": np.mean(
                        np.sign(joined["log_odds_effect_primary"])
                        == np.sign(joined["log_odds_effect_sensitivity"])
                    ),
                }
            )
    sensitivity_summary = pd.DataFrame(sensitivity_summary)

    method = pd.DataFrame(
        {
            "parameter": [
                "input_h5", "input_sha256", "reference", "umi_min", "spatial_radius_um",
                "center_bin", "feature", "transform", "patients", "reference_per_patient",
                "reference_k", "anchor_proportion", "anchor_refinement", "primary_membership_k",
                "sensitivity_membership_k", "support_filter", "primary_model", "permutations",
                "multiple_testing", "supported_finding_rule", "seed",
            ],
            "value": [
                str(input_h5), sha256(input_h5), "six-dataset Reference18", 50, 50,
                "included (inherited frozen assay)", "18-state local spatial-neighborhood composition",
                "exact frozen-workflow VoltRon zero-safe CLR", 21, REFERENCE_PER_PATIENT,
                REFERENCE_K, ANCHOR_PROPORTION,
                "random disease-blind reference vertices -> median of reference-kNN -> nearest observed reference bin -> deduplicate",
                PRIMARY_K, ";".join(map(str, SENSITIVITY_K)),
                f">={MIN_BINS_PER_SUPPORTED_PATIENT} bins in >={MIN_SUPPORTED_PATIENTS} patients; >=2 patients in >=2 chips; max patient fraction <={MAX_PATIENT_FRACTION}",
                "HC3 OLS on corrected-logit patient fraction: chip + disease",
                f"{PERMUTATIONS} disease-label permutations within chip",
                "BH within contrast and globally across all three contrasts",
                "global HC3 FDR<0.05 AND global permutation FDR<0.05 AND depth-adjusted same direction AND LOO direction fraction>=0.80",
                SEED,
            ],
        }
    )

    method.to_csv(output / "00_METHOD_AND_INPUT_LOCK.tsv", sep="\t", index=False)
    patients.to_csv(output / "01_PATIENT_DESIGN_AND_DEPTH.tsv", sep="\t", index=False)
    support.to_csv(output / "02_NEIGHBORHOOD_SUPPORT_QC.tsv", sep="\t", index=False)
    compositions.to_csv(output / "03_PRIMARY_NEIGHBORHOOD_COMPOSITION.tsv", sep="\t", index=False)
    patient_counts.to_csv(
        output / "04_PATIENT_NEIGHBORHOOD_COUNTS.tsv.gz", sep="\t", index=False, compression="gzip"
    )
    results.to_csv(output / "05_DISEASE_TEST_RESULTS.tsv", sep="\t", index=False)
    sensitivity_summary.to_csv(output / "06_SCALE_SENSITIVITY.tsv", sep="\t", index=False)

    membership_counts_per_bin = np.bincount(
        primary_indices.ravel(), minlength=len(composition)
    )
    validation_rows = [
        ("input_shape_126449_by_18", composition.shape == (126449, 18), str(composition.shape)),
        ("composition_finite", np.isfinite(composition).all(), "all finite"),
        ("composition_rows_sum_to_one", np.allclose(composition.sum(1), 1, atol=2e-5), str(composition.sum(1).min())),
        ("twenty_one_patients", len(patients) == 21, str(len(patients))),
        ("balanced_chip_disease_design", observed.equals(expected), observed.to_json()),
        ("reference_size_10500", len(reference_indices) == 10500, str(len(reference_indices))),
        ("anchors_after_refinement", len(anchor_full_indices) > 100, str(len(anchor_full_indices))),
        ("primary_neighborhood_exact_size", np.all(scale_counts[PRIMARY_K].sum(1) == PRIMARY_K), str(scale_counts[PRIMARY_K].sum(1).min())),
        ("primary_tested_neighborhoods", len(primary) > 0, str(primary["neighborhood_id"].nunique())),
        ("permutation_p_resolution", primary["permutation_p"].min() >= 1 / (PERMUTATIONS + 1), str(primary["permutation_p"].min())),
        ("finite_primary_effects", np.isfinite(primary["log_odds_effect"]).all(), "all finite"),
        ("membership_has_overlap", membership_counts_per_bin.max() > 1, str(membership_counts_per_bin.max())),
    ]
    validation = pd.DataFrame(validation_rows, columns=["check", "pass", "observed"])
    validation.to_csv(output / "10_FINAL_VALIDATION.tsv", sep="\t", index=False)
    if not validation["pass"].all():
        failed = validation.loc[~validation["pass"], "check"].tolist()
        raise RuntimeError(f"Validation failed: {failed}")

    make_figures(output, patients, scale_counts[PRIMARY_K], results, compositions)

    summary_lines = [
        "# Reference18 UMI50 overlapping composition-neighborhood analysis",
        "",
        "## Scope",
        "",
        "This analysis reused the frozen 50-um, center-included Reference18 UMI>=50 assay. "
        "It did not rerun RCTD, rebuild physical neighborhoods, or choose a discrete K.",
        "",
        "## Construction",
        "",
        f"- Patient-balanced reference: {len(reference_indices):,} bins ({REFERENCE_PER_PATIENT} per patient).",
        f"- Initial/refined unique anchors: {n_initial:,}/{len(anchor_full_indices):,}.",
        f"- Primary overlapping composition-neighborhood size: {PRIMARY_K} nearest bins; sensitivity: {SENSITIVITY_K[0]} and {SENSITIVITY_K[1]}.",
        f"- Primary covered bins: {(membership_counts_per_bin > 0).sum():,}/{len(composition):,} ({(membership_counts_per_bin > 0).mean():.1%}); median memberships among covered bins: {np.median(membership_counts_per_bin[membership_counts_per_bin > 0]):.0f}.",
        "",
        "## Patient-level inference",
        "",
        f"- Primary tested composition neighborhoods: {primary['neighborhood_id'].nunique():,}.",
        f"- Tests: {len(primary):,} across three disease contrasts.",
        f"- Global HC3 FDR <0.05: {(primary['hc3_fdr_global'] < 0.05).sum():,}.",
        f"- Global within-chip permutation FDR <0.05: {(primary['permutation_fdr_global'] < 0.05).sum():,}.",
        f"- Fully supported findings: {primary['supported_finding'].sum():,}.",
        "",
        "Disease labels were never used to construct anchors or memberships. Interpretation of "
        "the Reference18 composition was attached only after testing.",
    ]
    (output / "11_SUMMARY.md").write_text("\n".join(summary_lines) + "\n", encoding="utf-8")
    session = {
        "python": sys.version,
        "platform": platform.platform(),
        "numpy": np.__version__,
        "pandas": pd.__version__,
    }
    (output / "SESSION_INFO.json").write_text(json.dumps(session, indent=2), encoding="utf-8")
    command = f'"{sys.executable}" "{Path(__file__).resolve()}" --input-h5 "{input_h5}" --output-dir "{output}"'
    (output / "RUN_LOG.txt").write_text(command + "\n", encoding="utf-8")
    (output / "COMPLETE.ok").write_text("PASS\n", encoding="utf-8")
    print(output)


if __name__ == "__main__":
    main()
