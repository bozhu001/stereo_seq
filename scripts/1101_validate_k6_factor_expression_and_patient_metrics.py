#!/usr/bin/env python3
"""Validate frozen per-chip Neighborhood NMF K=6 factors with raw RNA.

Purpose
-------
This bespoke downstream analysis reuses persisted K=6 loadings/usages and does
not refit RCTD or NMF. It (1) creates patient-level continuous/spatial factor
metrics and (2) validates each factor against center-Bin50 raw RNA using the
six-dataset Reference18 primary-marker panels. Expression validation is done
within patient, with technical residualization and high-usage versus matched
low-usage pseudobulk contrasts. It is descriptive and performs no disease test.

Inputs
------
- Frozen per-chip Neighborhood NMF K=6 NPZ and loading tables.
- Frozen all-21 assignment metadata and raw-count H5AD.
- Six-dataset conservative Reference18 primary-marker definitions.

Outputs
-------
Timestamped directory under --output-root containing TSV tables, English
figures, a Chinese method/results summary, run log, session information and a
completion sentinel.

Exact run command
-----------------
/home/bozhu/miniconda3/envs/banksy_stereo/bin/python \
  scripts/1101_validate_k6_factor_expression_and_patient_metrics.py \
  --output-root results/reanalysis
"""

from __future__ import annotations

import argparse
import hashlib
import json
import platform
import sys
import time
from datetime import datetime
from pathlib import Path

import anndata as ad
import matplotlib
import numpy as np
import pandas as pd
import scipy.sparse as sp
from matplotlib.colors import TwoSlopeNorm
from scipy.sparse import csr_matrix
from scipy.sparse.csgraph import connected_components
from scipy.spatial import cKDTree
from scipy.stats import spearmanr

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import seaborn as sns


ROOT = Path(__file__).resolve().parents[1]
MODEL_RUN = (
    ROOT
    / "results/reanalysis/REFERENCE18_PER_CHIP_DUAL_NMF_20261012_204500_attempt4_alpha0"
)
ASSIGNMENTS = (
    ROOT
    / "results/reanalysis/ALL21_SHARED_UNBIASED_DISCRETE_NICHE_ATLAS_20261011_191652_attempt4"
    / "02_ALL21_SHARED_NICHE_ASSIGNMENTS.tsv.gz"
)
RAW_H5AD = (
    ROOT
    / "results/visium way/bin50_g15_umi25"
    / "BIN50_joint_reanalysis_G15_UMI25_raw_counts_260924170658.h5ad"
)
GENE_ANNOTATION = (
    ROOT
    / "results/visium way/04_gene_biotype_audit/gene_biotype_audit_260924174116"
    / "01_gene_annotation_audit_260924174116.tsv"
)
MARKER_DEFINITIONS = (
    ROOT
    / "results/reanalysis/SIX_DATASET_LEVEL2_REFERENCE_AUDIT_20260919_004011"
    / "08_CONSERVATIVE_REFERENCE_PROPOSAL.tsv"
)
PRIOR_AUDIT = (
    ROOT
    / "results/reanalysis/REFERENCE18_NMF_COMPREHENSIVE_METHOD_AUDIT_20261012_233739"
    / "04_FACTOR_BIOLOGICAL_AUDIT.tsv"
)

CHIPS = ("K8", "J2", "L3")
RANK = 6
SEED = 20261013
HIGH_QUANTILE = 0.90
LOW_CANDIDATE_QUANTILE = 0.50
SPATIAL_EDGE_RADIUS_UM = 36.0
EXPRESSION_NEIGHBORHOOD_RADIUS_UM = 50.0
EPS = 1e-12


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--output-root",
        type=Path,
        default=ROOT / "results/reanalysis",
    )
    return parser.parse_args()


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def safe_spearman(first: np.ndarray, second: np.ndarray) -> float:
    if len(first) < 3 or np.nanstd(first) <= 0 or np.nanstd(second) <= 0:
        return float("nan")
    return float(spearmanr(first, second).statistic)


def residualize(values: np.ndarray, technical: np.ndarray) -> np.ndarray:
    values = np.asarray(values, dtype=float)
    technical = np.asarray(technical, dtype=float)
    means = np.nanmean(technical, axis=0)
    scales = np.nanstd(technical, axis=0)
    scales[scales <= 0] = 1.0
    design = np.column_stack(
        [np.ones(len(values)), (technical - means[None, :]) / scales[None, :]]
    )
    coefficients = np.linalg.lstsq(design, values, rcond=None)[0]
    return values - design @ coefficients


def moran_i(values: np.ndarray, edges: np.ndarray) -> float:
    if len(edges) == 0:
        return float("nan")
    centered = np.asarray(values, dtype=float) - np.nanmean(values)
    denominator = float(np.nansum(centered**2))
    if denominator <= 0:
        return float("nan")
    numerator = float(np.sum(centered[edges[:, 0]] * centered[edges[:, 1]]))
    return float(len(values) / len(edges) * numerator / denominator)


def local_high_enrichment(high: np.ndarray, edges: np.ndarray) -> float:
    prevalence = float(np.mean(high))
    if prevalence <= 0 or len(edges) == 0:
        return float("nan")
    degrees = np.bincount(edges.ravel(), minlength=len(high))
    high_stubs = int(degrees[high].sum())
    if high_stubs == 0:
        return float("nan")
    high_high = int(np.sum(high[edges[:, 0]] & high[edges[:, 1]]))
    return float((2.0 * high_high / high_stubs) / prevalence)


def component_metrics(high: np.ndarray, edges: np.ndarray) -> tuple[int, float, float]:
    high_indices = np.flatnonzero(high)
    if len(high_indices) == 0:
        return 0, float("nan"), float("nan")
    remap = np.full(len(high), -1, dtype=int)
    remap[high_indices] = np.arange(len(high_indices))
    keep = high[edges[:, 0]] & high[edges[:, 1]]
    high_edges = edges[keep]
    if len(high_edges) == 0:
        return len(high_indices), 1.0 / len(high_indices), 1.0
    rows = remap[high_edges[:, 0]]
    columns = remap[high_edges[:, 1]]
    adjacency = csr_matrix(
        (
            np.ones(len(rows) * 2, dtype=np.uint8),
            (np.r_[rows, columns], np.r_[columns, rows]),
        ),
        shape=(len(high_indices), len(high_indices)),
    )
    n_components, labels = connected_components(adjacency, directed=False)
    sizes = np.bincount(labels, minlength=n_components)
    return int(n_components), float(sizes.max() / len(high_indices)), float(np.median(sizes))


def spatial_edges(frame: pd.DataFrame) -> np.ndarray:
    coordinates = frame[["x_um", "y_um"]].to_numpy(dtype=float)
    pairs = cKDTree(coordinates).query_pairs(
        SPATIAL_EDGE_RADIUS_UM + 1e-8,
        output_type="ndarray",
    )
    return pairs.astype(np.int64)


def raw_expression_neighborhood(
    frame: pd.DataFrame,
    marker_counts: sp.csr_matrix,
) -> tuple[sp.csr_matrix, np.ndarray, np.ndarray, np.ndarray]:
    """Aggregate center-excluded raw RNA over the exact 50-um footprint."""
    coordinates = frame[["x_um", "y_um"]].to_numpy(dtype=float)
    pairs = cKDTree(coordinates).query_pairs(
        EXPRESSION_NEIGHBORHOOD_RADIUS_UM + 1e-8,
        output_type="ndarray",
    )
    if len(pairs) == 0:
        raise RuntimeError("No 50-um raw-expression neighborhood edges")
    rows = np.r_[pairs[:, 0], pairs[:, 1]]
    columns = np.r_[pairs[:, 1], pairs[:, 0]]
    adjacency = csr_matrix(
        (np.ones(len(rows), dtype=np.float32), (rows, columns)),
        shape=(len(frame), len(frame)),
    )
    neighbor_counts = (adjacency @ marker_counts).tocsr()
    neighbor_umi = np.asarray(
        adjacency @ frame["nUMI"].to_numpy(dtype=float)
    ).ravel()
    neighbor_gene_sum = np.asarray(
        adjacency @ frame["nGene"].to_numpy(dtype=float)
    ).ravel()
    neighbor_n = np.asarray(adjacency.sum(axis=1)).ravel()
    return neighbor_counts, neighbor_umi, neighbor_gene_sum, neighbor_n


def load_marker_panels(
    raw: ad.AnnData,
) -> tuple[dict[str, list[str]], dict[str, str], pd.DataFrame]:
    definitions = pd.read_csv(MARKER_DEFINITIONS, sep="\t")
    definitions = definitions.loc[
        definitions["scheme"].eq("conservative")
        & definitions["recommended_for_rctd_pilot"].eq("YES")
    ].copy()
    annotation = pd.read_csv(
        GENE_ANNOTATION,
        sep="\t",
        usecols=["matrix_ensembl_gene_id", "gene_symbol"],
    ).dropna()
    annotation = annotation.drop_duplicates("gene_symbol", keep="first")
    symbol_to_ensembl = dict(
        zip(annotation["gene_symbol"].astype(str), annotation["matrix_ensembl_gene_id"].astype(str))
    )
    raw_genes = set(raw.var_names.astype(str))
    panels: dict[str, list[str]] = {}
    audit_rows = []
    for row in definitions.itertuples(index=False):
        state = str(row.proposed_label)
        requested = [gene.strip() for gene in str(row.primary_markers).split(";") if gene.strip()]
        covered = [
            gene
            for gene in requested
            if symbol_to_ensembl.get(gene, "") in raw_genes
        ]
        panels[state] = covered
        audit_rows.append(
            {
                "reference18_state": state,
                "requested_marker_n": len(requested),
                "covered_marker_n": len(covered),
                "requested_markers": ";".join(requested),
                "covered_markers": ";".join(covered),
                "coverage_fraction": len(covered) / max(1, len(requested)),
            }
        )
    if len(panels) != 18 or any(len(genes) == 0 for genes in panels.values()):
        raise RuntimeError("The six-dataset conservative marker panels did not cover all 18 states")
    return panels, symbol_to_ensembl, pd.DataFrame(audit_rows)


def load_models_and_metadata() -> tuple[
    dict[str, pd.DataFrame],
    dict[str, np.ndarray],
    dict[str, pd.DataFrame],
]:
    assignments = pd.read_csv(ASSIGNMENTS, sep="\t")
    metadata: dict[str, pd.DataFrame] = {}
    usage: dict[str, np.ndarray] = {}
    loadings: dict[str, pd.DataFrame] = {}
    for chip in CHIPS:
        model_path = MODEL_RUN / "models" / f"neighborhood50_{chip}_K06_MODEL_AND_USAGE.npz"
        with np.load(model_path, allow_pickle=True) as handle:
            keys = handle["entity_key"].astype(str).tolist()
            usage[chip] = handle["relative_usage"].astype(float)
        lookup = assignments.loc[assignments["chip_label"].eq(chip)].set_index("entity_key")
        frame = lookup.loc[keys].reset_index()
        if frame["entity_key"].astype(str).tolist() != keys:
            raise RuntimeError(f"{chip}: entity order mismatch")
        if usage[chip].shape != (len(frame), RANK):
            raise RuntimeError(f"{chip}: unexpected usage shape {usage[chip].shape}")
        metadata[chip] = frame
        loadings[chip] = pd.read_csv(
            MODEL_RUN / "models" / f"neighborhood50_{chip}_K06_LOADINGS.tsv",
            sep="\t",
            index_col=0,
        )
    return metadata, usage, loadings


def build_raw_marker_matrix(
    raw: ad.AnnData,
    panels: dict[str, list[str]],
    symbol_to_ensembl: dict[str, str],
) -> tuple[sp.csr_matrix, list[str], dict[str, list[int]]]:
    symbols = sorted({gene for genes in panels.values() for gene in genes})
    ensembl = [symbol_to_ensembl[gene] for gene in symbols]
    positions = raw.var_names.get_indexer(ensembl)
    if np.any(positions < 0):
        raise RuntimeError("Covered marker missing from raw matrix")
    order = np.argsort(positions)
    positions_sorted = positions[order]
    symbols_sorted = [symbols[index] for index in order]
    matrix = raw[:, positions_sorted].to_memory().X
    if not sp.issparse(matrix):
        matrix = sp.csr_matrix(matrix)
    else:
        matrix = matrix.tocsr()
    symbol_index = {gene: index for index, gene in enumerate(symbols_sorted)}
    panel_indices = {
        state: [symbol_index[gene] for gene in genes]
        for state, genes in panels.items()
    }
    return matrix, symbols_sorted, panel_indices


def program_scores(
    counts: sp.csr_matrix,
    total_counts: np.ndarray,
    state_order: list[str],
    panel_indices: dict[str, list[int]],
) -> np.ndarray:
    dense = counts.toarray().astype(float)
    log_cpm = np.log1p(dense / np.maximum(total_counts[:, None], 1.0) * 1e4)
    scores = np.column_stack(
        [log_cpm[:, panel_indices[state]].mean(axis=1) for state in state_order]
    )
    return scores


def matched_low_indices(
    frame: pd.DataFrame,
    values: np.ndarray,
    rng: np.random.Generator,
) -> tuple[np.ndarray, np.ndarray]:
    high_cutoff = float(np.quantile(values, HIGH_QUANTILE))
    low_cutoff = float(np.quantile(values, LOW_CANDIDATE_QUANTILE))
    high = np.flatnonzero(values >= high_cutoff)
    low = np.flatnonzero(values <= low_cutoff)
    if len(high) == 0 or len(low) < len(high):
        return high, np.array([], dtype=int)
    technical = np.column_stack(
        [
            np.log1p(frame["nUMI"].to_numpy(dtype=float)),
            np.log1p(frame["nGene"].to_numpy(dtype=float)),
            frame["official_full_radius_neighbor_n"].to_numpy(dtype=float),
        ]
    )
    means = technical.mean(axis=0)
    scales = technical.std(axis=0)
    scales[scales <= 0] = 1.0
    z = (technical - means[None, :]) / scales[None, :]
    shuffled_low = low.copy()
    rng.shuffle(shuffled_low)
    tree = cKDTree(z[shuffled_low])
    k = min(64, len(shuffled_low))
    _, neighbors = tree.query(z[high], k=k)
    if k == 1:
        neighbors = neighbors[:, None]
    selected = []
    used: set[int] = set()
    for candidate_row in neighbors:
        choice = None
        for candidate in np.atleast_1d(candidate_row):
            local = int(candidate)
            if local not in used:
                choice = local
                break
        if choice is None:
            remaining = np.setdiff1d(np.arange(len(shuffled_low)), np.fromiter(used, dtype=int))
            if len(remaining) == 0:
                break
            choice = int(rng.choice(remaining))
        used.add(choice)
        selected.append(shuffled_low[choice])
    matched = np.asarray(selected, dtype=int)
    return high[: len(matched)], matched


def weighted_marker_log2fc(
    gene_log2fc: np.ndarray,
    factor_loading: np.ndarray,
    state_order: list[str],
    panel_indices: dict[str, list[int]],
) -> float:
    top = np.argsort(factor_loading)[::-1][:3]
    weights = factor_loading[top]
    weights = weights / np.maximum(weights.sum(), EPS)
    value = 0.0
    for weight, state_index in zip(weights, top):
        genes = panel_indices[state_order[state_index]]
        value += float(weight) * float(np.median(gene_log2fc[genes]))
    return float(value)


def main() -> None:
    started = time.time()
    args = parse_args()
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    output = args.output_root.resolve() / f"REFERENCE18_K6_FACTOR_DOWNSTREAM_{timestamp}"
    figures = output / "figures"
    output.mkdir(parents=True, exist_ok=False)
    figures.mkdir()
    rng = np.random.default_rng(SEED)

    metadata, usage, loadings = load_models_and_metadata()
    raw = ad.read_h5ad(RAW_H5AD, backed="r")
    panels, symbol_to_ensembl, panel_audit = load_marker_panels(raw)
    marker_matrix, marker_symbols, panel_indices = build_raw_marker_matrix(
        raw,
        panels,
        symbol_to_ensembl,
    )

    state_order = loadings[CHIPS[0]].columns.astype(str).tolist()
    if state_order != list(panels):
        panels = {state: panels[state] for state in state_order}
        panel_indices = {state: panel_indices[state] for state in state_order}

    factor_definition_rows = []
    patient_metric_rows = []
    expression_rows = []

    for chip in CHIPS:
        frame = metadata[chip].copy()
        raw_indices = raw.obs_names.get_indexer(frame["bin_id"].astype(str))
        if np.any(raw_indices < 0):
            raise RuntimeError(f"{chip}: NMF centers missing from raw H5AD")
        chip_marker_counts = marker_matrix[raw_indices]
        chip_program_scores = program_scores(
            chip_marker_counts,
            frame["nUMI"].to_numpy(dtype=float),
            state_order,
            panel_indices,
        )
        loading_matrix = loadings[chip].to_numpy(dtype=float)

        factor_cutoffs = []
        for factor_index in range(RANK):
            patient_q90 = []
            for _, patient_frame in frame.groupby("patient_id", sort=True):
                patient_q90.append(
                    float(np.quantile(usage[chip][patient_frame.index, factor_index], HIGH_QUANTILE))
                )
            factor_cutoffs.append(float(np.median(patient_q90)))

        for factor_index, factor in enumerate(loadings[chip].index):
            loading = loading_matrix[factor_index]
            order = np.argsort(loading)[::-1]
            factor_definition_rows.append(
                {
                    "chip": chip,
                    "factor": factor,
                    "top_state_1": state_order[order[0]],
                    "loading_1": float(loading[order[0]]),
                    "top_state_2": state_order[order[1]],
                    "loading_2": float(loading[order[1]]),
                    "top_state_3": state_order[order[2]],
                    "loading_3": float(loading[order[2]]),
                    "balanced_patient_q90_cutoff": factor_cutoffs[factor_index],
                }
            )

        for patient, patient_frame in frame.groupby("patient_id", sort=True):
            local_indices = patient_frame.index.to_numpy(dtype=int)
            local_frame = patient_frame.reset_index(drop=True)
            local_usage = usage[chip][local_indices]
            local_program = chip_program_scores[local_indices]
            local_counts = chip_marker_counts[local_indices]
            technical = np.column_stack(
                [
                    np.log1p(local_frame["nUMI"].to_numpy(dtype=float)),
                    np.log1p(local_frame["nGene"].to_numpy(dtype=float)),
                    local_frame["official_full_radius_neighbor_n"].to_numpy(dtype=float),
                ]
            )
            residual_program = np.column_stack(
                [residualize(local_program[:, index], technical) for index in range(len(state_order))]
            )
            (
                neighbor_marker_counts,
                neighbor_umi,
                neighbor_gene_sum,
                raw_neighbor_n,
            ) = raw_expression_neighborhood(local_frame, local_counts)
            neighbor_program = program_scores(
                neighbor_marker_counts,
                neighbor_umi,
                state_order,
                panel_indices,
            )
            neighbor_technical = np.column_stack(
                [
                    np.log1p(neighbor_umi),
                    np.log1p(neighbor_gene_sum),
                    raw_neighbor_n,
                ]
            )
            residual_neighbor_program = np.column_stack(
                [
                    residualize(neighbor_program[:, index], neighbor_technical)
                    for index in range(len(state_order))
                ]
            )
            edges = spatial_edges(local_frame)

            for factor_index, factor in enumerate(loadings[chip].index):
                values = local_usage[:, factor_index]
                loading = loading_matrix[factor_index]
                top = np.argsort(loading)[::-1][:3]
                top_weights = loading[top] / np.maximum(loading[top].sum(), EPS)
                composite_score = residual_program[:, top] @ top_weights
                top_state_score = residual_program[:, top[0]]
                neighbor_composite_score = residual_neighbor_program[:, top] @ top_weights
                neighbor_top_state_score = residual_neighbor_program[:, top[0]]
                usage_residual = residualize(values, technical)
                high = values >= factor_cutoffs[factor_index]
                component_n, largest_component, median_component = component_metrics(high, edges)
                patient_metric_rows.append(
                    {
                        "chip": chip,
                        "patient_id": str(patient),
                        "disease": str(local_frame["disease"].iloc[0]),
                        "factor": factor,
                        "n_neighborhoods": len(local_frame),
                        "mean_usage": float(values.mean()),
                        "median_usage": float(np.median(values)),
                        "usage_q90": float(np.quantile(values, HIGH_QUANTILE)),
                        "balanced_q90_cutoff": factor_cutoffs[factor_index],
                        "high_usage_fraction": float(high.mean()),
                        "dominant_factor_fraction": float(
                            np.mean(np.argmax(local_usage, axis=1) == factor_index)
                        ),
                        "usage_umi_spearman": safe_spearman(values, technical[:, 0]),
                        "usage_n_genes_spearman": safe_spearman(values, technical[:, 1]),
                        "usage_neighbor_n_spearman": safe_spearman(values, technical[:, 2]),
                        "raw_usage_moran_i": moran_i(values, edges),
                        "technical_residual_usage_moran_i": moran_i(usage_residual, edges),
                        "high_usage_local_enrichment": local_high_enrichment(high, edges),
                        "high_usage_component_n": component_n,
                        "largest_high_component_fraction": largest_component,
                        "median_high_component_size": median_component,
                    }
                )

                high_indices, control_indices = matched_low_indices(local_frame, values, rng)
                if len(control_indices) == 0:
                    marker_log2fc = np.full(len(marker_symbols), np.nan)
                else:
                    high_sum = np.asarray(local_counts[high_indices].sum(axis=0)).ravel()
                    control_sum = np.asarray(local_counts[control_indices].sum(axis=0)).ravel()
                    high_library = float(local_frame.iloc[high_indices]["nUMI"].sum())
                    control_library = float(local_frame.iloc[control_indices]["nUMI"].sum())
                    high_cpm = high_sum / max(high_library, 1.0) * 1e6
                    control_cpm = control_sum / max(control_library, 1.0) * 1e6
                    marker_log2fc = np.log2((high_cpm + 0.5) / (control_cpm + 0.5))
                top_state = state_order[top[0]]
                expression_rows.append(
                    {
                        "chip": chip,
                        "patient_id": str(patient),
                        "factor": factor,
                        "top_state": top_state,
                        "top3_states": ";".join(state_order[index] for index in top),
                        "top3_composite_marker_residual_spearman": safe_spearman(
                            values,
                            composite_score,
                        ),
                        "top_state_marker_residual_spearman": safe_spearman(
                            values,
                            top_state_score,
                        ),
                        "neighbor_top3_marker_residual_spearman": safe_spearman(
                            values,
                            neighbor_composite_score,
                        ),
                        "neighbor_top_state_marker_residual_spearman": safe_spearman(
                            values,
                            neighbor_top_state_score,
                        ),
                        "matched_high_bin_n": len(high_indices),
                        "matched_control_bin_n": len(control_indices),
                        "top3_composite_marker_log2fc": weighted_marker_log2fc(
                            marker_log2fc,
                            loading,
                            state_order,
                            panel_indices,
                        ),
                        "top_state_marker_median_log2fc": float(
                            np.nanmedian(marker_log2fc[panel_indices[top_state]])
                        ),
                    }
                )

    factor_definitions = pd.DataFrame(factor_definition_rows)
    patient_metrics = pd.DataFrame(patient_metric_rows)
    expression = pd.DataFrame(expression_rows)

    expression_summary_rows = []
    patient_summary_rows = []
    for (chip, factor), factor_frame in expression.groupby(["chip", "factor"], sort=True):
        rho = factor_frame["neighbor_top3_marker_residual_spearman"]
        center_rho = factor_frame["top3_composite_marker_residual_spearman"]
        log2fc = factor_frame["top3_composite_marker_log2fc"]
        positive_rho_n = int((rho > 0).sum())
        positive_log2fc_n = int((log2fc > 0).sum())
        median_rho = float(rho.median())
        median_log2fc = float(log2fc.median())
        neighbor_strong = median_rho >= 0.10 and positive_rho_n >= 5
        neighbor_support = median_rho > 0 and positive_rho_n >= 4
        focal_strong = median_log2fc >= 0.25 and positive_log2fc_n >= 5
        focal_support = median_log2fc > 0 and positive_log2fc_n >= 4
        if neighbor_strong and focal_strong:
            evidence = "Concordant strong raw-expression support"
        elif neighbor_support and focal_support:
            evidence = "Concordant raw-expression support"
        elif neighbor_support:
            evidence = "Neighborhood raw-expression support"
        elif focal_support:
            evidence = "Focal center-bin enrichment support"
        else:
            evidence = "Raw-expression support unresolved"
        expression_summary_rows.append(
            {
                "chip": chip,
                "factor": factor,
                "median_neighbor_composite_marker_residual_spearman": median_rho,
                "positive_neighbor_composite_marker_rho_patient_n": positive_rho_n,
                "median_center_composite_marker_residual_spearman": float(center_rho.median()),
                "positive_center_composite_marker_rho_patient_n": int((center_rho > 0).sum()),
                "median_composite_marker_log2fc": median_log2fc,
                "positive_composite_marker_log2fc_patient_n": positive_log2fc_n,
                "median_neighbor_top_state_marker_residual_spearman": float(
                    factor_frame["neighbor_top_state_marker_residual_spearman"].median()
                ),
                "median_center_top_state_marker_residual_spearman": float(
                    factor_frame["top_state_marker_residual_spearman"].median()
                ),
                "median_top_state_marker_log2fc": float(
                    factor_frame["top_state_marker_median_log2fc"].median()
                ),
                "raw_expression_evidence": evidence,
            }
        )
    for (chip, factor), factor_frame in patient_metrics.groupby(["chip", "factor"], sort=True):
        mean_usage = factor_frame["mean_usage"]
        patient_summary_rows.append(
            {
                "chip": chip,
                "factor": factor,
                "patient_n": factor_frame["patient_id"].nunique(),
                "patient_high_usage_coverage_n": int((factor_frame["high_usage_fraction"] > 0).sum()),
                "max_patient_mean_usage_contribution": float(mean_usage.max() / mean_usage.sum()),
                "median_mean_usage": float(mean_usage.median()),
                "median_high_usage_fraction": float(factor_frame["high_usage_fraction"].median()),
                "median_largest_high_component_fraction": float(
                    factor_frame["largest_high_component_fraction"].median()
                ),
                "median_technical_residual_usage_moran_i": float(
                    factor_frame["technical_residual_usage_moran_i"].median()
                ),
                "median_high_usage_local_enrichment": float(
                    factor_frame["high_usage_local_enrichment"].median()
                ),
                "median_usage_umi_spearman": float(factor_frame["usage_umi_spearman"].median()),
                "median_usage_n_genes_spearman": float(
                    factor_frame["usage_n_genes_spearman"].median()
                ),
                "median_usage_neighbor_n_spearman": float(
                    factor_frame["usage_neighbor_n_spearman"].median()
                ),
            }
        )

    expression_summary = pd.DataFrame(expression_summary_rows)
    patient_summary = pd.DataFrame(patient_summary_rows)
    prior = pd.read_csv(PRIOR_AUDIT, sep="\t")
    prior = prior.loc[prior["rank"].eq(RANK)].copy()
    keep_prior = [
        "chip",
        "factor",
        "final_evidence_level",
        "retain_as_future_spatial_variable",
        "retention_basis",
        "model_restart_min_cosine",
        "model_lopo_min_cosine",
    ]
    final = (
        factor_definitions.merge(patient_summary, on=["chip", "factor"], how="left")
        .merge(expression_summary, on=["chip", "factor"], how="left")
        .merge(prior[keep_prior], on=["chip", "factor"], how="left")
    )
    final["recommended_downstream_use"] = np.where(
        final["retain_as_future_spatial_variable"].fillna(False),
        "Retain continuous usage; carry evidence label",
        "Exploratory only pending additional evidence",
    )

    factor_definitions.to_csv(output / "01_FACTOR_DEFINITIONS.tsv", sep="\t", index=False)
    patient_metrics.to_csv(output / "02_PATIENT_FACTOR_METRICS.tsv", sep="\t", index=False)
    expression.to_csv(
        output / "03_PATIENT_FACTOR_RAW_EXPRESSION_VALIDATION.tsv",
        sep="\t",
        index=False,
    )
    expression_summary.to_csv(
        output / "04_FACTOR_RAW_EXPRESSION_SUMMARY.tsv",
        sep="\t",
        index=False,
    )
    panel_audit.to_csv(output / "05_MARKER_PANEL_AUDIT.tsv", sep="\t", index=False)
    final.to_csv(output / "06_FINAL_FACTOR_ANALYSIS_TABLE.tsv", sep="\t", index=False)

    labels = final.apply(
        lambda row: f"{row['chip']} {row['factor']} | {row['top_state_1']}",
        axis=1,
    )
    figure_frame = final.set_index(labels)
    fig, axes = plt.subplots(1, 2, figsize=(15, 9), constrained_layout=True)
    sns.heatmap(
        figure_frame[["median_neighbor_composite_marker_residual_spearman"]],
        cmap="vlag",
        center=0,
        annot=True,
        fmt=".2f",
        cbar_kws={"label": "Median within-patient Spearman rho"},
        ax=axes[0],
    )
    axes[0].set_title("Raw RNA marker-program association")
    axes[0].set_xlabel("")
    axes[0].set_ylabel("Factor | dominant Reference18 state")
    sns.heatmap(
        figure_frame[["median_composite_marker_log2fc"]],
        cmap="vlag",
        center=0,
        annot=True,
        fmt=".2f",
        cbar_kws={"label": "Median matched high-vs-low log2 fold change"},
        ax=axes[1],
    )
    axes[1].set_title("Matched raw-count pseudobulk support")
    axes[1].set_xlabel("")
    axes[1].set_ylabel("")
    fig.savefig(figures / "FIG01_FACTOR_RAW_EXPRESSION_SUPPORT.png", dpi=220)
    plt.close(fig)

    expression_plot = expression.copy()
    expression_plot["factor_label"] = expression_plot.apply(
        lambda row: f"{row['chip']} {row['factor']}",
        axis=1,
    )
    rho_matrix = expression_plot.pivot(
        index="factor_label",
        columns="patient_id",
        values="neighbor_top3_marker_residual_spearman",
    )
    ordered_labels = [f"{row.chip} {row.factor}" for row in factor_definitions.itertuples(index=False)]
    rho_matrix = rho_matrix.reindex(ordered_labels)
    fig, axis = plt.subplots(figsize=(16, 8), constrained_layout=True)
    sns.heatmap(
        rho_matrix,
        cmap="vlag",
        center=0,
        vmin=-0.5,
        vmax=0.5,
        annot=True,
        fmt=".2f",
        cbar_kws={"label": "Within-patient residual Spearman rho"},
        ax=axis,
    )
    axis.set_title("K=6 factor association with expected raw-RNA marker programs")
    axis.set_xlabel("Patient")
    axis.set_ylabel("Factor")
    fig.savefig(figures / "FIG02_EXPRESSION_SUPPORT_BY_PATIENT.png", dpi=220)
    plt.close(fig)

    metric_columns = [
        "median_mean_usage",
        "max_patient_mean_usage_contribution",
        "median_largest_high_component_fraction",
        "median_technical_residual_usage_moran_i",
        "median_high_usage_local_enrichment",
    ]
    metric_labels = [
        "Median usage",
        "Max patient\ncontribution",
        "Largest high-usage\ncomponent",
        "Residual\nMoran's I",
        "Local high-usage\nenrichment",
    ]
    metric_frame = final.set_index(labels)[metric_columns].copy()
    metric_frame.columns = metric_labels
    metric_z = (metric_frame - metric_frame.mean(axis=0)) / metric_frame.std(axis=0).replace(0, 1)
    fig, axis = plt.subplots(figsize=(12, 9), constrained_layout=True)
    sns.heatmap(
        metric_z,
        cmap="vlag",
        center=0,
        annot=metric_frame,
        fmt=".2f",
        cbar_kws={"label": "Column z-score (numbers show original values)"},
        ax=axis,
    )
    axis.set_title("Patient-level abundance and spatial metrics for frozen K=6 factors")
    axis.set_xlabel("")
    axis.set_ylabel("Factor | dominant Reference18 state")
    axis.tick_params(axis="x", labelrotation=0, labelsize=9)
    fig.savefig(figures / "FIG03_FACTOR_PATIENT_LEVEL_METRICS.png", dpi=220)
    plt.close(fig)

    tech_columns = [
        "median_usage_umi_spearman",
        "median_usage_n_genes_spearman",
        "median_usage_neighbor_n_spearman",
    ]
    tech = final.set_index(labels)[tech_columns]
    tech.columns = ["UMI", "n_genes", "Effective neighbors"]
    fig, axis = plt.subplots(figsize=(10, 9), constrained_layout=True)
    sns.heatmap(
        tech,
        cmap="vlag",
        center=0,
        norm=TwoSlopeNorm(vmin=-0.8, vcenter=0, vmax=0.8),
        annot=True,
        fmt=".2f",
        cbar_kws={"label": "Median within-patient Spearman rho"},
        ax=axis,
    )
    axis.set_title("Technical associations of frozen K=6 factor usage")
    axis.set_xlabel("Technical variable")
    axis.set_ylabel("Factor | dominant Reference18 state")
    fig.savefig(figures / "FIG04_FACTOR_TECHNICAL_ASSOCIATIONS.png", dpi=220)
    plt.close(fig)

    evidence_counts = final["raw_expression_evidence"].value_counts().to_dict()
    retained_n = int(final["retain_as_future_spatial_variable"].fillna(False).sum())
    summary_lines = [
        "# Reference18 K=6 factor 原始表达验证与患者层面指标",
        "",
        "## 范围",
        "",
        "- 复用冻结的逐芯片 Neighborhood NMF K=6；未重新训练 NMF 或 RCTD。",
        "- 共分析 3 张芯片、21 位患者、18 个 factor 和 126,199 个可评价 neighborhood。",
        "- 未开展疾病差异检验，也未读取或修改冻结的 238 个 lymphoid aggregates。",
        "",
        "## 原始表达验证",
        "",
        "- marker 来自六数据集 Reference18 conservative reference 的 primary markers。",
        "- 主验证将原始 RNA 按与 NMF 相同的中心排除 50 μm 邻域汇总；中心 Bin50 原始 RNA 保留为独立的局灶证据。",
        "- 每位患者内对邻域和中心 marker program 分别残差化 UMI、n_genes 和有效邻居数。",
        "- 另比较患者内 top-10% usage Bin 与技术匹配的低-usage Bin 原始计数 pseudobulk。",
        "- 证据标签仅描述 RNA 支持程度，不作为删除 factor 的硬门槛。",
        "",
        "## 患者层面指标",
        "",
        "- 保存连续 mean/median/q90 usage、固定 chip 内患者平衡阈值以上的面积比例和 dominant-factor 比例。",
        "- 保存高 usage 连通组分、最大组分比例、局部富集、原始及技术残差 Moran's I。",
        "- 保存患者内 usage 与 UMI、n_genes、有效邻居数的相关性。",
        "",
        "## 自动汇总",
        "",
        f"- Concordant strong raw-expression support：{evidence_counts.get('Concordant strong raw-expression support', 0)} 个 factor。",
        f"- Concordant raw-expression support：{evidence_counts.get('Concordant raw-expression support', 0)} 个 factor。",
        f"- Neighborhood raw-expression support：{evidence_counts.get('Neighborhood raw-expression support', 0)} 个 factor。",
        f"- Focal center-bin enrichment support：{evidence_counts.get('Focal center-bin enrichment support', 0)} 个 factor。",
        f"- Raw-expression support unresolved：{evidence_counts.get('Raw-expression support unresolved', 0)} 个 factor。",
        f"- 既往综合审计建议保留为连续空间变量：{retained_n}/18 个 factor。",
        "",
        "详细的逐 factor 结论必须结合 `06_FINAL_FACTOR_ANALYSIS_TABLE.tsv`，不能仅凭单个阈值命名 biological niche。",
    ]
    (output / "07_METHOD_AND_RESULTS_ZH.md").write_text("\n".join(summary_lines), encoding="utf-8")

    input_paths = [
        ASSIGNMENTS,
        RAW_H5AD,
        GENE_ANNOTATION,
        MARKER_DEFINITIONS,
        PRIOR_AUDIT,
    ] + [
        MODEL_RUN / "models" / f"neighborhood50_{chip}_K06_MODEL_AND_USAGE.npz"
        for chip in CHIPS
    ] + [
        MODEL_RUN / "models" / f"neighborhood50_{chip}_K06_LOADINGS.tsv"
        for chip in CHIPS
    ]
    pd.DataFrame(
        {"path": [str(path) for path in input_paths], "sha256": [sha256(path) for path in input_paths]}
    ).to_csv(output / "00_INPUT_MANIFEST.tsv", sep="\t", index=False)

    session = {
        "created_at": datetime.now().isoformat(),
        "elapsed_seconds": time.time() - started,
        "python": sys.version,
        "platform": platform.platform(),
        "seed": SEED,
        "high_usage_quantile": HIGH_QUANTILE,
        "low_usage_candidate_quantile": LOW_CANDIDATE_QUANTILE,
        "spatial_edge_radius_um": SPATIAL_EDGE_RADIUS_UM,
        "expression_neighborhood_radius_um": EXPRESSION_NEIGHBORHOOD_RADIUS_UM,
        "raw_h5ad_shape": list(raw.shape),
    }
    (output / "SESSION_INFO.json").write_text(json.dumps(session, indent=2), encoding="utf-8")
    (output / "RUN_LOG.txt").write_text(
        "\n".join(
            [
                "Script: scripts/1101_validate_k6_factor_expression_and_patient_metrics.py",
                "Command: /home/bozhu/miniconda3/envs/banksy_stereo/bin/python scripts/1101_validate_k6_factor_expression_and_patient_metrics.py --output-root results/reanalysis",
                "NMF refit: no",
                "RCTD refit: no",
                "Disease comparison: no",
                "Frozen lymphoid aggregates accessed: no",
                f"Output: {output}",
                "Status: PASS",
            ]
        ),
        encoding="utf-8",
    )
    (output / "COMPLETE.ok").write_text("PASS\n", encoding="utf-8")
    print(output)


if __name__ == "__main__":
    main()
