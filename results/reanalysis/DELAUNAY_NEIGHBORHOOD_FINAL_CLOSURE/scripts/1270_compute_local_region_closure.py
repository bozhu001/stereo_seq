"""Compute the final, fixed-label Delaunay local-region closure tables.

The script reuses frozen unified PN1-PN8 labels and saved connected-component
IDs. It does not fit clustering, select K, rerun RCTD/NMF, or perform spatial
co-occurrence. New work is limited to deterministic region summaries and
interior/boundary/surrounding descriptions.
"""

from __future__ import annotations

import hashlib
import importlib.util
import json
import logging
import math
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import scipy.sparse as sp
from pandas.api.types import is_numeric_dtype
from scipy import stats
from scipy.spatial import cKDTree
from sklearn.metrics import adjusted_rand_score, normalized_mutual_info_score


ROOT = Path(__file__).resolve().parents[4]
REAN = ROOT / "results/reanalysis"
OUT = REAN / "DELAUNAY_NEIGHBORHOOD_FINAL_CLOSURE"
GRAPH_ROOT = REAN / "REFERENCE18_SPATIAL_GRAPH_NEIGHBORHOODS_V1"
FORMAL = GRAPH_ROOT / "ALL21_EXTERNAL_PROJECTION/FORMAL_RUN_MODEL_RECOVERED_EXACT"
LABEL_PATH = FORMAL / "02_ALL21_PROJECTION_LABELS.tsv.gz"
COMPONENT_PATH = FORMAL / "06_ALL_CONNECTED_COMPONENTS.tsv.gz"
MARKER_PATH = GRAPH_ROOT / "ALL21_PN_MARKER_VALIDATION/04_ALL_PN_MARKER_MODULE_SCORES.tsv.gz"
MARKER_INTERPRETATION = GRAPH_ROOT / "ALL21_PN_MARKER_VALIDATION/08_PN_BIOLOGICAL_INTERPRETATION.tsv"
SSDNA_INDEX = GRAPH_ROOT / "SSDNA_OVERLAY_REVIEW/00_INPUT_FILE_INDEX.tsv"
SSDNA_AUDIT = GRAPH_ROOT / "SSDNA_OVERLAY_COORDINATE_AUDIT/04_PATIENT_ALIGNMENT_QC.tsv"
BATCH_ROOT = REAN / "REFERENCE18_SPATIAL_GRAPH_BATCH_STRATEGY_AND_K_SELECTION"
K_REVIEW = REAN / "REFERENCE18_SPATIAL_GRAPH_K2_15_STABILITY_REVIEW"
NMF_DIR = REAN / "REFERENCE18_NEIGHBORHOOD_NMF_K7_K12_EXTENSION_20261015_attempt5/models"
NMF_MAP = REAN / "K7_FROZEN_UPSTREAM_BASELINE/05_K7_FACTOR_FAMILY_MAPPING.tsv"
PILOT_SCRIPT = GRAPH_ROOT / "scripts/1233_reference18_spatial_graph_neighborhoods_pilot.py"
GRAPH_SCRIPT = REAN / "REFERENCE18_SPATIAL_GRAPH_CHIP_SPECIFIC_K_SELECTION/scripts/1250_chip_specific_k_selection.py"
FROZEN_DIR = REAN / "REFERENCE18_UMI50_VOLTRON_PER_CHIP_NICHE_20261005_153036"

PN_ORDER = [f"PN{i}" for i in range(1, 9)]
CHIPS = ["K8", "J2", "L3"]
SEED = 20261026
BIN_AREA_UM2 = 625.0
BIN_HALF_UM = 12.5
FULL_REGION_MIN_BINS = 10
SMALL_REGION_MIN_BINS = 5


def ensure_dirs() -> None:
    for name in ("tables", "figures", "scripts", "logs", "report"):
        (OUT / name).mkdir(parents=True, exist_ok=True)


def logger() -> logging.Logger:
    result = logging.getLogger("delaunay_final_closure")
    result.setLevel(logging.INFO)
    result.handlers.clear()
    formatter = logging.Formatter("%(asctime)s\t%(levelname)s\t%(message)s")
    for handler in (
        logging.FileHandler(OUT / "logs/1270_compute_local_region_closure.log", encoding="utf-8"),
        logging.StreamHandler(sys.stdout),
    ):
        handler.setFormatter(formatter)
        result.addHandler(handler)
    return result


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def import_module(path: Path, name: str):
    spec = importlib.util.spec_from_file_location(name, path)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"Cannot import {path}")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def bh(values: np.ndarray) -> np.ndarray:
    values = np.asarray(values, float)
    result = np.full(values.shape, np.nan)
    valid = np.isfinite(values)
    if not valid.any():
        return result
    p = values[valid]
    order = np.argsort(p)
    ranked = p[order]
    adjusted = ranked * len(ranked) / np.arange(1, len(ranked) + 1)
    adjusted = np.minimum.accumulate(adjusted[::-1])[::-1]
    restored = np.empty_like(adjusted)
    restored[order] = np.minimum(adjusted, 1.0)
    result[valid] = restored
    return result


def load_core(log: logging.Logger):
    pilot = import_module(PILOT_SCRIPT, "closure_pilot")
    graph_module = import_module(GRAPH_SCRIPT, "closure_graph")
    data, weight_columns, _ = pilot.load_frozen_inputs(ROOT)
    data = data.reset_index(drop=True).copy()
    graph_pairs, graph_qc = graph_module.build_graphs(data, pilot)
    graph_pairs = np.asarray(graph_pairs, dtype=np.int64)
    labels = pd.read_csv(LABEL_PATH, sep="\t")
    if not labels.entity_key.is_unique or not data.entity_key.is_unique:
        raise RuntimeError("Entity keys are not unique")
    order = data[["entity_key"]].merge(
        labels, on="entity_key", how="left", validate="one_to_one", sort=False,
    )
    if order.patient_id.isna().any():
        raise RuntimeError("Formal labels do not cover all frozen bins")
    labels = order
    for column in ("patient_id", "chip_label", "disease", "roi_id", "tissue_component_id", "x_um", "y_um", "nUMI", "nGene"):
        source = data[column].to_numpy()
        target = labels[column].to_numpy()
        if is_numeric_dtype(data[column]):
            if not np.allclose(source.astype(float), target.astype(float), equal_nan=True):
                raise RuntimeError(f"Mismatch in {column}")
        elif not np.array_equal(source.astype(str), target.astype(str)):
            raise RuntimeError(f"Mismatch in {column}")
    labels["global_row"] = np.arange(len(labels))
    degree = np.bincount(graph_pairs.ravel(), minlength=len(data))
    if not np.array_equal(degree, labels.delaunay_degree.to_numpy(int)):
        raise RuntimeError("Reconstructed graph degree differs from frozen labels")
    adjacency = sp.coo_matrix(
        (
            np.ones(2 * len(graph_pairs), dtype=np.int8),
            (
                np.r_[graph_pairs[:, 0], graph_pairs[:, 1]],
                np.r_[graph_pairs[:, 1], graph_pairs[:, 0]],
            ),
        ),
        shape=(len(data), len(data)),
    ).tocsr()
    log.info("Loaded %d bins, %d graph edges, %d states", len(data), len(graph_pairs), len(weight_columns))
    return data, labels, weight_columns, graph_pairs, graph_qc, adjacency


def input_audit(data, labels, weights, pairs, graph_qc) -> pd.DataFrame:
    components = pd.read_csv(COMPONENT_PATH, sep="\t")
    rows = []
    sources = [
        ("frozen_reference18_K8", FROZEN_DIR / "K8/01_FROZEN_REFERENCE18_INPUT.tsv.gz", len(data[data.chip_label.eq("K8")])),
        ("frozen_reference18_J2", FROZEN_DIR / "J2/01_FROZEN_REFERENCE18_INPUT.tsv.gz", len(data[data.chip_label.eq("J2")])),
        ("frozen_reference18_L3", FROZEN_DIR / "L3/01_FROZEN_REFERENCE18_INPUT.tsv.gz", len(data[data.chip_label.eq("L3")])),
        ("formal_unified_pn_labels", LABEL_PATH, len(labels)),
        ("formal_connected_components", COMPONENT_PATH, len(components)),
        ("marker_module_scores", MARKER_PATH, 126449 * 18),
        ("ssdna_input_index", SSDNA_INDEX, 9),
        ("ssdna_alignment_audit", SSDNA_AUDIT, 9),
        ("frozen_nmf_factor_mapping", NMF_MAP, 21),
        ("delaunay_pilot_script", PILOT_SCRIPT, np.nan),
        ("delaunay_graph_script", GRAPH_SCRIPT, np.nan),
    ]
    for name, path, expected in sources:
        rows.append(
            {
                "input_name": name, "path": str(path.relative_to(ROOT)),
                "exists": path.exists(), "bytes": path.stat().st_size if path.exists() else 0,
                "sha256": sha256(path) if path.exists() else "MISSING",
                "expected_records": expected, "status": "PASS" if path.exists() and path.stat().st_size > 0 else "FAIL",
                "read_only_use": True,
            }
        )
    checks = [
        ("qc_bin_count", len(data), 126449),
        ("patient_count", labels.patient_id.nunique(), 21),
        ("chip_count", labels.chip_label.nunique(), 3),
        ("reference_state_count", len(weights), 18),
        ("graph_edge_count", len(pairs), 357959),
        ("graph_max_edge_um", float(graph_qc.observed_max_edge_um.max()), 25 * math.sqrt(2)),
        ("formal_evaluable_node_count", int(labels.classification_status.ne("NOT_EVALUABLE").sum()), 125471),
        ("saved_connected_component_count", len(components), 29750),
    ]
    for name, observed, expected in checks:
        passed = np.isclose(observed, expected) if isinstance(observed, float) else observed == expected
        rows.append(
            {
                "input_name": name, "path": "computed audit", "exists": True,
                "bytes": np.nan, "sha256": "NOT_APPLICABLE", "expected_records": expected,
                "observed_records": observed, "status": "PASS" if passed else "FAIL", "read_only_use": True,
            }
        )
    result = pd.DataFrame(rows)
    result["analysis_seed"] = SEED
    result["bin_area_um2"] = BIN_AREA_UM2
    result["maximum_delaunay_edge_um"] = 25 * math.sqrt(2)
    result.to_csv(OUT / "tables/01_INPUT_AND_VERSION_AUDIT.tsv", sep="\t", index=False)
    return result


def completeness() -> pd.DataFrame:
    records = [
        ("A", "Delaunay graph construction", GRAPH_SCRIPT, [FORMAL / "14_ALL21_GRAPH_QC.tsv"], 21, "COMPLETE", True),
        ("B", "Graph QC", GRAPH_SCRIPT, [FORMAL / "14_ALL21_GRAPH_QC.tsv"], 21, "COMPLETE", True),
        ("C", "Unified PN1-PN8", FORMAL / "scripts/1241_project_pn_all21.py", [LABEL_PATH], 21, "COMPLETE", True),
        ("D", "Chip-specific K=2-15 stability", K_REVIEW / "scripts/1260_compute_k2_15_stability.py", [K_REVIEW / "tables/02_FULL_K2_15_METRICS_BY_CHIP.tsv"], 21, "EXPLORATORY", True),
        ("E", "PN Marker Gene Concordance", GRAPH_ROOT / "ALL21_PN_MARKER_VALIDATION/scripts/1242_validate_all21_pn_markers.py", [MARKER_INTERPRETATION], 21, "COMPLETE", True),
        ("F", "Connected Components", FORMAL / "scripts/1241_project_pn_all21.py", [COMPONENT_PATH], 21, "COMPLETE", True),
        ("G", "ssDNA Overlay", GRAPH_ROOT / "SSDNA_OVERLAY_REVIEW/scripts/1234_map_graph_neighborhoods_to_ssdna.py", [SSDNA_INDEX, SSDNA_AUDIT], 9, "COMPLETE", True),
        ("H", "Region Interior / Boundary / Surrounding", OUT / "scripts/1270_compute_local_region_closure.py", [OUT / "tables/06_REGION_INTERNAL_COMPOSITION.tsv", OUT / "tables/07_REGION_BOUNDARY_AND_SURROUNDINGS.tsv"], 21, "COMPLETE", False),
        ("I", "Cross-chip Neighborhood Correspondence", BATCH_ROOT / "scripts/1252_phase3_phase4_chip_specific.py", [BATCH_ROOT / "tables/10B_CROSS_CHIP_CORRESPONDENCE.tsv"], 21, "EXPLORATORY", True),
        ("J", "Patient-level Disease Descriptive Analysis", OUT / "scripts/1270_compute_local_region_closure.py", [OUT / "tables/10_CHIP_AWARE_PATIENT_SUMMARY.tsv"], 21, "EXPLORATORY", False),
    ]
    rows = []
    for code, analysis, script, outputs, patients, status, reused in records:
        exists = script.exists() and all(path.exists() for path in outputs)
        rows.append(
            {
                "item": code, "analysis": analysis,
                "script_path": str(script.relative_to(ROOT)),
                "output_paths": "|".join(str(path.relative_to(ROOT)) for path in outputs),
                "completed_patient_n": patients if exists else 0,
                "version": "current formal/frozen chain",
                "completion_status": status if exists else "NOT_DONE",
                "directly_reused": reused,
                "evidence_exists": exists,
            }
        )
    result = pd.DataFrame(rows)
    result.to_csv(OUT / "tables/02_EXISTING_ANALYSIS_COMPLETENESS.tsv", sep="\t", index=False)
    return result


def load_marker_arrays(entity_order: pd.Series, state_names: list[str], log: logging.Logger):
    usecols = ["entity_key", "reference18_state", "fixed_depth_residual_score", "heldout_depth_residual_score"]
    marker = pd.read_csv(MARKER_PATH, sep="\t", usecols=usecols)
    fixed = marker.pivot(index="entity_key", columns="reference18_state", values="fixed_depth_residual_score")
    heldout = marker.pivot(index="entity_key", columns="reference18_state", values="heldout_depth_residual_score")
    fixed = fixed.reindex(index=entity_order, columns=state_names)
    heldout = heldout.reindex(index=entity_order, columns=state_names)
    log.info("Loaded marker arrays %s and %s", fixed.shape, heldout.shape)
    return fixed.to_numpy(float), heldout.to_numpy(float)


def size_stratum(bin_n: pd.Series | np.ndarray) -> np.ndarray:
    values = np.asarray(bin_n)
    return np.select(
        [values == 1, (values >= 2) & (values <= 4), (values >= 5) & (values <= 9), values >= 10],
        ["Singleton", "2-4 Bin50", "5-9 Bin50", ">=10 Bin50"],
        default="UNKNOWN",
    )


def region_table(data, labels, weight_columns, ssdna_patients: set[str]) -> pd.DataFrame:
    saved = pd.read_csv(COMPONENT_PATH, sep="\t")
    eligible = labels[labels.component_id.notna()].copy()
    grouped = eligible.groupby("component_id", sort=False)
    basic = grouped.agg(
        centroid_x_um=("x_um", "mean"), centroid_y_um=("y_um", "mean"),
        median_umi=("nUMI", "median"), median_ngenes=("nGene", "median"),
        mean_graph_degree=("delaunay_degree", "mean"),
        low_confidence_bin_fraction=("classification_status", lambda x: float(x.eq("LOW_CONFIDENCE").mean())),
        ood_unassigned_bin_fraction=("classification_status", lambda x: float(x.eq("OOD_UNASSIGNED").mean())),
        reliable_bin_fraction=("classification_status", lambda x: float(x.eq("RELIABLE").mean())),
    ).reset_index()
    means = pd.DataFrame(
        data.loc[eligible.global_row, weight_columns].to_numpy(float),
        columns=[f"mean_{column}" for column in weight_columns], index=eligible.index,
    )
    means["component_id"] = eligible.component_id.to_numpy()
    means = means.groupby("component_id", sort=False).mean().reset_index()
    result = saved.merge(basic, on="component_id", how="left", validate="one_to_one")
    result = result.merge(means, on="component_id", how="left", validate="one_to_one")
    result = result.rename(columns={"component_id": "local_region_id", "predicted_pn": "neighborhood_label"})
    result["size_stratum"] = size_stratum(result.bin_n)
    result["ssdna_alignment_status"] = np.where(
        result.patient_id.isin(ssdna_patients), "VALIDATED_ALIGNMENT_AVAILABLE", "SSDNA_NOT_AVAILABLE"
    )
    result["region_interpretation"] = "graph-connected same-label region; not automatically an aggregate"
    result["area_definition"] = "Bin50 count x 625 um2; holes not filled"
    result.to_csv(OUT / "tables/04_ALL_LOCAL_CONNECTED_REGIONS.tsv.gz", sep="\t", index=False, compression="gzip")
    strata = (
        result.groupby(["patient_id", "chip_label", "disease", "neighborhood_label", "size_stratum"], observed=True, as_index=False)
        .agg(region_n=("local_region_id", "size"), bin_n=("bin_n", "sum"), total_area_um2=("approximate_area_um2", "sum"), median_region_bin_n=("bin_n", "median"), maximum_region_bin_n=("bin_n", "max"))
    )
    strata.to_csv(OUT / "tables/05_LOCAL_REGION_SIZE_STRATA.tsv", sep="\t", index=False)
    return result


def boundary_nodes(region_indices: np.ndarray, adjacency: sp.csr_matrix) -> tuple[np.ndarray, np.ndarray]:
    membership = np.zeros(adjacency.shape[0], dtype=bool)
    membership[region_indices] = True
    boundary = []
    neighbor_outside = []
    for node in region_indices:
        neighbors = adjacency.indices[adjacency.indptr[node] : adjacency.indptr[node + 1]]
        outside = neighbors[~membership[neighbors]]
        if len(outside):
            boundary.append(node)
            neighbor_outside.extend(outside.tolist())
    return np.asarray(boundary, dtype=int), np.unique(np.asarray(neighbor_outside, dtype=int))


def zone_summary_rows(
    region_row,
    zone: str,
    indices: np.ndarray,
    data: pd.DataFrame,
    labels: pd.DataFrame,
    weight_columns: list[str],
    fixed_marker: np.ndarray,
    heldout_marker: np.ndarray,
    patient_background: np.ndarray,
    overlap_fraction: float,
    status: str,
) -> list[dict]:
    states = [column.removeprefix("weight__") for column in weight_columns]
    weights = data.loc[indices, weight_columns].to_numpy(float) if len(indices) else np.empty((0, len(states)))
    mean_weight = np.nanmean(weights, axis=0) if len(indices) else np.full(len(states), np.nan)
    median_weight = np.nanmedian(weights, axis=0) if len(indices) else np.full(len(states), np.nan)
    enrichment = np.log2((mean_weight + 1e-8) / (patient_background + 1e-8))
    fixed_mean = np.nanmean(fixed_marker[indices], axis=0) if len(indices) else np.full(len(states), np.nan)
    heldout_mean = np.nanmean(heldout_marker[indices], axis=0) if len(indices) else np.full(len(states), np.nan)
    row = {
        "patient_id": region_row.patient_id, "chip": region_row.chip_label,
        "disease": region_row.disease, "roi_id": region_row.roi_id,
        "tissue_component_id": region_row.tissue_component_id,
        "neighborhood_label": region_row.neighborhood_label,
        "local_region_id": region_row.local_region_id,
        "region_bin_n": region_row.bin_n, "region_size_stratum": region_row.size_stratum,
        "zone": zone, "zone_bin_n": len(indices), "zone_status": status,
        "median_umi": float(labels.loc[indices, "nUMI"].median()) if len(indices) else np.nan,
        "median_ngenes": float(labels.loc[indices, "nGene"].median()) if len(indices) else np.nan,
        "mean_graph_degree": float(labels.loc[indices, "delaunay_degree"].mean()) if len(indices) else np.nan,
        "low_confidence_fraction": float(labels.loc[indices, "classification_status"].eq("LOW_CONFIDENCE").mean()) if len(indices) else np.nan,
        "ood_fraction": float(labels.loc[indices, "classification_status"].eq("OOD_UNASSIGNED").mean()) if len(indices) else np.nan,
        "overlap_with_other_region_ring_fraction": overlap_fraction,
        "composition_layout": "wide: one row per region-zone; state names are column suffixes",
    }
    for state_i, state in enumerate(states):
        row[f"mean_weight__{state}"] = mean_weight[state_i]
        row[f"median_weight__{state}"] = median_weight[state_i]
        row[f"log2_enrichment_vs_background__{state}"] = enrichment[state_i]
        row[f"fixed_marker_mean__{state}"] = fixed_mean[state_i]
        row[f"heldout_marker_mean__{state}"] = heldout_mean[state_i]
    return [row]


def compute_zones(
    data, labels, regions, weight_columns, adjacency, fixed_marker, heldout_marker, log
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    eligible_regions = regions[regions.bin_n.ge(SMALL_REGION_MIN_BINS)].copy()
    large_regions = regions[regions.bin_n.ge(FULL_REGION_MIN_BINS)].copy()
    region_members = {
        region_id: group.global_row.to_numpy(int)
        for region_id, group in labels[labels.component_id.notna()].groupby("component_id", sort=False)
    }
    component_members = {
        key: group.global_row.to_numpy(int)
        for key, group in labels.groupby(["patient_id", "roi_id", "tissue_component_id"], dropna=False, sort=False)
    }
    patient_backgrounds = {
        key: data.loc[indices, weight_columns].to_numpy(float).mean(axis=0)
        for key, indices in component_members.items()
    }
    plans: dict[str, dict[str, np.ndarray]] = {}
    ring_count = np.zeros(len(labels), dtype=np.int32)
    for counter, region in enumerate(large_regions.itertuples(index=False), start=1):
        region_indices = region_members[region.local_region_id]
        boundary, immediate = boundary_nodes(region_indices, adjacency)
        interior = np.setdiff1d(region_indices, boundary, assume_unique=False)
        key = (region.patient_id, region.roi_id, region.tissue_component_id)
        component = component_members[key]
        outside_component = np.setdiff1d(component, region_indices, assume_unique=False)
        if len(boundary) and len(outside_component):
            tree = cKDTree(data.loc[boundary, ["x_um", "y_um"]].to_numpy(float))
            points = data.loc[outside_component, ["x_um", "y_um"]].to_numpy(float)
            _, nearest = tree.query(points, k=1)
            closest = data.loc[boundary[nearest], ["x_um", "y_um"]].to_numpy(float)
            delta = np.maximum(np.abs(points - closest) - BIN_HALF_UM, 0.0)
            distance_to_tile_boundary = np.sqrt((delta**2).sum(axis=1))
        else:
            distance_to_tile_boundary = np.full(len(outside_component), np.inf)
        is_immediate = np.isin(outside_component, immediate)
        ring_0_50 = outside_component[(distance_to_tile_boundary <= 50.0) & ~is_immediate]
        ring_50_100 = outside_component[(distance_to_tile_boundary > 50.0) & (distance_to_tile_boundary <= 100.0)]
        excluded = np.unique(np.r_[region_indices, immediate, ring_0_50, ring_50_100])
        background = np.setdiff1d(component, excluded, assume_unique=False)
        plans[region.local_region_id] = {
            "REGION_ALL": region_indices,
            "INTERIOR": interior,
            "INNER_BOUNDARY": boundary,
            "IMMEDIATE_OUTSIDE": immediate,
            "OUTSIDE_RING_0_50": ring_0_50,
            "OUTSIDE_RING_50_100": ring_50_100,
            "PATIENT_ROI_BACKGROUND": background,
        }
        ring_count[np.unique(np.r_[immediate, ring_0_50, ring_50_100])] += 1
        if counter % 250 == 0:
            log.info("Prepared surrounding zones for %d/%d large regions", counter, len(large_regions))

    internal_rows = []
    surrounding_rows = []
    zone_qc = []
    region_lookup = regions.set_index("local_region_id")
    for counter, region in enumerate(eligible_regions.itertuples(index=False), start=1):
        region_indices = region_members[region.local_region_id]
        boundary, _ = boundary_nodes(region_indices, adjacency)
        interior = np.setdiff1d(region_indices, boundary, assume_unique=False)
        key = (region.patient_id, region.roi_id, region.tissue_component_id)
        background_mean = patient_backgrounds[key]
        internal_zones = {"REGION_ALL": region_indices, "INTERIOR": interior, "INNER_BOUNDARY": boundary}
        for zone, indices in internal_zones.items():
            minimum = 1 if zone == "REGION_ALL" else 3
            status = "EVALUABLE" if len(indices) >= minimum else "NOT_EVALUABLE_TOO_FEW_BINS"
            internal_rows.extend(
                zone_summary_rows(region, zone, indices, data, labels, weight_columns, fixed_marker, heldout_marker, background_mean, 0.0, status)
            )
        if region.bin_n >= FULL_REGION_MIN_BINS:
            plan = plans[region.local_region_id]
            for zone in ("IMMEDIATE_OUTSIDE", "OUTSIDE_RING_0_50", "OUTSIDE_RING_50_100", "PATIENT_ROI_BACKGROUND"):
                indices = plan[zone]
                minimum = 3 if zone != "PATIENT_ROI_BACKGROUND" else 10
                status = "EVALUABLE" if len(indices) >= minimum else "NOT_EVALUABLE_TOO_FEW_BINS"
                overlap_fraction = (
                    float((ring_count[indices] > 1).mean())
                    if len(indices) and zone != "PATIENT_ROI_BACKGROUND" else 0.0
                )
                surrounding_rows.extend(
                    zone_summary_rows(region, zone, indices, data, labels, weight_columns, fixed_marker, heldout_marker, background_mean, overlap_fraction, status)
                )
            zone_qc.append(
                {
                    "local_region_id": region.local_region_id,
                    "patient_id": region.patient_id, "neighborhood_label": region.neighborhood_label,
                    "region_bin_n": region.bin_n, "interior_bin_n": len(plan["INTERIOR"]),
                    "inner_boundary_bin_n": len(plan["INNER_BOUNDARY"]),
                    "immediate_outside_bin_n": len(plan["IMMEDIATE_OUTSIDE"]),
                    "outside_ring_0_50_bin_n": len(plan["OUTSIDE_RING_0_50"]),
                    "outside_ring_50_100_bin_n": len(plan["OUTSIDE_RING_50_100"]),
                    "background_bin_n": len(plan["PATIENT_ROI_BACKGROUND"]),
                    "full_zone_analysis_status": "EVALUABLE" if all(
                        len(plan[name]) >= threshold
                        for name, threshold in (("INTERIOR", 3), ("INNER_BOUNDARY", 3), ("IMMEDIATE_OUTSIDE", 3), ("OUTSIDE_RING_0_50", 3), ("OUTSIDE_RING_50_100", 3), ("PATIENT_ROI_BACKGROUND", 10))
                    ) else "PARTIAL_ZONE_COVERAGE",
                    "ring_overlap_present": bool(any((ring_count[plan[name]] > 1).any() for name in ("IMMEDIATE_OUTSIDE", "OUTSIDE_RING_0_50", "OUTSIDE_RING_50_100") if len(plan[name]))),
                    "distance_definition": "minimum approximate Euclidean distance from outside Bin50 center to nearest 25x25um region tile boundary",
                    "ring_partition": "immediate outside excluded from 0-50um ring; zones are mutually exclusive within a focal region",
                }
            )
        if counter % 500 == 0:
            log.info("Summarized %d/%d eligible regions", counter, len(eligible_regions))
    internal = pd.DataFrame(internal_rows)
    surroundings = pd.DataFrame(surrounding_rows)
    zone_qc_frame = pd.DataFrame(zone_qc)
    internal.to_csv(OUT / "tables/06_REGION_INTERNAL_COMPOSITION.tsv", sep="\t", index=False)
    surroundings.to_csv(OUT / "tables/07_REGION_BOUNDARY_AND_SURROUNDINGS.tsv", sep="\t", index=False)
    zone_qc_frame.to_csv(OUT / "tables/07B_REGION_ZONE_EVALUABILITY.tsv", sep="\t", index=False)
    return internal, surroundings, zone_qc_frame


def neighborhood_summary(data, labels, regions, weight_columns) -> pd.DataFrame:
    eligible = labels[labels.predicted_pn.isin(PN_ORDER)].copy()
    rows = []
    weights = data[weight_columns].to_numpy(float)
    for pn in PN_ORDER:
        subset = eligible[eligible.predicted_pn.eq(pn)]
        indices = subset.global_row.to_numpy(int)
        mean = weights[indices].mean(axis=0)
        top = np.argsort(mean)[::-1][:5]
        patient_centroids = []
        for _, block in subset.groupby("patient_id"):
            patient_centroids.append(weights[block.global_row.to_numpy(int)].mean(axis=0))
        similarities = []
        for i in range(len(patient_centroids)):
            for j in range(i + 1, len(patient_centroids)):
                a, b = patient_centroids[i], patient_centroids[j]
                similarities.append(float(a @ b / (np.linalg.norm(a) * np.linalg.norm(b))))
        reg = regions[regions.neighborhood_label.eq(pn)]
        rows.append(
            {
                "neighborhood_label": pn, "bin_n": len(subset),
                "patient_coverage": subset.patient_id.nunique(), "chip_coverage": subset.chip_label.nunique(),
                "reliable_fraction": float(subset.classification_status.eq("RELIABLE").mean()),
                "low_confidence_fraction": float(subset.classification_status.eq("LOW_CONFIDENCE").mean()),
                "ood_unassigned_fraction": float(subset.classification_status.eq("OOD_UNASSIGNED").mean()),
                "top_reference18_states": ";".join(weight_columns[i].removeprefix("weight__") for i in top),
                "top_reference18_mean_weights": ";".join(f"{mean[i]:.6g}" for i in top),
                "patient_centroid_cosine_median": float(np.median(similarities)),
                "patient_centroid_cosine_min": float(np.min(similarities)),
                "all_region_n": len(reg), "region_ge5_n": int(reg.bin_n.ge(5).sum()),
                "region_ge10_n": int(reg.bin_n.ge(10).sum()),
                "patient_coverage_ge10_region": reg.loc[reg.bin_n.ge(10), "patient_id"].nunique(),
                "chip_coverage_ge10_region": reg.loc[reg.bin_n.ge(10), "chip_label"].nunique(),
                "classification_scope": "exploratory spatial organization; not proof of distinct pathological niche",
            }
        )
    result = pd.DataFrame(rows)
    result.to_csv(OUT / "tables/03_NEIGHBORHOOD_CLASSIFICATION_SUMMARY.tsv", sep="\t", index=False)
    return result


def ssdna_review(regions: pd.DataFrame, ssdna_index: pd.DataFrame) -> pd.DataFrame:
    available = set(ssdna_index.patient_id)
    rows = regions[regions.bin_n.ge(10)].copy()
    rows["ssdna_alignment_status"] = np.where(rows.patient_id.isin(available), "VALIDATED_ALIGNMENT_AVAILABLE", "SSDNA_NOT_AVAILABLE")
    rows["eligible_for_ssdna_overlay"] = rows.ssdna_alignment_status.eq("VALIDATED_ALIGNMENT_AVAILABLE")
    rows["review_scope"] = np.where(rows.eligible_for_ssdna_overlay, "region outline and representative zoom", "true-coordinate map only")
    rows["histology_interpretation_status"] = np.where(rows.eligible_for_ssdna_overlay, "VISUAL_REVIEW_REQUIRED", "NOT_EVALUABLE")
    rows[["patient_id", "chip_label", "disease", "roi_id", "tissue_component_id", "neighborhood_label", "local_region_id", "bin_n", "approximate_area_um2", "centroid_x_um", "centroid_y_um", "ssdna_alignment_status", "eligible_for_ssdna_overlay", "review_scope", "histology_interpretation_status"]].to_csv(
        OUT / "tables/08_SSDNA_REGION_REVIEW.tsv", sep="\t", index=False
    )
    return rows


def load_nmf_bin_usage(data: pd.DataFrame):
    mapping = pd.read_csv(NMF_MAP, sep="\t")
    frames = []
    for chip in CHIPS:
        model_path = NMF_DIR / f"neighborhood50_{chip}_K07_MODEL_AND_USAGE.npz"
        model = np.load(model_path, allow_pickle=True)
        entity = model["entity_key"].astype(str)
        usage = model["relative_usage"].astype(float)
        families = [mapping.loc[mapping.chip_short.eq(chip) & mapping.factor.eq(f"F{i+1:02d}"), "factor_family"].iloc[0] for i in range(7)]
        frame = pd.DataFrame(usage, columns=families)
        frame.insert(0, "entity_key", entity)
        frame["nmf_dominant_family"] = frame[families].idxmax(axis=1)
        frame["chip_label"] = chip
        frames.append(frame)
    result = pd.concat(frames, ignore_index=True)
    if not result.entity_key.is_unique:
        raise RuntimeError("NMF entity keys not unique")
    return result


def nmf_comparison(labels: pd.DataFrame, regions: pd.DataFrame) -> pd.DataFrame:
    nmf = load_nmf_bin_usage(labels)
    factor_columns = [column for column in nmf if column not in ("entity_key", "nmf_dominant_family", "chip_label")]
    joined = labels.merge(nmf.drop(columns="chip_label"), on="entity_key", how="left", validate="one_to_one")
    rows = []
    for chip, block in joined[joined.predicted_pn.isin(PN_ORDER)].groupby("chip_label"):
        ari = adjusted_rand_score(block.predicted_pn, block.nmf_dominant_family)
        nmi = normalized_mutual_info_score(block.predicted_pn, block.nmf_dominant_family)
        rows.append(
            {
                "record_type": "BIN_LEVEL_DOMINANT_COMPARISON", "chip": chip,
                "neighborhood_label": "ALL", "nmf_factor_family": "ALL",
                "common_bin_n": len(block), "ari": ari, "nmi": nmi,
                "value": np.nan,
                "interpretation": "low/moderate agreement reflects method and resolution differences; not independent validation",
            }
        )
    grouped = joined[joined.predicted_pn.isin(PN_ORDER)].groupby(["chip_label", "predicted_pn"], sort=True)
    for (chip, pn), block in grouped:
        for family in factor_columns:
            rows.append(
                {
                    "record_type": "PN_MEAN_NMF_USAGE", "chip": chip,
                    "neighborhood_label": pn, "nmf_factor_family": family,
                    "common_bin_n": len(block), "ari": np.nan, "nmi": np.nan,
                    "value": float(block[family].mean()),
                    "interpretation": "continuous frozen NMF K7 usage within fixed PN bins",
                }
            )
    large = regions[regions.bin_n.ge(10)][["local_region_id", "neighborhood_label", "patient_id", "chip_label", "bin_n"]]
    memberships = labels[labels.component_id.isin(large.local_region_id)].merge(
        nmf.drop(columns="chip_label"), on="entity_key", how="left", validate="one_to_one"
    )
    region_usage = memberships.groupby(["component_id", "predicted_pn", "patient_id", "chip_label"], sort=False)[factor_columns].mean().reset_index()
    for _, region in region_usage.iterrows():
        values = region[factor_columns].to_numpy(float)
        top = int(np.nanargmax(values))
        rows.append(
            {
                "record_type": "LARGE_REGION_NMF_PROFILE", "chip": region["chip_label"],
                "patient_id": region["patient_id"], "local_region_id": region["component_id"],
                "neighborhood_label": region["predicted_pn"], "nmf_factor_family": factor_columns[top],
                "common_bin_n": int(large.loc[large.local_region_id.eq(region["component_id"]), "bin_n"].iloc[0]),
                "ari": np.nan, "nmi": np.nan, "value": values[top],
                "interpretation": "dominant mean frozen NMF K7 usage in the Delaunay region; descriptive only",
            }
        )
    result = pd.DataFrame(rows)
    result.to_csv(OUT / "tables/09_DELAUNAY_VS_NMF_SUMMARY.tsv", sep="\t", index=False)
    return result


def ols_patient_contrast(frame: pd.DataFrame, outcome: str, contrast: str) -> dict:
    disease = pd.Categorical(frame.disease, categories=["HC", "IPF", "SSc-ILD"])
    chip = pd.Categorical(frame.chip_label, categories=["K8", "J2", "L3"])
    x = np.column_stack([
        np.ones(len(frame)), (disease == "IPF").astype(float),
        (disease == "SSc-ILD").astype(float), (chip == "J2").astype(float),
        (chip == "L3").astype(float),
    ])
    y = frame[outcome].to_numpy(float)
    if len(frame) <= x.shape[1] or np.linalg.matrix_rank(x) < x.shape[1]:
        return {"estimate": np.nan, "se": np.nan, "p": np.nan, "ci_low": np.nan, "ci_high": np.nan, "status": "NOT_EVALUABLE"}
    beta = np.linalg.lstsq(x, y, rcond=None)[0]
    residual = y - x @ beta
    df = len(y) - x.shape[1]
    covariance = (residual @ residual / df) * np.linalg.inv(x.T @ x)
    vectors = {
        "SSc-ILD_vs_IPF": np.array([0, -1, 1, 0, 0]),
        "SSc-ILD_vs_HC": np.array([0, 0, 1, 0, 0]),
        "IPF_vs_HC": np.array([0, 1, 0, 0, 0]),
    }
    c = vectors[contrast]
    estimate = float(c @ beta)
    se = float(np.sqrt(c @ covariance @ c))
    statistic = estimate / se if se > 0 else np.nan
    p = float(2 * stats.t.sf(abs(statistic), df)) if np.isfinite(statistic) else np.nan
    critical = float(stats.t.ppf(0.975, df))
    return {"estimate": estimate, "se": se, "p": p, "ci_low": estimate - critical * se, "ci_high": estimate + critical * se, "status": "EXPLORATORY"}


def patient_summary(labels: pd.DataFrame, regions: pd.DataFrame) -> pd.DataFrame:
    patient_meta = labels[["patient_id", "chip_label", "disease"]].drop_duplicates()
    total_bins = labels.groupby("patient_id").size()
    rows = []
    for patient_row in patient_meta.itertuples(index=False):
        patient = patient_row.patient_id
        patient_labels = labels[labels.patient_id.eq(patient)]
        patient_regions = regions[regions.patient_id.eq(patient)]
        for pn in PN_ORDER:
            block = patient_labels[patient_labels.predicted_pn.eq(pn)]
            reg = patient_regions[patient_regions.neighborhood_label.eq(pn)]
            rows.append(
                {
                    "record_type": "PATIENT_PN_SUMMARY", "patient_id": patient,
                    "chip_label": patient_row.chip_label, "disease": patient_row.disease,
                    "neighborhood_label": pn, "bin_n": len(block),
                    "bin_fraction": len(block) / total_bins[patient],
                    "reliable_fraction": float(block.classification_status.eq("RELIABLE").mean()) if len(block) else np.nan,
                    "low_confidence_fraction": float(block.classification_status.eq("LOW_CONFIDENCE").mean()) if len(block) else np.nan,
                    "region_n": len(reg), "region_ge5_n": int(reg.bin_n.ge(5).sum()),
                    "region_ge10_n": int(reg.bin_n.ge(10).sum()),
                    "region_ge10_area_fraction": reg.loc[reg.bin_n.ge(10), "bin_n"].sum() / total_bins[patient],
                    "largest_region_bin_fraction": (reg.bin_n.max() / total_bins[patient]) if len(reg) else 0.0,
                    "median_region_bin_n": float(reg.bin_n.median()) if len(reg) else np.nan,
                    "analysis_scope": "patient-level first; chip-aware exploratory description",
                }
            )
    patient_table = pd.DataFrame(rows)
    model_rows = []
    contrasts = ["SSc-ILD_vs_IPF", "SSc-ILD_vs_HC", "IPF_vs_HC"]
    metrics = ["bin_fraction", "region_ge10_n", "region_ge10_area_fraction", "largest_region_bin_fraction"]
    for pn in PN_ORDER:
        block = patient_table[patient_table.neighborhood_label.eq(pn)]
        for metric in metrics:
            for contrast in contrasts:
                fit = ols_patient_contrast(block, metric, contrast)
                model_rows.append(
                    {
                        "record_type": "CHIP_ADJUSTED_DISEASE_DESCRIPTION",
                        "neighborhood_label": pn, "metric": metric, "contrast": contrast,
                        "estimate": fit["estimate"], "standard_error": fit["se"],
                        "ci_low": fit["ci_low"], "ci_high": fit["ci_high"], "raw_p": fit["p"],
                        "model_status": fit["status"], "analysis_unit": "patient",
                        "interpretation": "exploratory; not used to select or rename neighborhoods",
                    }
                )
    models = pd.DataFrame(model_rows)
    for (metric, contrast), indices in models.groupby(["metric", "contrast"]).groups.items():
        models.loc[indices, "bh_fdr"] = bh(models.loc[indices, "raw_p"].to_numpy(float))
    result = pd.concat([patient_table, models], ignore_index=True, sort=False)
    result.to_csv(OUT / "tables/10_CHIP_AWARE_PATIENT_SUMMARY.tsv", sep="\t", index=False)
    return result


def evidence_classification(summary, regions, zone_qc, ssdna_review_table, patient_table) -> pd.DataFrame:
    model = patient_table[patient_table.record_type.eq("CHIP_ADJUSTED_DISEASE_DESCRIPTION")]
    records = [
        ("Delaunay graph and QC", "Computationally supported", "126,449 bins; 357,959 edges; 21 patients; boundary-constrained graph reproduced exactly"),
        ("Unified PN1-PN8 labels", "Spatially descriptive", "Frozen labels and confidence/OOD flags retained; not all PN are independent pathological niches"),
        ("Reference18 PN composition", "Computationally supported", f"8 PN; all have {summary.patient_coverage.min()}-{summary.patient_coverage.max()} patient coverage; composition is internal to the Reference18 feature source"),
        ("Marker concordance", "Exploratory / low confidence", "Depth-adjusted internal concordance exists for selected PN; no held-out marker module passed the prespecified FDR naming gate"),
        ("All connected regions", "Computationally supported", f"{len(regions):,} graph-connected regions retained; {int(regions.bin_n.ge(10).sum()):,} have >=10 Bin50"),
        ("Interior/boundary/surrounding composition", "Spatially descriptive", f"{int(zone_qc.full_zone_analysis_status.eq('EVALUABLE').sum()):,}/{len(zone_qc):,} large regions have all zones evaluable; ring overlap explicitly flagged"),
        ("ssDNA correspondence", "Biologically suggestive", f"Validated coordinate overlays available for {ssdna_review_table.loc[ssdna_review_table.eligible_for_ssdna_overlay, 'patient_id'].nunique()}/21 patients; morphology requires visual review"),
        ("Chip-adjusted disease evidence", "Exploratory / low confidence", f"{int(model.bh_fdr.lt(0.05).sum())} exploratory patient-level PN-metric contrasts at FDR<0.05; chip partial R2 remains a major limitation"),
        ("Frozen NMF K7 agreement", "Spatially descriptive", "Descriptive composition/usage agreement only; both methods derive from the same Reference18 source"),
        ("Independent pathological niche validation", "Not supported", "No independent pathology or held-out molecular validation establishes all PN as distinct niches"),
        ("BANKSY readiness", "Computationally supported", "Frozen Delaunay outputs provide a complete benchmark and handoff; BANKSY must remain a separate method comparison"),
    ]
    result = pd.DataFrame(records, columns=["evidence_item", "evidence_class", "evidence_basis"])
    result["closure_decision"] = "DELAUNAY_ANALYSIS_CLOSURE_READY"
    result.to_csv(OUT / "tables/11_FINAL_EVIDENCE_CLASSIFICATION.tsv", sep="\t", index=False)
    return result


def main() -> None:
    ensure_dirs()
    log = logger()
    log.info("START fixed-label Delaunay closure")
    data, labels, weight_columns, pairs, graph_qc, adjacency = load_core(log)
    audit = input_audit(data, labels, weight_columns, pairs, graph_qc)
    if not audit.status.eq("PASS").all():
        raise RuntimeError("Input audit failed")
    ssdna_index = pd.read_csv(SSDNA_INDEX, sep="\t")
    ssdna_patients = set(ssdna_index.patient_id)
    regions = region_table(data, labels, weight_columns, ssdna_patients)
    summary = neighborhood_summary(data, labels, regions, weight_columns)
    fixed_marker, heldout_marker = load_marker_arrays(data.entity_key, [column.removeprefix("weight__") for column in weight_columns], log)
    internal, surroundings, zone_qc = compute_zones(
        data, labels, regions, weight_columns, adjacency, fixed_marker, heldout_marker, log
    )
    ssdna_review_table = ssdna_review(regions, ssdna_index)
    nmf = nmf_comparison(labels, regions)
    patient = patient_summary(labels, regions)
    evidence_classification(summary, regions, zone_qc, ssdna_review_table, patient)
    completeness()
    run = {
        "status": "COMPUTE_COMPLETE",
        "seed": SEED,
        "patient_n": 21,
        "region_n": len(regions),
        "region_ge5_n": int(regions.bin_n.ge(5).sum()),
        "region_ge10_n": int(regions.bin_n.ge(10).sum()),
        "internal_rows": len(internal),
        "surrounding_rows": len(surroundings),
        "spatial_cooccurrence_run": False,
        "clustering_or_k_selection_run": False,
        "upstream_model_run": False,
    }
    (OUT / "logs/COMPUTE_VALIDATION.json").write_text(json.dumps(run, indent=2), encoding="utf-8")
    log.info("COMPLETE %s", run)


if __name__ == "__main__":
    main()
