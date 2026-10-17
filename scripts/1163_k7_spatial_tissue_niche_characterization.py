#!/usr/bin/env python3
"""Characterize frozen K=7 Reference18 composition factors in physical space.

This workflow is descriptive and disease-blind. It never refits RCTD, NMF,
edgeR, or GSEA. All spatial graphs are constructed within one patient only.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
from pathlib import Path

import h5py
import numpy as np
import pandas as pd
from scipy.sparse import csr_matrix
from scipy.sparse.csgraph import connected_components
from scipy.spatial import cKDTree


ROOT = Path(__file__).resolve().parents[1]
REANALYSIS = ROOT / "results" / "reanalysis"
ATLAS = REANALYSIS / "ALL21_SHARED_UNBIASED_DISCRETE_NICHE_ATLAS_20261011_191652_attempt4"
ASSIGNMENTS = ATLAS / "02_ALL21_SHARED_NICHE_ASSIGNMENTS.tsv.gz"
ASSAY = ATLAS / "12_SHARED_NICHE_ASSAY.h5"
MODEL_RUN = REANALYSIS / "REFERENCE18_NEIGHBORHOOD_NMF_K7_K12_EXTENSION_20261015_attempt5"
MAP_RUN = REANALYSIS / "REFERENCE18_NEIGHBORHOOD_NMF_K7_STATE_SENSITIVITY_20261016_210302_attempt5"
COUNT_H5 = REANALYSIS / "RCTD_SIX_DATASET_LEVEL1_ALL21_20260916_223802/prepared_inputs/six_dataset_level1_rctd_inputs.h5"
JOINT = REANALYSIS / "REFERENCE18_K7_FACTOR_DISEASE_EDGER_GSEA_261017223852"
GENCODE = ROOT / "refs/gencode_v50_GRCh38p14/gencode_v50_gene_annotation.tsv.gz"

CHIPS = ("K8", "J2", "L3")
QUANTILES = (0.90, 0.95, 0.99)
GRID_UM = 25.0
ADJACENCY_UM = math.sqrt(2) * GRID_UM + 1e-6
FAMILIES = (
    "LEC-associated", "Blood-endothelial", "AT1-associated",
    "Perivascular-stromal", "B-enriched/mixed-immune",
    "Ciliated-associated", "Fibroblast-associated",
)
MARKERS = {
    "B-enriched/mixed-immune": ["MS4A1", "CD79A", "CD79B", "BANK1", "MZB1", "XBP1", "JCHAIN"],
    "Ciliated-associated": ["FOXJ1", "RSPH1", "TPPP3", "SCGB1A1", "DNAH5", "DNAH9", "SPAG16"],
    "AT1-associated": ["AGER", "PDPN", "CAV1", "SFTPC", "SFTPA1"],
    "Fibroblast-associated": ["DCN", "LUM", "COL1A1", "COL1A2", "POSTN"],
    "Blood-endothelial": ["PECAM1", "VWF", "KDR", "PLVAP"],
    "LEC-associated": ["PROX1", "FLT4", "CCL21", "LYVE1"],
    "Perivascular-stromal": ["RGS5", "PDGFRB", "CSPG4", "ACTA2"],
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output-dir", required=True, type=Path)
    return parser.parse_args()


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def decode(values: np.ndarray) -> list[str]:
    return [value.decode() if isinstance(value, bytes) else str(value) for value in values]


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def exact_top_mask(values: np.ndarray, quantile: float, keys: np.ndarray) -> tuple[np.ndarray, float]:
    count = max(1, int(math.ceil((1.0 - quantile) * len(values))))
    order = np.lexsort((keys.astype(str), -values))
    mask = np.zeros(len(values), dtype=bool)
    mask[order[:count]] = True
    return mask, float(values[order[count - 1]])


def load_spatial_frame() -> tuple[pd.DataFrame, pd.DataFrame, list[str], np.ndarray]:
    metadata = pd.read_csv(ASSIGNMENTS, sep="\t").set_index("entity_key")
    family_map = pd.read_csv(MAP_RUN / "07_CROSS_CHIP_FAMILY_MAP.tsv", sep="\t")
    family_map = family_map.loc[
        family_map["scenario"].eq("main18") & family_map["family_label"].isin(FAMILIES)
    ].copy()
    require(set(family_map["family_label"]) == set(FAMILIES), "Incomplete family map")
    family_map["family_order"] = family_map["family_label"].map({x: i for i, x in enumerate(FAMILIES)})
    family_map = family_map.sort_values("family_order")

    parts = []
    for chip in CHIPS:
        path = MODEL_RUN / "models" / f"neighborhood50_{chip}_K07_MODEL_AND_USAGE.npz"
        with np.load(path, allow_pickle=True) as archive:
            keys = archive["entity_key"].astype(str)
            usage = archive["relative_usage"].astype(float)
        current = metadata.loc[keys].reset_index()
        for row in family_map.itertuples(index=False):
            factor = str(getattr(row, f"{chip}_factor"))
            current[row.family_label] = usage[:, int(factor[1:]) - 1]
        parts.append(current)
    frame = pd.concat(parts, ignore_index=True)
    require(frame["patient_id"].nunique() == 21, "Expected 21 patients")
    require(not frame["entity_key"].duplicated().any(), "Duplicated entity key")

    with h5py.File(ASSAY, "r") as handle:
        assay_keys = decode(handle["entity_key"][:])
        state_names = decode(handle["state_names"][:])
        composition = handle["composition"][:].astype(float)
    assay_index = pd.Series(np.arange(len(assay_keys)), index=assay_keys)
    take = assay_index.loc[frame["entity_key"]].to_numpy(int)
    return frame, family_map, state_names, composition[take]


def input_audit(frame: pd.DataFrame, family_map: pd.DataFrame, output: Path) -> None:
    rows = []
    for row in family_map.itertuples(index=False):
        for chip in CHIPS:
            factor = str(getattr(row, f"{chip}_factor"))
            load_path = MODEL_RUN / "models" / f"neighborhood50_{chip}_K07_LOADINGS.tsv"
            loadings = pd.read_csv(load_path, sep="\t", index_col=0)
            ordered = loadings.loc[factor].sort_values(ascending=False)
            rows.append({
                "record_type": "factor_map", "chip": chip,
                "family_label": row.family_label, "chip_specific_factor": factor,
                "patient_n": frame.loc[frame["chip_label"].eq(chip), "patient_id"].nunique(),
                "bin_n": int(frame["chip_label"].eq(chip).sum()),
                "top_state_1": ordered.index[0], "top_loading_1": ordered.iloc[0],
                "top_state_2": ordered.index[1], "top_loading_2": ordered.iloc[1],
                "top_state_3": ordered.index[2], "top_loading_3": ordered.iloc[2],
                "loading_sha256": sha256(load_path),
            })
    for path, label in [
        (ASSIGNMENTS, "spatial_metadata"), (ASSAY, "reference18_neighborhood_assay"),
        (COUNT_H5, "raw_counts"), (MAP_RUN / "07_CROSS_CHIP_FAMILY_MAP.tsv", "family_map"),
    ]:
        rows.append({"record_type": "input", "input_name": label,
                     "path": str(path.resolve()), "sha256": sha256(path)})
    pd.DataFrame(rows).to_csv(output / "01_INPUT_AND_FACTOR_MAP_AUDIT.tsv", sep="\t", index=False)


def graph_for_coordinates(coordinates: np.ndarray) -> tuple[csr_matrix, np.ndarray, np.ndarray]:
    pairs = cKDTree(coordinates).query_pairs(ADJACENCY_UM, output_type="ndarray")
    if len(pairs) == 0:
        adjacency = csr_matrix((len(coordinates), len(coordinates)), dtype=np.int8)
    else:
        rows = np.r_[pairs[:, 0], pairs[:, 1]]
        cols = np.r_[pairs[:, 1], pairs[:, 0]]
        adjacency = csr_matrix((np.ones(len(rows), dtype=np.int8), (rows, cols)),
                               shape=(len(coordinates), len(coordinates)))
    _, roi = connected_components(adjacency, directed=False)
    degree = np.asarray(adjacency.sum(axis=1)).ravel()
    return adjacency, roi, degree


def shape_metrics(coords: np.ndarray) -> tuple[float, float, float]:
    cells = {(int(round(x / GRID_UM)), int(round(y / GRID_UM))) for x, y in coords}
    perimeter_edges = sum(
        sum((gx + dx, gy + dy) not in cells for dx, dy in ((1, 0), (-1, 0), (0, 1), (0, -1)))
        for gx, gy in cells
    )
    area = len(coords) * GRID_UM * GRID_UM
    perimeter = perimeter_edges * GRID_UM
    compactness = float(4 * math.pi * area / (perimeter * perimeter)) if perimeter > 0 else np.nan
    if len(coords) < 3:
        elongation = np.nan
    else:
        eig = np.linalg.eigvalsh(np.cov(coords.T))
        elongation = float(math.sqrt(max(eig[-1], 1e-12) / max(eig[0], 1e-12)))
    return area, perimeter, compactness, elongation


def hotspot_components(frame: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame, dict]:
    component_rows = []
    summary_rows = []
    masks: dict[tuple[str, str, str], np.ndarray] = {}
    for patient, group in frame.groupby("patient_id", sort=True):
        indices = group.index.to_numpy(int)
        coords = group[["x_um", "y_um"]].to_numpy(float)
        adjacency, roi, tissue_degree = graph_for_coordinates(coords)
        chip = str(group["chip_label"].iloc[0])
        disease = str(group["disease"].iloc[0])
        for family in FAMILIES:
            values = group[family].to_numpy(float)
            for quantile in QUANTILES:
                label = f"q{int(quantile * 100)}"
                high, cutoff = exact_top_mask(values, quantile, group["entity_key"].to_numpy())
                global_mask = np.zeros(len(frame), dtype=bool)
                global_mask[indices] = high
                masks[(patient, family, label)] = global_mask
                sub = adjacency[high][:, high]
                component_n, labels = connected_components(sub, directed=False)
                high_local = np.flatnonzero(high)
                sizes = np.bincount(labels, minlength=component_n)
                for component in range(component_n):
                    local_members = high_local[labels == component]
                    component_coords = coords[local_members]
                    area, perimeter, compactness, elongation = shape_metrics(component_coords)
                    member_frame = group.iloc[local_members]
                    component_rows.append({
                        "record_type": "component", "patient_id": patient,
                        "chip": chip, "disease": disease, "family_label": family,
                        "threshold": label, "usage_cutoff": cutoff,
                        "component_id": component + 1, "roi_id": int(pd.Series(roi[local_members]).mode()[0]) + 1,
                        "bin_n": len(local_members), "area_um2": area,
                        "perimeter_um": perimeter, "compactness": compactness,
                        "elongation": elongation,
                        "edge_bin_fraction": float(np.mean(tissue_degree[local_members] < 8)),
                        "mean_nUMI": member_frame["nUMI"].mean(),
                        "mean_nGene": member_frame["nGene"].mean(),
                        "mean_usage": member_frame[family].mean(),
                        "centroid_x_um": component_coords[:, 0].mean(),
                        "centroid_y_um": component_coords[:, 1].mean(),
                        "isolated_bin": len(local_members) == 1,
                    })
                summary_rows.append({
                    "record_type": "threshold_summary", "patient_id": patient,
                    "chip": chip, "disease": disease, "family_label": family,
                    "threshold": label, "usage_cutoff": cutoff,
                    "hotspot_bin_n": int(high.sum()), "component_n": int(component_n),
                    "largest_component_bin_n": int(sizes.max()),
                    "largest_component_fraction": float(sizes.max() / high.sum()),
                    "isolated_component_n": int(np.sum(sizes == 1)),
                    "q99_low_coverage_flag": bool(label == "q99" and high.sum() < 20),
                    "usage_nUMI_spearman": group[[family, "nUMI"]].corr(method="spearman").iloc[0, 1],
                    "usage_nGene_spearman": group[[family, "nGene"]].corr(method="spearman").iloc[0, 1],
                })
    metrics = pd.concat([pd.DataFrame(summary_rows), pd.DataFrame(component_rows)], ignore_index=True, sort=False)
    return metrics, pd.DataFrame(summary_rows), masks


def composition_support(
    frame: pd.DataFrame, family_map: pd.DataFrame, state_names: list[str], composition: np.ndarray,
    masks: dict, output: Path,
) -> pd.DataFrame:
    rows = []
    for map_row in family_map.itertuples(index=False):
        family = map_row.family_label
        for chip in CHIPS:
            factor = str(getattr(map_row, f"{chip}_factor"))
            loading = pd.read_csv(
                MODEL_RUN / "models" / f"neighborhood50_{chip}_K07_LOADINGS.tsv",
                sep="\t", index_col=0,
            ).loc[factor]
            for state in state_names:
                rows.append({"record_type": "factor_loading", "family_label": family,
                             "chip": chip, "chip_specific_factor": factor,
                             "state": state, "value": float(loading[state])})
    for patient, group in frame.groupby("patient_id", sort=True):
        indices = group.index.to_numpy(int)
        chip = str(group["chip_label"].iloc[0])
        disease = str(group["disease"].iloc[0])
        for family in FAMILIES:
            high_global = masks[(patient, family, "q95")]
            high = high_global[indices]
            for region, take in (("q95_hotspot", high), ("non_hotspot", ~high)):
                means = composition[indices][take].mean(axis=0)
                for state, value in zip(state_names, means):
                    rows.append({"record_type": "observed_composition", "family_label": family,
                                 "chip": chip, "patient_id": patient, "disease": disease,
                                 "region": region, "state": state, "value": float(value)})
    table = pd.DataFrame(rows)
    table.to_csv(output / "04_FACTOR_REFERENCE18_COMPOSITION.tsv", sep="\t", index=False)
    return table


def marker_support(frame: pd.DataFrame, masks: dict, output: Path) -> pd.DataFrame:
    wanted = sorted({gene for genes in MARKERS.values() for gene in genes})
    rows = []
    with h5py.File(COUNT_H5, "r") as handle:
        genes = decode(handle["genes"][:])
        gene_index = {gene: index for index, gene in enumerate(genes)}
        available = [gene for gene in wanted if gene in gene_index]
        for patient, group in frame.groupby("patient_id", sort=True):
            sample_safe = patient.replace("/", "_")
            h5_group = handle[f"spatial/{sample_safe}"]
            bins = decode(h5_group["bins"][:])
            bin_index = {value: index for index, value in enumerate(bins)}
            columns = np.asarray([bin_index[value] for value in group["bin_id"].astype(str)], dtype=int)
            sparse = h5_group["counts"]
            from scipy.sparse import csc_matrix
            matrix = csc_matrix((sparse["data"][:], sparse["indices"][:], sparse["indptr"][:]),
                                shape=tuple(sparse["shape"][:]))
            marker_counts = matrix[[gene_index[g] for g in available]][:, columns].toarray().astype(float)
            library = np.asarray(matrix[:, columns].sum(axis=0)).ravel().astype(float)
            log_cpm = np.log2(marker_counts / np.maximum(library[None, :], 1) * 1e6 + 1)
            indices = group.index.to_numpy(int)
            for family in FAMILIES:
                high = masks[(patient, family, "q95")][indices]
                for gene, raw, expression in zip(available, marker_counts, log_cpm):
                    rows.append({
                        "record_type": "marker_gene", "patient_id": patient,
                        "chip": group["chip_label"].iloc[0], "disease": group["disease"].iloc[0],
                        "family_label": family, "gene_symbol": gene,
                        "expected_for_family": gene in MARKERS[family],
                        "hotspot_detection_fraction": float(np.mean(raw[high] > 0)),
                        "background_detection_fraction": float(np.mean(raw[~high] > 0)),
                        "hotspot_mean_log2cpm": float(expression[high].mean()),
                        "background_mean_log2cpm": float(expression[~high].mean()),
                        "hotspot_minus_background_log2cpm": float(expression[high].mean() - expression[~high].mean()),
                        "hotspot_mean_library_size": float(library[high].mean()),
                        "background_mean_library_size": float(library[~high].mean()),
                    })
    marker = pd.DataFrame(rows)

    edge = pd.read_csv(JOINT / "09_K7_FAMILY_DISEASE_EDGER_FULL.tsv.gz", sep="\t")
    edge = edge.loc[edge["gene_symbol"].isin(wanted)].copy()
    edge_rows = pd.DataFrame({
        "record_type": "existing_joint_edger", "family_label": edge["family_label"],
        "gene_symbol": edge["gene_symbol"], "endpoint": edge["endpoint"],
        "comparison": edge["comparison"], "effect": edge["log2_fc"],
        "p_value": edge["p_value"], "fdr": edge["global_fdr_within_endpoint"],
    })
    gsea = pd.read_csv(JOINT / "14_K7_FAMILY_DISEASE_PRERANKED_GSEA.tsv.gz", sep="\t")
    pattern = "B Cell|Plasma|Cili|Airway|Alveol|Extracellular Matrix|Collagen|Angiogen|Endothel|Vascular|Translation|ERAD|Interferon"
    gsea = gsea.loc[gsea["pathway"].str.contains(pattern, case=False, regex=True)].copy()
    gsea_rows = pd.DataFrame({
        "record_type": "existing_joint_gsea", "family_label": gsea["family_label"],
        "endpoint": gsea["endpoint"], "comparison": gsea["comparison"],
        "program": gsea["pathway"], "collection": gsea["collection"],
        "effect": gsea["NES"], "p_value": gsea["pval"],
        "fdr": gsea["global_fdr_within_endpoint_collection"],
    })
    result = pd.concat([marker, edge_rows, gsea_rows], ignore_index=True, sort=False)
    result.to_csv(output / "05_FACTOR_MARKER_AND_PROGRAM_SUPPORT.tsv", sep="\t", index=False)
    return result


def spatial_relationships(frame: pd.DataFrame, masks: dict, output: Path) -> pd.DataFrame:
    rows = []
    for patient, group in frame.groupby("patient_id", sort=True):
        indices = group.index.to_numpy(int)
        coords = group[["x_um", "y_um"]].to_numpy(float)
        tree_all = cKDTree(coords)
        total_cumulative = {r: tree_all.count_neighbors(tree_all, r) for r in (50.0, 100.0, 200.0)}
        total_bands = {
            "0-50": total_cumulative[50.0],
            "50-100": total_cumulative[100.0] - total_cumulative[50.0],
            "100-200": total_cumulative[200.0] - total_cumulative[100.0],
        }
        grid = {(int(round(x / GRID_UM)), int(round(y / GRID_UM))): i for i, (x, y) in enumerate(coords)}
        for left_index, left in enumerate(FAMILIES):
            left_mask = masks[(patient, left, "q95")][indices]
            left_coords = coords[left_mask]
            for right in FAMILIES[left_index + 1:]:
                right_mask = masks[(patient, right, "q95")][indices]
                right_coords = coords[right_mask]
                overlap = int(np.sum(left_mask & right_mask))
                union = int(np.sum(left_mask | right_mask))
                left_tree, right_tree = cKDTree(left_coords), cKDTree(right_coords)
                cumulative = {r: left_tree.count_neighbors(right_tree, r) for r in (50.0, 100.0, 200.0)}
                observed = {
                    "0-50": cumulative[50.0],
                    "50-100": cumulative[100.0] - cumulative[50.0],
                    "100-200": cumulative[200.0] - cumulative[100.0],
                }
                shared = 0
                left_only = np.flatnonzero(left_mask & ~right_mask)
                right_only_set = set(np.flatnonzero(right_mask & ~left_mask))
                for idx in left_only:
                    gx, gy = int(round(coords[idx, 0] / GRID_UM)), int(round(coords[idx, 1] / GRID_UM))
                    for dx, dy in ((1, 0), (-1, 0), (0, 1), (0, -1)):
                        neighbor = grid.get((gx + dx, gy + dy))
                        shared += int(neighbor in right_only_set)
                base = {
                    "patient_id": patient, "chip": group["chip_label"].iloc[0],
                    "disease": group["disease"].iloc[0], "family_a": left, "family_b": right,
                    "hotspot_definition": "patient-specific exact top 5% direct K=7 usage",
                    "overlap_bin_n": overlap, "jaccard_overlap": overlap / max(union, 1),
                    "shared_boundary_edges": shared,
                }
                for band in ("0-50", "50-100", "100-200"):
                    expected = total_bands[band] * left_mask.mean() * right_mask.mean()
                    rows.append({**base, "distance_band_um": band,
                                 "observed_pair_n": int(observed[band]),
                                 "random_mixing_expected_pair_n": float(expected),
                                 "enrichment_ratio": float(observed[band] / expected) if expected > 0 else np.nan,
                                 "log2_enrichment": float(np.log2((observed[band] + 0.5) / (expected + 0.5)))})
    result = pd.DataFrame(rows)
    result.to_csv(output / "06_FACTOR_SPATIAL_RELATIONSHIPS.tsv", sep="\t", index=False)
    return result


def reproducibility_and_classification(
    summaries: pd.DataFrame, composition_table: pd.DataFrame,
    marker_table: pd.DataFrame, relationships: pd.DataFrame, output: Path,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    observed = composition_table.loc[composition_table["record_type"].eq("observed_composition")]
    pivot = observed.pivot_table(index=["family_label", "chip", "patient_id", "region"],
                                 columns="state", values="value").reset_index()
    comp_rows = []
    for family, family_frame in pivot.loc[pivot["region"].eq("q95_hotspot")].groupby("family_label"):
        matrix = family_frame.drop(columns=["family_label", "chip", "patient_id", "region"]).to_numpy(float)
        corr = np.corrcoef(matrix)
        comp_rows.append((family, float(np.nanmedian(corr[np.triu_indices_from(corr, 1)]))))
    comp_consistency = dict(comp_rows)
    expected_markers = marker_table.loc[
        marker_table["record_type"].eq("marker_gene") & marker_table["expected_for_family"].eq(True)
    ]
    patient_marker = expected_markers.groupby(["family_label", "patient_id"])["hotspot_minus_background_log2cpm"].mean()
    marker_by_family = patient_marker.groupby("family_label").agg(["median", lambda x: np.mean(x > 0)])
    marker_by_family.columns = ["median_marker_delta", "positive_marker_patient_fraction"]
    q95 = summaries.loc[summaries["threshold"].eq("q95")]
    rows = []
    for family in FAMILIES:
        current = q95.loc[q95["family_label"].eq(family)]
        structured = (current["largest_component_bin_n"].ge(20) &
                      current["largest_component_fraction"].ge(0.20) &
                      current["isolated_component_n"].lt(current["component_n"] * 0.8))
        rows.append({
            "family_label": family, "patient_n": current["patient_id"].nunique(),
            "structured_patient_n": int(structured.sum()),
            "structured_patient_fraction": float(structured.mean()),
            "median_largest_component_fraction": current["largest_component_fraction"].median(),
            "median_component_n": current["component_n"].median(),
            "median_abs_usage_nUMI_spearman": current["usage_nUMI_spearman"].abs().median(),
            "median_abs_usage_nGene_spearman": current["usage_nGene_spearman"].abs().median(),
            "median_patient_composition_correlation": comp_consistency.get(family, np.nan),
            "median_expected_marker_delta_log2cpm": marker_by_family.loc[family, "median_marker_delta"],
            "expected_marker_positive_patient_fraction": marker_by_family.loc[family, "positive_marker_patient_fraction"],
        })
    reproducibility = pd.DataFrame(rows)
    reproducibility.to_csv(output / "07_PATIENT_REPRODUCIBILITY.tsv", sep="\t", index=False)

    classifications = []
    for row in reproducibility.itertuples(index=False):
        if (row.structured_patient_fraction >= 2 / 3 and
                row.median_patient_composition_correlation >= 0.75 and
                row.expected_marker_positive_patient_fraction >= 0.60 and
                row.median_abs_usage_nUMI_spearman < 0.50):
            category = "A. Biologically interpretable candidate tissue niche"
            rationale = "recurrent connected hotspots with composition and expected-marker support"
        elif (row.structured_patient_fraction >= 0.50 and
              row.median_patient_composition_correlation >= 0.70):
            category = "B. Composition-associated spatial domain"
            rationale = "recurrent composition-linked spatial domain; independent marker/shape support incomplete"
        elif row.structured_patient_n >= 3:
            category = "C. Patient-specific candidate"
            rationale = "coherent connected structure in a subset of patients"
        else:
            category = "D. Unresolved / technically confounded"
            rationale = "limited recurrence, fragmentation, or technical association"
        classifications.append({"family_label": row.family_label, "classification": category,
                                "rationale": rationale,
                                "disease_specificity_claim": "not tested in this disease-blind characterization",
                                "independent_validation_status": "not independently validated; composition is model-derived"})
    classification = pd.DataFrame(classifications)
    classification.to_csv(output / "08_CANDIDATE_NICHE_CLASSIFICATION.tsv", sep="\t", index=False)
    return reproducibility, classification


def main() -> None:
    args = parse_args()
    output = args.output_dir.resolve()
    if output.exists():
        raise FileExistsError(f"Refusing to overwrite {output}")
    output.mkdir(parents=True)
    (output / "figures").mkdir()
    (output / "work").mkdir()

    frame, family_map, states, composition = load_spatial_frame()
    input_audit(frame, family_map, output)
    metrics, summaries, masks = hotspot_components(frame)
    metrics.to_csv(output / "03_FACTOR_HOTSPOT_AND_COMPONENT_METRICS.tsv", sep="\t", index=False)
    comp = composition_support(frame, family_map, states, composition, masks, output)
    marker = marker_support(frame, masks, output)
    relationships = spatial_relationships(frame, masks, output)
    reproducibility, classification = reproducibility_and_classification(
        summaries, comp, marker, relationships, output
    )

    np.savez_compressed(
        output / "work/K7_SPATIAL_FRAME.npz",
        entity_key=frame["entity_key"].to_numpy(str),
        patient_id=frame["patient_id"].to_numpy(str),
        chip=frame["chip_label"].to_numpy(str),
        disease=frame["disease"].to_numpy(str),
        x_um=frame["x_um"].to_numpy(float), y_um=frame["y_um"].to_numpy(float),
        nUMI=frame["nUMI"].to_numpy(float), nGene=frame["nGene"].to_numpy(float),
        usage=frame[list(FAMILIES)].to_numpy(np.float32),
        composition=composition.astype(np.float32), state_names=np.asarray(states),
        family_names=np.asarray(FAMILIES),
    )
    validation = pd.DataFrame([
        ("patient_n", 21, frame["patient_id"].nunique()),
        ("family_n", 7, len(FAMILIES)),
        ("chip_n", 3, frame["chip_label"].nunique()),
        ("threshold_n", 3, summaries["threshold"].nunique()),
        ("patient_family_threshold_n", 441, len(summaries)),
        ("relationship_patient_pair_band_n", 1323, len(relationships)),
        ("classification_n", 7, len(classification)),
    ], columns=["check", "expected", "observed"])
    validation["passed"] = validation["expected"].astype(str).eq(validation["observed"].astype(str))
    validation.to_csv(output / "99_VALIDATION.tsv", sep="\t", index=False)
    require(validation["passed"].all(), "Validation failed")
    (output / "COMPUTATION_COMPLETE.ok").write_text("COMPLETE\n", encoding="utf-8")
    print(output)


if __name__ == "__main__":
    main()
