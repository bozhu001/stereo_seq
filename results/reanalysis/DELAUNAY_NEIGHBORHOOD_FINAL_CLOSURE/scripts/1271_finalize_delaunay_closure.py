#!/usr/bin/env python3
"""Create closure contrasts, figures, atlas, and Chinese reports.

This script only reads frozen Reference18, unified PN1-PN8, saved component,
marker, ssDNA-registration, and frozen NMF K7 outputs. It never fits a model.
"""

from __future__ import annotations

import importlib.util
import json
import logging
import math
import re
import sys
from pathlib import Path

import matplotlib
import numpy as np
import pandas as pd
import scipy.sparse as sp
from PIL import Image, ImageOps
from matplotlib.collections import PolyCollection
from matplotlib.lines import Line2D
from pypdf import PdfReader, PdfWriter
from scipy.spatial import cKDTree

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import seaborn as sns


ROOT = Path(__file__).resolve().parents[4]
REAN = ROOT / "results/reanalysis"
OUT = REAN / "DELAUNAY_NEIGHBORHOOD_FINAL_CLOSURE"
TABLES = OUT / "tables"
FIGURES = OUT / "figures"
GRAPH_ROOT = REAN / "REFERENCE18_SPATIAL_GRAPH_NEIGHBORHOODS_V1"
FORMAL = GRAPH_ROOT / "ALL21_EXTERNAL_PROJECTION/FORMAL_RUN_MODEL_RECOVERED_EXACT"
COMPUTE_SCRIPT = OUT / "scripts/1270_compute_local_region_closure.py"
COMPOSITION_PATH = FORMAL / "04_ALL21_PN_COMPOSITION.tsv"
MARKER_PATH = GRAPH_ROOT / "ALL21_PN_MARKER_VALIDATION/05_PN_MARKER_CONCORDANCE.tsv"
MARKER_INTERPRETATION = GRAPH_ROOT / "ALL21_PN_MARKER_VALIDATION/08_PN_BIOLOGICAL_INTERPRETATION.tsv"
SSDNA_INDEX = GRAPH_ROOT / "SSDNA_OVERLAY_REVIEW/00_INPUT_FILE_INDEX.tsv"
COORDINATE_AUDIT = ROOT / "results/ppt_qc/ssDNA_coordinate_audit.json"

PN_ORDER = [f"PN{i}" for i in range(1, 9)]
CHIP_ORDER = ["K8", "J2", "L3"]
DISEASE_ORDER = ["HC", "IPF", "SSc-ILD"]
PN_COLORS = {
    "PN1": "#1f77b4", "PN2": "#9ecae1", "PN3": "#ff7f0e", "PN4": "#fdbf6f",
    "PN5": "#2ca02c", "PN6": "#98df8a", "PN7": "#d62728", "PN8": "#ff9896",
    "UNASSIGNED": "#666666", "NOT_EVALUABLE": "#d0d0d0",
}
ZONE_COLORS = {
    "INTERIOR": "#238b45", "INNER_BOUNDARY": "#ffd92f",
    "IMMEDIATE_OUTSIDE": "#e7298a", "OUTSIDE_RING_0_50": "#41b6c4",
    "OUTSIDE_RING_50_100": "#225ea8", "OTHER": "#bdbdbd",
}
STATE_SHORT = {
    "Alveolar/resident macrophage": "Resident Mac", "Inflammatory monocyte/macrophage": "Inflam Mono/Mac",
    "Dendritic cell": "DC", "Mast/Basophil": "Mast/Baso", "Ciliated epithelial": "Ciliated",
    "AT2/alveolar transitional": "AT2/Transitional", "Secretory/other airway epithelial": "Secretory",
    "Blood endothelial": "Blood Endo", "Perivascular/mesothelial stromal": "Perivasc/Stromal",
    "Lymphatic endothelial": "Lymphatic Endo", "Basal-like epithelial": "Basal-like",
}


def setup() -> logging.Logger:
    for directory in (FIGURES, OUT / "logs", OUT / "report", FIGURES / "ssdna_regions"):
        directory.mkdir(parents=True, exist_ok=True)
    log = logging.getLogger("finalize_delaunay_closure")
    log.setLevel(logging.INFO)
    log.handlers.clear()
    formatter = logging.Formatter("%(asctime)s\t%(levelname)s\t%(message)s")
    for handler in (
        logging.FileHandler(OUT / "logs/1271_finalize_delaunay_closure.log", encoding="utf-8"),
        logging.StreamHandler(sys.stdout),
    ):
        handler.setFormatter(formatter)
        log.addHandler(handler)
    sns.set_theme(style="white", context="paper", font_scale=0.9)
    plt.rcParams.update({"pdf.fonttype": 42, "ps.fonttype": 42, "savefig.facecolor": "white"})
    return log


def import_module(path: Path, name: str):
    spec = importlib.util.spec_from_file_location(name, path)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"Cannot import {path}")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def save_figure(figure: plt.Figure, stem: str) -> Path:
    figure.savefig(FIGURES / f"{stem}.png", dpi=300, bbox_inches="tight", facecolor="white")
    pdf = FIGURES / f"{stem}.pdf"
    figure.savefig(pdf, bbox_inches="tight", facecolor="white")
    plt.close(figure)
    return pdf


def state_columns(frame: pd.DataFrame) -> list[str]:
    return [column.removeprefix("mean_weight__") for column in frame.columns if column.startswith("mean_weight__")]


def build_zone_contrasts(log: logging.Logger) -> pd.DataFrame:
    internal = pd.read_csv(TABLES / "06_REGION_INTERNAL_COMPOSITION.tsv", sep="\t")
    surrounding = pd.read_csv(TABLES / "07_REGION_BOUNDARY_AND_SURROUNDINGS.tsv", sep="\t")
    states = state_columns(internal)
    identity = ["patient_id", "chip", "disease", "neighborhood_label", "local_region_id"]
    all_zones = pd.concat([internal, surrounding], ignore_index=True, sort=False)
    values = all_zones[all_zones.zone_status.eq("EVALUABLE")].pivot_table(
        index=identity, columns="zone", values=[f"mean_weight__{state}" for state in states], aggfunc="first"
    )
    contrasts = {
        "INTERIOR_MINUS_INNER_BOUNDARY": ("INTERIOR", "INNER_BOUNDARY"),
        "REGION_MINUS_IMMEDIATE_OUTSIDE": ("REGION_ALL", "IMMEDIATE_OUTSIDE"),
        "REGION_MINUS_RING_0_50": ("REGION_ALL", "OUTSIDE_RING_0_50"),
        "REGION_MINUS_RING_50_100": ("REGION_ALL", "OUTSIDE_RING_50_100"),
        "REGION_MINUS_COMPONENT_BACKGROUND": ("REGION_ALL", "PATIENT_ROI_BACKGROUND"),
    }
    region_rows = []
    for contrast, (left, right) in contrasts.items():
        for state in states:
            key = f"mean_weight__{state}"
            if (key, left) not in values or (key, right) not in values:
                continue
            difference = (values[(key, left)] - values[(key, right)]).dropna()
            for index, value in difference.items():
                patient, chip, disease, pn, region_id = index
                region_rows.append({
                    "patient_id": patient, "chip": chip, "disease": disease,
                    "neighborhood_label": pn, "local_region_id": region_id,
                    "contrast": contrast, "reference18_state": state,
                    "weight_difference": float(value),
                })
    region_level = pd.DataFrame(region_rows)
    patient_level = (
        region_level.groupby(
            ["patient_id", "chip", "disease", "neighborhood_label", "contrast", "reference18_state"],
            as_index=False, observed=True,
        ).agg(patient_median_difference=("weight_difference", "median"), region_n=("local_region_id", "nunique"))
    )
    summary = (
        patient_level.groupby(["neighborhood_label", "contrast", "reference18_state"], as_index=False, observed=True)
        .agg(
            patient_n=("patient_id", "nunique"), chip_n=("chip", "nunique"),
            median_patient_difference=("patient_median_difference", "median"),
            q1_patient_difference=("patient_median_difference", lambda x: x.quantile(0.25)),
            q3_patient_difference=("patient_median_difference", lambda x: x.quantile(0.75)),
            positive_patient_fraction=("patient_median_difference", lambda x: float((x > 0).mean())),
            total_region_n=("region_n", "sum"),
        )
    )
    summary["direction_consistency"] = np.maximum(
        summary.positive_patient_fraction, 1.0 - summary.positive_patient_fraction
    )
    summary["evidence_scope"] = "descriptive patient-first summary; no region-level inferential test"
    summary.to_csv(TABLES / "07C_REGION_ZONE_CONTRAST_SUMMARY.tsv", sep="\t", index=False)
    region_level.to_csv(TABLES / "07D_REGION_ZONE_CONTRASTS.tsv.gz", sep="\t", index=False, compression="gzip")
    log.info("Wrote %d patient-first zone contrast summaries", len(summary))
    return summary


def figure1_graph_qc(labels: pd.DataFrame, regions: pd.DataFrame) -> Path:
    figure, axes = plt.subplots(2, 2, figsize=(11, 8))
    patient = labels.groupby(["patient_id", "chip_label"], as_index=False).agg(
        bins=("entity_key", "size"), evaluable=("analysis_pn", lambda x: x.isin(PN_ORDER).sum()),
        mean_degree=("delaunay_degree", "mean"),
    )
    sns.barplot(data=patient, x="chip_label", y="bins", order=CHIP_ORDER, color="#4c78a8", ax=axes[0, 0])
    axes[0, 0].set(title="QC Bin50 nodes by chip", xlabel="Chip", ylabel="Bins per patient")
    sns.boxplot(data=labels, x="chip_label", y="delaunay_degree", order=CHIP_ORDER, color="#72b7b2", showfliers=False, ax=axes[0, 1])
    axes[0, 1].set(title="Delaunay degree", xlabel="Chip", ylabel="Degree")
    status = labels.classification_status.value_counts().reindex(["RELIABLE", "LOW_CONFIDENCE", "OOD_UNASSIGNED", "NOT_EVALUABLE"], fill_value=0)
    axes[1, 0].bar(status.index, status.values, color=["#2ca02c", "#ffbf00", "#666666", "#bdbdbd"])
    axes[1, 0].tick_params(axis="x", rotation=25)
    axes[1, 0].set(title="Unified PN classification status", ylabel="Bin50 count")
    strata = regions.size_stratum.value_counts().reindex(["Singleton", "2-4 Bin50", "5-9 Bin50", ">=10 Bin50"], fill_value=0)
    axes[1, 1].bar(strata.index, strata.values, color=["#d9d9d9", "#bdbdbd", "#74a9cf", "#0570b0"])
    axes[1, 1].tick_params(axis="x", rotation=20)
    axes[1, 1].set(title="Connected-region size strata", ylabel="Region count")
    figure.suptitle("Figure 1. Spatial graph and input QC", fontsize=15, fontweight="bold")
    figure.tight_layout()
    return save_figure(figure, "FIGURE_1_SPATIAL_GRAPH_AND_INPUT_QC")


def figure2_composition_marker() -> Path:
    composition = pd.read_csv(COMPOSITION_PATH, sep="\t")
    composition = composition[composition.record_type.eq("global_pn")]
    comp = composition.pivot(index="predicted_pn", columns="reference18_state", values="mean_ego_weight").reindex(PN_ORDER)
    marker = pd.read_csv(MARKER_PATH, sep="\t")
    mark = marker.pivot(index="pn", columns="reference18_state", values="depth_residual_mean_effect").reindex(PN_ORDER)
    states = comp.columns.tolist()
    figure, axes = plt.subplots(2, 1, figsize=(15, 7.5), gridspec_kw={"height_ratios": [1, 1]})
    sns.heatmap(comp[states], cmap="viridis", ax=axes[0], cbar_kws={"label": "Mean ego weight"})
    axes[0].set(title="Reference18 ego composition", xlabel="", ylabel="Unified neighborhood")
    sns.heatmap(mark.reindex(columns=states), cmap="vlag", center=0, ax=axes[1], cbar_kws={"label": "Depth-residual marker effect"})
    axes[1].set(title="Marker concordance (internal RNA evidence)", xlabel="Reference18 state", ylabel="Unified neighborhood")
    for axis in axes:
        axis.set_xticklabels([STATE_SHORT.get(x.get_text(), x.get_text()) for x in axis.get_xticklabels()], rotation=40, ha="right")
    figure.suptitle("Figure 2. Neighborhood composition and marker concordance", fontsize=15, fontweight="bold")
    figure.tight_layout()
    return save_figure(figure, "FIGURE_2_NEIGHBORHOOD_COMPOSITION_AND_MARKER_CONCORDANCE")


def figure3_spatial_maps(labels: pd.DataFrame) -> Path:
    order = labels[["chip_label", "disease", "patient_id"]].drop_duplicates().copy()
    order["chip_i"] = order.chip_label.map({v: i for i, v in enumerate(CHIP_ORDER)})
    order["disease_i"] = order.disease.map({v: i for i, v in enumerate(DISEASE_ORDER)})
    order = order.sort_values(["chip_i", "disease_i", "patient_id"])
    figure, axes = plt.subplots(3, 7, figsize=(18, 8.5))
    for axis, row in zip(axes.ravel(), order.itertuples(index=False)):
        frame = labels[labels.patient_id.eq(row.patient_id)]
        colors = frame.analysis_pn.map(PN_COLORS).fillna(PN_COLORS["UNASSIGNED"])
        axis.scatter(frame.x_um, frame.y_um, c=colors, s=2.2, marker="s", linewidths=0, rasterized=True)
        axis.set_aspect("equal"); axis.invert_yaxis(); axis.set_xticks([]); axis.set_yticks([])
        axis.set_title(f"{row.patient_id}\n{row.chip_label}, {row.disease}", fontsize=7)
    handles = [Line2D([0], [0], marker="s", linestyle="", color=PN_COLORS[pn], label=pn) for pn in PN_ORDER]
    figure.legend(handles=handles, loc="lower center", ncol=8, frameon=False)
    figure.suptitle("Figure 3. Unified spatial neighborhoods across 21 patients", fontsize=15, fontweight="bold")
    figure.tight_layout(rect=(0, 0.05, 1, 0.96))
    return save_figure(figure, "FIGURE_3_SPATIAL_NEIGHBORHOOD_MAPS_ACROSS_PATIENTS")


def figure4_regions(regions: pd.DataFrame) -> Path:
    figure, axes = plt.subplots(1, 3, figsize=(14, 4.6))
    sns.boxplot(data=regions, x="neighborhood_label", y="bin_n", order=PN_ORDER, showfliers=False, color="#80b1d3", ax=axes[0])
    axes[0].set_yscale("log"); axes[0].set(title="Region size distribution", xlabel="Neighborhood", ylabel="Bin50 count (log scale)")
    strata_order = ["Singleton", "2-4 Bin50", "5-9 Bin50", ">=10 Bin50"]
    counts = regions.groupby(["neighborhood_label", "size_stratum"]).size().unstack(fill_value=0).reindex(PN_ORDER)[strata_order]
    proportions = counts.div(counts.sum(axis=1), axis=0)
    proportions.plot(kind="bar", stacked=True, color=["#eeeeee", "#cccccc", "#74a9cf", "#0570b0"], ax=axes[1])
    axes[1].set(title="Region-count fragmentation", xlabel="Neighborhood", ylabel="Fraction of regions")
    axes[1].legend(title="Size", fontsize=7)
    large = regions.assign(in_large=regions.bin_n.ge(10)).groupby("neighborhood_label").apply(
        lambda x: pd.Series({"fraction": x.loc[x.in_large, "bin_n"].sum() / x.bin_n.sum()}), include_groups=False
    ).reindex(PN_ORDER)
    axes[2].bar(large.index, large["fraction"], color="#2b8cbe")
    axes[2].set(title="Bins in regions >=10 Bin50", xlabel="Neighborhood", ylabel="Fraction of PN bins", ylim=(0, 1))
    figure.suptitle("Figure 4. Local connected-region size and fragmentation", fontsize=15, fontweight="bold")
    figure.tight_layout()
    return save_figure(figure, "FIGURE_4_LOCAL_CONNECTED_REGION_SIZE_AND_FRAGMENTATION")


def load_ssdna_transforms(index: pd.DataFrame) -> tuple[dict[str, np.ndarray], dict[str, Path]]:
    audit = json.loads(COORDINATE_AUDIT.read_text(encoding="utf-8"))
    transforms = {}
    image_paths = {}
    for chip, chip_id in {"K8": "Y40102K8", "J2": "Y40105J2", "L3": "Y40105L3"}.items():
        raw_to_registered = np.asarray(audit["chips"][chip_id]["ipr_transform"]["raw_to_registered_bin1_matrix"], float)
        transforms[chip] = np.linalg.inv(raw_to_registered)
        image_paths[chip] = Path(index.loc[index.chip.eq(chip), "ssDNA_source_file"].iloc[0])
    return transforms, image_paths


def transform_points(points: np.ndarray, matrix: np.ndarray) -> np.ndarray:
    hom = np.column_stack([points, np.ones(len(points))])
    return (matrix @ hom.T).T[:, :2]


def polygons(frame: pd.DataFrame, matrix: np.ndarray) -> np.ndarray:
    offsets = np.asarray([[-25, -25], [25, -25], [25, 25], [-25, 25]], float)
    raw = frame[["x_raw", "y_raw"]].to_numpy(float)[:, None, :] + offsets[None, :, :]
    return transform_points(raw.reshape(-1, 2), matrix).reshape(len(frame), 4, 2)


def boundary_and_zones(region_indices: np.ndarray, labels: pd.DataFrame, adjacency: sp.csr_matrix) -> dict[str, np.ndarray]:
    membership = np.zeros(len(labels), dtype=bool); membership[region_indices] = True
    boundary, immediate = [], []
    for node in region_indices:
        neighbors = adjacency.indices[adjacency.indptr[node]:adjacency.indptr[node + 1]]
        outside = neighbors[~membership[neighbors]]
        if len(outside):
            boundary.append(node); immediate.extend(outside.tolist())
    boundary = np.asarray(boundary, int); immediate = np.unique(np.asarray(immediate, int))
    interior = np.setdiff1d(region_indices, boundary)
    first = labels.loc[region_indices[0]]
    same_component = labels[
        labels.patient_id.eq(first.patient_id) & labels.roi_id.eq(first.roi_id)
        & labels.tissue_component_id.eq(first.tissue_component_id)
    ].index.to_numpy(int)
    outside = np.setdiff1d(same_component, region_indices)
    if len(boundary) and len(outside):
        tree = cKDTree(labels.loc[boundary, ["x_um", "y_um"]].to_numpy(float))
        points = labels.loc[outside, ["x_um", "y_um"]].to_numpy(float)
        _, nearest = tree.query(points)
        closest = labels.loc[boundary[nearest], ["x_um", "y_um"]].to_numpy(float)
        delta = np.maximum(np.abs(points - closest) - 12.5, 0)
        distance = np.sqrt((delta ** 2).sum(axis=1))
    else:
        distance = np.full(len(outside), np.inf)
    is_immediate = np.isin(outside, immediate)
    ring0 = outside[(distance <= 50) & ~is_immediate]
    ring50 = outside[(distance > 50) & (distance <= 100)]
    return {"INTERIOR": interior, "INNER_BOUNDARY": boundary, "IMMEDIATE_OUTSIDE": immediate,
            "OUTSIDE_RING_0_50": ring0, "OUTSIDE_RING_50_100": ring50}


def figure5_ssdna(labels: pd.DataFrame, adjacency: sp.csr_matrix, regions: pd.DataFrame, log: logging.Logger) -> tuple[Path, pd.DataFrame]:
    index = pd.read_csv(SSDNA_INDEX, sep="\t")
    transforms, image_paths = load_ssdna_transforms(index)
    zone_qc = pd.read_csv(TABLES / "07B_REGION_ZONE_EVALUABILITY.tsv", sep="\t")
    eligible = regions.merge(zone_qc[["local_region_id", "full_zone_analysis_status"]], on="local_region_id", how="left")
    eligible = eligible[eligible.patient_id.isin(index.patient_id) & eligible.bin_n.ge(10)].copy()
    eligible["full_score"] = eligible.full_zone_analysis_status.eq("EVALUABLE").astype(int)
    eligible["quality_score"] = eligible.full_score * 1_000_000 + eligible.reliable_bin_fraction * 100_000 + np.minimum(eligible.bin_n, 5000)
    selected = eligible.sort_values(["patient_id", "quality_score"], ascending=[True, False]).groupby("patient_id", as_index=False).head(1)
    selected.to_csv(TABLES / "08B_REPRESENTATIVE_SSDNA_REGION_SELECTION.tsv", sep="\t", index=False)
    images = {}
    Image.MAX_IMAGE_PIXELS = None
    for chip, path in image_paths.items():
        images[chip] = Image.open(path)
    figure, axes = plt.subplots(3, 3, figsize=(13, 12))
    rows = selected.merge(index[["patient_id", "chip", "disease"]], left_on=["patient_id", "chip_label", "disease"], right_on=["patient_id", "chip", "disease"], how="left")
    rows["chip_i"] = rows.chip_label.map({v: i for i, v in enumerate(CHIP_ORDER)})
    rows["disease_i"] = rows.disease.map({v: i for i, v in enumerate(DISEASE_ORDER)})
    rows = rows.sort_values(["chip_i", "disease_i", "patient_id"])
    for axis, row in zip(axes.ravel(), rows.itertuples(index=False)):
        rid = row.local_region_id
        region_indices = labels.index[labels.component_id.eq(rid)].to_numpy(int)
        zones = boundary_and_zones(region_indices, labels, adjacency)
        all_indices = np.unique(np.concatenate(list(zones.values())))
        frame = labels.loc[all_indices].copy()
        matrix = transforms[row.chip_label]
        centers = transform_points(frame[["x_raw", "y_raw"]].to_numpy(float), matrix)
        xmin, ymin = np.floor(centers.min(axis=0) - 150).astype(int)
        xmax, ymax = np.ceil(centers.max(axis=0) + 150).astype(int)
        image = images[row.chip_label]
        xmin, ymin = max(0, xmin), max(0, ymin); xmax, ymax = min(image.width, xmax), min(image.height, ymax)
        background = np.asarray(ImageOps.autocontrast(image.crop((xmin, ymin, xmax, ymax)).convert("L"), cutoff=0.5))
        axis.imshow(background, cmap="gray", vmin=0, vmax=255, origin="upper", extent=(xmin, xmax, ymax, ymin), interpolation="nearest")
        for zone in ["OUTSIDE_RING_50_100", "OUTSIDE_RING_0_50", "IMMEDIATE_OUTSIDE", "INNER_BOUNDARY", "INTERIOR"]:
            idx = zones[zone]
            if not len(idx):
                continue
            collection = PolyCollection(polygons(labels.loc[idx], matrix), facecolors=ZONE_COLORS[zone], edgecolors="none", alpha=0.55, rasterized=True)
            axis.add_collection(collection)
        region_poly = polygons(labels.loc[region_indices], matrix)
        axis.add_collection(PolyCollection(region_poly, facecolors="none", edgecolors="white", linewidths=0.4, alpha=0.8, rasterized=True))
        axis.set_xlim(xmin, xmax); axis.set_ylim(ymax, ymin); axis.set_aspect("equal")
        axis.set_xticks([]); axis.set_yticks([])
        axis.set_title(f"{row.patient_id} | {row.chip_label} | {row.disease}\n{row.neighborhood_label}, {rid}, n={row.bin_n}", fontsize=8)
    handles = [Line2D([0], [0], marker="s", linestyle="", color=ZONE_COLORS[z], label=z.replace("_", " ").title()) for z in ZONE_COLORS if z != "OTHER"]
    figure.legend(handles=handles, loc="lower center", ncol=5, frameon=False)
    figure.suptitle("Figure 5. Connected regions overlaid on validated ssDNA", fontsize=15, fontweight="bold")
    figure.tight_layout(rect=(0, 0.04, 1, 0.96))
    pdf = save_figure(figure, "FIGURE_5_CONNECTED_REGIONS_OVERLAID_ON_SSDNA")
    for image in images.values():
        image.close()
    log.info("Rendered %d representative ssDNA regions", len(rows))
    return pdf, selected


def figure6_zone_composition(contrast: pd.DataFrame) -> Path:
    subset = contrast[contrast.contrast.eq("REGION_MINUS_COMPONENT_BACKGROUND")].copy()
    matrix = subset.pivot(index="neighborhood_label", columns="reference18_state", values="median_patient_difference").reindex(PN_ORDER)
    figure, axes = plt.subplots(2, 1, figsize=(15, 8), gridspec_kw={"height_ratios": [1.0, 1.15]})
    sns.heatmap(matrix, cmap="vlag", center=0, ax=axes[0], cbar_kws={"label": "Median patient-level weight difference"})
    axes[0].set(title="Region composition minus component background", xlabel="Reference18 state", ylabel="Neighborhood")
    ring = contrast[contrast.contrast.isin(["REGION_MINUS_IMMEDIATE_OUTSIDE", "REGION_MINUS_RING_0_50", "REGION_MINUS_RING_50_100"])].copy()
    ring["row"] = ring.neighborhood_label + " | " + ring.contrast.str.replace("REGION_MINUS_", "", regex=False)
    ring_matrix = ring.pivot(index="row", columns="reference18_state", values="median_patient_difference")
    row_order = [f"{pn} | {suffix}" for pn in PN_ORDER for suffix in ["IMMEDIATE_OUTSIDE", "RING_0_50", "RING_50_100"]]
    ring_matrix = ring_matrix.reindex(row_order)
    sns.heatmap(ring_matrix, cmap="vlag", center=0, ax=axes[1], cbar_kws={"label": "Median patient-level weight difference"})
    axes[1].set(title="Region composition relative to surrounding zones", xlabel="Reference18 state", ylabel="Region contrast")
    for axis in axes:
        axis.set_xticklabels([STATE_SHORT.get(x.get_text(), x.get_text()) for x in axis.get_xticklabels()], rotation=40, ha="right")
    figure.suptitle("Figure 6. Region interior, boundary, and surrounding composition", fontsize=15, fontweight="bold")
    figure.tight_layout()
    return save_figure(figure, "FIGURE_6_REGION_INTERIOR_BOUNDARY_AND_SURROUNDING_COMPOSITION")


def figure7_nmf() -> Path:
    nmf = pd.read_csv(TABLES / "09_DELAUNAY_VS_NMF_SUMMARY.tsv", sep="\t")
    agreement = nmf[nmf.record_type.eq("BIN_LEVEL_DOMINANT_COMPARISON")]
    usage = nmf[nmf.record_type.eq("PN_MEAN_NMF_USAGE")]
    matrix = usage.pivot_table(index=["chip", "neighborhood_label"], columns="nmf_factor_family", values="value").sort_index()
    figure, axes = plt.subplots(1, 2, figsize=(14, 5.5), gridspec_kw={"width_ratios": [1, 3]})
    positions = np.arange(len(agreement)); width = 0.34
    axes[0].bar(positions - width / 2, agreement.ari, width, label="ARI", color="#4c78a8")
    axes[0].bar(positions + width / 2, agreement.nmi, width, label="NMI", color="#f58518")
    axes[0].set_xticks(positions, agreement.chip); axes[0].set_ylim(0, 0.08)
    axes[0].set(title="Dominant-label agreement", xlabel="Chip", ylabel="Score"); axes[0].legend(frameon=False)
    sns.heatmap(matrix, cmap="mako", ax=axes[1], cbar_kws={"label": "Mean frozen NMF usage"})
    axes[1].set(title="Continuous NMF K7 usage within unified PN bins", xlabel="Frozen NMF factor family", ylabel="Chip and PN")
    axes[1].tick_params(axis="x", rotation=35)
    figure.suptitle("Figure 7. Delaunay and frozen NMF K7 comparison", fontsize=15, fontweight="bold")
    figure.tight_layout()
    return save_figure(figure, "FIGURE_7_DELAUNAY_AND_FROZEN_NMF_K7_COMPARISON")


def figure8_patient_summary() -> Path:
    summary = pd.read_csv(TABLES / "10_CHIP_AWARE_PATIENT_SUMMARY.tsv", sep="\t")
    patient = summary[summary.record_type.eq("PATIENT_PN_SUMMARY")].copy()
    matrix = patient.pivot(index="patient_id", columns="neighborhood_label", values="bin_fraction")
    order = patient[["patient_id", "chip_label", "disease"]].drop_duplicates().copy()
    order["chip_i"] = order.chip_label.map({v: i for i, v in enumerate(CHIP_ORDER)})
    order["disease_i"] = order.disease.map({v: i for i, v in enumerate(DISEASE_ORDER)})
    order = order.sort_values(["chip_i", "disease_i", "patient_id"])
    matrix = matrix.reindex(order.patient_id).reindex(columns=PN_ORDER)
    figure, axes = plt.subplots(1, 2, figsize=(12.5, 7), gridspec_kw={"width_ratios": [1.6, 1]})
    sns.heatmap(matrix, cmap="Blues", ax=axes[0], cbar_kws={"label": "Fraction of patient bins"})
    axes[0].set(title="Patient-level unified PN composition", xlabel="Neighborhood", ylabel="Patient")
    plot_data = patient.copy()
    sns.boxplot(data=plot_data, x="neighborhood_label", y="region_ge10_area_fraction", hue="chip_label", hue_order=CHIP_ORDER, order=PN_ORDER, showfliers=False, ax=axes[1])
    axes[1].set(title="Large-region area fraction by chip", xlabel="Neighborhood", ylabel="Fraction of patient bins")
    axes[1].legend(title="Chip", frameon=False)
    figure.suptitle("Figure 8. Chip-aware patient-level summary", fontsize=15, fontweight="bold")
    figure.tight_layout()
    return save_figure(figure, "FIGURE_8_CHIP_AWARE_PATIENT_LEVEL_SUMMARY")


def top_composition_text(regions: pd.DataFrame) -> list[str]:
    states = state_columns(regions)
    rows = []
    for pn in PN_ORDER:
        frame = regions[regions.neighborhood_label.eq(pn)]
        means = {state: float(np.average(frame[f"mean_weight__{state}"], weights=frame.bin_n)) for state in states}
        top = sorted(means, key=means.get, reverse=True)[:4]
        rows.append(f"- {pn}: " + ", ".join(f"{state} ({means[state]:.3f})" for state in top))
    return rows


def write_reports(regions: pd.DataFrame, labels: pd.DataFrame, contrast: pd.DataFrame, selected: pd.DataFrame) -> None:
    zone_qc = pd.read_csv(TABLES / "07B_REGION_ZONE_EVALUABILITY.tsv", sep="\t")
    classification = pd.read_csv(TABLES / "03_NEIGHBORHOOD_CLASSIFICATION_SUMMARY.tsv", sep="\t")
    patient = pd.read_csv(TABLES / "10_CHIP_AWARE_PATIENT_SUMMARY.tsv", sep="\t")
    disease_models = patient[patient.record_type.eq("CHIP_ADJUSTED_DISEASE_DESCRIPTION")]
    nmf = pd.read_csv(TABLES / "09_DELAUNAY_VS_NMF_SUMMARY.tsv", sep="\t")
    agreement = nmf[nmf.record_type.eq("BIN_LEVEL_DOMINANT_COMPARISON")]
    top_contrasts = contrast[
        contrast.patient_n.ge(9) & contrast.direction_consistency.ge(0.75)
    ].assign(abs_effect=lambda x: x.median_patient_difference.abs()).sort_values("abs_effect", ascending=False).head(12)
    contrast_lines = [
        f"- {r.neighborhood_label} | {r.contrast} | {r.reference18_state}: median patient-level difference={r.median_patient_difference:+.3f}, direction consistency={r.direction_consistency:.0%}, n={int(r.patient_n)} patients."
        for r in top_contrasts.itertuples(index=False)
    ] or ["- 没有达到预设描述性筛选条件（至少9位患者且方向一致性≥75%）的区域-周围模式。"]
    main = f"""# Delaunay Spatial Neighborhood 最终收尾结论

## 最终状态

**DELAUNAY_ANALYSIS_CLOSURE_READY**

冻结Delaunay图、统一PN1–PN8标签、技术标志和既有K选择结果均未改变。本轮没有运行RCTD、NMF、K-means、edgeR、GSEA、LR、空间共现或新批次校正。

## 完成度

- 126,449个QC Bin50全部纳入版本审计；125,471个具有统一PN标签并进入连通区域分析。
- 冻结图包含357,959条无向边，按patient、ROI和tissue component限制，最大边长35.36 μm。
- 共保留{len(regions):,}个connected regions：Singleton {int((regions.bin_n==1).sum()):,}，2–4 Bin50 {int(regions.bin_n.between(2,4).sum()):,}，5–9 Bin50 {int(regions.bin_n.between(5,9).sum()):,}，≥10 Bin50 {int(regions.bin_n.ge(10).sum()):,}。
- {int(zone_qc.full_zone_analysis_status.eq('EVALUABLE').sum()):,}/{len(zone_qc):,}个≥10 Bin50区域具有全部预设zone；其余保留为PARTIAL_ZONE_COVERAGE。
- ssDNA配准仅对9/21位患者有已审核映射。本轮实际生成9个代表区域overlay；其余患者仅使用真实Bin50坐标，不伪造组织学配准。

## 可保留的PN组成描述

Reference18为混合Bin50的连续权重。以下为按区域Bin数加权的主要组成，不是单细胞计数，也不是独立病理验证：

{chr(10).join(top_composition_text(regions))}

全部PN在3张芯片均出现。具有≥10 Bin50区域的患者覆盖为：{'; '.join(f'{r.neighborhood_label} {int(r.patient_coverage_ge10_region)}/21' for r in classification.itertuples(index=False))}。PN标签仍应保持中性命名。

## Local region与周围组成

区域严格由同一patient、ROI、tissue component、PN及真实Delaunay边定义。面积为Bin50数×625 μm²，不填补空洞。以下模式通过患者内region中位数后再跨患者汇总，属于描述性证据：

{chr(10).join(contrast_lines)}

大型区域的外环高度重叠：{zone_qc.ring_overlap_present.mean():.1%}的区域至少一个环带与其他区域环带共享Bin。因此region-level环带不能作为独立生物学重复；本报告不对环带做伪重复显著性检验。

## ssDNA、芯片与疾病解释

- 9位已有坐标审计的患者可以在ssDNA上定位connected regions。图中保持25×25 μm tile，不做平滑。
- 其余12位没有可靠ssDNA映射，因此不能声称完成21例组织学验证。
- 既有特征空间审计显示：疾病调整后chip partial R²约0.598；进一步调整UMI、nGenes和graph degree后仍约0.408。该差异可能同时包含技术批次、取材和患者组成异质性。
- 患者级、chip调整的探索性疾病比较共有{len(disease_models)}个预设contrast；FDR<0.05为{int((disease_models.bh_fdr<0.05).sum())}。因此本轮没有可正式宣称的疾病相关PN/region差异。

## 与冻结NMF K=7的关系

按dominant label比较，ARI为{', '.join(f'{r.chip}={r.ari:.3f}' for r in agreement.itertuples(index=False))}，NMI为{', '.join(f'{r.chip}={r.nmi:.3f}' for r in agreement.itertuples(index=False))}。两种方法使用同一Reference18来源，低一致性说明分辨率与目标不同，不构成独立验证。Delaunay增加的是同标签真实邻接region及其局部边界/周围描述；NMF保留连续factor usage。

## 论文可用范围

- 可正式保留：图构建/QC、固定标签的可复现映射、全部connected-component清单、区域大小/碎片化、患者级覆盖、9例已验证ssDNA定位。
- 可作为空间描述：PN组成、≥10 Bin50区域、内部/边界/周围组成及Delaunay-NMF对应。
- 仅探索性：跨疾病比较、跨芯片PN比例、特定生物学niche命名、外环差异。
- 不支持：把所有component称为aggregate、把PN称为独立病理niche、把环带当独立重复、宣称21例均完成组织学验证。

## 收尾判断

Delaunay路线的预定收尾工作已完成，没有必须通过再次选K或重新聚类解决的缺口。剩余问题是独立组织学验证、未配准的12例ssDNA、批次/组成解耦及更稳健的空间表征，适合作为BANKSY阶段的明确交接问题。本任务到此停止，不自动启动BANKSY。
"""
    methods = f"""# 方法与局限

## 冻结输入

- G15/UMI50 Bin50，21位患者，Reference18连续RCTD权重。
- 每位patient、ROI、tissue component独立构建distance-limited Delaunay图；最大边长35.36 μm。
- 使用既有统一PN1–PN8、classification confidence、OOD_UNASSIGNED和NOT_EVALUABLE。没有重新聚类或选K。

## Connected region定义

同一PN Bin50仅在存在真实Delaunay边时连通。全部{len(regions):,}个component均保留。Singleton、2–4、5–9和≥10 Bin50分别报告；“connected region”不等于“aggregate”。面积定义为Bin50数×625 μm²，组织空洞不填补。

## Interior、boundary与surrounding

- Inner boundary：region内部且至少有一个图邻居位于region外。
- Interior：region中非inner-boundary的Bin50。
- Immediate outside：通过一条Delaunay边直接连接到region的外部Bin50。
- 0–50 μm与50–100 μm ring：从region boundary tile边缘计算最近近似欧氏距离；immediate outside从0–50 μm ring中排除。
- 所有zone限制在同一patient、ROI和tissue component；不同focal region的ring允许重叠但明确标记。
- ≥5 Bin50区域计算region/interior/boundary；≥10 Bin50区域再计算surrounding。zone少于3 Bin（background少于10 Bin）标记NOT_EVALUABLE。

## 汇总与统计单位

组成差异先在每个region计算，再在患者内取region中位数，最后跨患者描述。没有把Bin50或重叠ring作为独立患者。疾病比较仅使用患者级指标，并以chip作为协变量；BH-FDR覆盖预设全部比较。

## 组织学与NMF

仅对既有审计通过的9位患者复用芯片特异仿射变换。每个tile保持25×25 μm。不对其余12位构造伪配准。冻结NMF K=7仅加载既有usage进行描述性对照，没有重新拟合。

## 关键局限

1. PN与marker、NMF均源于同一RNA/Reference18体系，不能作为独立生物学验证。
2. chip-associated variation较强，统一PN的跨芯片可比性仅适合带chip信息的描述。
3. {zone_qc.ring_overlap_present.mean():.1%}的大型区域存在ring重叠，区域周围结果高度相关。
4. 只有{selected.patient_id.nunique()}/21位患者具有可复用的已验证ssDNA配准。
5. 连通性受固定图和离散PN标签影响；碎片化不应通过事后平滑隐藏。
6. 本轮未做新的空间共现、LR、edgeR、GSEA、批次校正或疾病发现分析。
"""
    handoff = """# BANKSY交接说明

## 可直接继承的对象

- 冻结Reference18连续权重和126,449个QC Bin50。
- patient、chip、disease、ROI、tissue component、真实坐标及UMI/nGenes。
- 357,959条经审核Delaunay边及degree/boundary QC。
- 统一PN1–PN8、confidence、OOD/NOT_EVALUABLE标志。
- 全部29,750个connected regions及固定的area/fragmentation定义。
- 9位患者经审核的ssDNA仿射配准。

## BANKSY应优先回答的问题

1. 在不依赖当前PN硬标签的情况下，空间平滑表征是否能减少碎片化且保留真实边界？
2. chip variation在空间表征中是否仍占主导？必须按患者留出与chip敏感性报告，不能用训练内稳定性代替验证。
3. BANKSY识别的区域是否与Delaunay ≥10 Bin50 regions在位置、组成和尺度上对应？
4. 当前仅505个大型region具备全部zone，BANKSY是否能定义更稳健、非重叠或层级化的boundary/surrounding结构？
5. 9例ssDNA支持的形态位置能否复现；其余12例必须先获得可靠配准才能做组织学结论。

## 建议的比较基线

- 以patient为评估单位；保持相同QC、Reference18、物理坐标和组织边界。
- 与Delaunay比较：患者覆盖、composition separation、同标签邻接、singleton/2–4比例、≥10区域Bin比例、UMI/nGenes/degree敏感性和chip依赖。
- 保留原始BANKSY标签，不根据疾病P值或预期免疫故事选参数。
- 不要求BANKSY复现PN1–PN8编号；使用组成与空间特征建立对应，并允许无可靠对应。

## 启动门槛

本次Delaunay收尾状态为 **DELAUNAY_ANALYSIS_CLOSURE_READY**。这表示Delaunay结果、限制和比较基线已整理完毕，并不表示BANKSY已获自动授权启动。本任务结束后停止。
"""
    (OUT / "12_MAIN_FINDINGS_ZH.md").write_text(main, encoding="utf-8")
    (OUT / "13_METHODS_AND_LIMITATIONS_ZH.md").write_text(methods, encoding="utf-8")
    (OUT / "14_BANKSY_HANDOFF_ZH.md").write_text(handoff, encoding="utf-8")
    for name, content in [("12_MAIN_FINDINGS_ZH.md", main), ("13_METHODS_AND_LIMITATIONS_ZH.md", methods), ("14_BANKSY_HANDOFF_ZH.md", handoff)]:
        (OUT / "report" / name).write_text(content, encoding="utf-8")


def merge_atlas(pdf_paths: list[Path]) -> Path:
    atlas = FIGURES / "FINAL_DELAUNAY_ATLAS.pdf"
    writer = PdfWriter()
    for pdf in pdf_paths:
        reader = PdfReader(pdf)
        for page in reader.pages:
            writer.add_page(page)
    writer.add_metadata({"/Title": "Final Delaunay Spatial Neighborhood Atlas", "/Author": "Reproducible analysis pipeline"})
    with atlas.open("wb") as handle:
        writer.write(handle)
    return atlas


def main() -> None:
    log = setup()
    compute = import_module(COMPUTE_SCRIPT, "closure_compute")
    data, labels, weight_columns, graph_pairs, graph_qc, adjacency = compute.load_core(log)
    labels = labels.copy()
    for column in ("x_raw", "y_raw"):
        labels[column] = data[column].to_numpy()
    regions = pd.read_csv(TABLES / "04_ALL_LOCAL_CONNECTED_REGIONS.tsv.gz", sep="\t")
    contrast = build_zone_contrasts(log)
    pdfs = [
        figure1_graph_qc(labels, regions),
        figure2_composition_marker(),
        figure3_spatial_maps(labels),
        figure4_regions(regions),
    ]
    figure5, selected = figure5_ssdna(labels, adjacency, regions, log)
    pdfs.extend([figure5, figure6_zone_composition(contrast), figure7_nmf(), figure8_patient_summary()])
    write_reports(regions, labels, contrast, selected)
    atlas = merge_atlas(pdfs)
    result = {
        "status": "DELAUNAY_ANALYSIS_CLOSURE_READY", "figure_pdf_n": len(pdfs),
        "atlas_pages": len(PdfReader(atlas).pages), "representative_ssdna_regions": len(selected),
        "region_n": len(regions), "region_ge10_n": int(regions.bin_n.ge(10).sum()),
    }
    (OUT / "logs/FINALIZATION_SUMMARY.json").write_text(json.dumps(result, indent=2), encoding="utf-8")
    log.info("Finalization complete: %s", result)


if __name__ == "__main__":
    main()
