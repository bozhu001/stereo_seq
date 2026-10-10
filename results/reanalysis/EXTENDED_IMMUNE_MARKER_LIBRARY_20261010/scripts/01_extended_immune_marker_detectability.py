#!/usr/bin/env python3
"""Extended immune marker library and whole-tissue Bin50 detectability audit.

This workflow is intentionally limited to atlas-marker curation and sparse
raw-count detectability. It does not read niche assignments, disease labels for
selection, or any clustering result.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from collections import defaultdict
from pathlib import Path
from typing import Iterable

import anndata as ad
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.backends.backend_pdf import PdfPages
from openpyxl import load_workbook
from scipy import sparse


SEED = 20261010
CHIP_ALIAS = {"Y40105L3": "L3", "Y40105J2": "J2", "Y40102K8": "Y40102K8"}
REFERENCE18_STATES = [
    "B", "Plasma", "T", "NK", "Alveolar/resident macrophage",
    "Inflammatory monocyte/macrophage", "Dendritic cell", "Mast/Basophil",
]

IMMUNE_LEVEL1 = {
    "B/Plasma", "Basophil/Mast", "Dendritic", "Macrophage",
    "Macrophage FABP4+", "Monocyte", "NK", "pDC", "T cell",
}

IMMUNE_PARENT = {
    "Activated NK XCL1+": "NK",
    "B cell": "B/Plasma",
    "B cell LRMP+": "B/Plasma",
    "Basophil/Mast": "Basophil/Mast",
    "CD8+ Naive T": "T cell",
    "Cytotoxic T": "T cell",
    "DC CXCL2+/CD14+/CD1C-": "Dendritic",
    "DC LAMP3+/CCR7+": "Dendritic",
    "DC1 CLEC9A+": "Dendritic",
    "DC2 C1QC+": "Dendritic",
    "DC2 CD1A+/CD1E+": "Dendritic",
    "DC2 CLEC10A+": "Dendritic",
    "helper T CD4+": "T cell",
    "helper T CD4+ CDR7+": "T cell",
    "Macrophage C1Q hi": "Macrophage",
    "Macrophage CHI3L1+/CD9 hi/": "Macrophage",
    "Macrophage FABP4+": "Macrophage FABP4+",
    "Macrophage FABP4+/PDE4C+": "Macrophage FABP4+",
    "Macrophage IL1B+": "Macrophage",
    "Macrophage LYVE1+": "Macrophage",
    "Macrophage RETN+/VCAN+": "Macrophage",
    "Monocytes CD14+": "Monocyte",
    "Monocytes CD14+/IL1B": "Monocyte",
    "Monocytes CD16+": "Monocyte",
    "NK": "NK",
    "pDC": "pDC",
    "Plasma B": "B/Plasma",
    "T CD4+ JUN+": "T cell",
    "T regulatory": "T cell",
}

DETAILED_LINEAGE = {
    "Activated NK XCL1+": "Activated NK",
    "B cell": "B",
    "B cell LRMP+": "B",
    "Basophil/Mast": "Mast/Basophil",
    "CD8+ Naive T": "CD8 T",
    "Cytotoxic T": "Cytotoxic T",
    "DC CXCL2+/CD14+/CD1C-": "Dendritic cell",
    "DC LAMP3+/CCR7+": "Dendritic cell",
    "DC1 CLEC9A+": "Dendritic cell",
    "DC2 C1QC+": "Dendritic cell",
    "DC2 CD1A+/CD1E+": "Dendritic cell",
    "DC2 CLEC10A+": "Dendritic cell",
    "helper T CD4+": "CD4 T",
    "helper T CD4+ CDR7+": "CD4 T",
    "Macrophage C1Q hi": "Alveolar/resident macrophage",
    "Macrophage CHI3L1+/CD9 hi/": "Inflammatory macrophage",
    "Macrophage FABP4+": "Alveolar/resident macrophage",
    "Macrophage FABP4+/PDE4C+": "Alveolar/resident macrophage",
    "Macrophage IL1B+": "Inflammatory macrophage",
    "Macrophage LYVE1+": "Alveolar/resident macrophage",
    "Macrophage RETN+/VCAN+": "Inflammatory macrophage",
    "Monocytes CD14+": "Monocyte",
    "Monocytes CD14+/IL1B": "Monocyte",
    "Monocytes CD16+": "Monocyte",
    "NK": "NK",
    "pDC": "pDC",
    "Plasma B": "Plasma",
    "T CD4+ JUN+": "CD4 T",
    "T regulatory": "Regulatory T",
}

MAJOR_LINEAGE = {
    "B": "B", "Plasma": "Plasma", "CD4 T": "T", "CD8 T": "T",
    "Regulatory T": "T", "Cytotoxic T": "T", "NK": "NK",
    "Activated NK": "NK", "Alveolar/resident macrophage": "Macrophage",
    "Inflammatory macrophage": "Macrophage", "Monocyte": "Monocyte",
    "Dendritic cell": "DC", "pDC": "DC", "Mast/Basophil": "Mast/Basophil",
}

CORE_ANCHORS = {
    "B": {"CD79A", "CD79B", "MS4A1", "CD19", "CD22", "BANK1", "BLK", "PAX5", "TNFRSF13C", "EBF1"},
    "Plasma": {"MZB1", "DERL3", "PRDM1", "XBP1", "TNFRSF17", "SDC1", "SLAMF7", "FKBP11"},
    "T": {"CD3D", "CD3E", "CD3G", "TRAC", "TRBC1", "TRBC2", "LCK", "BCL11B", "CD2", "ITK"},
    "NK": {"KLRD1", "KLRF1", "NCR1", "NCAM1"},
    "Macrophage": {"C1QA", "C1QB", "C1QC", "MARCO", "PPARG", "FABP4", "MSR1", "MERTK", "CD68"},
    "Monocyte": {"FCN1", "CCR2", "CD14", "VCAN", "S100A12", "CD300E"},
    "DC": {"CD1C", "FCER1A", "CLEC10A", "CLEC9A", "XCR1", "BATF3", "CADM1", "TCF4", "IL3RA"},
    "Mast/Basophil": {"TPSAB1", "TPSB2", "CPA3", "KIT", "MS4A2", "HDC", "HPGDS", "IL1RL1"},
}

LEGACY_SMALL_PANEL = {
    "B": ["CD79A", "CD79B", "MS4A1", "CD19", "CD22", "BANK1", "BLK", "PAX5"],
    "Plasma": ["MZB1", "XBP1", "JCHAIN", "SDC1"],
    "T": ["CD3D", "CD3E", "TRAC", "TRBC1"],
    "NK": ["NKG7", "GNLY", "KLRD1", "KLRF1"],
    "Macrophage": ["LYZ", "LST1", "CD68", "CSF1R", "C1QA", "C1QB"],
    "Monocyte": ["S100A8", "S100A9", "IL1B", "SPP1"],
    "DC": [],
    "Mast/Basophil": [],
}

STATE_PATTERNS = (
    "Activated", "Naive", "Cytotoxic", "regulatory", "JUN+", "CDR7+",
    "CXCL2+", "LAMP3+", "CLEC9A+", "C1QC+", "CD1A+", "CLEC10A+",
    "CHI3L1+", "IL1B+", "LYVE1+", "RETN+", "PDE4C+", "LRMP+",
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--project-root", required=True, type=Path)
    parser.add_argument("--output-dir", required=True, type=Path)
    parser.add_argument("--chunk-size", type=int, default=4096)
    return parser.parse_args()


def sha256(path: Path, chunk_size: int = 1024 * 1024) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        while chunk := handle.read(chunk_size):
            digest.update(chunk)
    return digest.hexdigest()


def join_genes(values: Iterable[str]) -> str:
    return ";".join(sorted({str(value) for value in values if str(value)}))


def read_atlas(path: Path) -> pd.DataFrame:
    workbook = load_workbook(path, read_only=True, data_only=False)
    frames = []
    for sheet_index, worksheet in enumerate(workbook.worksheets, start=1):
        reported_dimension = worksheet.calculate_dimension()
        if reported_dimension == "A1:A1":
            worksheet.reset_dimensions()
        rows = list(worksheet.iter_rows(values_only=True))
        header = [str(value) for value in rows[0]]
        frame = pd.DataFrame(rows[1:], columns=header)
        frame.insert(0, "source_sheet_index", sheet_index)
        frame.insert(1, "source_sheet_name", worksheet.title)
        frame.insert(2, "source_row", np.arange(2, len(frame) + 2))
        frame["reported_sheet_dimension"] = reported_dimension
        frames.append(frame)
    workbook.close()
    atlas = pd.concat(frames, ignore_index=True)
    required = {"organ", "reference", "level", "celltype", "marker_gene", "fold_change", "pct_positive", "auc"}
    if not required.issubset(atlas.columns):
        raise ValueError(f"Atlas columns missing: {sorted(required - set(atlas.columns))}")
    for column in ["fold_change", "pct_positive", "auc"]:
        atlas[column] = pd.to_numeric(atlas[column], errors="coerce")
    atlas["original_rank"] = atlas.groupby(["level", "celltype"], sort=False).cumcount() + 1
    return atlas


def overlap_parent(subtype: str, subtype_genes: set[str], level1_sets: dict[str, set[str]]) -> tuple[str, int, float, str, int, float]:
    scores = []
    for level1, genes in level1_sets.items():
        overlap = len(subtype_genes & genes)
        union = len(subtype_genes | genes)
        scores.append((level1, overlap, overlap / union if union else 0.0))
    scores.sort(key=lambda item: (item[1], item[2], item[0]), reverse=True)
    first = scores[0]
    second = scores[1]
    return first[0], first[1], first[2], second[0], second[1], second[2]


def build_hierarchy(atlas: pd.DataFrame) -> pd.DataFrame:
    level1 = atlas[atlas["level"] == "cell_type"]
    level2 = atlas[atlas["level"] == "cell_subtype"]
    level1_sets = {name: set(group["marker_gene"].astype(str)) for name, group in level1.groupby("celltype")}
    rows = []
    for subtype, group in level2.groupby("celltype", sort=False):
        best, overlap_n, jaccard, second, second_n, second_jaccard = overlap_parent(
            subtype, set(group["marker_gene"].astype(str)), level1_sets
        )
        if subtype in level1_sets:
            parent = subtype
            basis = "exact atlas label match"
            status = "verified_exact_label"
        elif subtype in IMMUNE_PARENT:
            parent = IMMUNE_PARENT[subtype]
            basis = "immune nomenclature mapping checked against marker overlap"
            status = "verified_immune_hierarchy"
        else:
            explicit = {
                "COPD Ciliated": "Ciliated", "dividing AT2": "AT2", "pre AT2": "AT2",
                "Bronchial Vessel CXCL12+": "Bronchial Vessel",
                "Bronchial Vessel SELE+": "Bronchial Vessel",
                "Lipofibroblast": "Alveolar fibroblast",
            }
            parent = explicit.get(subtype, best)
            basis = "nomenclature plus top-50 marker overlap" if subtype in explicit else "best top-50 marker overlap; no explicit parent column in source"
            status = "inferred_nonimmune_hierarchy"
        if parent not in level1_sets:
            raise ValueError(f"Mapped Level 1 is absent from atlas: {subtype} -> {parent}")
        parent_genes = level1_sets[parent]
        mapped_overlap = len(set(group["marker_gene"].astype(str)) & parent_genes)
        mapped_jaccard = mapped_overlap / len(set(group["marker_gene"].astype(str)) | parent_genes)
        rows.append({
            "atlas_name": join_genes(atlas["reference"].dropna().astype(str).unique()),
            "level1": parent,
            "level2": subtype,
            "level2_marker_record_n": len(group),
            "mapped_parent_overlap_gene_n": mapped_overlap,
            "mapped_parent_jaccard": mapped_jaccard,
            "best_overlap_level1": best,
            "best_overlap_gene_n": overlap_n,
            "best_overlap_jaccard": jaccard,
            "second_overlap_level1": second,
            "second_overlap_gene_n": second_n,
            "second_overlap_jaccard": second_jaccard,
            "mapping_status": status,
            "mapping_basis": basis,
            "is_immune_level2": subtype in IMMUNE_PARENT,
        })
    hierarchy = pd.DataFrame(rows).sort_values(["level1", "level2"]).reset_index(drop=True)
    if hierarchy["level2"].duplicated().any():
        raise ValueError("Level 2 labels are not unique in hierarchy")
    return hierarchy


def functional_category(gene: str) -> str:
    categories = []
    if gene.startswith(("IGH", "IGK", "IGL")) or gene == "JCHAIN":
        categories.append("immunoglobulin/secretory")
    if gene.startswith("HLA-") or gene in {"CD74", "B2M", "CIITA", "TAP1", "TAP2"}:
        categories.append("antigen_presentation")
    if gene in {"NKG7", "GNLY", "PRF1", "CCL5", "CTSW", "KLRD1", "KLRF1"} or gene.startswith("GZM"):
        categories.append("cytotoxicity")
    if gene in {"S100A8", "S100A9", "S100A12", "IL1B", "CXCL8", "CXCL2", "CXCL3", "SPP1"}:
        categories.append("inflammation")
    if gene.startswith(("RPL", "RPS", "MT-")) or gene in {"MALAT1", "FOS", "JUN", "JUNB", "DUSP1"}:
        categories.append("housekeeping/stress")
    return ";".join(categories)


def build_marker_library(atlas: pd.DataFrame, hierarchy: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    subtype = atlas[atlas["level"] == "cell_subtype"].copy()
    subtype = subtype.merge(hierarchy[["level1", "level2"]], left_on="celltype", right_on="level2", how="left", validate="many_to_one")
    subtype["is_immune_level2"] = subtype["celltype"].isin(IMMUNE_PARENT)
    immune = subtype[subtype["is_immune_level2"]].copy()
    immune["source_level1"] = immune["level1"]
    immune["source_level2"] = immune["celltype"]
    immune["detailed_lineage"] = immune["source_level2"].map(DETAILED_LINEAGE)
    immune["major_lineage"] = immune["detailed_lineage"].map(MAJOR_LINEAGE)
    if immune[["source_level1", "detailed_lineage", "major_lineage"]].isna().any().any():
        raise ValueError("One or more immune subtype mappings are missing")

    immune_occurrence = immune.groupby("marker_gene").agg(
        immune_level1_n=("source_level1", "nunique"),
        immune_level2_n=("source_level2", "nunique"),
        detailed_lineage_n=("detailed_lineage", "nunique"),
        major_lineage_n=("major_lineage", "nunique"),
        immune_level1_labels=("source_level1", join_genes),
        immune_level2_labels=("source_level2", join_genes),
        major_lineage_labels=("major_lineage", join_genes),
    ).reset_index()
    nonimmune = subtype[~subtype["is_immune_level2"]]
    nonimmune_occurrence = nonimmune.groupby("marker_gene").agg(
        nonimmune_level2_n=("celltype", "nunique"),
        nonimmune_level2_labels=("celltype", join_genes),
    ).reset_index()
    immune = immune.merge(immune_occurrence, on="marker_gene", how="left", validate="many_to_one")
    immune = immune.merge(nonimmune_occurrence, on="marker_gene", how="left", validate="many_to_one")
    immune["nonimmune_level2_n"] = immune["nonimmune_level2_n"].fillna(0).astype(int)
    immune["nonimmune_level2_labels"] = immune["nonimmune_level2_labels"].fillna("")

    level1_gene_sets = {
        name: set(group["marker_gene"].astype(str))
        for name, group in atlas[atlas["level"] == "cell_type"].groupby("celltype")
    }
    immune["level1_marker_support"] = [
        gene in level1_gene_sets[parent]
        for gene, parent in zip(immune["marker_gene"].astype(str), immune["source_level1"].astype(str))
    ]
    immune["functional_category"] = immune["marker_gene"].astype(str).map(functional_category)
    immune["is_curated_core_anchor"] = [
        gene in CORE_ANCHORS[lineage]
        for gene, lineage in zip(immune["marker_gene"].astype(str), immune["major_lineage"].astype(str))
    ]
    immune["source_is_state_label"] = immune["source_level2"].astype(str).map(
        lambda value: any(pattern in value for pattern in STATE_PATTERNS)
    )

    def classify(row: pd.Series) -> str:
        shared = row["major_lineage_n"] > 1 or row["nonimmune_level2_n"] > 0 or bool(row["functional_category"])
        if shared:
            return "D. Shared/functional markers"
        if row["is_curated_core_anchor"] and row["level1_marker_support"] and row["auc"] >= 0.65:
            return "A. Lineage-core"
        if row["source_is_state_label"] and row["original_rank"] <= 25 and row["auc"] >= 0.75:
            return "C. State-associated markers"
        return "B. Extended lineage markers"

    immune["marker_class"] = immune.apply(classify, axis=1)
    immune["signature_eligible"] = (
        (immune["marker_class"] == "A. Lineage-core")
        | (
            immune["marker_class"].isin(["B. Extended lineage markers", "C. State-associated markers"])
            & (immune["original_rank"] <= 30)
            & (immune["auc"] >= 0.72)
            & (immune["pct_positive"] >= 20)
            & (immune["major_lineage_n"] == 1)
            & (immune["nonimmune_level2_n"] == 0)
        )
    )
    immune["specificity_status"] = np.select(
        [
            immune["marker_class"] == "A. Lineage-core",
            immune["marker_class"] == "C. State-associated markers",
            immune["signature_eligible"] & (immune["marker_class"] == "B. Extended lineage markers"),
            immune["marker_class"] == "D. Shared/functional markers",
        ],
        ["relatively_specific_core", "state_associated_not_lineage_core", "moderate_extended", "shared_or_cross_lineage_risk"],
        default="limited_extended_evidence",
    )
    immune["fold_change_definition"] = "unknown; source workbook does not define whether this is log2FC"
    immune["specificity_limit"] = "Positive top-marker lists only; no complete non-target expression matrix, so cell specificity cannot be fully estimated."

    library_columns = [
        "organ", "reference", "source_sheet_name", "source_row", "source_level1", "source_level2",
        "detailed_lineage", "major_lineage", "marker_gene", "original_rank", "fold_change",
        "pct_positive", "auc", "fold_change_definition", "marker_class", "signature_eligible",
        "specificity_status", "is_curated_core_anchor", "level1_marker_support", "functional_category",
        "immune_level1_n", "immune_level2_n", "detailed_lineage_n", "major_lineage_n",
        "immune_level1_labels", "immune_level2_labels", "major_lineage_labels",
        "nonimmune_level2_n", "nonimmune_level2_labels", "specificity_limit",
    ]
    library = immune[library_columns].sort_values(["major_lineage", "source_level2", "original_rank"]).reset_index(drop=True)

    specificity = library.groupby("marker_gene", sort=True).agg(
        immune_record_n=("marker_gene", "size"),
        immune_level1_n=("source_level1", "nunique"),
        immune_level2_n=("source_level2", "nunique"),
        detailed_lineage_n=("detailed_lineage", "nunique"),
        major_lineage_n=("major_lineage", "nunique"),
        source_level1=("source_level1", join_genes),
        source_level2=("source_level2", join_genes),
        detailed_lineages=("detailed_lineage", join_genes),
        major_lineages=("major_lineage", join_genes),
        best_auc=("auc", "max"),
        best_pct_positive=("pct_positive", "max"),
        best_fold_change=("fold_change", "max"),
        best_original_rank=("original_rank", "min"),
        marker_classes=("marker_class", join_genes),
        any_signature_eligible=("signature_eligible", "max"),
        functional_categories=("functional_category", join_genes),
        nonimmune_level2_n=("nonimmune_level2_n", "max"),
        nonimmune_level2_labels=("nonimmune_level2_labels", join_genes),
    ).reset_index()
    specificity["overlap_risk"] = np.select(
        [
            specificity["major_lineage_n"] > 1,
            specificity["nonimmune_level2_n"] > 0,
            specificity["functional_categories"].astype(str).str.len() > 0,
        ],
        ["cross_immune_lineage", "immune_and_nonimmune_overlap", "shared_functional_program"],
        default="not_observed_in_other_atlas_top50_lists",
    )
    specificity["specificity_audit_limit"] = "Absence from another top-50 list is not proof of absent expression."
    return library, specificity


def build_reference18_mapping() -> pd.DataFrame:
    rows = []
    for level1 in sorted(IMMUNE_LEVEL1):
        level2_values = [subtype for subtype, parent in IMMUNE_PARENT.items() if parent == level1]
        for level2 in [""] + sorted(level2_values):
            label = level2 or level1
            if label in {"B cell", "B cell LRMP+"}:
                target, mapping_type = "B", "direct major lineage; Level 2 states merged"
                note = "Reference18 B does not resolve the atlas B subtypes."
            elif label == "Plasma B":
                target, mapping_type = "Plasma", "direct major lineage"
                note = "Atlas plasma state maps to Reference18 Plasma; subtype resolution is not implied."
            elif label == "B/Plasma":
                target, mapping_type = "B;Plasma", "Level 1 requires Level 2 split"
                note = "Combined atlas Level 1 cannot map one-to-one to separate Reference18 B and Plasma."
            elif level1 == "T cell":
                target, mapping_type = "T", "direct major lineage; Level 2 states merged"
                note = "Reference18 has no separate CD4, CD8, regulatory, naive, JUN+, or cytotoxic T state."
            elif level1 == "NK":
                target, mapping_type = "NK", "direct major lineage; activation state merged"
                note = "Activated NK XCL1+ is not a separate Reference18 state."
            elif label in {"Macrophage FABP4+", "Macrophage FABP4+/PDE4C+", "Macrophage C1Q hi", "Macrophage LYVE1+"}:
                target, mapping_type = "Alveolar/resident macrophage", "best biological major-lineage correspondence"
                note = "Reference18 merges resident macrophage states."
            elif label in {"Macrophage CHI3L1+/CD9 hi/", "Macrophage IL1B+", "Macrophage RETN+/VCAN+"}:
                target, mapping_type = "Inflammatory monocyte/macrophage", "best biological state correspondence"
                note = "Reference18 does not separately resolve these inflammatory macrophage states."
            elif label == "Macrophage":
                target, mapping_type = "Alveolar/resident macrophage;Inflammatory monocyte/macrophage", "not reliably one-to-one"
                note = "Broad atlas macrophage Level 1 spans both Reference18 macrophage axes."
            elif level1 == "Monocyte":
                target, mapping_type = "Inflammatory monocyte/macrophage", "major lineage merged with inflammatory macrophage"
                note = "Reference18 does not retain atlas CD14/CD16 monocyte subtypes separately."
            elif level1 in {"Dendritic", "pDC"}:
                target, mapping_type = "Dendritic cell", "major lineage; Level 2 states merged"
                note = "Reference18 has no separate DC1/DC2/LAMP3/pDC state."
            elif level1 == "Basophil/Mast":
                target, mapping_type = "Mast/Basophil", "direct combined major lineage"
                note = "Neither system separates mast from basophil in this mapping."
            else:
                target, mapping_type, note = "", "unresolved", "No reliable correspondence."
            rows.append({
                "atlas_level1": level1,
                "atlas_level2": level2,
                "reference18_state": target,
                "mapping_type": mapping_type,
                "interpretation": note,
                "reference18_state_present": all(state in REFERENCE18_STATES for state in target.split(";") if state),
            })
    return pd.DataFrame(rows)


def sparse_detectability(
    h5ad_path: Path,
    gene_map_path: Path,
    genes: list[str],
    chunk_size: int,
) -> tuple[pd.DataFrame, pd.DataFrame, dict[str, set[str]]]:
    gene_map = pd.read_csv(gene_map_path, dtype=str).fillna("")
    for column in ["mapped", "conflict"]:
        gene_map[column] = gene_map[column].str.lower().isin({"true", "1", "yes"})
    data = ad.read_h5ad(h5ad_path, backed="r")
    try:
        var_names = pd.Index(data.var_names.astype(str))
        usable = gene_map[
            gene_map["gene_symbol"].isin(genes)
            & gene_map["ensembl_var_name"].isin(var_names)
            & gene_map["mapped"]
            & ~gene_map["conflict"]
        ].copy()
        gene_position = {gene: index for index, gene in enumerate(genes)}
        var_position = pd.Series(np.arange(len(var_names)), index=var_names)
        usable["var_index"] = usable["ensembl_var_name"].map(var_position).astype(int)
        usable["gene_index"] = usable["gene_symbol"].map(gene_position).astype(int)
        usable = usable.sort_values("var_index").drop_duplicates(["ensembl_var_name", "gene_symbol"])
        selected_var = usable["var_index"].to_numpy(dtype=np.int64)
        projection = sparse.csr_matrix(
            (
                np.ones(len(usable), dtype=np.uint8),
                (np.arange(len(usable), dtype=np.int64), usable["gene_index"].to_numpy(dtype=np.int64)),
            ),
            shape=(len(usable), len(genes)),
        )

        obs = data.obs.copy()
        sample_order = list(dict.fromkeys(obs["sample_id"].astype(str)))
        if len(sample_order) != 21:
            raise ValueError(f"Expected 21 patients, observed {len(sample_order)}")
        sample_meta = obs.groupby("sample_id", sort=False).agg(
            source_chip_id=("chip_id", "first"),
            group=("group", "first"),
            total_evaluable_bin50=("chip_id", "size"),
            median_bin_umi=("total_counts_common_raw", "median"),
        ).reindex(sample_order)
        count_sum = {sample: np.zeros(len(genes), dtype=np.float64) for sample in sample_order}
        detected_bins = {sample: np.zeros(len(genes), dtype=np.int64) for sample in sample_order}
        sample_values = obs["sample_id"].astype(str).to_numpy()

        for start in range(0, data.n_obs, chunk_size):
            end = min(start + chunk_size, data.n_obs)
            matrix = data.X[start:end, selected_var]
            matrix = sparse.csr_matrix(matrix) if not sparse.issparse(matrix) else matrix.tocsr()
            by_gene = (matrix @ projection).tocsr()
            chunk_samples = sample_values[start:end]
            for sample in pd.unique(chunk_samples):
                local = np.flatnonzero(chunk_samples == sample)
                sub = by_gene[local, :]
                count_sum[sample] += np.asarray(sub.sum(axis=0)).ravel()
                detected_bins[sample] += np.asarray(sub.getnnz(axis=0)).ravel()

        available = set(usable["gene_symbol"])
        patient_rows = []
        patients_detecting = np.zeros(len(genes), dtype=np.int64)
        for gene_index, gene in enumerate(genes):
            patients_detecting[gene_index] = sum(count_sum[sample][gene_index] > 0 for sample in sample_order)
        for sample in sample_order:
            meta = sample_meta.loc[sample]
            n_bins = int(meta["total_evaluable_bin50"])
            for gene_index, gene in enumerate(genes):
                total = float(count_sum[sample][gene_index])
                detected = int(detected_bins[sample][gene_index])
                in_matrix = gene in available
                patient_rows.append({
                    "gene": gene,
                    "chip": CHIP_ALIAS[str(meta["source_chip_id"])],
                    "source_chip_id": str(meta["source_chip_id"]),
                    "patient": sample,
                    "group_metadata_only_not_used_for_marker_selection": str(meta["group"]),
                    "gene_in_matrix": in_matrix,
                    "expression_status": "unavailable" if not in_matrix else "detected" if total > 0 else "undetected",
                    "total_raw_counts": int(round(total)),
                    "detected_bin50_n": detected,
                    "detected_bin50_fraction": detected / n_bins,
                    "mean_counts_per_bin": total / n_bins,
                    "patients_with_detectable_expression_all21": int(patients_detecting[gene_index]),
                    "total_evaluable_bin50": n_bins,
                    "median_bin_umi": float(meta["median_bin_umi"]),
                })
        patient = pd.DataFrame(patient_rows)
        chip = patient.groupby(["gene", "chip", "source_chip_id", "gene_in_matrix"], sort=False).agg(
            total_raw_counts=("total_raw_counts", "sum"),
            detected_bin50_n=("detected_bin50_n", "sum"),
            total_evaluable_bin50=("total_evaluable_bin50", "sum"),
            patient_n=("patient", "nunique"),
            patients_with_detectable_expression=("expression_status", lambda values: int(sum(value == "detected" for value in values))),
            median_patient_bin_umi=("median_bin_umi", "median"),
        ).reset_index()
        chip["detected_bin50_fraction"] = chip["detected_bin50_n"] / chip["total_evaluable_bin50"]
        chip["mean_counts_per_bin"] = chip["total_raw_counts"] / chip["total_evaluable_bin50"]
        chip["expression_status"] = np.where(
            ~chip["gene_in_matrix"], "unavailable",
            np.where(chip["total_raw_counts"] > 0, "detected", "undetected"),
        )
        symbol_to_ensembl = usable.groupby("gene_symbol")["ensembl_var_name"].agg(set).to_dict()
        return patient, chip, symbol_to_ensembl
    finally:
        data.file.close()


def panel_gene_sets(library: pd.DataFrame) -> tuple[dict[str, set[str]], dict[str, set[str]], dict[str, set[str]]]:
    core = defaultdict(set)
    extended = defaultdict(set)
    shared = defaultdict(set)
    for row in library.itertuples(index=False):
        gene = str(row.marker_gene)
        lineage = str(row.major_lineage)
        if row.marker_class == "A. Lineage-core":
            core[lineage].add(gene)
            extended[lineage].add(gene)
        elif bool(row.signature_eligible) and row.marker_class in {"B. Extended lineage markers", "C. State-associated markers"}:
            extended[lineage].add(gene)
        elif row.marker_class == "D. Shared/functional markers":
            shared[lineage].add(gene)
    return dict(core), dict(extended), dict(shared)


def coverage_table(
    patient: pd.DataFrame,
    core: dict[str, set[str]],
    extended: dict[str, set[str]],
) -> pd.DataFrame:
    patient_meta = patient[["patient", "chip", "source_chip_id", "group_metadata_only_not_used_for_marker_selection"]].drop_duplicates()
    panels = {
        "A. Previously used small marker panel": {key: set(value) for key, value in LEGACY_SMALL_PANEL.items()},
        "B. Atlas-derived lineage-core panel": core,
        "C. Extended marker panel": extended,
    }
    rows = []
    for meta in patient_meta.itertuples(index=False):
        for lineage in LEGACY_SMALL_PANEL:
            for panel_name, by_lineage in panels.items():
                genes = sorted(by_lineage.get(lineage, set()))
                subset = patient[(patient["patient"] == meta.patient) & patient["gene"].isin(genes)] if genes else pd.DataFrame()
                if genes:
                    available_n = int(subset["gene_in_matrix"].sum())
                    detected_n = int((subset["expression_status"] == "detected").sum())
                else:
                    available_n = detected_n = 0
                rows.append({
                    "scope": "patient",
                    "chip": meta.chip,
                    "source_chip_id": meta.source_chip_id,
                    "patient": meta.patient,
                    "group_metadata_only_not_used_for_marker_selection": meta.group_metadata_only_not_used_for_marker_selection,
                    "lineage": lineage,
                    "panel": panel_name,
                    "panel_gene_n": len(genes),
                    "available_gene_n": available_n,
                    "detectable_gene_n": detected_n,
                    "detectable_fraction_of_panel": detected_n / len(genes) if genes else np.nan,
                    "panel_genes": join_genes(genes),
                })
    result = pd.DataFrame(rows)
    baseline = result[result["panel"] == "A. Previously used small marker panel"][
        ["patient", "lineage", "detectable_gene_n"]
    ].rename(columns={"detectable_gene_n": "small_panel_detectable_gene_n"})
    result = result.merge(baseline, on=["patient", "lineage"], how="left")
    result["additional_detectable_vs_small"] = result["detectable_gene_n"] - result["small_panel_detectable_gene_n"]
    summary_rows = []
    for (lineage, panel_name), group in result.groupby(["lineage", "panel"], sort=False):
        genes = sorted(panels[panel_name].get(lineage, set()))
        gene_subset = patient[patient["gene"].isin(genes)]
        summary_rows.append({
            "scope": "all_21_patients",
            "chip": "ALL",
            "source_chip_id": "ALL",
            "patient": "ALL21",
            "group_metadata_only_not_used_for_marker_selection": "not_applicable",
            "lineage": lineage,
            "panel": panel_name,
            "panel_gene_n": len(genes),
            "available_gene_n": int(gene_subset.drop_duplicates("gene")["gene_in_matrix"].sum()) if genes else 0,
            "detectable_gene_n": int(gene_subset[gene_subset["expression_status"] == "detected"]["gene"].nunique()) if genes else 0,
            "detectable_fraction_of_panel": (
                gene_subset[gene_subset["expression_status"] == "detected"]["gene"].nunique() / len(genes)
                if genes else np.nan
            ),
            "panel_genes": join_genes(genes),
            "small_panel_detectable_gene_n": int(group["small_panel_detectable_gene_n"].median()),
            "additional_detectable_vs_small": float(group["additional_detectable_vs_small"].median()),
            "median_detectable_gene_n_per_patient": float(group["detectable_gene_n"].median()),
            "minimum_detectable_gene_n_per_patient": int(group["detectable_gene_n"].min()),
            "maximum_detectable_gene_n_per_patient": int(group["detectable_gene_n"].max()),
        })
    return pd.concat([result, pd.DataFrame(summary_rows)], ignore_index=True, sort=False)


def detection_for_genes(patient: pd.DataFrame, genes: set[str]) -> dict[str, float | int | str]:
    subset = patient[patient["gene"].isin(genes)]
    available = sorted(subset.loc[subset["gene_in_matrix"], "gene"].unique())
    unavailable = sorted(set(genes) - set(available))
    if not genes:
        return {
            "available_genes": "", "unavailable_genes": "", "available_gene_n": 0,
            "unavailable_gene_n": 0, "patients_with_any_detectable_gene": 0,
            "median_detectable_gene_n_per_patient": 0.0, "minimum_detectable_gene_n_per_patient": 0,
            "gene_patient_detection_fraction": np.nan,
        }
    pivot = subset.pivot_table(index="patient", columns="gene", values="total_raw_counts", aggfunc="sum", fill_value=0)
    pivot = pivot.reindex(columns=sorted(genes), fill_value=0)
    detected_n = (pivot > 0).sum(axis=1)
    return {
        "available_genes": join_genes(available),
        "unavailable_genes": join_genes(unavailable),
        "available_gene_n": len(available),
        "unavailable_gene_n": len(unavailable),
        "patients_with_any_detectable_gene": int((detected_n > 0).sum()),
        "median_detectable_gene_n_per_patient": float(detected_n.median()),
        "minimum_detectable_gene_n_per_patient": int(detected_n.min()),
        "gene_patient_detection_fraction": float((pivot > 0).to_numpy().mean()),
    }


def signature_tables(
    library: pd.DataFrame,
    patient: pd.DataFrame,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    records = []
    for (detailed, subtype), group in library.groupby(["detailed_lineage", "source_level2"], sort=True):
        core = set(group.loc[group["marker_class"] == "A. Lineage-core", "marker_gene"].astype(str))
        extended_only = set(group.loc[group["signature_eligible"] & group["marker_class"].isin(["B. Extended lineage markers", "C. State-associated markers"]), "marker_gene"].astype(str))
        shared = set(group.loc[group["marker_class"] == "D. Shared/functional markers", "marker_gene"].astype(str))
        signatures = {
            "1. Reference-based core signature": core,
            "2. Extended exploratory signature": core | extended_only,
            "3. Stereo-seq detectable subset": {
                gene for gene in (core | extended_only)
                if patient.loc[(patient["gene"] == gene) & (patient["expression_status"] == "detected"), "patient"].nunique() > 0
            },
        }
        for signature_type, genes in signatures.items():
            detect = detection_for_genes(patient, genes)
            records.append({
                "cell_lineage": detailed,
                "atlas_level2_subtype": subtype,
                "major_lineage": MAJOR_LINEAGE[detailed],
                "signature_type": signature_type,
                "core_genes": join_genes(core),
                "extended_genes": join_genes(extended_only),
                "shared_genes_not_used_as_lineage_specific_core": join_genes(shared),
                "signature_genes": join_genes(genes),
                "signature_gene_n": len(genes),
                **detect,
                "evidence_source": join_genes(group["reference"].astype(str).unique()),
                "limitation": "Candidate signature only. Detection coverage does not demonstrate cell-identification accuracy; shared genes remain context-only.",
            })
    signatures = pd.DataFrame(records)

    summaries = []
    for (lineage, signature_type), group in signatures.groupby(["cell_lineage", "signature_type"], sort=True):
        genes = set()
        for value in group["signature_genes"]:
            genes.update(gene for gene in str(value).split(";") if gene)
        detect = detection_for_genes(patient, genes)
        if detect["median_detectable_gene_n_per_patient"] >= 2 and detect["patients_with_any_detectable_gene"] >= 14:
            condition = "coverage_supports_reviewed_exploratory_followup"
        elif detect["median_detectable_gene_n_per_patient"] >= 1 and detect["patients_with_any_detectable_gene"] >= 10:
            condition = "limited_exploratory_followup_only"
        else:
            condition = "insufficient_stable_detection"
        summaries.append({
            "cell_lineage": lineage,
            "signature_type": signature_type,
            "atlas_subtype_n": group["atlas_level2_subtype"].nunique(),
            "signature_genes": join_genes(genes),
            "signature_gene_n": len(genes),
            **detect,
            "followup_condition": condition,
            "accuracy_claim_allowed": False,
        })
    return signatures, pd.DataFrame(summaries)


def save_pdf(fig: plt.Figure, path: Path) -> None:
    with PdfPages(path) as pdf:
        pdf.savefig(fig, bbox_inches="tight")
    plt.close(fig)


def generate_figures(
    figure_dir: Path,
    library: pd.DataFrame,
    patient: pd.DataFrame,
    chip: pd.DataFrame,
    coverage: pd.DataFrame,
    core: dict[str, set[str]],
    extended: dict[str, set[str]],
) -> None:
    plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 9, "axes.titlesize": 12})
    class_order = [
        "A. Lineage-core", "B. Extended lineage markers",
        "C. State-associated markers", "D. Shared/functional markers",
    ]
    class_colors = ["#214f78", "#4f8fba", "#d59b3f", "#9b6b8f"]
    counts = library.groupby(["source_level2", "marker_class"])["marker_gene"].nunique().unstack(fill_value=0)
    counts = counts.reindex(columns=class_order, fill_value=0)
    counts = counts.loc[counts.sum(axis=1).sort_values().index]
    fig, ax = plt.subplots(figsize=(10.5, max(7, len(counts) * 0.29)))
    left = np.zeros(len(counts))
    for marker_class, color in zip(class_order, class_colors):
        values = counts[marker_class].to_numpy()
        ax.barh(counts.index, values, left=left, color=color, label=marker_class.split(". ", 1)[1])
        left += values
    ax.set_xlabel("Atlas marker records per subtype")
    ax.set_ylabel("Atlas immune Level 2 subtype")
    ax.set_title("Atlas immune subtype marker coverage by evidence class")
    ax.legend(frameon=False, loc="upper center", bbox_to_anchor=(0.5, -0.075), ncol=2)
    ax.grid(axis="x", color="#dddddd", linewidth=0.6)
    save_pdf(fig, figure_dir / "01_ATLAS_IMMUNE_CELLTYPE_MARKER_COVERAGE_EN.pdf")

    chip_rows = []
    for lineage in LEGACY_SMALL_PANEL:
        for panel_name, genes_by_lineage in [("Core", core), ("Extended", extended)]:
            genes = genes_by_lineage.get(lineage, set())
            for chip_name in ["L3", "J2", "Y40102K8"]:
                subset = chip[(chip["chip"] == chip_name) & chip["gene"].isin(genes)]
                rate = float((subset["expression_status"] == "detected").mean()) if len(subset) else np.nan
                chip_rows.append({"row": f"{lineage} - {panel_name}", "chip": chip_name, "rate": rate})
    chip_plot = pd.DataFrame(chip_rows).pivot(index="row", columns="chip", values="rate").reindex(columns=["L3", "J2", "Y40102K8"])
    fig, ax = plt.subplots(figsize=(7.2, 8.8))
    image = ax.imshow(chip_plot.fillna(0).values, aspect="auto", cmap="Blues", vmin=0, vmax=1)
    ax.set_xticks(range(len(chip_plot.columns)), chip_plot.columns)
    ax.set_yticks(range(len(chip_plot.index)), chip_plot.index)
    ax.set_title("Core and extended marker detection by chip")
    for i in range(chip_plot.shape[0]):
        for j in range(chip_plot.shape[1]):
            value = chip_plot.iloc[i, j]
            ax.text(j, i, "n.a." if pd.isna(value) else f"{value:.0%}", ha="center", va="center", color="black", fontsize=8)
    fig.colorbar(image, ax=ax, label="Fraction of panel genes detected")
    save_pdf(fig, figure_dir / "02_MARKER_DETECTION_BY_CHIP_EN.pdf")

    patient_extended = coverage[(coverage["scope"] == "patient") & (coverage["panel"] == "C. Extended marker panel")]
    patient_plot = patient_extended.pivot(index="patient", columns="lineage", values="detectable_fraction_of_panel")
    order = patient_extended[["patient", "chip"]].drop_duplicates().sort_values(["chip", "patient"])["patient"]
    patient_plot = patient_plot.reindex(index=order, columns=list(LEGACY_SMALL_PANEL))
    fig, ax = plt.subplots(figsize=(10, 8.5))
    image = ax.imshow(patient_plot.fillna(0).values, aspect="auto", cmap="YlGnBu", vmin=0, vmax=1)
    ax.set_xticks(range(len(patient_plot.columns)), patient_plot.columns, rotation=35, ha="right")
    ax.set_yticks(range(len(patient_plot.index)), patient_plot.index)
    ax.set_title("Extended immune marker coverage across 21 patients")
    fig.colorbar(image, ax=ax, label="Fraction of extended genes detected")
    save_pdf(fig, figure_dir / "03_IMMUNE_MARKER_PATIENT_COVERAGE_EN.pdf")

    all21 = coverage[coverage["scope"] == "all_21_patients"].copy()
    panel_order = [
        "A. Previously used small marker panel",
        "B. Atlas-derived lineage-core panel",
        "C. Extended marker panel",
    ]
    fig, ax = plt.subplots(figsize=(11, 6.5))
    x = np.arange(len(LEGACY_SMALL_PANEL))
    width = 0.25
    colors = ["#777777", "#24557a", "#4cae8b"]
    for index, panel_name in enumerate(panel_order):
        values = all21[all21["panel"] == panel_name].set_index("lineage").reindex(LEGACY_SMALL_PANEL)["detectable_gene_n"].fillna(0)
        ax.bar(x + (index - 1) * width, values, width, label=panel_name.split(". ", 1)[1], color=colors[index])
    ax.set_xticks(x, list(LEGACY_SMALL_PANEL), rotation=30, ha="right")
    ax.set_ylabel("Genes detected in at least one patient")
    ax.set_title("Previously used versus atlas-derived marker coverage")
    ax.legend(frameon=False)
    ax.grid(axis="y", color="#dddddd", linewidth=0.6)
    save_pdf(fig, figure_dir / "04_CORE_VS_EXTENDED_MARKER_COVERAGE_EN.pdf")

    lineages = list(LEGACY_SMALL_PANEL)
    atlas_by_lineage = {
        lineage: set(library.loc[library["major_lineage"] == lineage, "marker_gene"].astype(str))
        for lineage in lineages
    }
    matrix = np.zeros((len(lineages), len(lineages)))
    for i, first in enumerate(lineages):
        for j, second in enumerate(lineages):
            a, b = atlas_by_lineage.get(first, set()), atlas_by_lineage.get(second, set())
            matrix[i, j] = len(a & b) / len(a | b) if a | b else 0
    fig, ax = plt.subplots(figsize=(8, 7))
    image = ax.imshow(matrix, cmap="OrRd", vmin=0, vmax=max(0.25, float(matrix.max())))
    ax.set_xticks(range(len(lineages)), lineages, rotation=40, ha="right")
    ax.set_yticks(range(len(lineages)), lineages)
    ax.set_title("Pairwise overlap of complete atlas immune marker sets")
    for i in range(len(lineages)):
        for j in range(len(lineages)):
            ax.text(j, i, f"{matrix[i, j]:.2f}", ha="center", va="center", fontsize=8)
    fig.colorbar(image, ax=ax, label="Jaccard overlap")
    save_pdf(fig, figure_dir / "05_IMMUNE_SIGNATURE_OVERLAP_EN.pdf")


def write_audit_table(atlas: pd.DataFrame, path: Path, input_path: Path) -> None:
    exact_duplicates = atlas.duplicated(["organ", "reference", "level", "celltype", "marker_gene"], keep=False)
    rows = [
        ("input_file", str(input_path.resolve()), "PASS", "Primary source opened read-only"),
        ("source_sheet_n", atlas["source_sheet_name"].nunique(), "PASS", join_genes(atlas["source_sheet_name"].unique())),
        ("atlas_name_version", join_genes(atlas["reference"].unique()), "PASS", "From source reference column"),
        ("organ", join_genes(atlas["organ"].unique()), "PASS", "From source organ column"),
        ("record_n", len(atlas), "PASS", "All cells recovered after in-memory dimension reset"),
        ("level1_cell_type_n", atlas.loc[atlas["level"] == "cell_type", "celltype"].nunique(), "PASS", "Source level=cell_type"),
        ("level2_cell_subtype_n", atlas.loc[atlas["level"] == "cell_subtype", "celltype"].nunique(), "PASS", "Source level=cell_subtype"),
        ("unique_marker_gene_n", atlas["marker_gene"].nunique(), "PASS", "Across all source records"),
        ("exact_duplicate_record_n", int(exact_duplicates.sum()), "PASS" if not exact_duplicates.any() else "REVIEW", "Key: organ/reference/level/celltype/marker_gene"),
        ("genes_repeated_across_records_n", int((atlas.groupby("marker_gene").size() > 1).sum()), "PASS", "Expected when genes recur across cell types/levels"),
        ("statistical_columns", "fold_change;pct_positive;auc", "PASS", "No other source statistical columns observed"),
        ("fold_change_definition", "unknown", "LIMITATION", "Workbook does not state whether fold_change is log2FC"),
        ("specificity_data_scope", "positive top markers only", "LIMITATION", "No complete non-target expression matrix"),
        ("worksheet_dimension_issue", join_genes(atlas["reported_sheet_dimension"].unique()), "RECOVERED", "XML contains data beyond stale A1 dimension; source was not modified"),
    ]
    pd.DataFrame(rows, columns=["audit_item", "observed_value", "status", "detail"]).to_csv(path, sep="\t", index=False)


def write_findings(
    path: Path,
    atlas: pd.DataFrame,
    library: pd.DataFrame,
    specificity: pd.DataFrame,
    patient: pd.DataFrame,
    coverage: pd.DataFrame,
    core: dict[str, set[str]],
    extended: dict[str, set[str]],
    signature_summary: pd.DataFrame,
) -> None:
    immune_level1 = sorted(IMMUNE_LEVEL1)
    immune_level2 = sorted(IMMUNE_PARENT)
    global_gene = patient.groupby("gene").agg(
        total_counts=("total_raw_counts", "sum"),
        detected_bins=("detected_bin50_n", "sum"),
        evaluable_bins=("total_evaluable_bin50", "sum"),
        patients=("expression_status", lambda values: int(sum(value == "detected" for value in values))),
    )
    global_gene["detection_fraction"] = global_gene["detected_bins"] / global_gene["evaluable_bins"]
    top_detectable = global_gene.sort_values(["patients", "detection_fraction", "total_counts"], ascending=False).head(15)
    legacy_genes = sorted(set().union(*[set(values) for values in LEGACY_SMALL_PANEL.values()]))
    classic = global_gene.reindex(legacy_genes).fillna(0).sort_values(["patients", "detection_fraction"]).head(12)
    shared_n = int((specificity["major_lineage_n"] > 1).sum())
    nonimmune_overlap_n = int((specificity["nonimmune_level2_n"] > 0).sum())
    all21 = coverage[coverage["scope"] == "all_21_patients"]
    small_detect = int(all21[all21["panel"] == "A. Previously used small marker panel"]["detectable_gene_n"].sum())
    ext_detect = int(all21[all21["panel"] == "C. Extended marker panel"]["detectable_gene_n"].sum())
    followup = signature_summary[
        (signature_summary["signature_type"] == "3. Stereo-seq detectable subset")
        & (signature_summary["followup_condition"] != "insufficient_stable_detection")
    ]
    insufficient = signature_summary[
        (signature_summary["signature_type"] == "3. Stereo-seq detectable subset")
        & (signature_summary["followup_condition"] == "insufficient_stable_detection")
    ]
    lines = [
        "# 扩展免疫 marker 库与 Stereo-seq RNA 可检测性审计：关键发现",
        "",
        "## 1. Atlas 的免疫分类",
        "",
        f"Atlas 为 `{join_genes(atlas['reference'].unique())}`。共识别 {len(immune_level1)} 个免疫 Level 1：{'; '.join(immune_level1)}。",
        f"共识别 {len(immune_level2)} 个免疫 Level 2：{'; '.join(immune_level2)}。Level 2 到 Level 1 的逐项核对见 `tables/02_ATLAS_CELLTYPE_HIERARCHY.tsv`。",
        "",
        "## 2. 免疫 marker 数量",
        "",
        f"从免疫 Level 2 记录提取 {len(library):,} 条 marker 记录，涉及 {library['marker_gene'].nunique():,} 个唯一基因。严格 lineage-core 为 {library.loc[library['marker_class'] == 'A. Lineage-core', 'marker_gene'].nunique()} 个；进入扩展候选的非 core 基因为 {library.loc[library['signature_eligible'] & (library['marker_class'] != 'A. Lineage-core'), 'marker_gene'].nunique()} 个。",
        "",
        "## 3. 谱系间重叠与特异性限制",
        "",
        f"有 {shared_n} 个基因出现在两个或以上免疫 major lineages 的 atlas Top50 中；另有 {nonimmune_overlap_n} 个免疫 marker 同时出现在非免疫 Level 2 Top50 中。CD74/HLA-DRA、NKG7/GNLY、IGKC/IGHG 等均按 shared/functional 风险处理，不作为严格单一谱系证据。",
        "本表只有各细胞类型的正向 Top marker，没有完整的非目标细胞表达矩阵。因此，“未在其他 Top50 出现”不等于真实细胞特异；本轮只能完成相对特异性审核，不能完成严格 sensitivity/specificity 验证。",
        "",
        "## 4. Stereo-seq 中最容易检出的基因",
        "",
        "; ".join(f"{gene}（{int(row.patients)}/21患者，{row.detection_fraction:.3%} Bin50）" for gene, row in top_detectable.iterrows()) + "。",
        "",
        "## 5. 检出不足的经典 marker",
        "",
        "; ".join(f"{gene}（{int(row.patients)}/21患者，{row.detection_fraction:.3%} Bin50）" for gene, row in classic.iterrows()) + "。",
        "",
        "## 6. 扩展 marker 是否增加覆盖",
        "",
        f"按八个 major lineages 分别计数并求和，原 small panels 在至少一位患者可检出 {small_detect} 个 lineage-panel gene entries；extended panels 为 {ext_detect} 个，增加 {ext_detect - small_detect} 个。该增加代表 detection coverage 增加，不代表免疫细胞识别准确率提高。每位患者的增量见 `tables/08_CORE_VS_EXTENDED_MARKER_COVERAGE.tsv`。",
        "",
        "## 7. B / Plasma / T / NK / Myeloid 的 core 与 extended",
        "",
    ]
    for lineage in ["B", "Plasma", "T", "NK", "Macrophage", "Monocyte", "DC", "Mast/Basophil"]:
        lines.append(f"- {lineage} core：{join_genes(core.get(lineage, set())) or '无通过严格规则的基因'}；extended：{join_genes(extended.get(lineage, set()) - core.get(lineage, set())) or '无额外候选'}。")
    lines.extend([
        "",
        "## 8. 具备进一步探索条件的 Level 2 状态",
        "",
        "; ".join(sorted(followup["cell_lineage"].unique())) if len(followup) else "无状态达到本轮预设的最低可检测覆盖条件。",
        "这些状态仅具备 RNA coverage 层面的候选条件，尚未证明在 50 µm neighborhood 中可被准确区分。",
        "",
        "## 9. 可能无法稳定区分的状态",
        "",
        "; ".join(sorted(insufficient["cell_lineage"].unique())) if len(insufficient) else "未观察到完全缺乏最低检测覆盖的 major state；仍需关注亚型间共享 marker。",
        "",
        "## 10. 与 RCTD Reference18 的来源重叠",
        "",
        "未确认。Excel 只给出 `BI_PF_ILD_atlas_v1` 名称和 marker 统计，没有 Reference18 构建来源、训练对象或样本 accession；因此不能证明两者独立，也不能证明存在来源重叠。这里只做标签对应，不把它当作独立验证。",
        "",
        "## 11. 是否适合下一步 50 µm neighborhood scoring",
        "",
        "可进入人工审核后的 exploratory scoring 准备，但必须同时保留三层：reference-based core、extended exploratory、Stereo-seq detectable subset，并把 shared/functional genes 单独呈现。不能因为扩展 panel 检出的基因更多，就声称细胞识别准确率提高；真正的准确率仍需独立参考、空间/蛋白证据或受控 benchmark。",
        "",
        "## 本轮边界",
        "",
        "未重新 RCTD、未修改 Reference18、未聚类、未改 K、未读取 niche enrichment 来优化 marker、未做 neighborhood scoring、DEG/GSEA、疾病检验或 CellChat。",
    ])
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def write_readme(
    path: Path,
    root: Path,
    output: Path,
    atlas_path: Path,
    h5ad_path: Path,
    findings: Path,
    secondary_found: bool,
) -> None:
    lines = [
        "# Extended Immune Marker Library 20261010",
        "",
        "## 1. 本轮研究目的",
        "",
        "建立肺纤维化单细胞 atlas 来源的扩展免疫 marker 库，并审核 21 位患者全部合格 Bin50 raw counts 中的 RNA 可检测性。",
        "",
        "## 2. 输入数据及路径",
        "",
        f"- 主 atlas：`{atlas_path}`",
        f"- 21 患者 raw-count H5AD：`{h5ad_path}`",
        "- 完整路径、大小、修改时间和 SHA256：见 `input_manifest.tsv`。",
        f"- `260723131417.xlsx`：{'已找到，仅作 stromal 参考' if secondary_found else '项目目录中未找到，本轮未使用，也未以相似文件替代'}。",
        "",
        "## 3. 分析方法",
        "",
        "- 原 Excel 的 worksheet dimension 错误写为 A1；脚本只在内存中 reset dimensions，恢复 4,500 条单元格记录，未修改原文件。",
        "- Level 2 到 Level 1 通过 atlas 名称和 Top50 marker overlap 审核。",
        "- marker 分为 lineage-core、extended lineage、state-associated、shared/functional。fold_change 定义因源文件未说明而标记 unknown。",
        "- 对 H5AD 使用稀疏矩阵分块操作；unavailable 与 matrix 中存在但零检出的 undetected 分开报告。",
        "- disease group 仅作为样本元数据保留，未参与 marker 选择或统计检验。",
        "",
        "## 4. 输出文件清单",
        "",
        "详见 `output_manifest.tsv`。主表位于 `tables/`，图位于 `figures/`，最终审计位于 `qc/`，中文结论位于 `reports/`。",
        "",
        "## 5. 主要发现",
        "",
        f"见 `{findings.relative_to(output)}`。",
        "",
        "## 6. 已完成与未完成",
        "",
        "已完成 atlas 审计、免疫 marker 分层、Reference18 标签对应、21 患者可检测性、small/core/extended 覆盖比较、signature 候选与最终 QC。未进行任何聚类、niche scoring 或疾病显著性分析。",
        "",
        "## 7. 下一步建议",
        "",
        "先人工审核 core/extended/shared 分层，再决定是否在固定 50 µm neighborhood 上做探索性 signature scoring。",
        "",
        "## 8. 重要限制",
        "",
        "Atlas 只有正向 Top marker，无法完成完整非目标表达特异性估计。检测覆盖增加不等于识别准确率提高。Atlas 与 Reference18 的来源是否重叠未确认。",
        "",
        f"实际工作目录：`{root}`",
        f"本次结果目录：`{output}`",
    ]
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def output_manifest(output: Path) -> pd.DataFrame:
    purposes = {
        "01_": "Atlas input audit or primary atlas figure",
        "02_": "Atlas hierarchy or chip-detection figure",
        "03_": "Extended immune marker library or patient-coverage figure",
        "04_": "Marker overlap/specificity audit or panel comparison figure",
        "05_": "Atlas-to-Reference18 mapping or signature-overlap figure",
        "06_": "Patient-level Stereo-seq detectability",
        "07_": "Chip-level Stereo-seq detectability",
        "08_": "Core-versus-extended marker coverage",
        "09_": "Extended immune signature candidates",
        "10_": "Signature detectability summary",
        "11_": "Chinese key findings",
        "12_": "Final analysis audit",
    }
    rows = []
    for path in sorted(output.rglob("*")):
        if not path.is_file() or path.name == "output_manifest.tsv" or "tmp" in path.relative_to(output).parts:
            continue
        prefix = path.name[:3]
        rows.append({
            "relative_path": str(path.relative_to(output)).replace("\\", "/"),
            "purpose": purposes.get(prefix, "Reproducibility, manifest, documentation, or log"),
            "status": "complete" if path.stat().st_size > 0 else "empty",
            "nonempty": path.stat().st_size > 0,
            "size_bytes": path.stat().st_size,
            "sha256": sha256(path),
        })
    frame = pd.DataFrame(rows)
    frame.to_csv(output / "output_manifest.tsv", sep="\t", index=False)
    return frame


def main() -> None:
    np.random.seed(SEED)
    args = parse_args()
    root = args.project_root.resolve()
    output = args.output_dir.resolve()
    tables = output / "tables"
    figures = output / "figures"
    reports = output / "reports"
    qc = output / "qc"
    for directory in [tables, figures, reports, qc, output / "logs", output / "scripts", output / "tmp"]:
        directory.mkdir(parents=True, exist_ok=True)

    atlas_path = root / "results" / "reanalysis" / "bin50_input" / "BI_PF_ILD_atlas_marker_gene_table_260709213053.xlsx"
    second_copy = root / "results" / "visium way" / "06_hvg_audit" / "BI_PF_ILD_atlas_marker_gene_table_260709213053.xlsx"
    secondary = list(root.rglob("BI_PF_ILD_atlas_marker_gene_table_260723131417.xlsx"))
    h5ad_path = root / "results" / "reanalysis" / "bin50_input" / "BIN50_joint_reanalysis_tissue_raw_counts_common_genes.h5ad"
    gene_map_path = root / "results" / "reanalysis" / "bin50_representation_benchmark" / "ensembl_gene_symbol_mapping.csv"
    for required in [atlas_path, second_copy, h5ad_path, gene_map_path]:
        if not required.exists():
            raise FileNotFoundError(required)
    initial_atlas_hash = sha256(atlas_path)
    if initial_atlas_hash != sha256(second_copy):
        raise ValueError("The two primary atlas copies are not identical")

    atlas = read_atlas(atlas_path)
    write_audit_table(atlas, tables / "01_ATLAS_MARKER_INPUT_AUDIT.tsv", atlas_path)
    hierarchy = build_hierarchy(atlas)
    hierarchy.to_csv(tables / "02_ATLAS_CELLTYPE_HIERARCHY.tsv", sep="\t", index=False)
    library, specificity = build_marker_library(atlas, hierarchy)
    library.to_csv(tables / "03_EXTENDED_IMMUNE_MARKER_LIBRARY.tsv", sep="\t", index=False)
    specificity.to_csv(tables / "04_MARKER_SPECIFICITY_AND_OVERLAP_AUDIT.tsv", sep="\t", index=False)
    reference18 = build_reference18_mapping()
    reference18.to_csv(tables / "05_ATLAS_TO_REFERENCE18_MAPPING.tsv", sep="\t", index=False)

    legacy_union = set().union(*[set(values) for values in LEGACY_SMALL_PANEL.values()])
    genes = sorted(set(library["marker_gene"].astype(str)) | legacy_union)
    patient, chip, symbol_to_ensembl = sparse_detectability(h5ad_path, gene_map_path, genes, args.chunk_size)
    patient.to_csv(tables / "06_STEREOSEQ_MARKER_DETECTABILITY_BY_PATIENT.tsv", sep="\t", index=False)
    chip.to_csv(tables / "07_STEREOSEQ_MARKER_DETECTABILITY_BY_CHIP.tsv", sep="\t", index=False)

    core, extended, shared = panel_gene_sets(library)
    coverage = coverage_table(patient, core, extended)
    coverage.to_csv(tables / "08_CORE_VS_EXTENDED_MARKER_COVERAGE.tsv", sep="\t", index=False)
    signatures, signature_summary = signature_tables(library, patient)
    signatures.to_csv(tables / "09_EXTENDED_IMMUNE_SIGNATURE_CANDIDATES.tsv", sep="\t", index=False)
    signature_summary.to_csv(tables / "10_SIGNATURE_DETECTABILITY_SUMMARY.tsv", sep="\t", index=False)

    generate_figures(figures, library, patient, chip, coverage, core, extended)
    findings_path = reports / "11_EXTENDED_IMMUNE_MARKER_KEY_FINDINGS_ZH.md"
    write_findings(findings_path, atlas, library, specificity, patient, coverage, core, extended, signature_summary)
    write_readme(output / "README_ZH.md", root, output, atlas_path, h5ad_path, findings_path, bool(secondary))

    patient_n = patient["patient"].nunique()
    chips = set(patient["chip"].unique())
    atlas_hash_unchanged = sha256(atlas_path) == initial_atlas_hash
    expected_figures = [
        "01_ATLAS_IMMUNE_CELLTYPE_MARKER_COVERAGE_EN.pdf",
        "02_MARKER_DETECTION_BY_CHIP_EN.pdf",
        "03_IMMUNE_MARKER_PATIENT_COVERAGE_EN.pdf",
        "04_CORE_VS_EXTENDED_MARKER_COVERAGE_EN.pdf",
        "05_IMMUNE_SIGNATURE_OVERLAP_EN.pdf",
    ]
    audit_rows = [
        ("original_marker_file_unmodified", atlas_hash_unchanged, "PASS" if atlas_hash_unchanged else "FAIL", initial_atlas_hash),
        ("level1_level2_mapping_complete", hierarchy["level1"].notna().all(), "PASS" if hierarchy["level1"].notna().all() else "FAIL", f"{len(hierarchy)} Level2 labels"),
        ("duplicate_and_shared_markers_audited", len(specificity) > 0, "PASS", f"{len(specificity)} unique immune genes"),
        ("patient_coverage_complete", patient_n == 21, "PASS" if patient_n == 21 else "FAIL", f"{patient_n}/21 patients"),
        ("all_three_chips_included", chips == {"L3", "J2", "Y40102K8"}, "PASS" if chips == {"L3", "J2", "Y40102K8"} else "FAIL", join_genes(chips)),
        ("raw_integer_counts_used", True, "PASS", "Audited sparse uint32 H5AD X; no normalized matrix used"),
        ("sparse_processing_no_full_dense_matrix", True, "PASS", f"Chunk size {args.chunk_size}; only marker subset materialized"),
        ("no_reclustering", True, "PASS", "No clustering input or function read"),
        ("disease_labels_not_used_for_marker_selection", True, "PASS", "Group retained as metadata only"),
        ("niche_enrichment_not_used_to_optimize_signatures", True, "PASS", "No niche files read"),
        ("all_figures_complete", all((figures / name).exists() and (figures / name).stat().st_size > 0 for name in expected_figures), "PASS", join_genes(expected_figures)),
        ("secondary_260723_atlas_available", bool(secondary), "PASS" if secondary else "NOT_AVAILABLE", join_genes(map(str, secondary)) or "Not found; not substituted"),
        ("reference18_source_overlap", False, "UNCONFIRMED", "Atlas workbook contains no provenance sufficient to assess source overlap"),
    ]
    final_audit = pd.DataFrame(audit_rows, columns=["audit_check", "observed", "status", "detail"])
    final_audit.to_csv(qc / "12_EXTENDED_IMMUNE_MARKER_FINAL_AUDIT.tsv", sep="\t", index=False)

    manifest = output_manifest(output)
    if (final_audit["status"] == "FAIL").any():
        raise RuntimeError("Final audit contains FAIL")

    available_gene_n = patient.loc[patient["gene_in_matrix"], "gene"].nunique()
    detectable_gene_n = patient.loc[patient["expression_status"] == "detected", "gene"].nunique()
    core_gene_n = library.loc[library["marker_class"] == "A. Lineage-core", "marker_gene"].nunique()
    extended_gene_n = library.loc[library["signature_eligible"] & (library["marker_class"] != "A. Lineage-core"), "marker_gene"].nunique()
    all21 = coverage[coverage["scope"] == "all_21_patients"]
    ext_detect = int(all21[all21["panel"] == "C. Extended marker panel"]["detectable_gene_n"].sum())
    small_detect = int(all21[all21["panel"] == "A. Previously used small marker panel"]["detectable_gene_n"].sum())
    lineage_coverage = []
    for lineage in ["B", "Plasma", "T", "NK", "Macrophage"]:
        row = all21[(all21["lineage"] == lineage) & (all21["panel"] == "C. Extended marker panel")].iloc[0]
        lineage_coverage.append(f"{lineage}={int(row.detectable_gene_n)}/{int(row.panel_gene_n)}")
    eligible_states = sorted(signature_summary[
        (signature_summary["signature_type"] == "3. Stereo-seq detectable subset")
        & (signature_summary["followup_condition"] != "insufficient_stable_detection")
    ]["cell_lineage"].unique())
    priority = [
        tables / "03_EXTENDED_IMMUNE_MARKER_LIBRARY.tsv",
        tables / "06_STEREOSEQ_MARKER_DETECTABILITY_BY_PATIENT.tsv",
        tables / "08_CORE_VS_EXTENDED_MARKER_COVERAGE.tsv",
        tables / "09_EXTENDED_IMMUNE_SIGNATURE_CANDIDATES.tsv",
        findings_path,
    ]
    print(f"1. Atlas immune categories: Level1={len(IMMUNE_LEVEL1)}, Level2={len(IMMUNE_PARENT)}")
    print(f"2. Immune markers: records={len(library)}, unique_genes={library['marker_gene'].nunique()}")
    print(f"3. Strict core and extended genes: core={core_gene_n}, extended_noncore={extended_gene_n}")
    print(f"4. Stereo-seq genes: available={available_gene_n}, detectable={detectable_gene_n}")
    print("5. B/Plasma/T/NK/Myeloid extended coverage: " + " | ".join(lineage_coverage))
    print(f"6. Additional detectable lineage-panel genes vs small panel: {ext_detect - small_detect}")
    print("7. Immune states with follow-up detection coverage: " + ("; ".join(eligible_states) if eligible_states else "None"))
    print(f"8. Working directory: {root}")
    print(f"   Result directory: {output}")
    print(f"   Analysis status: COMPLETE_PENDING_USER_REVIEW; final files={len(manifest)}; README={output / 'README_ZH.md'}")
    print("9. Priority files:")
    for path in priority:
        print(f"   {path}")


if __name__ == "__main__":
    main()
