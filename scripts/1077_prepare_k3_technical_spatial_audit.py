#!/usr/bin/env python3
"""Prepare patient-niche QC covariates and unsmoothed spatial counts.

Representative genes are selected from the already frozen robust DEG table for
visual auditing only. No additional discovery or threshold optimization occurs.
"""

from __future__ import annotations

from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import scipy.sparse as sp


ROOT = Path(__file__).resolve().parents[1]
OUT = (
    ROOT
    / "results/reanalysis/ALL21_K3_NICHE_PSEUDOBULK_DEG_20261011_201008"
)
ASSIGNMENTS = (
    ROOT
    / "results/reanalysis/ALL21_SHARED_UNBIASED_DISCRETE_NICHE_ATLAS_20261011_191652_attempt4"
    / "02_ALL21_SHARED_NICHE_ASSIGNMENTS.tsv.gz"
)
H5AD = (
    ROOT
    / "results/reanalysis/bin50_input"
    / "BIN50_joint_reanalysis_tissue_raw_counts_common_genes.h5ad"
)
GENES = ["SCGB1A1", "DRC1", "MS4A8", "DNAI1", "NWD1", "NEK5"]


def main() -> None:
    assignments = pd.read_csv(ASSIGNMENTS, sep="\t")
    assignments = assignments.loc[
        assignments["niche"].isin(["CN1", "CN2", "CN3"])
    ].copy()
    gene_meta = pd.read_csv(OUT / "02_GENE_METADATA.tsv", sep="\t")
    symbol_rows = gene_meta.loc[gene_meta["gene_symbol"].isin(GENES)].copy()
    duplicates = symbol_rows.groupby("gene_symbol").size()
    if not duplicates.reindex(GENES).fillna(0).eq(1).all():
        raise RuntimeError(f"Representative symbols are not one-to-one: {duplicates}")
    symbol_rows = symbol_rows.set_index("gene_symbol").loc[GENES].reset_index()

    adata = ad.read_h5ad(H5AD, backed="r")
    obs_names = pd.Index(adata.obs_names.astype(str))
    bin_ids = pd.Index(assignments["bin_id"].astype(str))
    positions = obs_names.get_indexer(bin_ids)
    if np.any(positions < 0):
        raise RuntimeError("Some K3 bins are absent from the raw H5AD")

    qc_columns = [
        "total_counts_common_raw",
        "n_genes_common_raw",
        "pct_counts_mt",
        "total_counts",
        "n_genes_by_counts",
    ]
    obs = adata.obs.iloc[positions][qc_columns].reset_index(drop=True)
    meta = assignments[
        [
            "bin_id",
            "patient_id",
            "chip_label",
            "disease",
            "niche",
            "x_um",
            "y_um",
        ]
    ].reset_index(drop=True)
    bin_qc = pd.concat([meta, obs], axis=1)
    qc_summary = (
        bin_qc.groupby(
            ["patient_id", "chip_label", "disease", "niche"], observed=True
        )
        .agg(
            bin_n=("bin_id", "size"),
            median_common_umi=("total_counts_common_raw", "median"),
            mean_common_umi=("total_counts_common_raw", "mean"),
            median_common_ngene=("n_genes_common_raw", "median"),
            mean_common_ngene=("n_genes_common_raw", "mean"),
            median_pct_mt=("pct_counts_mt", "median"),
            mean_pct_mt=("pct_counts_mt", "mean"),
            median_input_total_counts=("total_counts", "median"),
            median_input_ngene=("n_genes_by_counts", "median"),
        )
        .reset_index()
    )
    qc_summary.insert(
        0,
        "pseudobulk_id",
        qc_summary["patient_id"].astype(str) + "::" + qc_summary["niche"],
    )
    qc_summary["log10_bin_n"] = np.log10(qc_summary["bin_n"])
    qc_summary["log10_median_common_umi_plus1"] = np.log10(
        qc_summary["median_common_umi"] + 1
    )
    qc_summary["log10_median_common_ngene_plus1"] = np.log10(
        qc_summary["median_common_ngene"] + 1
    )
    qc_summary.to_csv(OUT / "30_PATIENT_NICHE_QC_COVARIATES.tsv", sep="\t", index=False)

    var_index = pd.Index(adata.var_names.astype(str))
    gene_positions = var_index.get_indexer(symbol_rows["ensembl_id"].astype(str))
    if np.any(gene_positions < 0):
        raise RuntimeError("Representative genes are absent from H5AD")
    raw = adata.X[:, gene_positions]
    if sp.issparse(raw):
        raw = raw.toarray()
    else:
        raw = np.asarray(raw)
    selected_raw = raw[positions]
    if hasattr(adata, "file"):
        adata.file.close()
    if np.any(selected_raw < 0) or np.any(selected_raw != np.floor(selected_raw)):
        raise RuntimeError("Representative spatial values are not raw integer counts")

    spatial = meta.copy()
    for index, symbol in enumerate(GENES):
        spatial[symbol] = selected_raw[:, index].astype(np.uint32)
    spatial.to_csv(
        OUT / "31_REPRESENTATIVE_GENE_UNSMOOTHED_SPATIAL_COUNTS.tsv.gz",
        sep="\t",
        index=False,
        compression="gzip",
    )
    validation = pd.DataFrame(
        [
            ("evaluable_bins", len(spatial), 126_199, "PASS"),
            ("patient_niche_qc_rows", len(qc_summary), 63, "PASS"),
            ("representative_genes", len(GENES), 6, "PASS"),
            ("negative_values", int((selected_raw < 0).sum()), 0, "PASS"),
            (
                "noninteger_values",
                int((selected_raw != np.floor(selected_raw)).sum()),
                0,
                "PASS",
            ),
        ],
        columns=["check", "observed", "expected", "status"],
    )
    validation.to_csv(OUT / "32_TECHNICAL_SPATIAL_PREPARATION_VALIDATION.tsv", sep="\t", index=False)
    print(OUT)


if __name__ == "__main__":
    main()
