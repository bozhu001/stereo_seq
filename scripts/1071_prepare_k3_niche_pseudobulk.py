#!/usr/bin/env python3
"""Prepare raw-count patient-by-K3-niche pseudobulks for step 7.

This script is deliberately limited to data preparation and audit. It reads the
frozen disease-blind K=3 niche assignments and the joint Bin50 raw integer-count
H5AD, excludes only NOT_EVALUABLE bins, and sums counts once per patient/niche.
It does not rerun RCTD, neighborhood construction, or clustering.

Exact command
-------------
envs/stereo_seq_py/Scripts/python.exe scripts/1071_prepare_k3_niche_pseudobulk.py
"""

from __future__ import annotations

import argparse
import gzip
import hashlib
import json
import shutil
from datetime import datetime
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import scipy.io
import scipy.sparse as sp


ROOT = Path(__file__).resolve().parents[1]
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
GENE_MAP = (
    ROOT
    / "results/reanalysis/bin50_representation_benchmark"
    / "ensembl_gene_symbol_mapping.csv"
)
OUT_PARENT = ROOT / "results/reanalysis"
NICHES = ["CN1", "CN2", "CN3"]


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def write_tsv(frame: pd.DataFrame, path: Path) -> None:
    frame.to_csv(path, sep="\t", index=False)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path)
    parser.add_argument("--chunk-size", type=int, default=2048)
    args = parser.parse_args()

    stamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    output = args.output or (
        OUT_PARENT / f"ALL21_K3_NICHE_PSEUDOBULK_DEG_{stamp}"
    )
    output.mkdir(parents=True, exist_ok=False)
    (output / "work").mkdir()
    (output / "logs").mkdir()

    assignments = pd.read_csv(ASSIGNMENTS, sep="\t")
    required = {
        "bin_id",
        "patient_id",
        "chip_label",
        "disease",
        "niche",
        "nUMI",
    }
    require(required.issubset(assignments.columns), "Assignment columns missing")
    require(assignments["bin_id"].is_unique, "Assignment bin_id is not unique")
    analysis_bins = assignments.loc[assignments["niche"].isin(NICHES)].copy()
    require(len(analysis_bins) == 126_199, "Unexpected evaluable K3 bin count")
    require(analysis_bins["patient_id"].nunique() == 21, "Expected 21 patients")
    require(set(analysis_bins["niche"]) == set(NICHES), "Expected CN1-CN3")

    sample_meta = (
        analysis_bins.groupby(
            ["patient_id", "chip_label", "disease", "niche"], observed=True
        )
        .agg(bin_n=("bin_id", "size"), assignment_n_umi=("nUMI", "sum"))
        .reset_index()
        .sort_values(["niche", "chip_label", "disease", "patient_id"])
        .reset_index(drop=True)
    )
    require(len(sample_meta) == 63, "Expected 21 patients x 3 niches")
    sample_meta.insert(
        0,
        "pseudobulk_id",
        sample_meta["patient_id"].astype(str) + "::" + sample_meta["niche"],
    )
    sample_meta["low_bin_coverage_lt20"] = sample_meta["bin_n"] < 20
    sample_index = dict(zip(sample_meta["pseudobulk_id"], sample_meta.index))
    bin_to_sample = dict(
        zip(
            analysis_bins["bin_id"].astype(str),
            (
                analysis_bins["patient_id"].astype(str)
                + "::"
                + analysis_bins["niche"].astype(str)
            ).map(sample_index),
        )
    )

    adata = ad.read_h5ad(H5AD, backed="r")
    require(adata.X is not None, "H5AD X is missing")
    require(adata.n_vars == 27_618, "Unexpected common-gene count")
    h5_names = pd.Index(adata.obs_names.astype(str))
    missing = pd.Index(analysis_bins["bin_id"].astype(str)).difference(h5_names)
    require(len(missing) == 0, f"{len(missing)} K3 bins absent from H5AD")

    pseudobulk = np.zeros((len(sample_meta), adata.n_vars), dtype=np.int64)
    h5ad_common_library = np.zeros(len(sample_meta), dtype=np.int64)
    selected_rows = 0
    noninteger_values = 0
    negative_values = 0
    for start in range(0, adata.n_obs, args.chunk_size):
        end = min(start + args.chunk_size, adata.n_obs)
        names = h5_names[start:end]
        local = np.array([bin_to_sample.get(name, -1) for name in names], dtype=int)
        keep = np.flatnonzero(local >= 0)
        if not len(keep):
            continue
        matrix = adata.X[start:end]
        if not sp.issparse(matrix):
            matrix = sp.csr_matrix(np.asarray(matrix))
        else:
            matrix = matrix.tocsr()
        matrix = matrix[keep]
        if matrix.nnz:
            negative_values += int(np.sum(matrix.data < 0))
            noninteger_values += int(np.sum(matrix.data != np.floor(matrix.data)))
        selected_rows += len(keep)
        local_samples = local[keep]
        common_totals = np.asarray(
            adata.obs.iloc[start:end]["total_counts_common_raw"], dtype=np.int64
        )[keep]
        h5ad_common_library += np.bincount(
            local_samples,
            weights=common_totals,
            minlength=len(sample_meta),
        ).astype(np.int64)
        membership = sp.csr_matrix(
            (
                np.ones(len(keep), dtype=np.int8),
                (np.arange(len(keep)), local_samples),
            ),
            shape=(len(keep), len(sample_meta)),
        )
        pseudobulk += (membership.T @ matrix).toarray().astype(np.int64)
    if hasattr(adata, "file"):
        adata.file.close()

    require(selected_rows == len(analysis_bins), "Not every K3 bin was summed once")
    require(negative_values == 0, "Negative raw counts detected")
    require(noninteger_values == 0, "Non-integer raw counts detected")
    sample_meta["raw_library_size"] = pseudobulk.sum(axis=1)
    sample_meta["h5ad_total_counts_common_raw"] = h5ad_common_library
    sample_meta["raw_minus_h5ad_common_total"] = (
        sample_meta["raw_library_size"]
        - sample_meta["h5ad_total_counts_common_raw"]
    )
    sample_meta["library_minus_assignment_n_umi"] = (
        sample_meta["raw_library_size"] - sample_meta["assignment_n_umi"]
    )

    gene_map = pd.read_csv(GENE_MAP)
    require(len(gene_map) == adata.n_vars, "Gene map length mismatch")
    require(
        gene_map["ensembl_var_name"].astype(str).tolist()
        == list(map(str, adata.var_names)),
        "Gene-map order differs from H5AD var order",
    )
    gene_meta = gene_map.rename(
        columns={"ensembl_var_name": "ensembl_id"}
    ).copy()
    gene_meta.insert(0, "gene_index_1based", np.arange(1, len(gene_meta) + 1))

    write_tsv(sample_meta, output / "01_PSEUDOBULK_SAMPLE_METADATA.tsv")
    write_tsv(gene_meta, output / "02_GENE_METADATA.tsv")
    matrix_path = output / "03_PSEUDOBULK_RAW_COUNTS.mtx"
    scipy.io.mmwrite(matrix_path, sp.coo_matrix(pseudobulk.T), field="integer")
    with matrix_path.open("rb") as source, gzip.open(
        str(matrix_path) + ".gz", "wb", compresslevel=6
    ) as target:
        shutil.copyfileobj(source, target)
    matrix_path.unlink()

    input_manifest = pd.DataFrame(
        [
            {
                "role": "frozen_k3_assignments",
                "path": str(ASSIGNMENTS.relative_to(ROOT)),
                "sha256": sha256(ASSIGNMENTS),
            },
            {
                "role": "raw_common_gene_h5ad",
                "path": str(H5AD.relative_to(ROOT)),
                "sha256": sha256(H5AD),
            },
            {
                "role": "ensembl_symbol_map",
                "path": str(GENE_MAP.relative_to(ROOT)),
                "sha256": sha256(GENE_MAP),
            },
        ]
    )
    write_tsv(input_manifest, output / "00_INPUT_MANIFEST.tsv")
    validation = pd.DataFrame(
        [
            ("evaluable_k3_bins", selected_rows, 126_199, "PASS"),
            ("patients", sample_meta["patient_id"].nunique(), 21, "PASS"),
            ("patient_niche_libraries", len(sample_meta), 63, "PASS"),
            ("genes", pseudobulk.shape[1], 27_618, "PASS"),
            ("negative_count_values", negative_values, 0, "PASS"),
            ("noninteger_count_values", noninteger_values, 0, "PASS"),
            (
                "total_counts_reconcile_h5ad_common_raw",
                int(pseudobulk.sum()),
                int(h5ad_common_library.sum()),
                "PASS"
                if int(pseudobulk.sum()) == int(h5ad_common_library.sum())
                else "FAIL",
            ),
            (
                "rctd_nUMI_total_for_context_only",
                int(analysis_bins["nUMI"].sum()),
                "not expected to equal common-gene raw H5AD total",
                "INFO",
            ),
            (
                "low_coverage_patient_niche_lt20_bins",
                int(sample_meta["low_bin_coverage_lt20"].sum()),
                "reported_not_hidden",
                "PASS",
            ),
        ],
        columns=["check", "observed", "expected", "status"],
    )
    write_tsv(validation, output / "04_PREPARATION_VALIDATION.tsv")
    run = {
        "created_at": datetime.now().astimezone().isoformat(timespec="seconds"),
        "script": str(Path(__file__).relative_to(ROOT)),
        "output": str(output.relative_to(ROOT)),
        "chunk_size": args.chunk_size,
        "scope": "frozen K3; raw Bin50 counts; patient-by-niche pseudobulk",
        "excluded": "NOT_EVALUABLE only",
    }
    (output / "RUN_PARAMETERS.json").write_text(
        json.dumps(run, indent=2), encoding="utf-8"
    )
    (output / "work/PREP_COMPLETE.ok").write_text("PASS\n", encoding="utf-8")
    print(output)


if __name__ == "__main__":
    main()
