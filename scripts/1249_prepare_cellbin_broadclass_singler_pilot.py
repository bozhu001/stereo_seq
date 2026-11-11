"""Prepare frozen Reference18 pseudobulks and a balanced CellBin pilot.

Reference cells are aggregated by donor and unchanged Reference18 state. Query
pilot cells are deterministically sampled across the nGenes range within each
of all 21 samples (500 per sample). Full matrices are always sparse.
"""

from __future__ import annotations

import json
from pathlib import Path

import h5py
import numpy as np
import pandas as pd
import pyarrow.parquet as pq
from scipy import sparse


ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "results/cell bin/CellBin_BroadClass_pilot_21samples"
WORK = OUT / "work"
REF_ROOT = (
    ROOT
    / "results/reanalysis/RCTD_SIX_DATASET_LEVEL1_ALL21_20260916_223802"
)
REF_H5 = REF_ROOT / "prepared_inputs/six_dataset_level1_rctd_inputs.h5"
REF_META = REF_ROOT / "prepared_inputs/selected_reference_cells.tsv.gz"
REF_MAPPING = (
    ROOT
    / "results/reanalysis/RCTD_LEVEL2_REFERENCE18_ALL21_20260919_150957/"
    "prepared_inputs/REFERENCE_18_CELL_MAPPING.tsv"
)
JOINT = (
    ROOT
    / "results/cell bin/CellBin_joint_raw_21samples/"
    "Joint_CellBin_raw_with_sampleID.h5ad"
)
QC = (
    ROOT
    / "results/cell bin/CellBin_QC_diagnostics_21samples/"
    "CellBin_QC_metrics.parquet"
)
MAP_FILE = OUT / "Reference18_to_BroadClass_mapping.tsv"
N_PER_SAMPLE = 500


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def decode(values: np.ndarray) -> list[str]:
    return [
        item.decode() if isinstance(item, (bytes, np.bytes_)) else str(item)
        for item in values
    ]


def read_h5ad_frame(handle: h5py.File, group_name: str) -> pd.DataFrame:
    group = handle[group_name]
    index_key = group.attrs.get("_index", "_index")
    if isinstance(index_key, bytes):
        index_key = index_key.decode()
    columns: dict[str, object] = {}
    for key, obj in group.items():
        if key == index_key:
            continue
        if isinstance(obj, h5py.Dataset):
            values = obj[:]
            columns[key] = decode(values) if values.dtype.kind in "SOU" else values
        elif isinstance(obj, h5py.Group) and {
            "categories",
            "codes",
        }.issubset(obj.keys()):
            categories = decode(obj["categories"][:])
            columns[key] = pd.Categorical.from_codes(
                obj["codes"][:].astype(int), categories
            )
    return pd.DataFrame(columns, index=decode(group[index_key][:]))


def read_custom_csc(group: h5py.Group) -> sparse.csc_matrix:
    shape = tuple(int(x) for x in group["shape"][:])
    return sparse.csc_matrix(
        (group["data"][:], group["indices"][:], group["indptr"][:]),
        shape=shape,
    )


def read_h5ad_csr(group: h5py.Group) -> sparse.csr_matrix:
    shape = tuple(int(x) for x in group.attrs["shape"])
    return sparse.csr_matrix(
        (group["data"][:], group["indices"][:], group["indptr"][:]),
        shape=shape,
    )


def write_csc(path: Path, matrix: sparse.spmatrix, row_names: list[str],
              col_names: list[str]) -> None:
    matrix = matrix.tocsc()
    with h5py.File(path, "w") as handle:
        group = handle.create_group("matrix")
        group.create_dataset("data", data=matrix.data, compression="gzip")
        group.create_dataset("indices", data=matrix.indices, compression="gzip")
        group.create_dataset("indptr", data=matrix.indptr, compression="gzip")
        group.create_dataset("shape", data=np.asarray(matrix.shape, dtype=np.int64))
        strings = h5py.string_dtype("utf-8")
        handle.create_dataset("row_names", data=np.asarray(row_names, dtype=object), dtype=strings)
        handle.create_dataset("col_names", data=np.asarray(col_names, dtype=object), dtype=strings)


def choose_pilot(qc: pd.DataFrame) -> np.ndarray:
    selected: list[int] = []
    for _, frame in qc.groupby("sample_id", sort=True):
        frame = frame.sort_values(["n_genes", "n_umi", "global_cellbin_id"])
        require(len(frame) >= N_PER_SAMPLE, "A sample has too few passed cells")
        positions = np.linspace(0, len(frame) - 1, N_PER_SAMPLE).round().astype(int)
        selected.extend(frame.iloc[positions]["row_index"].astype(int).tolist())
    selected_array = np.asarray(selected, dtype=np.int64)
    require(len(selected_array) == 21 * N_PER_SAMPLE, "Pilot size mismatch")
    require(len(np.unique(selected_array)) == len(selected_array), "Pilot IDs repeat")
    return selected_array


def main() -> None:
    WORK.mkdir(parents=True, exist_ok=True)
    broad_map = pd.read_csv(MAP_FILE, sep="\t").set_index("state18")["broad_class"]

    ref_meta = pd.read_csv(REF_META, sep="\t")
    mapping = pd.read_csv(REF_MAPPING, sep="\t")
    with h5py.File(REF_H5, "r") as handle:
        ref_cells = decode(handle["reference_cells"][:])
        genes = decode(handle["genes"][:])
        ref_counts = read_custom_csc(handle["reference"])
    require(ref_counts.shape == (18_966, 57_492), "Reference shape mismatch")
    require(ref_cells == mapping["cell_id"].tolist(), "Reference mapping order mismatch")
    ref_meta = ref_meta.set_index("cell_id").loc[ref_cells].reset_index()
    ref_meta["reference18_state"] = mapping["rctd_state"].to_numpy()
    ref_meta["broad_class"] = ref_meta["reference18_state"].map(broad_map)
    require(ref_meta["broad_class"].notna().all(), "Unmapped Reference18 state")

    group_table = (
        ref_meta[["donor_id", "study", "reference18_state", "broad_class"]]
        .drop_duplicates()
        .sort_values(["donor_id", "reference18_state"])
        .reset_index(drop=True)
    )
    group_table["profile_id"] = [f"reference_profile_{i + 1:04d}" for i in range(len(group_table))]
    key_to_group = {
        (row.donor_id, row.reference18_state): i
        for i, row in group_table.iterrows()
    }
    group_index = np.asarray(
        [key_to_group[(row.donor_id, row.reference18_state)] for row in ref_meta.itertuples()],
        dtype=np.int32,
    )
    membership = sparse.csc_matrix(
        (
            np.ones(len(ref_meta), dtype=np.uint8),
            (np.arange(len(ref_meta), dtype=np.int32), group_index),
        ),
        shape=(len(ref_meta), len(group_table)),
    )
    profile_counts = (ref_counts @ membership).tocsc()
    cell_counts = np.asarray(membership.sum(axis=0)).ravel().astype(int)
    group_table["n_reference_cells"] = cell_counts
    group_table["total_counts"] = np.asarray(profile_counts.sum(axis=0)).ravel()
    require(int(cell_counts.sum()) == 57_492, "Reference cells lost in aggregation")
    group_table.to_csv(WORK / "reference_profiles.tsv", sep="\t", index=False)
    write_csc(
        WORK / "reference_profiles_counts.h5",
        profile_counts,
        genes,
        group_table["profile_id"].tolist(),
    )
    del ref_counts, membership, profile_counts

    schema = pq.read_schema(QC).names
    id_col = "global_cellbin_id" if "global_cellbin_id" in schema else "cellbin_id"
    qc = pq.read_table(
        QC,
        columns=[id_col, "cellbin_id", "sample_id", "chip_id", "disease_group", "x", "y", "n_genes", "n_umi"],
    ).to_pandas()
    if id_col != "global_cellbin_id":
        qc = qc.rename(columns={id_col: "global_cellbin_id"})
    qc["row_index"] = np.arange(len(qc), dtype=np.int64)
    passed = (
        ~qc["sample_id"].isin(["unassigned", "Unknown", ""])
        & qc["n_genes"].ge(10)
        & qc["n_umi"].ge(20)
    )
    passed_qc = qc.loc[passed].copy()
    require(len(passed_qc) == 359_329, "G10/U20 passed count mismatch")
    selected = choose_pilot(passed_qc)
    pilot_meta = qc.iloc[selected].copy()

    with h5py.File(JOINT, "r") as handle:
        obs = read_h5ad_frame(handle, "obs")
        var = read_h5ad_frame(handle, "var")
        require(obs.index.tolist() == qc["global_cellbin_id"].tolist(), "Joint/QC order mismatch")
        joint_symbols = var["gene_symbol"].astype(str).str.upper()
        ref_lookup = {gene.upper(): i for i, gene in enumerate(genes)}
        feature_rows: list[int] = []
        feature_cols: list[int] = []
        for feature_index, symbol in enumerate(joint_symbols):
            if symbol in ref_lookup:
                feature_rows.append(feature_index)
                feature_cols.append(ref_lookup[symbol])
        require(len(set(feature_cols)) == len(genes), "Not all reference genes map to Joint")
        mapper = sparse.csr_matrix(
            (
                np.ones(len(feature_rows), dtype=np.uint8),
                (np.asarray(feature_rows), np.asarray(feature_cols)),
            ),
            shape=(len(var), len(genes)),
        )
        joint_x = read_h5ad_csr(handle["X"])
        pilot_full = joint_x[selected, :]
        pilot_common = (pilot_full @ mapper).tocsr()

    full_umi = np.asarray(pilot_full.sum(axis=1)).ravel()
    common_umi = np.asarray(pilot_common.sum(axis=1)).ravel()
    pilot_meta["common_gene_umi"] = common_umi
    pilot_meta["common_gene_umi_fraction"] = np.divide(
        common_umi,
        full_umi,
        out=np.zeros_like(common_umi, dtype=float),
        where=full_umi > 0,
    )
    require(np.all(full_umi == pilot_meta["n_umi"].to_numpy()), "Pilot UMI mismatch")
    pilot_meta.to_csv(WORK / "pilot_query_metadata.tsv", sep="\t", index=False)
    write_csc(
        WORK / "pilot_query_counts.h5",
        pilot_common.T,
        genes,
        pilot_meta["global_cellbin_id"].tolist(),
    )

    summary = {
        "reference_profiles": int(len(group_table)),
        "reference_donors": int(group_table["donor_id"].nunique()),
        "reference_cells_retained": int(cell_counts.sum()),
        "pilot_query_cells": int(len(pilot_meta)),
        "pilot_samples": int(pilot_meta["sample_id"].nunique()),
        "pilot_chips": int(pilot_meta["chip_id"].nunique()),
        "pilot_median_common_umi_fraction": float(pilot_meta["common_gene_umi_fraction"].median()),
        "full_dense_matrix_created": False,
    }
    (WORK / "pilot_preparation_summary.json").write_text(
        json.dumps(summary, indent=2), encoding="utf-8"
    )
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
