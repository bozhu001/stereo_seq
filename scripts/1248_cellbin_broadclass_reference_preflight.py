"""Preflight audit for the 21-sample CellBin broad-class annotation pilot.

This script is read-only with respect to Joint V1 and all frozen Reference18
inputs. It writes only small audit/mapping/method-selection artifacts into the
new pilot result directory. No prediction is performed here.
"""

from __future__ import annotations

import hashlib
import json
from pathlib import Path

import h5py
import numpy as np
import pandas as pd
import pyarrow.parquet as pq


ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "results/cell bin/CellBin_BroadClass_pilot_21samples"
REF_ROOT = (
    ROOT
    / "results/reanalysis/RCTD_SIX_DATASET_LEVEL1_ALL21_20260916_223802"
)
REF_H5 = REF_ROOT / "prepared_inputs/six_dataset_level1_rctd_inputs.h5"
REF_META = REF_ROOT / "prepared_inputs/selected_reference_cells.tsv.gz"
REF18_ROOT = (
    ROOT
    / "results/reanalysis/RCTD_LEVEL2_REFERENCE18_ALL21_20260919_150957"
)
REF18_MAPPING = REF18_ROOT / "prepared_inputs/REFERENCE_18_CELL_MAPPING.tsv"
REF18_STATES = REF18_ROOT / "02_REFERENCE18_STATE_DEFINITION.tsv"
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

BROAD_MAP = {
    "Alveolar/resident macrophage": "Immune",
    "Inflammatory monocyte/macrophage": "Immune",
    "Dendritic cell": "Immune",
    "T": "Immune",
    "NK": "Immune",
    "Mast/Basophil": "Immune",
    "B": "Immune",
    "Plasma": "Immune",
    "Ciliated epithelial": "Epithelial",
    "AT2/alveolar transitional": "Epithelial",
    "Secretory/other airway epithelial": "Epithelial",
    "AT1": "Epithelial",
    "Basal-like epithelial": "Epithelial",
    "Blood endothelial": "Endothelial / vascular",
    "Lymphatic endothelial": "Endothelial / vascular",
    "Perivascular/mesothelial stromal": "Stromal",
    "Fibroblast": "Stromal",
    "Myofibroblast": "Stromal",
}

RATIONALE = {
    "Alveolar/resident macrophage": "Myeloid immune lineage.",
    "Inflammatory monocyte/macrophage": "Monocyte/macrophage immune lineage.",
    "Dendritic cell": "Antigen-presenting immune lineage.",
    "T": "Lymphoid immune lineage.",
    "NK": "Lymphoid immune lineage.",
    "Mast/Basophil": "Granulocytic immune lineage.",
    "B": "Lymphoid immune lineage.",
    "Plasma": "Antibody-secreting immune lineage.",
    "Ciliated epithelial": "Airway epithelial lineage.",
    "AT2/alveolar transitional": "Alveolar epithelial lineage.",
    "Secretory/other airway epithelial": "Airway epithelial lineage.",
    "AT1": "Alveolar epithelial lineage.",
    "Basal-like epithelial": "Basal epithelial lineage.",
    "Blood endothelial": "Blood-vessel endothelial lineage.",
    "Lymphatic endothelial": "Lymphatic endothelial lineage.",
    "Perivascular/mesothelial stromal": (
        "Pericyte, smooth-muscle and mesothelial stromal lineages."
    ),
    "Fibroblast": "Mesenchymal stromal lineage.",
    "Myofibroblast": "Contractile mesenchymal stromal lineage.",
}


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def decode(values: np.ndarray) -> list[str]:
    return [
        item.decode() if isinstance(item, (bytes, np.bytes_)) else str(item)
        for item in values
    ]


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        while block := handle.read(8 * 1024 * 1024):
            digest.update(block)
    return digest.hexdigest()


def h5ad_dataframe(handle: h5py.File, group_name: str) -> pd.DataFrame:
    group = handle[group_name]
    index_key = group.attrs.get("_index", "_index")
    if isinstance(index_key, bytes):
        index_key = index_key.decode()
    result: dict[str, object] = {}
    for key in group.keys():
        if key == index_key:
            continue
        obj = group[key]
        if isinstance(obj, h5py.Dataset):
            values = obj[:]
            if values.dtype.kind in {"S", "O", "U"}:
                result[key] = decode(values)
            else:
                result[key] = values
        elif isinstance(obj, h5py.Group) and {
            "categories",
            "codes",
        }.issubset(obj.keys()):
            categories = np.asarray(decode(obj["categories"][:]), dtype=object)
            codes = obj["codes"][:].astype(int)
            result[key] = pd.Categorical.from_codes(codes, categories)
    index = decode(group[index_key][:])
    return pd.DataFrame(result, index=index)


def inspect_sparse(group: h5py.Group) -> dict[str, object]:
    if "shape" in group:
        shape = tuple(int(x) for x in group["shape"][:])
    else:
        shape = tuple(int(x) for x in group.attrs["shape"])
    data = group["data"]
    sample_n = min(1_000_000, len(data))
    values = data[:sample_n]
    return {
        "shape": shape,
        "nnz": int(len(data)),
        "sample_n": sample_n,
        "sample_min": float(values.min()) if sample_n else 0.0,
        "sample_max": float(values.max()) if sample_n else 0.0,
        "sample_nonnegative": bool(np.all(values >= 0)),
        "sample_integer_like": bool(np.all(values == np.round(values))),
        "encoding_type": str(
            group.attrs.get(
                "encoding-type",
                "csc_matrix" if len(group["indptr"]) == shape[1] + 1 else "csr_matrix",
            )
        ),
    }


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    required = [REF_H5, REF_META, REF18_MAPPING, REF18_STATES, JOINT, QC]
    for path in required:
        require(path.exists(), f"Missing required input: {path}")

    states = pd.read_csv(REF18_STATES, sep="\t").sort_values("state_order")
    mapping = pd.read_csv(REF18_MAPPING, sep="\t")
    require(len(states) == 18, "Reference18 state definition does not have 18 rows")
    require(set(states["state18"]) == set(BROAD_MAP), "Broad map/state mismatch")
    require(mapping["cell_id"].is_unique, "Reference18 mapping cell IDs not unique")
    require(mapping["rctd_state"].nunique() == 18, "Reference mapping lacks states")

    map_out = states[
        ["state_order", "state18", "n_reference_cells", "original_fine_labels"]
    ].copy()
    map_out["broad_class"] = map_out["state18"].map(BROAD_MAP)
    map_out["biological_basis"] = map_out["state18"].map(RATIONALE)
    map_out["mapping_status"] = "Mapped without changing Reference18 label"
    map_out.to_csv(
        OUT / "Reference18_to_BroadClass_mapping.tsv", sep="\t", index=False
    )

    ref_meta = pd.read_csv(REF_META, sep="\t")
    with h5py.File(REF_H5, "r") as handle:
        h5_keys = sorted(handle.keys())
        ref_cells = decode(handle["reference_cells"][:])
        genes = decode(handle["genes"][:])
        ref_group = handle["reference"]
        ref_sparse = inspect_sparse(ref_group)
    require(len(ref_cells) == 57_492, "Unexpected reference-cell count")
    require(len(genes) == 18_966, "Unexpected reference-gene count")
    require(ref_cells == mapping["cell_id"].tolist(), "Reference mapping order mismatch")
    require(set(ref_cells) == set(ref_meta["cell_id"]), "Reference metadata IDs mismatch")

    with h5py.File(JOINT, "r") as handle:
        obs = h5ad_dataframe(handle, "obs")
        var = h5ad_dataframe(handle, "var")
        joint_sparse = inspect_sparse(handle["X"])
    require(len(obs) == 799_172, "Unexpected Joint V1 cell count")
    require(len(var) == 32_089, "Unexpected Joint V1 gene count")
    require(obs.index.is_unique, "Joint obs_names are not unique")
    require("gene_symbol" in var.columns, "Joint var lacks gene_symbol")

    qc_cols = pq.read_schema(QC).names
    id_col = "global_cellbin_id" if "global_cellbin_id" in qc_cols else "cellbin_id"
    qc = pq.read_table(
        QC,
        columns=[id_col, "sample_id", "chip_id", "n_genes", "n_umi"],
    ).to_pandas()
    require(len(qc) == len(obs), "QC/Joint row-count mismatch")
    require(qc[id_col].is_unique, "QC cell IDs are not unique")
    require(qc[id_col].tolist() == obs.index.tolist(), "QC/Joint ID order mismatch")
    assigned = ~qc["sample_id"].isin(["unassigned", "Unknown", "", None])
    passed = assigned & qc["n_genes"].ge(10) & qc["n_umi"].ge(20)
    require(int(passed.sum()) == 359_329, "Assigned G10/U20 count mismatch")
    require(qc.loc[passed, "sample_id"].nunique() == 21, "Not all samples represented")

    ref_gene_upper = pd.Index(pd.Series(genes, dtype=str).str.upper().unique())
    joint_symbols = var["gene_symbol"].astype(str).str.upper()
    joint_gene_upper = pd.Index(joint_symbols.unique())
    common = ref_gene_upper.intersection(joint_gene_upper)

    metadata_columns = pd.DataFrame(
        {
            "column": ref_meta.columns,
            "nonmissing": [int(ref_meta[c].notna().sum()) for c in ref_meta.columns],
            "unique_values": [int(ref_meta[c].nunique(dropna=True)) for c in ref_meta.columns],
        }
    )
    metadata_columns.to_csv(
        OUT / "Reference18_metadata_columns.tsv", sep="\t", index=False
    )
    candidate_donor_cols = [
        c
        for c in ref_meta.columns
        if any(token in c.lower() for token in ["donor", "patient", "sample", "dataset"])
    ]

    audit_rows = [
        ("reference_cells", len(ref_cells), "PASS"),
        ("reference_states", mapping["rctd_state"].nunique(), "PASS"),
        ("reference_genes", len(genes), "PASS"),
        ("reference_sparse_encoding", ref_sparse["encoding_type"], "PASS"),
        ("reference_sparse_nnz", ref_sparse["nnz"], "PASS"),
        ("reference_sample_nonnegative", ref_sparse["sample_nonnegative"], "PASS"),
        ("reference_sample_integer_like", ref_sparse["sample_integer_like"], "PASS"),
        ("reference_mapping_id_order", True, "PASS"),
        ("joint_shape", str(joint_sparse["shape"]), "PASS"),
        ("joint_sparse_encoding", joint_sparse["encoding_type"], "PASS"),
        ("joint_sample_nonnegative", joint_sparse["sample_nonnegative"], "PASS"),
        ("joint_sample_integer_like", joint_sparse["sample_integer_like"], "PASS"),
        ("qc_joint_id_order", True, "PASS"),
        ("assigned_g10_u20_cellbins", int(passed.sum()), "PASS"),
        ("samples", qc.loc[passed, "sample_id"].nunique(), "PASS"),
        ("reference_unique_gene_symbols", len(ref_gene_upper), "PASS"),
        ("joint_unique_gene_symbols", len(joint_gene_upper), "PASS"),
        ("common_unique_gene_symbols", len(common), "PASS"),
        ("candidate_donor_metadata_columns", ";".join(candidate_donor_cols), "REVIEW"),
        ("reference_h5_top_keys", ";".join(h5_keys), "PASS"),
    ]
    pd.DataFrame(audit_rows, columns=["check", "value", "status"]).to_csv(
        OUT / "Reference18_input_audit.tsv", sep="\t", index=False
    )

    state_counts = (
        mapping.assign(broad_class=mapping["rctd_state"].map(BROAD_MAP))
        .groupby(["rctd_state", "broad_class"], observed=True)
        .size()
        .rename("n_reference_cells")
        .reset_index()
    )
    state_counts.to_csv(
        OUT / "Reference18_state_cell_counts.tsv", sep="\t", index=False
    )

    method = f"""# Broad-class method selection (pre-execution)

## Selected primary method

**SingleR 2.12.0 with donor/state-aggregated frozen Reference18 profiles**, followed by SingleR pruning. Reference18 labels remain unchanged; only the documented 18-state-to-four-broad-class mapping is used for the pilot.

## Why this method

- The CellBin query is FFPE and very sparse. SingleR uses rank/correlation evidence across genes rather than requiring a dense full matrix or de novo clustering.
- Aggregating reference cells within donor × Reference18 state reduces reference noise and memory while preserving independent donor profiles when the frozen metadata supports them.
- `pruned.labels`, `delta.next`, and per-class scores provide explicit uncertainty outputs. Pruned predictions are retained as `Unknown`; no spatial or nearest-neighbor filling is allowed.
- Scores/deltas are method-specific evidence, **not calibrated probabilities or accuracy estimates**.
- Query raw counts remain unchanged. Library-size log normalization is used only as a transient input representation required for comparison.
- The main risks are FFPE sparsity, ambient RNA, reference/query technology mismatch, and chip-depth effects. These are evaluated in the pilot and full validation.

## Alternatives reviewed

- **scANVI:** not selected because `scvi-tools` is not installed in the project environment, joint training would be substantially more resource intensive, and it would add a learned latent model whose depth/batch behavior requires a separate validation program.
- **RCTD/cell2location:** not selected for the primary method because their continuous mixture outputs should not be converted directly into single-CellBin identities; this task requests provisional broad labels with Unknown.
- **Canonical-marker scoring alone:** retained only as orthogonal validation, not as the classifier.

## Frozen inputs and audit results

- Reference18: {len(ref_cells):,} cells, 18 unchanged states, {len(genes):,} raw-count genes.
- Query: 359,329 assigned G10/U20 CellBins from all 21 samples.
- Common unique gene symbols: {len(common):,}.
- Reference donor/sample candidate columns: {', '.join(candidate_donor_cols) if candidate_donor_cols else 'none detected'}.

## Pilot gate

The small pilot will include balanced cells from all three chips and the six prespecified caution samples where possible. Expansion is allowed only if:

1. all IDs and score rows match the pilot input exactly;
2. reference self-validation is materially above chance at broad-class level;
3. pruning produces a nonzero Unknown set rather than forcing every label;
4. predictions are not nearly mono-class within a chip/sample without marker support;
5. score/delta behavior is not overwhelmingly explained by nGenes/nUMI;
6. sparse processing completes without constructing the 359,329 × 32,089 matrix densely.

## Planned outputs

Each provisional query row will retain the unique Joint V1 ID, original cellbin ID, chip/sample/disease fields, x/y coordinates, predicted broad class, all broad-class scores, delta/uncertainty fields, `unknown_flag`, and the explicit label `Exploratory / Provisional annotation`.
"""
    (OUT / "BroadClass_method_selection.md").write_text(method, encoding="utf-8")

    provenance = pd.DataFrame(
        [
            {
                "input": str(path),
                "size_bytes": path.stat().st_size,
                "sha256": sha256(path),
            }
            for path in required
        ]
    )
    provenance.to_csv(OUT / "Input_provenance.tsv", sep="\t", index=False)
    payload = {
        "reference_cells": len(ref_cells),
        "reference_states": 18,
        "reference_genes": len(genes),
        "query_g10_u20_assigned": int(passed.sum()),
        "query_samples": int(qc.loc[passed, "sample_id"].nunique()),
        "common_unique_gene_symbols": len(common),
        "candidate_donor_columns": candidate_donor_cols,
        "joint_size_bytes": JOINT.stat().st_size,
        "joint_mtime_ns": JOINT.stat().st_mtime_ns,
    }
    (OUT / "preflight.json").write_text(
        json.dumps(payload, indent=2), encoding="utf-8"
    )
    print(json.dumps(payload, indent=2))


if __name__ == "__main__":
    main()
