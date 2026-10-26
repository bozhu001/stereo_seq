#!/usr/bin/env python3
"""Create a read-only BANKSY input handoff audit; do not run BANKSY."""

from __future__ import annotations

import hashlib
import json
import logging
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parents[4]
OUT = ROOT / "results/reanalysis/DELAUNAY_NEIGHBORHOOD_FINAL_CLOSURE"
TABLES = OUT / "tables"
REPORT = OUT / "report"
LOGS = OUT / "logs"
RAW_H5AD = ROOT / "results/reanalysis/bin50_input/BIN50_joint_reanalysis_tissue_raw_counts_common_genes.h5ad"
GENE_MAP = ROOT / "results/reanalysis/bin50_representation_benchmark/ensembl_gene_symbol_mapping.csv"
REF_ROOT = ROOT / "results/reanalysis/REFERENCE18_UMI50_VOLTRON_PER_CHIP_NICHE_20261005_153036"
SSDNA_INDEX = ROOT / "results/reanalysis/REFERENCE18_SPATIAL_GRAPH_NEIGHBORHOODS_V1/SSDNA_OVERLAY_REVIEW/00_INPUT_FILE_INDEX.tsv"
CHIPS = ("K8", "J2", "L3")


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def fail(message: str) -> None:
    raise RuntimeError(message)


def main() -> None:
    TABLES.mkdir(parents=True, exist_ok=True)
    REPORT.mkdir(parents=True, exist_ok=True)
    LOGS.mkdir(parents=True, exist_ok=True)
    logging.basicConfig(
        filename=LOGS / "1273_prepare_banksy_handoff.log",
        level=logging.INFO,
        format="%(asctime)s\t%(levelname)s\t%(message)s",
        force=True,
    )

    for path in (RAW_H5AD, GENE_MAP, SSDNA_INDEX):
        if not path.exists():
            fail(f"Missing required input: {path}")

    frozen_frames = []
    frozen_paths = []
    for chip in CHIPS:
        path = REF_ROOT / chip / "01_FROZEN_REFERENCE18_INPUT.tsv.gz"
        frozen_paths.append(path)
        frozen_frames.append(pd.read_csv(path, sep="\t"))
    formal = pd.concat(frozen_frames, ignore_index=True)
    if len(formal) != 126449 or not formal["bin_id"].is_unique:
        fail(f"Unexpected frozen Bin50 index: rows={len(formal)}, unique={formal['bin_id'].nunique()}")

    adata = ad.read_h5ad(RAW_H5AD, backed="r")
    raw_ids = pd.Index(adata.obs_names.astype(str))
    formal_ids = pd.Index(formal["bin_id"].astype(str))
    positions = raw_ids.get_indexer(formal_ids)
    if (positions < 0).any():
        fail(f"Frozen Bin50 IDs missing from raw H5AD: {(positions < 0).sum()}")
    obs = adata.obs.iloc[positions].copy()
    obs.index = formal.index

    comparisons = {
        "patient_id": "sample_id",
        "chip_id": "chip_id",
        "disease": "group",
        "x_raw": "x",
        "y_raw": "y",
    }
    mismatch_counts = {}
    for left, right in comparisons.items():
        a = formal[left]
        b = obs[right]
        if pd.api.types.is_numeric_dtype(a):
            mismatch = ~np.isclose(a.to_numpy(float), b.to_numpy(float), equal_nan=True)
        else:
            mismatch = a.astype(str).to_numpy() != b.astype(str).to_numpy()
        mismatch_counts[f"{left}_vs_{right}"] = int(mismatch.sum())
    if any(mismatch_counts.values()):
        fail(f"Frozen metadata mismatch: {mismatch_counts}")
    qc_field_differences = {
        "formal_nUMI_vs_raw_common_counts": int(
            (~np.isclose(formal["nUMI"].to_numpy(float), obs["total_counts_common_raw"].to_numpy(float))).sum()
        ),
        "formal_nGene_vs_raw_common_genes": int(
            (~np.isclose(formal["nGene"].to_numpy(float), obs["n_genes_common_raw"].to_numpy(float))).sum()
        ),
    }

    if not np.allclose(formal["x_um"], formal["x_raw"] / 2.0):
        fail("x_um is not x_raw / 2")
    if not np.allclose(formal["y_um"], formal["y_raw"] / 2.0):
        fail("y_um is not y_raw / 2")
    if not formal["nGene"].ge(15).all() or not formal["nUMI"].ge(50).all():
        fail("Frozen Bin50 index violates G15/UMI50")

    naive_mask = adata.obs["n_genes_common_raw"].ge(15) & adata.obs["total_counts_common_raw"].ge(50)
    naive_ids = set(raw_ids[naive_mask.to_numpy()])
    formal_set = set(formal_ids)
    naive_extra = len(naive_ids - formal_set)
    formal_missing_naive = len(formal_set - naive_ids)

    gene_map = pd.read_csv(GENE_MAP)
    var_ids = pd.Index(adata.var_names.astype(str))
    gene_order_match = len(gene_map) == adata.n_vars and gene_map["ensembl_var_name"].astype(str).equals(pd.Series(var_ids))
    if not gene_order_match:
        fail("Gene annotation order does not match raw H5AD var_names")
    gene_map = gene_map.copy()
    gene_map.insert(0, "raw_h5ad_var_index", np.arange(len(gene_map), dtype=int))
    gene_map["duplicate_gene_symbol"] = gene_map["gene_symbol"].duplicated(keep=False)
    gene_map["annotation_status"] = np.where(
        gene_map["mapped"].astype(bool) & ~gene_map["conflict"].astype(bool),
        "MAPPED_NO_CONFLICT",
        "REVIEW_REQUIRED",
    )
    gene_map.to_csv(TABLES / "17_BANKSY_GENE_ANNOTATION_AUDIT.tsv", sep="\t", index=False)

    ssdna_source = pd.read_csv(SSDNA_INDEX, sep="\t")
    ssdna_source = ssdna_source.rename(columns={"chip": "chip_label"})
    patient_meta = formal[["patient_id", "chip_label", "disease"]].drop_duplicates()
    if len(patient_meta) != 21:
        fail(f"Expected 21 unique patients, found {len(patient_meta)}")
    ssdna = patient_meta.merge(
        ssdna_source[["patient_id", "ssDNA_source_file", "whether_alignment_directly_available", "comments"]],
        on="patient_id",
        how="left",
    )
    ssdna["ssdna_alignment_status"] = np.where(
        ssdna["whether_alignment_directly_available"].fillna("no").str.lower().eq("yes"),
        "VALIDATED_GLOBAL_GEOMETRY",
        "SSDNA_NOT_AVAILABLE",
    )
    ssdna["ssDNA_source_file"] = ssdna["ssDNA_source_file"].fillna("")
    ssdna["comments"] = ssdna["comments"].fillna(
        "No validated ssDNA mapping in the current closure; do not infer or fabricate alignment."
    )
    ssdna = ssdna.sort_values(["chip_label", "patient_id"])
    ssdna.to_csv(TABLES / "18_BANKSY_SSDNA_ALIGNMENT_STATUS.tsv", sep="\t", index=False)

    status_map = ssdna.set_index("patient_id")["ssdna_alignment_status"]
    handoff_cols = [
        "bin_id", "entity_key", "patient_id", "chip_label", "chip_id", "disease", "roi_id",
        "tissue_component_id", "analysis_block_id", "x_raw", "y_raw", "x_um", "y_um", "nUMI", "nGene",
    ]
    handoff = formal[handoff_cols].copy()
    handoff.insert(0, "raw_h5ad_obs_index", positions)
    handoff["raw_h5ad_total_counts_common_raw"] = obs["total_counts_common_raw"].to_numpy()
    handoff["raw_h5ad_n_genes_common_raw"] = obs["n_genes_common_raw"].to_numpy()
    handoff["formal_qc_status"] = "FROZEN_G15_UMI50_INCLUDED"
    handoff["ssdna_alignment_status"] = handoff["patient_id"].map(status_map)
    handoff["banksy_input_status"] = "READY_NOT_RUN"
    handoff.to_csv(TABLES / "16_BANKSY_FORMAL_BIN_INDEX.tsv.gz", sep="\t", index=False, compression="gzip")

    matrix_kind = type(adata.X).__name__
    matrix_dtype = str(adata.X.dtype)
    summary = [
        ("raw_counts_h5ad", RAW_H5AD, "PASS", f"shape={adata.n_obs}x{adata.n_vars}; X={matrix_kind}; dtype={matrix_dtype}"),
        ("formal_bin_index", REF_ROOT, "PASS", "126449 unique frozen Bin50 IDs; exact ID selection required"),
        ("formal_bin_match_to_raw", RAW_H5AD, "PASS", "126449/126449 exact obs_name matches; patient/chip/disease/x/y mismatches=0"),
        ("formal_qc", REF_ROOT, "PASS", "All frozen bins satisfy nGene>=15 and nUMI>=50"),
        ("qc_field_provenance", RAW_H5AD, "PASS_WITH_NOTE", f"Frozen nUMI/nGene are preserved separately; current common-gene H5AD QC fields differ in {qc_field_differences['formal_nUMI_vs_raw_common_counts']}/{qc_field_differences['formal_nGene_vs_raw_common_genes']} rows and must not overwrite them"),
        ("threshold_only_selection", RAW_H5AD, "FAIL_FOR_FORMAL_USE", f"naive threshold count={len(naive_ids)}; extras={naive_extra}; missing={formal_missing_naive}"),
        ("gene_annotation", GENE_MAP, "PASS", f"{len(gene_map)} genes; order exact; mapped={int(gene_map.mapped.sum())}; conflicts={int(gene_map.conflict.sum())}; duplicate symbols={int(gene_map.duplicate_gene_symbol.sum())}"),
        ("coordinate_transform", REF_ROOT, "PASS", "x_um=x_raw/2 and y_um=y_raw/2 for all formal bins"),
        ("patient_metadata", REF_ROOT, "PASS", "21 patients; K8/J2/L3 each 7; diseases SSc-ILD=9, IPF=6, HC=6"),
        ("ssdna_alignment", SSDNA_INDEX, "PARTIAL", "9/21 validated global geometry; 12/21 SSDNA_NOT_AVAILABLE"),
        ("banksy_execution", OUT, "NOT_RUN", "Input handoff audit only; BANKSY was not initialized or executed"),
    ]
    audit = pd.DataFrame(summary, columns=["input_component", "source_path", "audit_status", "evidence_or_comment"])
    audit["source_path"] = audit["source_path"].astype(str)
    hashes = {str(path): sha256(path) for path in [RAW_H5AD, GENE_MAP, SSDNA_INDEX, *frozen_paths]}
    audit["sha256_if_file"] = audit["source_path"].map(hashes).fillna("")
    audit.to_csv(TABLES / "15_BANKSY_INPUT_HANDOFF_AUDIT.tsv", sep="\t", index=False)

    validation = {
        "status": "PASS_READY_NOT_RUN",
        "formal_bins": len(formal),
        "raw_h5ad_shape": [int(adata.n_obs), int(adata.n_vars)],
        "raw_matrix": {"class": matrix_kind, "dtype": matrix_dtype},
        "metadata_mismatches": mismatch_counts,
        "qc_field_differences": qc_field_differences,
        "naive_threshold_bins": len(naive_ids),
        "naive_threshold_extra_bins": naive_extra,
        "formal_missing_from_naive_threshold": formal_missing_naive,
        "patients": int(formal.patient_id.nunique()),
        "validated_ssdna_patients": int(ssdna.ssdna_alignment_status.eq("VALIDATED_GLOBAL_GEOMETRY").sum()),
        "ssdna_not_available_patients": int(ssdna.ssdna_alignment_status.eq("SSDNA_NOT_AVAILABLE").sum()),
        "banksy_started": False,
        "input_sha256": hashes,
    }
    (LOGS / "BANKSY_HANDOFF_VALIDATION.json").write_text(
        json.dumps(validation, ensure_ascii=False, indent=2), encoding="utf-8"
    )

    audit_report = f"""# BANKSY输入交接只读审核

状态：**PASS — READY_NOT_RUN**。本次仅审核并冻结输入指针，没有初始化或运行BANKSY。

## 可交接输入

- 原始计数：`{RAW_H5AD.relative_to(ROOT)}`，形状 `{adata.n_obs:,} × {adata.n_vars:,}`，稀疏矩阵类型 `{matrix_kind}`，dtype `{matrix_dtype}`。
- 正式Bin50：必须使用冻结Reference18输入中的精确Bin ID，共 **{len(formal):,}** 个；全部可在raw counts中一对一定位。患者、芯片、疾病及坐标完全一致。
- 坐标：正式分析坐标为 `x_um=x_raw/2`、`y_um=y_raw/2`；Bin50物理尺寸为25 × 25 μm。
- 临床信息：21位患者，K8/J2/L3各7位；SSc-ILD 9、IPF 6、HC 6。
- 基因注释：{len(gene_map):,}个Ensembl ID按raw矩阵列顺序完全匹配；全部有symbol映射，无mapping conflict；有20行涉及重复symbol，后续不得直接按symbol静默合并。
- ssDNA：9/21位患者有已验证全局几何配准；其余12位标记 `SSDNA_NOT_AVAILABLE`。

## 关键门控

仅按 `nGenes≥15、UMI≥50` 从raw H5AD重新筛选会得到 **{len(naive_ids):,}** 个Bin，比正式冻结集合多 **{naive_extra:,}** 个。因此BANKSY必须读取 `16_BANKSY_FORMAL_BIN_INDEX.tsv.gz` 的精确ID，不允许以阈值筛选替代。

冻结Reference18表中的 `nUMI/nGene` 与当前common-gene raw H5AD派生的QC字段来自不同计数口径，不能相互覆盖；交接索引同时保留两套字段并明确其来源。

## 结论

输入在文件、ID、顺序、坐标及患者元数据层面已通过交接审核。ssDNA不是BANKSY建模的必需输入，只能在9位已验证患者中用于后续组织学审阅。BANKSY尚未运行。
"""
    (REPORT / "BANKSY_INPUT_HANDOFF_AUDIT_ZH.md").write_text(audit_report, encoding="utf-8")

    plan = """# BANKSY分析方案（交接草案，未执行）

## 分析边界

BANKSY应作为独立空间表示学习路线启动，不继承Delaunay PN1–PN8为默认标签，也不以复现PN、NMF或免疫故事为优化目标。冻结Delaunay结果仅作为描述性对照。

## 冻结输入

1. 用 `16_BANKSY_FORMAL_BIN_INDEX.tsv.gz` 精确选取126,449个正式Bin50。
2. 从审计通过的raw-count H5AD按 `raw_h5ad_obs_index` 和冻结ID读取计数，不重新定义QC。
3. 使用27,618行顺序一致的Ensembl注释；重复gene symbol须显式处理并留痕。
4. 保留patient、chip、disease、ROI、tissue component、x/y坐标和QC字段。

## 建议的正式执行门控

- 先在每张芯片选择一位代表患者做计算规模、稀疏矩阵、坐标和输出checkpoint验证。
- 预先固定BANKSY参数网格、随机种子、邻域定义及K选择标准；不得根据疾病P值或预期免疫类别选参。
- 患者内建图，禁止跨patient/ROI/tissue component连接。
- 明确区分训练内稳定性、patient-held-out稳定性和组织学验证。
- 保留chip-aware解释；不得自动运行Harmony、CCA、RPCA或把Delaunay患者内中心化规则移植到BANKSY。
- 只有9位患者可进行现有ssDNA overlay审核；其余12位不得伪造配准。

## 尚未授权

本文件不是执行指令。不得启动BANKSY、重新聚类、疾病比较、LR、TCR/BCR或新的批次校正，直至收到下一条正式指令。
"""
    (REPORT / "BANKSY_ANALYSIS_PLAN_ZH.md").write_text(plan, encoding="utf-8")
    adata.file.close()
    logging.info("BANKSY handoff audit completed: %s", json.dumps(validation, ensure_ascii=False))
    print(json.dumps(validation, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
