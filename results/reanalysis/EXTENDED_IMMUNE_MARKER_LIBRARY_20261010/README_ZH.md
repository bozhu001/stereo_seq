# Extended Immune Marker Library 20261010

## 1. 本轮研究目的

建立肺纤维化单细胞 atlas 来源的扩展免疫 marker 库，并审核 21 位患者全部合格 Bin50 raw counts 中的 RNA 可检测性。

## 2. 输入数据及路径

- 主 atlas：`D:\Lung ST\codex agent\stereo_seq\stereo_seq\results\reanalysis\bin50_input\BI_PF_ILD_atlas_marker_gene_table_260709213053.xlsx`
- 21 患者 raw-count H5AD：`D:\Lung ST\codex agent\stereo_seq\stereo_seq\results\reanalysis\bin50_input\BIN50_joint_reanalysis_tissue_raw_counts_common_genes.h5ad`
- 完整路径、大小、修改时间和 SHA256：见 `input_manifest.tsv`。
- `260723131417.xlsx`：项目目录中未找到，本轮未使用，也未以相似文件替代。

## 3. 分析方法

- 原 Excel 的 worksheet dimension 错误写为 A1；脚本只在内存中 reset dimensions，恢复 4,500 条单元格记录，未修改原文件。
- Level 2 到 Level 1 通过 atlas 名称和 Top50 marker overlap 审核。
- marker 分为 lineage-core、extended lineage、state-associated、shared/functional。fold_change 定义因源文件未说明而标记 unknown。
- 对 H5AD 使用稀疏矩阵分块操作；unavailable 与 matrix 中存在但零检出的 undetected 分开报告。
- disease group 仅作为样本元数据保留，未参与 marker 选择或统计检验。

## 4. 输出文件清单

详见 `output_manifest.tsv`。主表位于 `tables/`，图位于 `figures/`，最终审计位于 `qc/`，中文结论位于 `reports/`。

## 5. 主要发现

见 `reports\11_EXTENDED_IMMUNE_MARKER_KEY_FINDINGS_ZH.md`。

## 6. 已完成与未完成

已完成 atlas 审计、免疫 marker 分层、Reference18 标签对应、21 患者可检测性、small/core/extended 覆盖比较、signature 候选与最终 QC。未进行任何聚类、niche scoring 或疾病显著性分析。

## 7. 下一步建议

先人工审核 core/extended/shared 分层，再决定是否在固定 50 µm neighborhood 上做探索性 signature scoring。

## 8. 重要限制

Atlas 只有正向 Top marker，无法完成完整非目标表达特异性估计。检测覆盖增加不等于识别准确率提高。Atlas 与 Reference18 的来源是否重叠未确认。

实际工作目录：`D:\Lung ST\codex agent\stereo_seq\stereo_seq`
本次结果目录：`D:\Lung ST\codex agent\stereo_seq\stereo_seq\results\reanalysis\EXTENDED_IMMUNE_MARKER_LIBRARY_20261010`
