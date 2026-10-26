# BANKSY输入交接只读审核

状态：**PASS — READY_NOT_RUN**。本次仅审核并冻结输入指针，没有初始化或运行BANKSY。

## 可交接输入

- 原始计数：`results\reanalysis\bin50_input\BIN50_joint_reanalysis_tissue_raw_counts_common_genes.h5ad`，形状 `178,044 × 27,618`，稀疏矩阵类型 `_CSRDataset`，dtype `uint32`。
- 正式Bin50：必须使用冻结Reference18输入中的精确Bin ID，共 **126,449** 个；全部可在raw counts中一对一定位。患者、芯片、疾病及坐标完全一致。
- 坐标：正式分析坐标为 `x_um=x_raw/2`、`y_um=y_raw/2`；Bin50物理尺寸为25 × 25 μm。
- 临床信息：21位患者，K8/J2/L3各7位；SSc-ILD 9、IPF 6、HC 6。
- 基因注释：27,618个Ensembl ID按raw矩阵列顺序完全匹配；全部有symbol映射，无mapping conflict；有20行涉及重复symbol，后续不得直接按symbol静默合并。
- ssDNA：9/21位患者有已验证全局几何配准；其余12位标记 `SSDNA_NOT_AVAILABLE`。

## 关键门控

仅按 `nGenes≥15、UMI≥50` 从raw H5AD重新筛选会得到 **133,026** 个Bin，比正式冻结集合多 **6,577** 个。因此BANKSY必须读取 `16_BANKSY_FORMAL_BIN_INDEX.tsv.gz` 的精确ID，不允许以阈值筛选替代。

冻结Reference18表中的 `nUMI/nGene` 与当前common-gene raw H5AD派生的QC字段来自不同计数口径，不能相互覆盖；交接索引同时保留两套字段并明确其来源。

## 结论

输入在文件、ID、顺序、坐标及患者元数据层面已通过交接审核。ssDNA不是BANKSY建模的必需输入，只能在9位已验证患者中用于后续组织学审阅。BANKSY尚未运行。
