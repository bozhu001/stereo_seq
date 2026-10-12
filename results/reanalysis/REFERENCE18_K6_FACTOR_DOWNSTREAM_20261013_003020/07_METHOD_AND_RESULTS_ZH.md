# Reference18 K=6 factor 原始表达验证与患者层面指标

## 范围

- 复用冻结的逐芯片 Neighborhood NMF K=6；未重新训练 NMF 或 RCTD。
- 共分析 3 张芯片、21 位患者、18 个 factor 和 126,199 个可评价 neighborhood。
- 未开展疾病差异检验，也未读取或修改冻结的 238 个 lymphoid aggregates。

## 原始表达验证

- marker 来自六数据集 Reference18 conservative reference 的 primary markers。
- 主验证将原始 RNA 按与 NMF 相同的中心排除 50 μm 邻域汇总；中心 Bin50 原始 RNA 保留为独立的局灶证据。
- 每位患者内对邻域和中心 marker program 分别残差化 UMI、n_genes 和有效邻居数。
- 另比较患者内 top-10% usage Bin 与技术匹配的低-usage Bin 原始计数 pseudobulk。
- 证据标签仅描述 RNA 支持程度，不作为删除 factor 的硬门槛。

## 患者层面指标

- 保存连续 mean/median/q90 usage、固定 chip 内患者平衡阈值以上的面积比例和 dominant-factor 比例。
- 保存高 usage 连通组分、最大组分比例、局部富集、原始及技术残差 Moran's I。
- 保存患者内 usage 与 UMI、n_genes、有效邻居数的相关性。

## 自动汇总

- Concordant strong raw-expression support：1 个 factor。
- Concordant raw-expression support：9 个 factor。
- Neighborhood raw-expression support：6 个 factor。
- Focal center-bin enrichment support：1 个 factor。
- Raw-expression support unresolved：1 个 factor。
- 既往综合审计建议保留为连续空间变量：9/18 个 factor。

详细的逐 factor 结论必须结合 `06_FINAL_FACTOR_ANALYSIS_TABLE.tsv`，不能仅凭单个阈值命名 biological niche。