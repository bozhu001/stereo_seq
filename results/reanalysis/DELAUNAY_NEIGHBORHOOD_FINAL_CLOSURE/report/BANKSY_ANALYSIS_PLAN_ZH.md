# BANKSY分析方案（交接草案，未执行）

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
