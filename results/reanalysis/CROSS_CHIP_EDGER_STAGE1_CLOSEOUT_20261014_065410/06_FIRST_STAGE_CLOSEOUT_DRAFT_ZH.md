# 第一阶段 edgeR 跨芯片证据收口（分析草案）

## 范围与主终点

- 仅复用冻结的 K=6 factor-family 映射、family edgeR、单芯片 edgeR 和既有患者 LOO；未重跑 RCTD、NMF、pseudobulk 或 edgeR。
- 主终点固定为 `paired_high_minus_control`，即同一患者 factor-high 区域减 matched control 区域后的疾病差异。
- 联合 family edgeR（控制芯片）仍为主分析；芯片级 meta-analysis 仅用于复现与异质性敏感性评估。

## 新增方法

- 由单芯片 edgeR 的 1-df QL-F 统计量近似恢复 effect-scale SE：`abs(log2FC) / sqrt(QL-F)`。
- 计算 inverse-variance fixed-effect meta、DerSimonian-Laird random-effect meta、Cochran Q、I² 与 tau²。
- 因每个 family 只有 2–3 张芯片，random-effect 和 I² 均不作为单独淘汰标准。
- Meta FDR 在全部可估计 protein-coding family × comparison × gene 检验中统一校正；同时保留 family/contrast 内 FDR。

## 数据概览

- 既有候选信号：20。
- Joint global coding FDR < 0.05：7。
- Fixed meta global FDR < 0.05：12。
- 所有可估计芯片方向与 Joint 一致或效应中性：19/20。
- 全 meta 检验中 fixed global FDR < 0.05：23。

## 描述性证据等级

- Cross-chip supported: 8
- Heterogeneous or patient-sensitive: 5
- Unresolved / within-model only: 3
- Cross-chip robust: 2
- Directionally consistent: 2

证据等级不是新的筛除门槛。联合模型、单芯片结果、meta、异质性和患者 LOO 应并列阅读；B-associated 等患者敏感信号继续保留，但不得表述为三芯片独立显著。

## 关键限制

- 芯片级 SE 来自已保存 QL-F 的近似反推，并非重新拟合 edgeR 后直接导出的协方差矩阵。
- 仅 2–3 张芯片时，I²、tau² 和 random-effect 估计很不稳定。
- Meta 与 joint family 模型使用相同患者数据，是互补敏感性分析，不是独立验证队列。
- 候选表来自既有 edgeR 候选规则；完整 meta 结果另存，避免只报告被 Joint 预筛选的基因。
