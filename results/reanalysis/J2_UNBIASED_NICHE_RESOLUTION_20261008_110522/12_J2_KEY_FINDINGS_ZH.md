# J2 unbiased neighborhood analysis：主要发现

## 范围与输入

J2包含7位患者和58,840个UMI≥50的Bin50；疾病构成为HC=2、IPF=2、SSc-ILD=3。官方50 µm radius图中58,772个Bin有至少一个邻居，68个空邻域Bin被标记为not evaluable并排除聚类分母。图包含325,776条边，无跨患者或跨ROI边。中心Bin未额外纳入邻域组成。

输入为冻结Reference18连续权重；邻域组成经VoltRon原生CLR处理。frozen aggregate未读取，疾病标签仅作为metadata保留，未参与邻域、聚类、K选择或命名。

## K分辨率结论

推荐主分辨率为K=5，欠分辨对照为K=3，高分辨率敏感性为K=7。K=5在20/20次seed中收敛，ARI中位数0.988、NMI中位数0.979，五个niche均覆盖7/7患者，最小niche占15.7%，患者内UMI eta-squared为0.040。

K=3的自动综合分数最高且技术稳定性最好，但合并了多个可解释组成轴，因此保留为欠分辨对照。K=7仍无<1%小簇且覆盖全部患者，但重复拟合ARI中位数降至0.797，仅作为高分辨率敏感性。K选择未使用疾病分离或与L3的相似程度。

## K=5组成与跨患者重复性

- **N1：alveolar-associated mixed niche**。核心正向状态：AT1；核心状态UMI稳健=TRUE；核心状态单患者驱动=FALSE。
- **N2：B-enriched mixed niche**。核心正向状态：B; Ciliated epithelial; Fibroblast; Plasma；核心状态UMI稳健=TRUE；核心状态单患者驱动=FALSE。
- **N3：airway-associated mixed niche**。核心正向状态：Basal-like epithelial; Perivascular/mesothelial stromal; Ciliated epithelial; T; AT2/alveolar transitional; Lymphatic endothelial; Mast/Basophil; NK; Dendritic cell; Myofibroblast; Secretory/other airway epithelial; Inflammatory monocyte/macrophage; Alveolar/resident macrophage；核心状态UMI稳健=TRUE；核心状态单患者驱动=FALSE。
- **N4：endothelial-associated mixed niche**。核心正向状态：Lymphatic endothelial；核心状态UMI稳健=TRUE；核心状态单患者驱动=FALSE。
- **N5：endothelial-associated mixed niche**。核心正向状态：Blood endothelial; Fibroblast; Perivascular/mesothelial stromal; T; Dendritic cell; Alveolar/resident macrophage; Inflammatory monocyte/macrophage；核心状态UMI稳健=TRUE；核心状态单患者驱动=FALSE。

最明确的跨患者轴包括N1的AT1富集、N2的B富集、N4的淋巴内皮富集及N5的血管内皮富集；这些核心状态均为7/7患者同方向、留一患者方向稳定，并在log1p(Bin UMI)调整及患者内UMI分位分层后保持方向。N2因此可保守称为B-enriched mixed niche，但不等于纯B细胞区域。

部分niche存在`any_state_single_patient_driven=TRUE`，涉及的是非核心或较弱状态；五个niche的`core_state_single_patient_driven`均为FALSE。该区分防止用任一边缘状态否定一个跨患者重复的核心组成轴。

## UMI敏感性与限制

K=5整体患者内UMI关联较低，但这不表示UMI影响被完全消除。UMI模型与分层仅用于敏感性检查；niche仍是50 µm邻域的Reference18推断组成模式，不是单细胞类型或真实细胞数量。

## 与L3的盲性描述

J2独立出现了B富集、肺泡、内皮、气道/基质及其他混合组成主题，这些主题在宽泛生物学层面与L3曾见主题相似。但本轮未进行正式跨芯片对齐，不能声明J2的任何Nx等同于L3的某个Ny。完成第三张芯片的独立发现后，才适合进入正式cross-chip alignment。

## 本轮停止范围

- frozen aggregate：未使用。
- disease label参与clustering：否。
- disease significance：未做。
- DEG/GSEA/pathway：未做。
- CellChat/ligand-receptor：未做。
- L3/J2 consensus matching：未做。

当前J2结果通过输入、邻域、K分辨率、跨患者方向和UMI敏感性审核，适合进入下一阶段；后续仍应保留K=3和K=7作为分辨率敏感性，而不能把K=5编号直接映射到其他芯片。
