# J2 K=5 spatial architecture audit：主要发现

## 范围与方法

本轮复用冻结的J2 K=5 assignment，对58,772个neighborhood-evaluable Bin50进行对称空间审核。Queen邻接仅在同一患者、同一ROI内建立，要求|dx|和|dy|均不超过25 µm且每条物理边只计一次。
共纳入218614条唯一直接邻接边；五个niche和全部10个heterotypic niche pair均进入分析。observed/expected参照来自每位患者/ROI内999次标签置换。

## 空间连续性

按跨患者最大连通分量占比中位数排序，连续性最高的三个niche为：N2 (median largest-component fraction 0.169)；N4 (median largest-component fraction 0.109)；N1 (median largest-component fraction 0.100)。
连通分量结果描述K=5标签场在真实坐标中的连续或碎片化程度；它不等于对应细胞类型形成了纯组织区室。

## niche–niche直接邻接

共有9个heterotypic pair在至少6/7患者中保持同一方向。方向最一致的关系包括：N1-N3 [Negative in >=6/7 patients, median log2(O/E)=-2.066]；N4-N5 [Negative in >=6/7 patients, median log2(O/E)=-0.643]；N2-N4 [Negative in >=6/7 patients, median log2(O/E)=-0.571]；N2-N5 [Negative in >=6/7 patients, median log2(O/E)=-0.533]；N1-N2 [Negative in >=6/7 patients, median log2(O/E)=-0.348]。
高同类邻接不是正文的主要生物学证据。K=5标签来自重叠的50 µm邻域组成，置换参照又破坏内生空间自相关，因此这里仅支持空间连续性或界面偏好，不能视为独立验证、细胞互作或因果关系。

## 解释边界

- 未读取或使用frozen aggregate标签、成员或边界。
- disease仅用于患者图注，未参与空间定义，未做疾病显著性检验。
- 未做DEG/GSEA、pathway、CellChat或ligand-receptor分析。
- 未进行跨芯片consensus或将J2 niche编号映射到L3。
