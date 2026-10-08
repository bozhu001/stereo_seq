# L3 K=8 niche内部空间架构：主要发现

## 范围与解释边界

本分析复用冻结K=8标签，覆盖39982个neighborhood-evaluable Bin50、7位患者和7个ROI。
Queen邻接仅在同一患者、同一ROI内建立；中心距离在x、y方向均不超过25 um，且每条物理边只计一次。
置换零模型在每个patient/ROI内部固定坐标和niche数量，随机打乱标签999次。
K8标签本身来自50 um邻域Reference18组成。因此，标签的局部相似不是独立验证；本轮新增信息是连通分量、真实Bin-to-Bin界面及跨患者接触方向的一致性。

## 1. 空间连续性

按患者最大连通分量占比的中位数排序，连续性最高的三个niche为：N7 (0.154); N4 (0.135); N5 (0.112)。
碎片化指数中位数最高的三个niche为：N3 (0.935); N6 (0.928); N1 (0.922)。
N3在7/7位患者均可评估；其最大连通分量占比中位数为0.065，每位患者的连通分量数范围为25-275。这说明N3可形成真实空间patch，同时也可能包含多个分离component，而不是单一连续区块。

## 2. 全部28个heterotypic pair的无偏比较

N1-N8均以完全相同的方法进入空间连续性、直接邻接、patient/ROI内置换及跨患者一致性分析；没有使用N3或既往B-cell问题筛选pair。
综合跨患者方向一致性、非单患者驱动、observed/expected效应量、两端niche连续性、Reference18可解释性及UMI稳健性后，排序最前的三个pair为：N4-N7 (7/7 same direction; median |log2 O/E|=4.607); N2-N7 (7/7 same direction; median |log2 O/E|=3.673); N7-N8 (7/7 same direction; median |log2 O/E|=3.117)。完整28个pair及全部判据见07A表。
该排序仅用于选择后续生物学表征对象，不改变任何空间结果，也不是新的显著性门槛。

## 3. N3 targeted follow-up

按未校正的直接接触边数，N3最常接触的三个niche为：N1 (7127 unique edges); N7 (5927 unique edges); N5 (5369 unique edges)。
校正各patient/ROI内niche丰度后，N3方向最一致的三个空间关系为：N4 (7/7 same direction; median log2 O/E=-2.813); N6 (7/7 same direction; median log2 O/E=-2.413); N2 (7/7 same direction; median log2 O/E=-2.350)。负值表示低于随机标签期望的直接界面，而不是“没有空间结构”。
N3-N2、N3-N4、N3-N7的结果分别为：N2 0/7正向、7/7负向; N4 0/7正向、7/7负向; N7 3/7正向、4/7负向。应结合06表中的效应量与置换P值逐项解释。
直接接触次数受niche丰度影响；正文判断优先使用patient/ROI内置换后的observed/expected，而不只看raw edge count。
该置换保持每个ROI内niche数量，但会破坏标签的空间自相关，也不保留相邻中心Bin之间50 um输入邻域重叠带来的相关性。因此，广泛的heterotypic负向富集可反映标签patch化与构建尺度，不能全部解释为特异性生物学排斥。

## 4. 跨患者重复性与单患者驱动

共有11/28个heterotypic pair在至少6/7患者中方向一致。
这些pair包括：N4-N7; N2-N7; N7-N8; N3-N4; N3-N6; N2-N3; N1-N4; N2-N5; N6-N8; N5-N6; N6-N7。
按单患者贡献>=50%规则，可能单患者驱动的pair为：N1-N7。
方向一致不等同于每位患者置换检验显著；两者在05表中分别保留。

## 5. 是否存在有组织的spatial architecture

K8标签场相对于可交换标签零模型呈现明显空间组织：多个niche形成连通component，部分heterotypic界面在多个患者中保持同向observed/expected偏移。但K8由重叠的50 um邻域组成构建，置换又破坏了这种内生空间自相关；因此这里只能确认有组织的niche label field，不能把全部偏移视为独立生物学验证、细胞互作或因果关系。
全局绝对效应排序优先的pair为：N4-N7; N2-N7; N7-N8；当前三者均为负向observed/expected，主要代表重复的空间分离/少接触。正向界面偏好应单独查看05与07A表，不能与负向分离混为同一故事。

## 6. 明确未做事项

- 未读取或使用frozen aggregate标签、成员或边界。
- disease仅作为患者图注，未参与空间定义，也未进行疾病显著性检验。
- 未做DEG、GSEA、pathway、CellChat、ligand-receptor或cell-cell communication。
- adjacency只表示直接Bin50空间接触或界面偏好，不等于communication。
- 同一Bin50可包含多个细胞；niche label不是单细胞类型。
