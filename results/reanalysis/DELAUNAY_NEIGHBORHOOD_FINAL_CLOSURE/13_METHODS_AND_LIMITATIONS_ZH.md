# 方法与局限

## 冻结输入

- G15/UMI50 Bin50，21位患者，Reference18连续RCTD权重。
- 每位patient、ROI、tissue component独立构建distance-limited Delaunay图；最大边长35.36 μm。
- 使用既有统一PN1–PN8、classification confidence、OOD_UNASSIGNED和NOT_EVALUABLE。没有重新聚类或选K。

## Connected region定义

同一PN Bin50仅在存在真实Delaunay边时连通。全部29,750个component均保留。Singleton、2–4、5–9和≥10 Bin50分别报告；“connected region”不等于“aggregate”。面积定义为Bin50数×625 μm²，组织空洞不填补。

## Interior、boundary与surrounding

- Inner boundary：region内部且至少有一个图邻居位于region外。
- Interior：region中非inner-boundary的Bin50。
- Immediate outside：通过一条Delaunay边直接连接到region的外部Bin50。
- 0–50 μm与50–100 μm ring：从region boundary tile边缘计算最近近似欧氏距离；immediate outside从0–50 μm ring中排除。
- 所有zone限制在同一patient、ROI和tissue component；不同focal region的ring允许重叠但明确标记。
- ≥5 Bin50区域计算region/interior/boundary；≥10 Bin50区域再计算surrounding。zone少于3 Bin（background少于10 Bin）标记NOT_EVALUABLE。

## 汇总与统计单位

组成差异先在每个region计算，再在患者内取region中位数，最后跨患者描述。没有把Bin50或重叠ring作为独立患者。疾病比较仅使用患者级指标，并以chip作为协变量；BH-FDR覆盖预设全部比较。

## 组织学与NMF

仅对既有审计通过的9位患者复用芯片特异仿射变换。每个tile保持25×25 μm。不对其余12位构造伪配准。冻结NMF K=7仅加载既有usage进行描述性对照，没有重新拟合。

## 关键局限

1. PN与marker、NMF均源于同一RNA/Reference18体系，不能作为独立生物学验证。
2. chip-associated variation较强，统一PN的跨芯片可比性仅适合带chip信息的描述。
3. 99.7%的大型区域存在ring重叠，区域周围结果高度相关。
4. 只有9/21位患者具有可复用的已验证ssDNA配准。
5. 连通性受固定图和离散PN标签影响；碎片化不应通过事后平滑隐藏。
6. 本轮未做新的空间共现、LR、edgeR、GSEA、批次校正或疾病发现分析。
