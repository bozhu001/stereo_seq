# Delaunay Spatial Neighborhood 最终收尾结论

## 最终状态

**DELAUNAY_ANALYSIS_CLOSURE_READY**

冻结Delaunay图、统一PN1–PN8标签、技术标志和既有K选择结果均未改变。本轮没有运行RCTD、NMF、K-means、edgeR、GSEA、LR、空间共现或新批次校正。

## 完成度

- 126,449个QC Bin50全部纳入版本审计；125,471个具有统一PN标签并进入连通区域分析。
- 冻结图包含357,959条无向边，按patient、ROI和tissue component限制，最大边长35.36 μm。
- 共保留29,750个connected regions：Singleton 13,148，2–4 Bin50 9,655，5–9 Bin50 4,665，≥10 Bin50 2,282。
- 505/2,282个≥10 Bin50区域具有全部预设zone；其余保留为PARTIAL_ZONE_COVERAGE。
- ssDNA配准仅对9/21位患者有已审核映射。本轮实际生成9个代表区域overlay；其余患者仅使用真实Bin50坐标，不伪造组织学配准。

## 可保留的PN组成描述

Reference18为混合Bin50的连续权重。以下为按区域Bin数加权的主要组成，不是单细胞计数，也不是独立病理验证：

- PN1: AT1 (0.192), B (0.191), Lymphatic endothelial (0.131), Ciliated epithelial (0.112)
- PN2: AT1 (0.216), B (0.148), Lymphatic endothelial (0.116), Fibroblast (0.094)
- PN3: AT1 (0.197), B (0.173), Fibroblast (0.144), Lymphatic endothelial (0.139)
- PN4: AT1 (0.223), B (0.130), Lymphatic endothelial (0.120), Fibroblast (0.092)
- PN5: AT1 (0.210), B (0.159), Lymphatic endothelial (0.119), Fibroblast (0.106)
- PN6: B (0.146), AT1 (0.144), Lymphatic endothelial (0.128), Fibroblast (0.084)
- PN7: AT1 (0.230), B (0.157), Blood endothelial (0.112), Fibroblast (0.086)
- PN8: AT1 (0.208), Lymphatic endothelial (0.153), B (0.136), Perivascular/mesothelial stromal (0.082)

全部PN在3张芯片均出现。具有≥10 Bin50区域的患者覆盖为：PN1 21/21; PN2 20/21; PN3 21/21; PN4 20/21; PN5 18/21; PN6 20/21; PN7 18/21; PN8 21/21。PN标签仍应保持中性命名。

## Local region与周围组成

区域严格由同一patient、ROI、tissue component、PN及真实Delaunay边定义。面积为Bin50数×625 μm²，不填补空洞。以下模式通过患者内region中位数后再跨患者汇总，属于描述性证据：

- PN2 | INTERIOR_MINUS_INNER_BOUNDARY | Inflammatory monocyte/macrophage: median patient-level difference=+0.117, direction consistency=100%, n=17 patients.
- PN6 | REGION_MINUS_RING_0_50 | AT1: median patient-level difference=-0.092, direction consistency=100%, n=20 patients.
- PN7 | REGION_MINUS_RING_0_50 | Lymphatic endothelial: median patient-level difference=-0.089, direction consistency=100%, n=18 patients.
- PN6 | REGION_MINUS_COMPONENT_BACKGROUND | AT1: median patient-level difference=-0.073, direction consistency=100%, n=20 patients.
- PN7 | REGION_MINUS_COMPONENT_BACKGROUND | Lymphatic endothelial: median patient-level difference=-0.065, direction consistency=100%, n=18 patients.
- PN7 | REGION_MINUS_RING_50_100 | Lymphatic endothelial: median patient-level difference=-0.065, direction consistency=100%, n=18 patients.
- PN6 | REGION_MINUS_RING_50_100 | AT1: median patient-level difference=-0.064, direction consistency=100%, n=20 patients.
- PN6 | INTERIOR_MINUS_INNER_BOUNDARY | AT1: median patient-level difference=-0.052, direction consistency=93%, n=15 patients.
- PN7 | REGION_MINUS_RING_0_50 | Blood endothelial: median patient-level difference=+0.050, direction consistency=100%, n=18 patients.
- PN4 | INTERIOR_MINUS_INNER_BOUNDARY | Alveolar/resident macrophage: median patient-level difference=+0.049, direction consistency=100%, n=15 patients.
- PN1 | REGION_MINUS_RING_0_50 | Ciliated epithelial: median patient-level difference=+0.043, direction consistency=100%, n=21 patients.
- PN8 | REGION_MINUS_RING_0_50 | Fibroblast: median patient-level difference=-0.043, direction consistency=100%, n=21 patients.

大型区域的外环高度重叠：99.7%的区域至少一个环带与其他区域环带共享Bin。因此region-level环带不能作为独立生物学重复；本报告不对环带做伪重复显著性检验。

## ssDNA、芯片与疾病解释

- 9位已有坐标审计的患者可以在ssDNA上定位connected regions。图中保持25×25 μm tile，不做平滑。
- 其余12位没有可靠ssDNA映射，因此不能声称完成21例组织学验证。
- 既有特征空间审计显示：疾病调整后chip partial R²约0.598；进一步调整UMI、nGenes和graph degree后仍约0.408。该差异可能同时包含技术批次、取材和患者组成异质性。
- 患者级、chip调整的探索性疾病比较共有96个预设contrast；FDR<0.05为0。因此本轮没有可正式宣称的疾病相关PN/region差异。

## 与冻结NMF K=7的关系

按dominant label比较，ARI为J2=0.018, K8=0.052, L3=0.011，NMI为J2=0.037, K8=0.036, L3=0.025。两种方法使用同一Reference18来源，低一致性说明分辨率与目标不同，不构成独立验证。Delaunay增加的是同标签真实邻接region及其局部边界/周围描述；NMF保留连续factor usage。

## 论文可用范围

- 可正式保留：图构建/QC、固定标签的可复现映射、全部connected-component清单、区域大小/碎片化、患者级覆盖、9例已验证ssDNA定位。
- 可作为空间描述：PN组成、≥10 Bin50区域、内部/边界/周围组成及Delaunay-NMF对应。
- 仅探索性：跨疾病比较、跨芯片PN比例、特定生物学niche命名、外环差异。
- 不支持：把所有component称为aggregate、把PN称为独立病理niche、把环带当独立重复、宣称21例均完成组织学验证。

## 收尾判断

Delaunay路线的预定收尾工作已完成，没有必须通过再次选K或重新聚类解决的缺口。剩余问题是独立组织学验证、未配准的12例ssDNA、批次/组成解耦及更稳健的空间表征，适合作为BANKSY阶段的明确交接问题。本任务到此停止，不自动启动BANKSY。
