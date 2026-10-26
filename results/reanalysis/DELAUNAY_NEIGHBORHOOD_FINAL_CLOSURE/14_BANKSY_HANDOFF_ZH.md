# BANKSY交接说明

## 可直接继承的对象

- 冻结Reference18连续权重和126,449个QC Bin50。
- patient、chip、disease、ROI、tissue component、真实坐标及UMI/nGenes。
- 357,959条经审核Delaunay边及degree/boundary QC。
- 统一PN1–PN8、confidence、OOD/NOT_EVALUABLE标志。
- 全部29,750个connected regions及固定的area/fragmentation定义。
- 9位患者经审核的ssDNA仿射配准。

## BANKSY应优先回答的问题

1. 在不依赖当前PN硬标签的情况下，空间平滑表征是否能减少碎片化且保留真实边界？
2. chip variation在空间表征中是否仍占主导？必须按患者留出与chip敏感性报告，不能用训练内稳定性代替验证。
3. BANKSY识别的区域是否与Delaunay ≥10 Bin50 regions在位置、组成和尺度上对应？
4. 当前仅505个大型region具备全部zone，BANKSY是否能定义更稳健、非重叠或层级化的boundary/surrounding结构？
5. 9例ssDNA支持的形态位置能否复现；其余12例必须先获得可靠配准才能做组织学结论。

## 建议的比较基线

- 以patient为评估单位；保持相同QC、Reference18、物理坐标和组织边界。
- 与Delaunay比较：患者覆盖、composition separation、同标签邻接、singleton/2–4比例、≥10区域Bin比例、UMI/nGenes/degree敏感性和chip依赖。
- 保留原始BANKSY标签，不根据疾病P值或预期免疫故事选参数。
- 不要求BANKSY复现PN1–PN8编号；使用组成与空间特征建立对应，并允许无可靠对应。

## 启动门槛

本次Delaunay收尾状态为 **DELAUNAY_ANALYSIS_CLOSURE_READY**。这表示Delaunay结果、限制和比较基线已整理完毕，并不表示BANKSY已获自动授权启动。本任务结束后停止。
