# 三类非B NMF factor外部微环境图谱（疾病盲法）

## 范围

本分析只复用冻结的逐芯片K=6 Neighborhood NMF连续usage与逐Bin50 Reference18比例。
未重新运行RCTD或NMF，未开展疾病比较，未分析B-associated factor，也未进行factor-factor接触分析。

## 主方法

- 主外环为100–200 μm；0–50 μm仅作内部正对照，50–100及200–300 μm为径向敏感性。
- 外环只在同一patient、ROI及tissue_component内汇总，避免跨组织间隙错误连接。
- 每个外环至少20个Bin50并覆盖8个角扇区中的至少4个，才进入100–200 μm主分析。
- 先按150/200/250 μm空间块聚合，再在患者内计算factor usage与外环细胞状态比例的相关。
- 主结果为200 μm块的患者级技术调整Spearman rho；调整nUMI、nGene、原50 μm有效邻居数和外环Bin数。
- 21位患者是复现单位；疾病标签仅保留作溯源，未进入发现、筛选或检验。
- 正向关系解释为外环富集/连续性；负向关系仅表示组成比例相对降低，受总和为1的闭合约束影响，不能直接解释为细胞排斥。

## 数据完整性

- 100–200 μm主环有效中心中位覆盖率：0.976。
- 主分析患者×factor×state相关记录：1,134。
- BH校正后跨芯片可复现关系：9；部分复现关系：7。

## 证据分级

### 跨芯片可复现

- Ciliated → Ciliated epithelial: median adjusted rho=0.472, 患者同向比例=1.00, 3个芯片同向数=3, P=9.537e-07, BH q=1.717e-05.
- Ciliated → Blood endothelial: median adjusted rho=-0.121, 患者同向比例=0.90, 3个芯片同向数=3, P=0.0001612, BH q=0.002176.
- Ciliated → Fibroblast: median adjusted rho=-0.176, 患者同向比例=0.86, 3个芯片同向数=3, P=0.000721, BH q=0.006489.
- Ciliated → NK: median adjusted rho=0.107, 患者同向比例=0.76, 3个芯片同向数=3, P=0.003753, BH q=0.02252.
- Fibroblast → Fibroblast: median adjusted rho=0.393, 患者同向比例=1.00, 3个芯片同向数=3, P=9.537e-07, BH q=1.717e-05.
- Fibroblast → Perivascular/mesothelial stromal: median adjusted rho=-0.137, 患者同向比例=0.86, 3个芯片同向数=3, P=0.0004263, BH q=0.004604.
- Fibroblast → Lymphatic endothelial: median adjusted rho=-0.165, 患者同向比例=0.81, 3个芯片同向数=3, P=0.001859, BH q=0.01434.
- Fibroblast → Ciliated epithelial: median adjusted rho=-0.102, 患者同向比例=0.81, 3个芯片同向数=3, P=0.002857, BH q=0.01929.
- Lymphatic → Lymphatic endothelial: median adjusted rho=0.475, 患者同向比例=1.00, 3个芯片同向数=3, P=9.537e-07, BH q=1.717e-05.

### 部分复现

- Ciliated → Inflammatory monocyte/macrophage: median adjusted rho=-0.079, 患者同向比例=0.76, 3个芯片同向数=3, P=0.01135, BH q=0.06127.
- Ciliated → Secretory/other airway epithelial: median adjusted rho=0.060, 患者同向比例=0.76, 3个芯片同向数=3, P=0.01755, BH q=0.08614.
- Ciliated → Alveolar/resident macrophage: median adjusted rho=-0.053, 患者同向比例=0.62, 3个芯片同向数=3, P=0.03192, BH q=0.1326.
- Fibroblast → Secretory/other airway epithelial: median adjusted rho=-0.094, 患者同向比例=0.67, 3个芯片同向数=3, P=0.03192, BH q=0.1326.
- Fibroblast → AT2/alveolar transitional: median adjusted rho=0.133, 患者同向比例=0.76, 3个芯片同向数=3, P=0.03844, BH q=0.1384.
- Lymphatic → Blood endothelial: median adjusted rho=-0.135, 患者同向比例=0.71, 3个芯片同向数=3, P=0.03506, BH q=0.1352.
- Lymphatic → Fibroblast: median adjusted rho=-0.099, 患者同向比例=0.71, 3个芯片同向数=3, P=0.04208, BH q=0.142.

## 解释边界

这些结果描述的是factor连续usage较高处外围Reference18组成的重复关系，不等同于细胞直接接触、配体受体通讯或疾病效应。
空间分块降低了重叠邻域造成的伪重复，但不是严格的组织形状保持空间零模型；因此未通过跨患者/跨芯片标准的关系保留为unresolved。
