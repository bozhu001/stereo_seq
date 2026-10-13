# 四个primary factor-family连通块统一敏感性核查

## 范围

本轮将已批准的B-associated核查规则不加选择地应用于B、Ciliated、Fibroblast和Lymphatic四个三芯片primary family。未重跑NMF/RCTD，未修改factor-family映射，未开展factor–factor分析。

## 固定方法

- 阈值：芯片内患者q85/q90/q95中位数，以及患者自身q90。
- 指标：最大块绝对Bin数、最大块/全部有效Bin、最大块/全部high Bin。
- 邻接：患者内36 μm queen-like邻接；不跨患者或坐标图组织分量。
- 推断：患者为单位；指标在family×chip内转换并标准化；疾病标签仅在芯片内置换100,000次。
- 完整报告IPF-HC、SSc-ILD-HC和SSc-ILD-IPF，探索性BH同时按全部144项和每个contrast内48项计算。

## patient-specific q90下IPF vs HC集中度

| Family | Effect (within-chip SD) | P | Global q | 3-chip same direction | LOO direction fraction |
|---|---:|---:|---:|---|---:|
| B-associated | 1.429 | 0.0074 | 0.446 | True | 1.000 |
| Ciliated-associated | -0.267 | 0.6590 | 0.869 | False | 0.833 |
| Fibroblast-associated | 0.695 | 0.2377 | 0.760 | True | 1.000 |
| Lymphatic-associated | 1.053 | 0.0668 | 0.656 | True | 1.000 |

## 统一多重校正与保留原则

全部144项探索性敏感性检验中，通过global BH q<0.05的结果为0项。因此没有任何疾病关联被升级为确认性结论。patient-specific q90的相对集中度可作为后续患者级空间变量，但其疾病方向按`13_FINAL_VARIABLE_RECOMMENDATION.tsv`分为nominal cross-chip exploratory、cross-chip directional、secondary descriptive或unresolved。

## 统一阈值覆盖异常

原chip-median patient q90规则下，高usage覆盖超过50%的family-patient组合共有6个，完整列于`14_EXTREME_BASELINE_COVERAGE.tsv`。这类异常不限于B-associated，因此芯片统一阈值不宜作为跨患者面积比较的唯一规则；患者自身q90应保留为主要拓扑敏感性。

## 解释原则

`largest_component_over_high`衡量top集合内部的集中程度；`largest_component_over_valid`衡量患者内相对组织占据；`largest_component_bin_n`才对应绝对物理块大小。患者自身q90固定约10% high，因此前两个比例近似成比例，不能视为两项独立证据。

原1104结果原样保存在`09_ORIGINAL_1104_PRIMARY_COMPONENT_RESULTS.tsv`，本轮敏感性结果不替代其P值和FDR。只有阈值敏感、三芯片方向、LOO和三种分母共同支持时，才建议冻结为后续空间变量；其余保留为描述性或unresolved。

原1104基准行数：12。完整144项结果见`04_ALL_DISEASE_COMPONENT_SENSITIVITY.tsv`。

最终推荐行数：12（4 families × 3 contrasts）；所有`confirmed_disease_association`均为False。
