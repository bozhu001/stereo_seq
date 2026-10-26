# Pilot K-means模型恢复结论

## 结论

**MODEL_RECOVERED_EXACT**

在最终Pilot脚本、冻结Reference18输入、原9位患者、原Delaunay构图规则、原特征转换、K=8及随机种子20261025下，仅重新拟合K=8一次。没有重选K，没有使用21位患者，也没有运行RCTD、NMF、edgeR、GSEA、LR或Marker分析。

## 可重复性结果

- 原Pilot可评估Bin50：50,972
- 新旧标签完全一致：50,972/50,972（100.000000000000%）
- 不一致Bin：0
- ARI：1.000000000000
- NMI：1.000000000000
- Delaunay边集合完全一致：True
- 标准化特征矩阵与原函数逐元素完全一致：True；最大绝对差=0
- 未进行标签置换或事后PN重命名；cluster 0–7直接映射PN1–PN8。

## 对21例External Projection的判断

现在可以严格支持下一阶段21例External Projection。恢复包已保存训练中心、scale、训练Bin索引、Reference18顺序和PN映射；下一阶段仍须明确区分9位in-sample与12位out-of-sample。
