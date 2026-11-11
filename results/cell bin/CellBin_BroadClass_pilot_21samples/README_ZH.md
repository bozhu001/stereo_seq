# CellBin Broad-class pilot（已在试运行门控停止）

本目录记录全部21例进入候选试运行后的方法审计。没有排除患者，也没有修改Joint V1或G10/U20标准。

## 当前结论

SingleR供体×Reference18状态聚合方案在Reference18内部供体分组自验证中表现良好，但对10,500个稀疏CellBin的默认pruning没有产生任何Unknown，并出现明显的参考—查询分数域偏移。因此未扩展到359,329个CellBin。

## 已生成

- `Reference18_to_BroadClass_mapping.tsv`：18个冻结状态到四大类的可追溯映射。
- `BroadClass_method_selection.md`：执行前方法选择、风险与门控规则。
- `Reference18_input_audit.tsv`：reference/counts/gene/ID审计。
- `BroadClass_pilot_gate_summary.tsv`：试运行门控结果。
- `Pilot_reference_CV_summary.tsv`、`Pilot_summary_by_sample.tsv`、`Pilot_prediction_composition.tsv`：试运行诊断。
- `PILOT_STOP_REPORT.md`：停止原因与可执行恢复路径。
- `pilot_figures/`：英文试运行诊断图。
- `work/`：可复核的试运行输入、预测和R环境信息；其中标签不得视为正式结果。

## 有意未生成

由于门控失败，以下正式输出不存在：`CellBin_BroadClass_predictions.parquet`、全队列样本汇总/不确定性/marker验证/可行性表、21张正式空间预测图和三芯片正式比较图。

## 输入保护

Joint V1运行前记录的大小为183428028 bytes、mtime_ns=1794388172328082700；本分析仅读取该文件。
