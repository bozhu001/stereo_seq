# B-lineage覆盖与immune subniche映射：关键结果

## 分析边界

本轮只交叉映射冻结的29,645-bin Strategy C免疫候选池、K=3/4/5 assignments和最终2,669-bin B-lineage覆盖mask。未重新运行RCTD、邻域构建、聚类、K选择、marker筛选或疾病比较。患者是重复单位；Bin不作为独立生物学重复。

## 患者重复的B富集候选

- J2 K=3 I3: B覆盖=249/3111 (8.0%); 捕获该芯片B-mask=45.7%；患者同向=6/7；患者中位富集比=1.91；core佐证=45.8%；extended gain=74.7%；分类=recurrent_B_enriched_patient_imbalanced。
- J2 K=4 I1: B覆盖=209/2385 (8.8%); 捕获该芯片B-mask=38.3%；患者同向=6/7；患者中位富集比=1.99；core佐证=48.3%；extended gain=70.3%；分类=recurrent_B_enriched_patient_imbalanced。
- J2 K=5 I3: B覆盖=177/1721 (10.3%); 捕获该芯片B-mask=32.5%；患者同向=6/7；患者中位富集比=1.92；core佐证=52.5%；extended gain=65.5%；分类=recurrent_B_enriched_patient_imbalanced。
- L3 K=3 I1: B覆盖=1032/6344 (16.3%); 捕获该芯片B-mask=74.5%；患者同向=6/7；患者中位富集比=3.08；core佐证=73.5%；extended gain=58.9%；分类=recurrent_B_enriched_patient_imbalanced。
- L3 K=4 I4: B覆盖=947/4937 (19.2%); 捕获该芯片B-mask=68.3%；患者同向=6/7；患者中位富集比=3.60；core佐证=77.0%；extended gain=56.1%；分类=recurrent_B_enriched_patient_imbalanced。
- L3 K=5 I3: B覆盖=909/4577 (19.9%); 捕获该芯片B-mask=65.6%；患者同向=6/7；患者中位富集比=3.31；core佐证=77.4%；extended gain=55.6%；分类=recurrent_B_enriched_patient_imbalanced。
- Y40102K8 K=3 I1: B覆盖=734/2080 (35.3%); 捕获该芯片B-mask=99.5%；患者同向=5/7；患者中位富集比=4.90；core佐证=97.3%；extended gain=11.9%；分类=recurrent_B_enriched_patient_imbalanced。
- Y40102K8 K=5 I3: B覆盖=724/1797 (40.3%); 捕获该芯片B-mask=98.1%；患者同向=5/7；患者中位富集比=8.22；core佐证=97.4%；extended gain=10.9%；分类=recurrent_B_enriched_patient_imbalanced。

## K3→K4→K5主分支

- L3: K3 I1→K4 I4 (B=943, overlap=4910); K4 I4→K5 I3 (B=906, overlap=4525)。
- J2: K3 I3→K4 I1 (B=207, overlap=2314); K4 I1→K5 I3 (B=171, overlap=1630)。
- Y40102K8: K3 I1→K4 I4 (B=723, overlap=1803); K4 I4→K5 I3 (B=722, overlap=1751)。

## 解释限制

这里的分类是映射审计标签，不是最终immune subniche命名。B-lineage mask表示RNA覆盖而非纯B细胞判定；K3/K4/K5来自同一候选池且相互嵌套。最终全免疫注释仍需在T、NK、Plasma、Myeloid和Mast/Basophil覆盖扩展后联合冻结。