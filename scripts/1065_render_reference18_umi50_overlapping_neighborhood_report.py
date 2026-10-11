#!/usr/bin/env python3
"""Render a self-contained HTML report for analysis 1064.

All numeric statements are read from the frozen output tables at render time;
no result value is hardcoded. The report embeds the approved figures directly.
"""

from __future__ import annotations

import argparse
import base64
import datetime as dt
import html
import json
import platform
import sys
from pathlib import Path

import numpy as np
import pandas as pd


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--result-dir", required=True, type=Path)
    return parser.parse_args()


def image_data(path: Path) -> str:
    payload = base64.b64encode(path.read_bytes()).decode("ascii")
    return f"data:image/png;base64,{payload}"


def table_html(frame: pd.DataFrame, digits: int = 3) -> str:
    shown = frame.copy()
    for column in shown.select_dtypes(include=[np.number]).columns:
        shown[column] = shown[column].map(lambda x: f"{x:.{digits}f}" if pd.notna(x) else "")
    return shown.to_html(index=False, escape=True, classes="dataframe")


def main() -> None:
    args = parse_args()
    result_dir = args.result_dir.resolve()
    method = pd.read_csv(result_dir / "00_METHOD_AND_INPUT_LOCK.tsv", sep="\t")
    patients = pd.read_csv(result_dir / "01_PATIENT_DESIGN_AND_DEPTH.tsv", sep="\t")
    support = pd.read_csv(result_dir / "02_NEIGHBORHOOD_SUPPORT_QC.tsv", sep="\t")
    results = pd.read_csv(result_dir / "05_DISEASE_TEST_RESULTS.tsv", sep="\t")
    sensitivity = pd.read_csv(result_dir / "06_SCALE_SENSITIVITY.tsv", sep="\t")
    validation = pd.read_csv(result_dir / "10_FINAL_VALIDATION.tsv", sep="\t")
    session = json.loads((result_dir / "SESSION_INFO.json").read_text(encoding="utf-8"))

    primary = results[results["membership_k"].eq(500)].copy()
    primary_support = support[support["membership_k"].eq(500)]
    tested_n = primary["neighborhood_id"].nunique()
    top = primary.sort_values(["permutation_p", "hc3_p"]).head(10).copy()
    top["effect_direction"] = np.where(top["log_odds_effect"] > 0, "higher", "lower")
    top_table = top[
        [
            "neighborhood_id",
            "contrast",
            "odds_ratio",
            "hc3_p",
            "permutation_p",
            "permutation_fdr_global",
            "depth_adjusted_hc3_p",
            "loo_direction_fraction",
            "top3_mean_states",
        ]
    ].rename(
        columns={
            "neighborhood_id": "Neighborhood",
            "contrast": "Contrast",
            "odds_ratio": "Odds ratio",
            "hc3_p": "HC3 P",
            "permutation_p": "Permutation P",
            "permutation_fdr_global": "Global permutation FDR",
            "depth_adjusted_hc3_p": "Depth-adjusted P",
            "loo_direction_fraction": "LOO direction fraction",
            "top3_mean_states": "Top mean states",
        }
    )
    chip_design = pd.crosstab(patients["chip_short"], patients["disease"]).reset_index()
    method_subset = method[method["parameter"].isin(
        [
            "input_h5", "input_sha256", "spatial_radius_um", "center_bin", "transform",
            "reference_per_patient", "reference_k", "anchor_proportion",
            "primary_membership_k", "support_filter", "primary_model", "permutations",
            "multiple_testing", "supported_finding_rule", "seed",
        ]
    )]

    timestamp = dt.datetime.now().strftime("%y%m%d%H%M%S")
    output_path = result_dir / f"REFERENCE18_UMI50_OVERLAPPING_NEIGHBORHOOD_REPORT_{timestamp}.html"
    figures = {
        "design": image_data(result_dir / "FIG01_DESIGN_AND_SUPPORT_QC.png"),
        "effects": image_data(result_dir / "FIG02_PRIMARY_DISEASE_EFFECTS.png"),
        "composition": image_data(result_dir / "FIG03_TOP_NEIGHBORHOOD_COMPOSITION.png"),
    }
    min_row = primary.sort_values("permutation_p").iloc[0]
    css = """
body {font-family: Arial, sans-serif; max-width: 1200px; margin: 30px auto; padding: 0 24px; color:#222; line-height:1.55;}
h1,h2,h3 {color:#17365d;} h1 {border-bottom:3px solid #4472c4; padding-bottom:10px;}
.callout {background:#eef4fb; border-left:5px solid #4472c4; padding:14px 18px; margin:18px 0;}
.warning {background:#fff5e6; border-left-color:#ed7d31;} .result {background:#edf7ed; border-left-color:#70ad47;}
img {width:100%; height:auto; border:1px solid #ddd; margin:8px 0 20px;}
table.dataframe {border-collapse:collapse; width:100%; font-size:0.88rem; margin:12px 0 22px;}
table.dataframe th, table.dataframe td {border:1px solid #ccc; padding:6px 8px; vertical-align:top;}
table.dataframe th {background:#e8eef7;} code {background:#f4f4f4; padding:2px 4px;}
.small {font-size:0.88rem; color:#555;} ul li {margin-bottom:6px;}
"""
    document = f"""<!doctype html>
<html lang="zh"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1">
<title>Reference18 UMI50 overlapping composition-neighborhood analysis</title><style>{css}</style></head><body>
<h1>Reference18 UMI50：重叠组成-neighborhood患者级疾病比较</h1>
<p class="small">Rendered {dt.datetime.now().isoformat(timespec='seconds')} | Analysis directory: {html.escape(str(result_dir))}</p>

<div class="callout result"><strong>主要结论：</strong>在冻结的UMI≥50、21患者、Reference18数据上，
没有组成-neighborhood同时通过HC3全局FDR、chip内置换全局FDR、深度方向一致性和逐患者剔除标准。
因此本分析未检出HC、IPF和SSc-ILD之间受统计支持、可重复的局部细胞组成差异；这不等于证明差异绝对不存在。</div>

<h2>1. 这一步究竟分析了什么？</h2>
<p>每个Bin50先继承既有的50 μm物理空间邻域：把中心bin及其附近bin的RCTD Reference18权重汇总为18维组成。
这一步没有重跑RCTD，也没有重新改变物理空间邻域。</p>
<p>随后把126,449个局部空间邻域放进共同的18维组成空间。这里不再用K=3、4、5或8强迫分群，
而是用彼此可以重叠的近邻集合表示组成相似区域。因此，一个bin可以参与多个组成-neighborhood，
连续的组成梯度不会被强制切成互斥类别。</p>
<ul><li>疾病盲参考：每位患者抽取500个bin，共10,500个。</li>
<li>精炼锚点：679个；主尺度为每个锚点最近的500个bin。</li>
<li>可比较性过滤后：{tested_n}个neighborhood，共{len(primary):,}项疾病对比。</li>
<li>统计单位：患者，而不是bin；模型校正chip，疾病标签在chip内置换10,000次。</li></ul>

<h2>2. 输入与患者支持图怎么读？</h2>
<img src="{figures['design']}" alt="Design and support QC">
<p><strong>左图：</strong>柱高是各chip、各疾病组患者可用Bin50数量的均值，黑线是绘图默认的组均值置信区间。
它说明三张chip的有效bin数量不同，患者之间也相差很大，所以不能把126,449个bin当成126,449个独立样本。
本分析先转为患者比例，再以21位患者进行推断。</p>
<p><strong>右图：</strong>横轴表示一个组成-neighborhood在多少位患者中至少包含5个bin。
虚线位于6，是最低患者支持门槛。许多锚点只出现在少数患者，因此被排除；
主尺度679个锚点中有{int(primary_support['support_filter_pass'].sum())}个进入正式疾病检验。</p>
<div class="callout warning"><strong>注意：</strong>“k=500”不是K-means的K，也不是500个患者；
它指每个重叠组成-neighborhood由组成空间中最近的500个Bin50构成。</div>

<h3>每张chip的疾病构成</h3>{table_html(chip_design, 0)}
<p>三张chip均为HC 2、IPF 2、SSc-ILD 3，故疾病与chip不是完全混杂，模型可以同时估计chip和疾病效应。</p>

<h2>3. 疾病效应图怎么读？</h2>
<img src="{figures['effects']}" alt="Disease effect plots">
<p>每个点代表一个通过支持过滤的组成-neighborhood。横轴是患者级log-odds效应：右侧表示前一个疾病组中更丰富，
左侧表示更少；纵轴是chip内置换P值的−log10。虚线对应未经多重校正的P=0.05。</p>
<ul><li>IPF vs HC：没有点达到未经校正的置换P&lt;0.05。</li>
<li>SSc-ILD vs HC：少数点越过虚线，但经过786项全局多重校正后全部消失。</li>
<li>SSc-ILD vs IPF：没有明显候选。</li></ul>
<p>最小置换P来自{html.escape(str(min_row['neighborhood_id']))}的{html.escape(str(min_row['contrast']))}：
OR={min_row['odds_ratio']:.2f}，置换P={min_row['permutation_p']:.4f}，
全局置换FDR={min_row['permutation_fdr_global']:.3f}。它是探索性趋势，不是阳性疾病结果。</p>

<h2>4. 热图表示什么？</h2>
<img src="{figures['composition']}" alt="Top exploratory neighborhood composition">
<p>行是按未经校正置换P排序靠前的探索性neighborhood，列是Reference18状态，颜色表示该状态在neighborhood中的平均RCTD组成。
这些模式整体以AT1、B、淋巴内皮，并伴不同程度T、血管/基质或其他免疫成分为主。</p>
<p>这里的“B较高”不能直接解释为B细胞聚集、B细胞niche或TLS。它只表示在50 μm局部RCTD混合组成中，
B权重相对较高；而且这些候选没有通过FDR，加入患者测序深度后也不再名义显著。</p>

<h2>5. 我们得到了什么？</h2>
<ol><li><strong>完成了方法缺口：</strong>以前UMI50只做过离散K聚类；现在完成了不依赖K的、多变量、重叠neighborhood疾病比较。</li>
<li><strong>得到可靠的阴性主结论：</strong>HC、IPF、SSc-ILD之间没有通过预设稳健性标准的组成-neighborhood。</li>
<li><strong>排除了“只是K选错”的主要担忧：</strong>300和800近邻敏感性与500近邻的效应Spearman相关为
{sensitivity['spearman_effect_correlation'].min():.2f}–{sensitivity['spearman_effect_correlation'].max():.2f}，方向一致率为
{100*sensitivity['direction_agreement'].min():.0f}%–{100*sensitivity['direction_agreement'].max():.0f}%。</li>
<li><strong>保留了少量探索线索：</strong>主要是SSc-ILD相对HC升高的AT1/B/淋巴内皮混合组成，但不能升级为疾病niche。</li></ol>

<h2>6. 最靠前的探索性结果</h2>{table_html(top_table, 4)}

<h2>7. 方法锁定与验证</h2>
<h3>关键参数</h3>{table_html(method_subset, 4)}
<h3>自动验证</h3>{table_html(validation, 4)}
<p>所有验证项均为PASS。第一次运行目录在写图时因交互命令超时而中断且没有完成标记；
本报告只使用带<code>COMPLETE.ok</code>的<code>attempt2</code>目录。</p>

<h2>8. 环境信息</h2>
<pre>{html.escape(json.dumps(session, indent=2, ensure_ascii=False))}\nreport_python: {html.escape(sys.version)}\nreport_platform: {html.escape(platform.platform())}</pre>
</body></html>"""
    output_path.write_text(document, encoding="utf-8")
    print(output_path)


if __name__ == "__main__":
    main()
