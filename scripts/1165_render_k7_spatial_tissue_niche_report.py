#!/usr/bin/env python3
"""Build a timestamped, self-contained HTML report for the frozen K=7 atlas.

Purpose: document the bespoke spatial characterization performed by scripts
1163 and 1164 without rerunning RCTD, NMF, edgeR, GSEA, or disease testing.

Render command:
  envs/stereo_seq_py/Scripts/python.exe scripts/1165_render_k7_spatial_tissue_niche_report.py \
    --analysis-dir results/reanalysis/K7_SPATIAL_TISSUE_NICHE_CHARACTERIZATION
"""

from __future__ import annotations

import argparse
import base64
import contextlib
import datetime as dt
import html
import io
import platform
from pathlib import Path

import matplotlib
import numpy as np
import pandas as pd
import seaborn as sns
import session_info2


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--analysis-dir", required=True, type=Path)
    return parser.parse_args()


def image_uri(path: Path) -> str:
    encoded = base64.b64encode(path.read_bytes()).decode("ascii")
    return f"data:image/png;base64,{encoded}"


def markdown_to_html(text: str) -> str:
    blocks: list[str] = []
    in_list = False
    for raw_line in text.splitlines():
        line = raw_line.strip()
        if line.startswith("# "):
            if in_list:
                blocks.append("</ul>")
                in_list = False
            blocks.append(f"<h1>{html.escape(line[2:])}</h1>")
        elif line.startswith("## "):
            if in_list:
                blocks.append("</ul>")
                in_list = False
            blocks.append(f"<h2>{html.escape(line[3:])}</h2>")
        elif line.startswith("- "):
            if not in_list:
                blocks.append("<ul>")
                in_list = True
            blocks.append(f"<li>{html.escape(line[2:])}</li>")
        elif line:
            if in_list:
                blocks.append("</ul>")
                in_list = False
            blocks.append(f"<p>{html.escape(line)}</p>")
    if in_list:
        blocks.append("</ul>")
    return "\n".join(blocks)


def main() -> None:
    args = parse_args()
    analysis_dir = args.analysis_dir.resolve()
    figures_dir = analysis_dir / "figures"
    reports_dir = analysis_dir / "reports"
    reports_dir.mkdir(parents=True, exist_ok=True)

    required = [
        analysis_dir / "01_INPUT_AND_FACTOR_MAP_AUDIT.tsv",
        analysis_dir / "03_FACTOR_HOTSPOT_AND_COMPONENT_METRICS.tsv",
        analysis_dir / "04_FACTOR_REFERENCE18_COMPOSITION.tsv",
        analysis_dir / "05_FACTOR_MARKER_AND_PROGRAM_SUPPORT.tsv",
        analysis_dir / "06_FACTOR_SPATIAL_RELATIONSHIPS.tsv",
        analysis_dir / "07_PATIENT_REPRODUCIBILITY.tsv",
        analysis_dir / "08_CANDIDATE_NICHE_CLASSIFICATION.tsv",
        analysis_dir / "09_METHODS_AND_FINDINGS_ZH.md",
        analysis_dir / "10_NEXT_STEPS_ZH.md",
        analysis_dir / "99_VALIDATION.tsv",
        analysis_dir / "02_PATIENT_FACTOR_SPATIAL_ATLAS.pdf",
    ]
    missing = [str(path) for path in required if not path.exists()]
    if missing:
        raise FileNotFoundError(f"Missing report inputs: {missing}")

    audit = pd.read_csv(required[0], sep="\t")
    reproducibility = pd.read_csv(required[5], sep="\t")
    classification = pd.read_csv(required[6], sep="\t")
    validation = pd.read_csv(required[9], sep="\t")
    if not validation["passed"].astype(str).str.lower().eq("true").all():
        raise RuntimeError("Analysis validation contains failed checks")

    summary_figures = sorted(figures_dir.glob("SUMMARY_*.png"))
    if len(summary_figures) != 6:
        raise RuntimeError(f"Expected 6 summary figures, found {len(summary_figures)}")

    tstamp = dt.datetime.now().strftime("%y%m%d%H%M%S")
    output = reports_dir / f"K7_SPATIAL_TISSUE_NICHE_CHARACTERIZATION_REPORT_{tstamp}.html"

    figure_html = "\n".join(
        f'<section class="figure"><h3>{html.escape(path.stem.replace("_", " ").title())}</h3>'
        f'<img src="{image_uri(path)}" alt="{html.escape(path.stem)}"></section>'
        for path in summary_figures
    )
    methods = markdown_to_html(required[7].read_text(encoding="utf-8"))
    next_steps = markdown_to_html(required[8].read_text(encoding="utf-8"))

    session_buffer = io.StringIO()
    with contextlib.redirect_stdout(session_buffer):
        session_info2.session_info(dependencies=True, gpu=True)
    session_text = session_buffer.getvalue().strip()
    if not session_text:
        session_text = (
            f"Python {platform.python_version()} | pandas {pd.__version__} | "
            f"numpy {np.__version__} | matplotlib {matplotlib.__version__} | "
            f"seaborn {sns.__version__}"
        )

    css = """
    body {font-family: Arial, sans-serif; color: #20242a; max-width: 1480px;
          margin: 0 auto; padding: 32px; line-height: 1.55;}
    h1, h2, h3 {color: #17365d;} h1 {border-bottom: 3px solid #5b9bd5; padding-bottom: 10px;}
    .notice {background: #fff3cd; border-left: 5px solid #d39e00; padding: 14px 18px;}
    .meta {color: #5a6573;} table {border-collapse: collapse; width: 100%; margin: 16px 0 28px;}
    th, td {border: 1px solid #c9d1d9; padding: 7px 9px; text-align: left; font-size: 13px;}
    th {background: #eaf2f8;} tr:nth-child(even) {background: #f7f9fb;}
    .figure {page-break-inside: avoid; margin: 30px 0;} .figure img {max-width: 100%; height: auto;}
    pre {white-space: pre-wrap; background: #f4f5f7; padding: 14px; overflow-wrap: anywhere;}
    code {font-family: Consolas, monospace;} a {color: #1f6fb2;}
    """
    document = f"""<!doctype html>
<html lang="zh-CN"><head><meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>K=7 Spatial Tissue Niche Characterization</title><style>{css}</style></head><body>
<h1>K=7 Spatial Tissue Niche Characterization</h1>
<p class="meta">Rendered {html.escape(tstamp)} from frozen K=7 inputs.</p>
<div class="notice"><strong>Interpretation boundary:</strong> composition-derived factors are not
automatically validated biological niches. This analysis is disease-blind and did not rerun RCTD, NMF,
edgeR, GSEA, or disease tests.</div>
<h2>Reproducibility summary</h2>
{reproducibility.to_html(index=False, border=0, float_format=lambda value: f"{value:.4g}")}
<h2>Candidate classification</h2>
{classification.to_html(index=False, border=0)}
<h2>Input and factor-map audit</h2>
{audit.to_html(index=False, border=0)}
<h2>Validation</h2>
{validation.to_html(index=False, border=0)}
<h2>Summary figures</h2>
{figure_html}
<h2>Methods and findings (Chinese)</h2>
{methods}
<h2>Next steps (Chinese)</h2>
{next_steps}
<h2>Reproducible commands</h2>
<pre><code>powershell.exe -NoProfile -ExecutionPolicy Bypass -File WORK/run_k7_spatial_tissue_niche_characterization.ps1
powershell.exe -NoProfile -ExecutionPolicy Bypass -File WORK/run_k7_spatial_tissue_niche_render.ps1
envs/stereo_seq_py/Scripts/python.exe scripts/1165_render_k7_spatial_tissue_niche_report.py --analysis-dir results/reanalysis/K7_SPATIAL_TISSUE_NICHE_CHARACTERIZATION</code></pre>
<h2>Session information</h2><pre>{html.escape(session_text)}</pre>
</body></html>"""
    output.write_text(document, encoding="utf-8")
    if output.stat().st_size < 100_000:
        raise RuntimeError("Rendered HTML is unexpectedly small")
    print(output)


if __name__ == "__main__":
    main()
