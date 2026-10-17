#!/usr/bin/env python3
"""Render the approved K=7-on-ssDNA analysis as a standalone HTML report.

Typical command:
    envs/stereo_seq_py/Scripts/python.exe \
        scripts/1158_render_reference18_k7_ssdna_report.py \
        --analysis-dir results/reanalysis/REFERENCE18_K7_SSDNA_MAPS_261017220004
"""

from __future__ import annotations

import argparse
import base64
import datetime as dt
import hashlib
import html
import json
import platform
import sys
from pathlib import Path

import pandas as pd


ROOT = Path(__file__).resolve().parents[1]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--analysis-dir", required=True, type=Path)
    return parser.parse_args()


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def embedded_png(path: Path, alt: str) -> str:
    encoded = base64.b64encode(path.read_bytes()).decode("ascii")
    return (
        f'<img src="data:image/png;base64,{encoded}" '
        f'alt="{html.escape(alt)}" loading="lazy">'
    )


def table_html(frame: pd.DataFrame) -> str:
    return frame.to_html(index=False, border=0, classes="data-table", escape=True)


def main() -> None:
    args = parse_args()
    analysis_dir = args.analysis_dir.resolve()
    required = {
        "validation": analysis_dir / "03_VALIDATION.json",
        "samples": analysis_dir / "01_SAMPLE_MAP_MANIFEST.tsv",
        "coordinates": analysis_dir / "02_COORDINATE_VALIDATION.tsv",
        "inputs": analysis_dir / "00_INPUT_MANIFEST.tsv",
        "dominant_figure": analysis_dir
        / "FIG01_ALL21_K7_DOMINANT_FACTORS_ON_SSDNA.png",
        "b_figure": analysis_dir / "FIG02_ALL21_B_MIXED_USAGE_ON_SSDNA.png",
        "run_log": analysis_dir / "RUN_LOG.txt",
        "session": analysis_dir / "SESSION_INFO.txt",
        "complete": analysis_dir / "COMPLETE.ok",
    }
    missing = [str(path) for path in required.values() if not path.exists()]
    if missing:
        raise FileNotFoundError(f"Missing required approved-analysis files: {missing}")

    validation = json.loads(required["validation"].read_text(encoding="utf-8"))
    samples = pd.read_csv(required["samples"], sep="\t")
    coordinates = pd.read_csv(required["coordinates"], sep="\t")
    inputs = pd.read_csv(required["inputs"], sep="\t")
    disease_counts = (
        samples[["disease", "patient_id"]]
        .drop_duplicates()
        .groupby("disease", observed=True)
        .size()
        .rename("patient_n")
        .reset_index()
    )
    family_labels = validation["families"]
    stamp = dt.datetime.now().strftime("%y%m%d%H%M%S")
    report_dir = analysis_dir / "reports"
    report_dir.mkdir(exist_ok=True)
    output = report_dir / f"REFERENCE18_K7_SSDNA_MAPS_REPORT_{stamp}.html"

    dominant = embedded_png(
        required["dominant_figure"],
        "Dominant harmonized K=7 factor for all 21 samples on ssDNA",
    )
    b_mixed = embedded_png(
        required["b_figure"],
        "Continuous B-enriched mixed-immune K=7 usage for all 21 samples on ssDNA",
    )
    family_list = "".join(f"<li>{html.escape(label)}</li>" for label in family_labels)
    run_log = html.escape(required["run_log"].read_text(encoding="utf-8"))
    session = html.escape(required["session"].read_text(encoding="utf-8"))

    document = f"""<!doctype html>
<html lang="en">
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>Harmonized K=7 neighborhood-factor maps on ssDNA</title>
<style>
body {{ margin: 0; color: #17212b; background: #f5f7fa; font-family: Arial, sans-serif; line-height: 1.55; }}
main {{ max-width: 1500px; margin: 0 auto; padding: 34px 42px 70px; background: white; }}
h1 {{ margin-bottom: 0.25rem; font-size: 2rem; }}
h2 {{ margin-top: 2.4rem; border-bottom: 2px solid #d9e2ec; padding-bottom: 0.35rem; }}
.subtitle {{ color: #52606d; margin-top: 0; }}
.callout {{ background: #eef6ff; border-left: 5px solid #2878b5; padding: 14px 18px; margin: 22px 0; }}
.figure {{ margin: 24px 0 36px; }}
.figure img {{ width: 100%; height: auto; border: 1px solid #d9e2ec; }}
.caption {{ color: #52606d; font-size: 0.94rem; margin-top: 8px; }}
.data-table {{ border-collapse: collapse; width: 100%; margin: 12px 0 24px; font-size: 0.92rem; }}
.data-table th, .data-table td {{ border: 1px solid #d9e2ec; padding: 7px 9px; text-align: left; }}
.data-table th {{ background: #eef2f6; }}
pre {{ overflow-x: auto; background: #f3f4f6; padding: 14px; border-radius: 4px; }}
code {{ font-family: Consolas, monospace; }}
ul {{ columns: 2; }}
</style>
</head>
<body><main>
<h1>Harmonized K=7 neighborhood-factor maps on ssDNA</h1>
<p class="subtitle">Approved visualization-only analysis · rendered {dt.datetime.now().isoformat(timespec='seconds')}</p>

<div class="callout"><strong>Scope.</strong> Frozen disease-blind K=7 Neighborhood-NMF results were mapped onto the chip-specific ssDNA images. No RCTD/NMF model was refitted and no disease hypothesis test was performed.</div>

<h2>Goal</h2>
<p>Visualize the spatial distribution of the seven harmonized K=7 neighborhood-factor families across all 21 lung samples on their corresponding ssDNA backgrounds, and separately display continuous B-enriched/mixed-immune factor usage.</p>

<h2>Methods</h2>
<p>Each frozen modelled Bin50 was linked to its stored patient and registered spatial coordinates. The previously validated inverse chip-specific affine transform was applied to map registered Bin50 centers into raw ssDNA pixel space. The categorical map assigns each Bin50 to the factor with the highest relative usage. The continuous map displays the chip-matched B-enriched/mixed-immune factor on a fixed 0–1 scale.</p>
<p>The seven frozen families were:</p><ul>{family_list}</ul>

<h2>Coverage and validation</h2>
<p><strong>{validation['patient_n']}</strong> patients and <strong>{validation['mapped_bin50_n']:,}</strong> modelled Bin50 were mapped. All transformed Bin50 centers fell within their chip-specific ssDNA image bounds.</p>
<h3>Patients by disease</h3>{table_html(disease_counts)}
<h3>Coordinate validation</h3>{table_html(coordinates)}

<h2>Dominant K=7 factor on ssDNA</h2>
<div class="figure">{dominant}<p class="caption">Figure 1. Each Bin50 is colored by its dominant frozen K=7 factor. Rows show HC, IPF, and SSc-ILD; panels are ordered by chip and patient. The grayscale background is ssDNA.</p></div>

<h2>B-enriched/mixed-immune continuous usage</h2>
<div class="figure">{b_mixed}<p class="caption">Figure 2. Continuous B-enriched/mixed-immune neighborhood usage on a fixed 0–1 scale. This quantity represents neighborhood exposure and must not be interpreted as absolute B-cell abundance.</p></div>

<h2>Interpretation boundaries</h2>
<p>The maps are descriptive spatial visualizations of frozen factors. Differences in apparent color coverage are not formal disease comparisons. The B-enriched/mixed-immune factor is an immune-focused neighborhood factor supported by patient-paired expression analysis, but it is not a pure B-cell count. Cross-chip family labels come from the frozen harmonization and are not re-estimated here.</p>

<h2>Input provenance</h2>{table_html(inputs)}
<h2>Execution log</h2><pre><code>{run_log}</code></pre>
<h2>Session information</h2><pre><code>{session}</code></pre>
<p><strong>Report SHA-256:</strong> written to the adjacent <code>.sha256.txt</code> file after rendering.</p>
</main></body></html>"""
    output.write_text(document, encoding="utf-8")
    digest = sha256(output)
    output.with_suffix(".sha256.txt").write_text(f"{digest}  {output.name}\n", encoding="utf-8")
    print(output)


if __name__ == "__main__":
    main()
