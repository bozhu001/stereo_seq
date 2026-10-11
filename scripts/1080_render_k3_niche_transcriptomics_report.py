#!/usr/bin/env python3
"""Render the approved K=3 niche transcriptomics audit as standalone HTML.

Purpose
-------
Integrate the patient-by-niche pseudobulk edgeR analysis, fixed-program tests,
low-bin sensitivity analysis, Reference18 composition adjustment, chip/LOPO/QC
audits, and unsmoothed spatial checks into one reproducible report.

Inputs
------
The frozen tabular and PNG outputs in
``results/reanalysis/ALL21_K3_NICHE_PSEUDOBULK_DEG_20261011_201008``.

Outputs
-------
A timestamped, self-contained HTML report and provenance JSON in ``results/``.

Exact render command
--------------------
``envs/stereo_seq_py/Scripts/python.exe scripts/1080_render_k3_niche_transcriptomics_report.py --analysis-dir results/reanalysis/ALL21_K3_NICHE_PSEUDOBULK_DEG_20261011_201008 --output-dir results``
Redirect stdout/stderr to a timestamped file under ``logs/`` when executed.
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

import numpy as np
import pandas as pd


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--analysis-dir", required=True, type=Path)
    parser.add_argument("--output-dir", required=True, type=Path)
    return parser.parse_args()


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def image_data(path: Path) -> str:
    encoded = base64.b64encode(path.read_bytes()).decode("ascii")
    return f"data:image/png;base64,{encoded}"


def table_html(frame: pd.DataFrame, digits: int = 3) -> str:
    display = frame.copy()
    numeric = display.select_dtypes(include=np.number).columns
    display[numeric] = display[numeric].round(digits)
    return display.to_html(index=False, border=0, classes="dataframe", escape=True)


def figure_block(path: Path, caption: str) -> str:
    return (
        "<figure>"
        f'<img src="{image_data(path)}" alt="{html.escape(caption)}">'
        f"<figcaption>{html.escape(caption)}</figcaption>"
        "</figure>"
    )


def metric(frame: pd.DataFrame, name: str) -> float:
    value = frame.loc[frame["metric"] == name, "value"]
    if value.empty:
        raise KeyError(f"Missing metric: {name}")
    return float(value.iloc[0])


def main() -> None:
    args = parse_args()
    analysis_dir = args.analysis_dir.resolve()
    output_dir = args.output_dir.resolve()
    output_dir.mkdir(parents=True, exist_ok=True)
    timestamp = dt.datetime.now().strftime("%y%m%d%H%M%S")
    output_path = output_dir / (
        f"ALL21_K3_NICHE_TRANSCRIPTOMICS_REPORT_{timestamp}.html"
    )

    required = {
        "input_manifest": analysis_dir / "00_INPUT_MANIFEST.tsv",
        "sample_metadata": analysis_dir / "01_PSEUDOBULK_SAMPLE_METADATA.tsv",
        "preparation_validation": analysis_dir / "04_PREPARATION_VALIDATION.tsv",
        "program_tests": analysis_dir / "09_PROGRAM_DISEASE_TESTS.tsv",
        "deg_summary": analysis_dir / "12_DEG_SUMMARY.tsv",
        "model_validation": analysis_dir / "14_FINAL_VALIDATION.tsv",
        "robustness": analysis_dir / "17_DEG_ROBUSTNESS_SUMMARY.tsv",
        "composition_tests": analysis_dir / "22_RCTD_COMPOSITION_DISEASE_TESTS.tsv",
        "composition_adjusted": analysis_dir / "23_COMPOSITION_ADJUSTED_ROBUST_DEG.tsv",
        "composition_summary": analysis_dir / "27_COMPOSITION_VS_STATE_SUMMARY.tsv",
        "composition_validation": analysis_dir / "28_COMPOSITION_VS_STATE_VALIDATION.tsv",
        "technical_summary": analysis_dir / "36_TECHNICAL_PATIENT_CONSISTENCY_SUMMARY.tsv",
        "technical_validation": analysis_dir / "37_TECHNICAL_PATIENT_CONSISTENCY_VALIDATION.tsv",
        "evidence_tiers": analysis_dir / "38_TECHNICAL_EVIDENCE_TIER_ASSIGNMENT.tsv",
        "program_figure": analysis_dir / "figures/FIG01_PROGRAM_DISEASE_EFFECT_HEATMAP.png",
        "volcano_figure": analysis_dir / "figures/FIG02_NICHE_DISEASE_DEG_VOLCANO.png",
        "sensitivity_figure": analysis_dir / "figures/FIG03_PRIMARY_VS_LOW_BIN_EXCLUSION_DEG.png",
        "composition_figure": analysis_dir / "figures/FIG04_CILIATED_COMPOSITION_BY_DISEASE.png",
        "adjusted_figure": analysis_dir / "figures/FIG05_ORIGINAL_VS_COMPOSITION_ADJUSTED_DEG.png",
        "chip_figure": analysis_dir / "figures/FIG06_CHIP_SPECIFIC_GENE_EFFECT_HEATMAP.png",
        "attenuation_figure": analysis_dir / "figures/FIG07_TECHNICAL_COVARIATE_ATTENUATION_HEATMAP.png",
        "spatial_ms4a8": analysis_dir / "figures/FIG08_SPATIAL_MS4A8_CN3.png",
        "spatial_nwd1": analysis_dir / "figures/FIG08_SPATIAL_NWD1_CN3.png",
        "session_info": analysis_dir / "SESSION_INFO.txt",
        "composition_session_info": analysis_dir / "SESSION_INFO_COMPOSITION_VS_STATE.txt",
    }
    missing = [str(path) for path in required.values() if not path.exists()]
    if missing:
        raise FileNotFoundError(f"Missing approved outputs: {missing}")

    samples = pd.read_csv(required["sample_metadata"], sep="\t")
    programs = pd.read_csv(required["program_tests"], sep="\t")
    deg_summary = pd.read_csv(required["deg_summary"], sep="\t")
    robustness = pd.read_csv(required["robustness"], sep="\t")
    composition_tests = pd.read_csv(required["composition_tests"], sep="\t")
    composition_summary = pd.read_csv(required["composition_summary"], sep="\t")
    tiers = pd.read_csv(required["evidence_tiers"], sep="\t")

    tier_a = tiers.loc[
        tiers["evidence_tier"] == "Tier_A_cross_chip_and_technical_stable"
    ].copy()
    tier_a_display = tier_a[
        [
            "niche",
            "comparison",
            "gene_symbol",
            "same_direction_chips",
            "same_direction_lopo",
            "maximum_technical_attenuation",
            "evidence_tier",
        ]
    ]
    tier_counts = (
        tiers["evidence_tier"]
        .value_counts()
        .rename_axis("evidence_tier")
        .reset_index(name="gene_contrast_rows")
    )

    robust_total = int(robustness["robust_genes"].sum())
    program_fdr = int((programs["global_fdr"] < 0.05).sum())
    composition_fdr_column = next(
        col for col in composition_tests.columns if "fdr" in col.lower() and "global" in col.lower()
    )
    composition_fdr = int((composition_tests[composition_fdr_column] < 0.05).sum())
    patient_count = int(samples["patient_id"].nunique())
    niche_count = int(samples["niche"].nunique())
    tiny_libraries = int((samples["bin_n"] < 20).sum())
    adjusted_nominal = int(
        metric(composition_summary, "primary_ciliated_adjusted_nominal_p_lt_005")
    )
    lopo_stable = int(tiers["all_21_lopo_same_direction"].sum())
    cross_chip_stable = int(tiers["all_three_chips_same_direction"].sum())

    validations = []
    for key, label in [
        ("preparation_validation", "pseudobulk preparation"),
        ("model_validation", "edgeR models"),
        ("composition_validation", "composition audit"),
        ("technical_validation", "technical/patient audit"),
    ]:
        frame = pd.read_csv(required[key], sep="\t")
        frame.insert(0, "section", label)
        validations.append(frame)
    validation = pd.concat(validations, ignore_index=True)

    css = """
    body {font-family: Arial, sans-serif; max-width: 1280px; margin: 32px auto; padding: 0 24px; color: #202124; line-height: 1.55;}
    h1, h2, h3 {color: #17365d;} h1 {border-bottom: 3px solid #4c78a8; padding-bottom: 10px;}
    .callout {background: #eef5fb; border-left: 5px solid #4c78a8; padding: 14px 18px; margin: 18px 0;}
    .warning {background: #fff4e5; border-left-color: #dd8452;}
    .strong {background: #edf7ed; border-left-color: #4c9f70;}
    table {border-collapse: collapse; width: 100%; margin: 16px 0 26px; font-size: 0.88rem;}
    th, td {border: 1px solid #d8dde3; padding: 7px; text-align: left;} th {background: #e9eef5;}
    figure {margin: 24px 0 36px;} img {max-width: 100%; height: auto; border: 1px solid #d8dde3;}
    figcaption {font-size: 0.9rem; color: #555; margin-top: 8px;}
    code {background: #f3f4f6; padding: 2px 4px;} .small {font-size: 0.84rem; color: #555;}
    """
    body = [
        "<h1>All21 K=3 niche transcriptomics: patient-level pseudobulk and robustness audit</h1>",
        f'<p class="small">Rendered {dt.datetime.now().isoformat(timespec="seconds")}; report ID {timestamp}</p>',
        '<div class="callout"><b>Approved scope.</b> Frozen Reference18 UMI50 RCTD weights and frozen disease-blind CN1/CN2/CN3 assignments were reused. No RCTD, neighborhood construction, K selection, or niche assignment was rerun.</div>',
        "<h2>1. Question and statistical unit</h2>",
        f"<p>The analysis asks whether the same shared niche has a different transcriptional state among HC, IPF and SSc-ILD. Raw integer counts were aggregated into {patient_count} patients × {niche_count} niches ({len(samples)} libraries). Patients—not bins—are the biological replicates. Separate edgeR models used <code>~ chip + disease</code>. The prespecified sensitivity excluded {tiny_libraries} libraries containing fewer than 20 member bins.</p>",
        "<h2>2. Fixed biological programs</h2>",
        f'<div class="callout warning"><b>Result:</b> {program_fdr}/{len(programs)} program-by-niche-by-contrast tests passed global BH-FDR &lt; 0.05. Thus the eight fixed programs do not support a broad disease-associated program shift at this resolution.</div>',
        figure_block(required["program_figure"], "Figure 1. Patient-level effects for eight prespecified gene programs."),
        "<h2>3. Exploratory niche-specific differential expression</h2>",
        f"<p>The primary models identified gene-level candidates, followed by a predeclared low-bin sensitivity analysis. Requiring direction agreement, sensitivity FDR &lt; 0.05 and |log2FC| ≥ 0.5 retained {robust_total} gene-by-contrast rows. These are candidates for auditing, not 21 independently confirmed discoveries.</p>",
        figure_block(required["volcano_figure"], "Figure 2. Niche-specific edgeR results for all disease contrasts."),
        figure_block(required["sensitivity_figure"], "Figure 3. Primary effects versus the low-bin exclusion sensitivity effects."),
        table_html(robustness),
        "<h2>4. Is the signal only altered Reference18 composition?</h2>",
        f"<p>No Reference18 state × niche × disease comparison passed global FDR ({composition_fdr}/{len(composition_tests)}). All {adjusted_nominal}/{robust_total} audited gene-by-contrast rows retained nominal P &lt; 0.05 after separately adjusting for the ciliated-epithelial fraction. This argues against broad ciliated composition as the sole explanation, but cannot exclude finer subtype shifts or state/composition coupling below Reference18 resolution.</p>",
        figure_block(required["composition_figure"], "Figure 4. Ciliated-epithelial composition across disease groups within each niche."),
        figure_block(required["adjusted_figure"], "Figure 5. Original and Reference18 composition-adjusted disease effects."),
        "<h2>5. Patient, chip, and technical robustness</h2>",
        f"<p>{lopo_stable}/{robust_total} rows preserved direction in every leave-one-patient-out refit; {cross_chip_stable}/{robust_total} preserved direction in all three chips. Independent adjustments for niche bin count, median raw UMI, median detected genes, and median mitochondrial fraction were used as diagnostics. Because these covariates may also carry genuine injury biology, attenuation indicates uncertainty rather than proving artifact.</p>",
        figure_block(required["chip_figure"], "Figure 6. Disease effects estimated separately in J2, K8, and L3."),
        figure_block(required["attenuation_figure"], "Figure 7. Effect attenuation after separate technical-covariate adjustments."),
        table_html(tier_counts),
        "<h2>6. Highest-confidence subset</h2>",
        f'<div class="callout strong"><b>Tier A:</b> {len(tier_a)} gene-by-contrast rows ({tier_a["gene_symbol"].nunique()} unique genes) pass all prespecified direction and attenuation checks. Every Tier-A row occurs in CN3 and represents lower expression in disease than HC. This is a within-CN3 state signal—not a new niche, a cell-type abundance result, or proof of causality.</div>',
        table_html(tier_a_display),
        "<h2>7. Direct unsmoothed spatial checks</h2>",
        "<p>Maps use direct raw Bin50 detections without smoothing: grey marks evaluable tissue, blue marks the target niche, and magma-colored points mark bins with detected transcript. Sparse single-bin detection is expected; inference remains patient-level pseudobulk.</p>",
        figure_block(required["spatial_ms4a8"], "Figure 8. Unsmoothed MS4A8 counts in CN3 across patients."),
        figure_block(required["spatial_nwd1"], "Figure 9. Unsmoothed NWD1 counts in CN3 across patients."),
        "<h2>8. Interpretation and limits</h2>",
        '<div class="callout"><ul><li>The original unbiased K=3 niches remain descriptive mixed-tissue neighborhoods; this analysis tests transcriptional state within them.</li><li>The strongest reproducible signal is CN3 disease-versus-HC downregulation, represented by the Tier-A set.</li><li>CN2 SCGB1A1 and DRC1 are cross-chip and patient-direction stable but QC-attenuated, so they remain provisional.</li><li>No fixed program passed global FDR, and no broad Reference18 composition change explains the candidate set.</li><li>Post-selection diagnostics refine confidence but are not independent replication. Biological naming requires unbiased ranked-pathway analysis and orthogonal marker/anatomical validation.</li></ul></div>',
        "<h2>9. Validation and provenance</h2>",
        table_html(validation),
        "<h3>Input fingerprints</h3>",
        table_html(pd.DataFrame({"artifact": list(required), "path": [str(path) for path in required.values()], "sha256": [sha256(path) for path in required.values()]}), digits=0),
        f'<h3>Session information</h3><pre>{html.escape(required["session_info"].read_text(encoding="utf-8", errors="replace"))}</pre>',
        f'<pre>{html.escape(required["composition_session_info"].read_text(encoding="utf-8", errors="replace"))}</pre>',
        f'<p class="small">Report renderer: Python {html.escape(sys.version)}<br>{html.escape(platform.platform())}<br>numpy {np.__version__}; pandas {pd.__version__}</p>',
    ]
    document = (
        "<!doctype html><html><head><meta charset='utf-8'>"
        "<meta name='viewport' content='width=device-width, initial-scale=1'>"
        f"<title>K3 niche transcriptomics {timestamp}</title><style>{css}</style></head>"
        f"<body>{''.join(body)}</body></html>"
    )
    output_path.write_text(document, encoding="utf-8")
    provenance = {
        "report": str(output_path),
        "timestamp": timestamp,
        "command": (
            "envs/stereo_seq_py/Scripts/python.exe "
            "scripts/1080_render_k3_niche_transcriptomics_report.py "
            "--analysis-dir results/reanalysis/"
            "ALL21_K3_NICHE_PSEUDOBULK_DEG_20261011_201008 "
            "--output-dir results"
        ),
        "inputs": {
            name: {"path": str(path), "sha256": sha256(path)}
            for name, path in required.items()
        },
        "python": sys.version,
        "platform": platform.platform(),
        "numpy": np.__version__,
        "pandas": pd.__version__,
    }
    output_path.with_suffix(".provenance.json").write_text(
        json.dumps(provenance, indent=2), encoding="utf-8"
    )
    print(output_path)


if __name__ == "__main__":
    main()
