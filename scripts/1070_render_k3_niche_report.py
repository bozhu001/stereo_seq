#!/usr/bin/env python3
"""Render the approved All21 K3 niche analysis as a self-contained HTML report.

This renderer is used because Quarto/nbconvert is unavailable in the project
environment. It reads only frozen tabular/figure outputs, embeds figures as
base64, and generates a timestamped HTML report plus provenance JSON.
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
    parser.add_argument("--atlas-dir", required=True, type=Path)
    parser.add_argument("--heatmap-dir", required=True, type=Path)
    parser.add_argument("--adjacency-dir", required=True, type=Path)
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


def main() -> None:
    args = parse_args()
    atlas = args.atlas_dir.resolve()
    heatmap = args.heatmap_dir.resolve()
    adjacency = args.adjacency_dir.resolve()
    output_dir = args.output_dir.resolve()
    output_dir.mkdir(parents=True, exist_ok=True)
    timestamp = dt.datetime.now().strftime("%y%m%d%H%M%S")
    output_path = output_dir / f"REFERENCE18_UMI50_ALL21_K3_NICHE_REPORT_{timestamp}.html"

    required = {
        "method": atlas / "00_METHOD_AND_INPUT_LOCK.tsv",
        "k_audit": atlas / "03_K_SELECTION_AUDIT.tsv",
        "composition": atlas / "05_SELECTED_K_NICHE_COMPOSITION.tsv",
        "patient_abundance": atlas / "06_PATIENT_NICHE_ABUNDANCE.tsv",
        "disease_tests": atlas / "09_NICHE_DISEASE_TESTS.tsv",
        "atlas_validation": atlas / "15_FINAL_VALIDATION.tsv",
        "celltype_heatmap": heatmap / "figures" / "FIG01_K3_NICHE_X_CELLTYPE_HEATMAP.png",
        "disease_distribution": heatmap / "figures" / "FIG02_PATIENT_NICHE_FRACTION_BY_DISEASE.png",
        "disease_effect": heatmap / "figures" / "FIG03_DISEASE_EFFECT_HEATMAP.png",
        "spatial_atlas": atlas / "figures" / "FIG03_ALL21_SPATIAL_NICHE_ATLAS.png",
        "adjacency_table": adjacency / "03_DIFFERENTIAL_NICHE_ADJACENCY.tsv",
        "adjacency_distribution": adjacency / "figures" / "FIG01_PATIENT_K3_NICHE_ADJACENCY_BY_DISEASE.png",
        "adjacency_effect": adjacency / "figures" / "FIG02_DIFFERENTIAL_K3_NICHE_ADJACENCY_HEATMAP.png",
        "adjacency_validation": adjacency / "04_FINAL_VALIDATION.tsv",
    }
    missing = [str(path) for path in required.values() if not path.exists()]
    if missing:
        raise FileNotFoundError(f"Missing approved outputs: {missing}")

    method = pd.read_csv(required["method"], sep="\t")
    k_audit = pd.read_csv(required["k_audit"], sep="\t")
    composition = pd.read_csv(required["composition"], sep="\t")
    disease_tests = pd.read_csv(required["disease_tests"], sep="\t")
    atlas_validation = pd.read_csv(required["atlas_validation"], sep="\t")
    adjacency_tests = pd.read_csv(required["adjacency_table"], sep="\t")
    adjacency_validation = pd.read_csv(required["adjacency_validation"], sep="\t")

    k_display = k_audit[
        [
            "k",
            "repeat_ari_median",
            "cross_chip_ari_median",
            "silhouette_mean",
            "spatial_edge_agreement",
            "min_patient_coverage_n",
            "chip_cramer_v",
            "selected",
            "selection_status",
        ]
    ]
    composition_display = composition[
        ["niche", "top3_states", "dominant_weight"]
    ].rename(columns={"dominant_weight": "largest_mean_rctd_fraction"})
    disease_display = disease_tests[
        [
            "niche",
            "contrast",
            "log_odds_effect",
            "odds_ratio",
            "permutation_p",
            "permutation_fdr_global",
            "supported_finding",
        ]
    ]
    adjacency_display = adjacency_tests[
        [
            "niche_pair",
            "contrast",
            "difference_in_log2_adjacency_enrichment",
            "within_chip_permutation_p",
            "permutation_fdr",
            "supported_finding",
        ]
    ]
    validation = pd.concat(
        [
            atlas_validation.assign(section="atlas"),
            adjacency_validation.assign(section="adjacency"),
        ],
        ignore_index=True,
    )[["section", "check", "pass", "observed"]]

    css = """
    body {font-family: Arial, sans-serif; max-width: 1280px; margin: 32px auto; padding: 0 24px; color: #202124; line-height: 1.5;}
    h1, h2, h3 {color: #17365d;} h1 {border-bottom: 3px solid #4c78a8; padding-bottom: 10px;}
    .callout {background: #eef5fb; border-left: 5px solid #4c78a8; padding: 14px 18px; margin: 18px 0;}
    .warning {background: #fff4e5; border-left-color: #dd8452;}
    table {border-collapse: collapse; width: 100%; margin: 16px 0 26px; font-size: 0.9rem;}
    th, td {border: 1px solid #d8dde3; padding: 7px; text-align: left;} th {background: #e9eef5;}
    figure {margin: 24px 0 36px;} img {max-width: 100%; height: auto; border: 1px solid #d8dde3;}
    figcaption {font-size: 0.9rem; color: #555; margin-top: 8px;}
    code {background: #f3f4f6; padding: 2px 4px;} .small {font-size: 0.85rem; color: #555;}
    """
    body = [
        "<h1>All21 Reference18 UMI≥50: disease-blind K=3 spatial niche analysis</h1>",
        f'<p class="small">Rendered {dt.datetime.now().isoformat(timespec="seconds")}; report ID {timestamp}</p>',
        '<div class="callout"><b>Approved scope.</b> One common disease-blind K=3 niche atlas across 21 patients, followed by patient-level HC, IPF and SSc-ILD comparisons. K=4 was an accidental exploratory run and is excluded from this report.</div>',
        "<h2>1. Goal and analysis design</h2>",
        "<p>Construct local spatial neighborhoods from frozen Reference18 RCTD weights, describe each niche by its cell-type mixture, and test disease-associated changes in niche abundance and niche–niche adjacency without treating bins as independent biological replicates.</p>",
        table_html(method),
        "<h2>2. K audit and interpretation boundary</h2>",
        "<p>K=3 was retained as the most parsimonious and most repeat-stable candidate. Chip-centering reduced chip association, but held-chip transfer remained below the predeclared validation threshold. K=3 is therefore an exploratory common resolution, not a fully validated universal atlas.</p>",
        table_html(k_display),
        "<h2>3. Niche × cell-type composition</h2>",
        "<p>The upper heatmap reports uncorrected mean RCTD fractions, preserving biological composition. The lower heatmap reports cell-type-wise z scores across niches to emphasize relative enrichment used for annotation.</p>",
        figure_block(required["celltype_heatmap"], "Figure 1. K=3 niche × Reference18 cell-type composition."),
        table_html(composition_display),
        "<h2>4. Spatial localization across patients</h2>",
        "<p>Assignments are displayed separately for all 21 patients. Visual recurrence and spatial coherence are necessary descriptive checks, while patient-level statistics below provide inference.</p>",
        figure_block(required["spatial_atlas"], "Figure 2. Spatial distribution of CN1–CN3 in all 21 patients."),
        "<h2>5. Disease comparison of niche abundance</h2>",
        "<p>Each point represents one patient. Disease effects were estimated from corrected-logit niche fractions with chip as a covariate and evaluated using 10,000 within-chip disease-label permutations.</p>",
        figure_block(required["disease_distribution"], "Figure 3. Patient-level niche fractions by disease."),
        figure_block(required["disease_effect"], "Figure 4. Chip-adjusted disease effects and permutation FDR."),
        table_html(disease_display),
        "<h2>6. Disease comparison of niche–niche adjacency</h2>",
        "<p>Within each patient, adjacency enrichment is log2(observed/expected) for bins separated by ≤36 µm. Expected counts preserve patient-specific niche abundance. Disease inference again uses patients as replicates, chip adjustment and within-chip permutations.</p>",
        figure_block(required["adjacency_distribution"], "Figure 5. Patient-level niche–niche adjacency enrichment."),
        figure_block(required["adjacency_effect"], "Figure 6. Differential niche–niche adjacency by disease."),
        table_html(adjacency_display),
        "<h2>7. Conclusions</h2>",
        '<div class="callout"><ul><li>The K=3 heatmap provides an interpretable descriptive summary of three mixed tissue-composition niches.</li><li>No niche-abundance contrast passed permutation FDR &lt;0.05.</li><li>No inter-niche adjacency contrast passed permutation FDR &lt;0.05.</li><li>These null results do not establish biological equivalence; the atlas is exploratory and between-patient heterogeneity is substantial.</li><li>The next justified layer is patient × niche pseudobulk analysis of raw gene counts and predeclared fibrosis, epithelial injury, vascular-remodeling and immune gene programs.</li></ul></div>',
        "<h2>8. Validation and provenance</h2>",
        table_html(validation),
        "<h3>Input fingerprints</h3>",
        table_html(pd.DataFrame({"artifact": list(required), "path": [str(path) for path in required.values()], "sha256": [sha256(path) for path in required.values()]}), digits=0),
        f'<p class="small">Python {html.escape(sys.version)}<br>{html.escape(platform.platform())}<br>numpy {np.__version__}; pandas {pd.__version__}</p>',
    ]
    document = (
        "<!doctype html><html><head><meta charset='utf-8'>"
        "<meta name='viewport' content='width=device-width, initial-scale=1'>"
        f"<title>All21 K3 niche report {timestamp}</title><style>{css}</style></head>"
        f"<body>{''.join(body)}</body></html>"
    )
    output_path.write_text(document, encoding="utf-8")
    provenance = {
        "report": str(output_path),
        "timestamp": timestamp,
        "inputs": {name: {"path": str(path), "sha256": sha256(path)} for name, path in required.items()},
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
