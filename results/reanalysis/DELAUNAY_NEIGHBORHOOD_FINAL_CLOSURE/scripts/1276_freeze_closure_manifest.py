#!/usr/bin/env python3
"""Freeze the approved closure with parameters, file hashes and a READY marker."""

from __future__ import annotations

import hashlib
import json
from datetime import datetime, timezone
from pathlib import Path

import pandas as pd


ROOT = Path(__file__).resolve().parents[4]
OUT = ROOT / "results/reanalysis/DELAUNAY_NEIGHBORHOOD_FINAL_CLOSURE"
MANIFEST = OUT / "logs/FROZEN_ANALYSIS_MANIFEST.tsv"
PARAMETERS = OUT / "logs/FROZEN_PARAMETERS.json"
MARKER = OUT / "DELAUNAY_ANALYSIS_CLOSURE_READY.ok"
HTML_VALIDATION = OUT / "logs/HTML_REVIEW_VALIDATION.json"


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def main() -> None:
    html_reports = sorted(
        (OUT / "report").glob("DELAUNAY_NEIGHBORHOOD_FINAL_CLOSURE_REVIEW_*.html"),
        key=lambda path: path.stat().st_mtime,
    )
    if not html_reports:
        raise RuntimeError("No timestamped HTML review report found")
    html_path = html_reports[-1]
    html_text = html_path.read_text(encoding="utf-8")
    html_checks = {
        "path": str(html_path),
        "bytes": html_path.stat().st_size,
        "has_title": "Delaunay Neighborhood Final Closure Review" in html_text,
        "has_ready_status": "DELAUNAY_ANALYSIS_CLOSURE_READY" in html_text,
        "has_banksy_not_run": "READY_NOT_RUN" in html_text,
        "embedded_png_count": html_text.count("data:image/png;base64"),
        "has_closing_html": html_text.rstrip().endswith("</html>"),
    }
    html_checks["status"] = "PASS" if (
        html_checks["bytes"] > 1_000_000
        and html_checks["has_title"]
        and html_checks["has_ready_status"]
        and html_checks["has_banksy_not_run"]
        and html_checks["embedded_png_count"] == 8
        and html_checks["has_closing_html"]
    ) else "FAIL"
    HTML_VALIDATION.write_text(json.dumps(html_checks, ensure_ascii=False, indent=2), encoding="utf-8")
    if html_checks["status"] != "PASS":
        raise RuntimeError(f"HTML validation failed: {html_checks}")

    parameters = {
        "status": "DELAUNAY_ANALYSIS_CLOSURE_READY",
        "scope": "project-specific frozen analysis; not a default workflow and not applied automatically to BANKSY",
        "formal_qc": {"nGene_min": 15, "nUMI_min": 50, "formal_bin_count": 126449},
        "bin50_physical_size_um": [25, 25],
        "delaunay_max_edge_um": 35.36,
        "graph_boundaries": ["patient", "ROI", "tissue_component"],
        "frozen_outputs": ["PN labels", "connected regions", "interior/boundary/surroundings"],
        "region_strata_bins": ["1", "2-4", "5-9", ">=10"],
        "outside_rings_um": ["immediate outside", "0-50", "50-100"],
        "retained_limitations": [
            "chip-associated variation", "region fragmentation", "overlapping rings",
            "RNA-depth sensitivity", "no disease contrasts at BH-FDR < 0.05",
            "9/21 patients with validated ssDNA alignment",
        ],
        "prohibited_in_this_freeze": [
            "K optimization", "reclustering", "spatial co-occurrence", "new batch correction",
            "BANKSY execution", "RCTD", "NMF", "edgeR", "GSEA", "LR",
        ],
        "banksy_status": "READY_NOT_RUN",
        "timestamped_html_review": str(html_path.relative_to(OUT)),
        "random_seed": "preserved in source scripts where stochastic operations occurred; no new stochastic analysis in wrap-up",
        "generated_utc": datetime.now(timezone.utc).isoformat(),
    }
    PARAMETERS.write_text(json.dumps(parameters, ensure_ascii=False, indent=2), encoding="utf-8")

    excluded = {MANIFEST.resolve(), MARKER.resolve()}
    rows = []
    for path in sorted(p for p in OUT.rglob("*") if p.is_file()):
        if path.resolve() in excluded or "__pycache__" in path.parts:
            continue
        rows.append({
            "scope": "closure_output",
            "relative_path": str(path.relative_to(OUT)),
            "bytes": path.stat().st_size,
            "sha256": sha256(path),
        })

    input_paths = [
        ROOT / "results/reanalysis/bin50_input/BIN50_joint_reanalysis_tissue_raw_counts_common_genes.h5ad",
        ROOT / "results/reanalysis/bin50_representation_benchmark/ensembl_gene_symbol_mapping.csv",
        ROOT / "results/reanalysis/REFERENCE18_SPATIAL_GRAPH_NEIGHBORHOODS_V1/ALL21_EXTERNAL_PROJECTION/FORMAL_RUN_MODEL_RECOVERED_EXACT/02_ALL21_PROJECTION_LABELS.tsv.gz",
        ROOT / "results/reanalysis/REFERENCE18_SPATIAL_GRAPH_NEIGHBORHOODS_V1/ALL21_EXTERNAL_PROJECTION/FORMAL_RUN_MODEL_RECOVERED_EXACT/06_ALL_CONNECTED_COMPONENTS.tsv.gz",
    ]
    for chip in ("K8", "J2", "L3"):
        input_paths.append(ROOT / f"results/reanalysis/REFERENCE18_UMI50_VOLTRON_PER_CHIP_NICHE_20261005_153036/{chip}/01_FROZEN_REFERENCE18_INPUT.tsv.gz")
    for path in input_paths:
        if path.exists():
            rows.append({
                "scope": "external_frozen_input",
                "relative_path": str(path.relative_to(ROOT)),
                "bytes": path.stat().st_size,
                "sha256": sha256(path),
            })
        else:
            rows.append({
                "scope": "external_frozen_input_missing",
                "relative_path": str(path.relative_to(ROOT)),
                "bytes": 0,
                "sha256": "",
            })
    manifest = pd.DataFrame(rows)
    manifest.to_csv(MANIFEST, sep="\t", index=False)
    missing = int(manifest.scope.eq("external_frozen_input_missing").sum())
    if missing:
        raise RuntimeError(f"Cannot freeze: {missing} required external inputs are missing")
    MARKER.write_text(
        "DELAUNAY_ANALYSIS_CLOSURE_READY\n"
        f"frozen_utc={datetime.now(timezone.utc).isoformat()}\n"
        "banksy_status=READY_NOT_RUN\n"
        "scope=project_specific_not_default\n",
        encoding="utf-8",
    )
    print(json.dumps({"status": "DELAUNAY_ANALYSIS_CLOSURE_READY", "manifest_rows": len(manifest)}))


if __name__ == "__main__":
    main()
