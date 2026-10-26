#!/usr/bin/env python3
"""Independent integrity checks for the final Delaunay closure outputs."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd
import pymupdf
from PIL import Image
from pypdf import PdfReader


ROOT = Path(__file__).resolve().parents[4]
OUT = ROOT / "results/reanalysis/DELAUNAY_NEIGHBORHOOD_FINAL_CLOSURE"
TABLES = OUT / "tables"
FIGURES = OUT / "figures"


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def add(rows: list[dict], check: str, passed: bool, observed, expected) -> None:
    rows.append({"check": check, "status": "PASS" if passed else "FAIL", "observed": observed, "expected": expected})


def main() -> None:
    rows: list[dict] = []
    regions = pd.read_csv(TABLES / "04_ALL_LOCAL_CONNECTED_REGIONS.tsv.gz", sep="\t")
    internal = pd.read_csv(TABLES / "06_REGION_INTERNAL_COMPOSITION.tsv", sep="\t")
    surroundings = pd.read_csv(TABLES / "07_REGION_BOUNDARY_AND_SURROUNDINGS.tsv", sep="\t")
    zone_qc = pd.read_csv(TABLES / "07B_REGION_ZONE_EVALUABILITY.tsv", sep="\t")
    labels = pd.read_csv(
        ROOT / "results/reanalysis/REFERENCE18_SPATIAL_GRAPH_NEIGHBORHOODS_V1/ALL21_EXTERNAL_PROJECTION/FORMAL_RUN_MODEL_RECOVERED_EXACT/02_ALL21_PROJECTION_LABELS.tsv.gz",
        sep="\t", usecols=["patient_id", "analysis_pn", "predicted_pn", "component_id"],
    )
    add(rows, "region row count", len(regions) == 29750, len(regions), 29750)
    add(rows, "region IDs unique", regions.local_region_id.is_unique, regions.local_region_id.nunique(), len(regions))
    add(rows, "region Bin50 total", int(regions.bin_n.sum()) == 125471, int(regions.bin_n.sum()), 125471)
    predicted_count = int(labels.predicted_pn.isin([f"PN{i}" for i in range(1, 9)]).sum())
    add(rows, "formal predicted PN Bin50 total", predicted_count == 125471, predicted_count, 125471)
    component_count = int(labels.component_id.ne("NOT_EVALUABLE").sum())
    add(rows, "component membership reconciles", component_count == 125471, component_count, 125471)
    add(rows, "area equals bin count x 625", np.allclose(regions.approximate_area_um2, regions.bin_n * 625), float((regions.approximate_area_um2 - regions.bin_n * 625).abs().max()), 0)
    add(rows, "regions >=5", int(regions.bin_n.ge(5).sum()) == 6947, int(regions.bin_n.ge(5).sum()), 6947)
    add(rows, "regions >=10", int(regions.bin_n.ge(10).sum()) == 2282, int(regions.bin_n.ge(10).sum()), 2282)
    add(rows, "internal zone rows", len(internal) == 6947 * 3, len(internal), 6947 * 3)
    add(rows, "surrounding zone rows", len(surroundings) == 2282 * 4, len(surroundings), 2282 * 4)
    add(rows, "unique internal region-zone", not internal.duplicated(["local_region_id", "zone"]).any(), int(internal.duplicated(["local_region_id", "zone"]).sum()), 0)
    add(rows, "unique surrounding region-zone", not surroundings.duplicated(["local_region_id", "zone"]).any(), int(surroundings.duplicated(["local_region_id", "zone"]).sum()), 0)
    add(rows, "full-zone QC rows", len(zone_qc) == 2282, len(zone_qc), 2282)
    add(rows, "full-zone evaluable count", int(zone_qc.full_zone_analysis_status.eq("EVALUABLE").sum()) == 505, int(zone_qc.full_zone_analysis_status.eq("EVALUABLE").sum()), 505)
    patient_summary = pd.read_csv(TABLES / "10_CHIP_AWARE_PATIENT_SUMMARY.tsv", sep="\t")
    patient_rows = patient_summary[patient_summary.record_type.eq("PATIENT_PN_SUMMARY")]
    models = patient_summary[patient_summary.record_type.eq("CHIP_ADJUSTED_DISEASE_DESCRIPTION")]
    add(rows, "patient-PN summary coverage", len(patient_rows) == 21 * 8, len(patient_rows), 168)
    add(rows, "disease models patient-level", set(models.analysis_unit.dropna()) == {"patient"}, sorted(models.analysis_unit.dropna().unique()), ["patient"])
    add(rows, "disease FDR significant count", int(models.bh_fdr.lt(0.05).sum()) == 0, int(models.bh_fdr.lt(0.05).sum()), 0)
    selected = pd.read_csv(TABLES / "08B_REPRESENTATIVE_SSDNA_REGION_SELECTION.tsv", sep="\t")
    add(rows, "validated ssDNA representative patients", selected.patient_id.nunique() == 9, selected.patient_id.nunique(), 9)
    add(rows, "ssDNA chip coverage", set(selected.chip_label) == {"K8", "J2", "L3"}, sorted(selected.chip_label.unique()), ["J2", "K8", "L3"])
    required_tables = [f"{i:02d}_" for i in range(1, 12)]
    names = [p.name for p in TABLES.iterdir() if p.is_file()]
    for prefix in required_tables:
        add(rows, f"required table prefix {prefix}", any(name.startswith(prefix) for name in names), [name for name in names if name.startswith(prefix)], "at least one file")
    reports = [OUT / f"{i}_{name}.md" for i, name in [(12, "MAIN_FINDINGS_ZH"), (13, "METHODS_AND_LIMITATIONS_ZH"), (14, "BANKSY_HANDOFF_ZH")]]
    for report in reports:
        text = report.read_text(encoding="utf-8")
        add(rows, f"report exists and UTF-8 {report.name}", len(text) > 500 and "�" not in text, len(text), ">500 chars, no replacement character")
    pdfs = sorted(FIGURES.glob("FIGURE_*.pdf"))
    pngs = sorted(FIGURES.glob("FIGURE_*.png"))
    add(rows, "figure PDF count", len(pdfs) == 8, len(pdfs), 8)
    add(rows, "figure PNG count", len(pngs) == 8, len(pngs), 8)
    for path in pdfs:
        add(rows, f"readable single-page PDF {path.name}", len(PdfReader(path).pages) == 1, len(PdfReader(path).pages), 1)
    for path in pngs:
        with Image.open(path) as image:
            passed = min(image.size) >= 1000 and image.mode in {"RGB", "RGBA"}
            add(rows, f"readable high-resolution PNG {path.name}", passed, f"{image.size[0]}x{image.size[1]} {image.mode}", "min dimension >=1000; RGB/RGBA")
    atlas = FIGURES / "FINAL_DELAUNAY_ATLAS.pdf"
    add(rows, "atlas page count", len(PdfReader(atlas).pages) == 8, len(PdfReader(atlas).pages), 8)
    render_dir = OUT / "logs/pdf_render_qa"
    render_dir.mkdir(parents=True, exist_ok=True)
    document = pymupdf.open(atlas)
    for page_index, page in enumerate(document):
        pixmap = page.get_pixmap(dpi=120, alpha=False)
        pixmap.save(render_dir / f"atlas_page_{page_index + 1:02d}.png")
    document.close()
    rendered = sorted(render_dir.glob("atlas_page_*.png"))
    add(rows, "atlas rendered page count", len(rendered) == 8, len(rendered), 8)
    validation = pd.DataFrame(rows)
    validation.to_csv(OUT / "logs/FINAL_VALIDATION.tsv", sep="\t", index=False)
    failures = validation[validation.status.eq("FAIL")]
    manifest_rows = []
    for path in sorted(p for p in OUT.rglob("*") if p.is_file() and "__pycache__" not in str(p)):
        manifest_rows.append({"relative_path": str(path.relative_to(OUT)), "bytes": path.stat().st_size, "sha256": sha256(path)})
    pd.DataFrame(manifest_rows).to_csv(OUT / "logs/ALL_OUTPUT_FILES.tsv", sep="\t", index=False)
    result = {"status": "PASS" if failures.empty else "FAIL", "checks": len(validation), "failures": failures.to_dict("records"), "files": len(manifest_rows)}
    (OUT / "logs/FINAL_VALIDATION.json").write_text(json.dumps(result, indent=2, ensure_ascii=False), encoding="utf-8")
    print(json.dumps(result, indent=2, ensure_ascii=False))
    if not failures.empty:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
