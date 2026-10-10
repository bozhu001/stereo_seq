#!/usr/bin/env python3
"""Build the approved Extended Immune Marker Library report source.

This is a reporting-only step. It reads the finalized tables, audits, figures,
and manifests produced by the approved analysis and creates a Quarto source.
It does not read the expression matrix or rerun scientific calculations.
"""

from __future__ import annotations

import argparse
import importlib.metadata
import os
import platform
from pathlib import Path

import pandas as pd
import pypdfium2 as pdfium


def markdown_table(frame: pd.DataFrame) -> str:
    """Return a compact GitHub-flavored Markdown table."""
    if frame.empty:
        return "_No rows._"
    columns = [str(column) for column in frame.columns]
    lines = [
        "| " + " | ".join(columns) + " |",
        "| " + " | ".join(["---"] * len(columns)) + " |",
    ]
    for row in frame.itertuples(index=False, name=None):
        values = [str(value).replace("|", "\\|").replace("\n", " ") for value in row]
        lines.append("| " + " | ".join(values) + " |")
    return "\n".join(lines)


def render_pdf_previews(figures: Path, assets: Path) -> list[tuple[str, Path]]:
    """Render page 1 of each final PDF for embedding in the HTML report."""
    assets.mkdir(parents=True, exist_ok=True)
    previews: list[tuple[str, Path]] = []
    for pdf_path in sorted(figures.glob("*.pdf")):
        document = pdfium.PdfDocument(str(pdf_path))
        page = document[0]
        image = page.render(scale=1.6).to_pil()
        preview_path = assets / f"{pdf_path.stem}.png"
        image.save(preview_path)
        previews.append((pdf_path.stem, preview_path))
        page.close()
        document.close()
    return previews


def package_version(name: str) -> str:
    try:
        return importlib.metadata.version(name)
    except importlib.metadata.PackageNotFoundError:
        return "not installed"


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output-dir", required=True, type=Path)
    args = parser.parse_args()

    output = args.output_dir.resolve()
    tables = output / "tables"
    qc = output / "qc"
    reports = output / "reports"
    figures = output / "figures"
    scripts = output / "scripts"
    assets = reports / "assets"

    reports.mkdir(parents=True, exist_ok=True)

    hierarchy = pd.read_csv(tables / "02_ATLAS_CELLTYPE_HIERARCHY.tsv", sep="\t")
    library = pd.read_csv(tables / "03_EXTENDED_IMMUNE_MARKER_LIBRARY.tsv", sep="\t")
    detect_patient = pd.read_csv(
        tables / "06_STEREOSEQ_MARKER_DETECTABILITY_BY_PATIENT.tsv", sep="\t"
    )
    coverage = pd.read_csv(tables / "08_CORE_VS_EXTENDED_MARKER_COVERAGE.tsv", sep="\t")
    signature = pd.read_csv(tables / "10_SIGNATURE_DETECTABILITY_SUMMARY.tsv", sep="\t")
    final_audit = pd.read_csv(qc / "12_EXTENDED_IMMUNE_MARKER_FINAL_AUDIT.tsv", sep="\t")
    validation = pd.read_csv(qc / "00_OUTPUT_VALIDATION.tsv", sep="\t")
    input_manifest = pd.read_csv(output / "input_manifest.tsv", sep="\t")
    findings = (reports / "11_EXTENDED_IMMUNE_MARKER_KEY_FINDINGS_ZH.md").read_text(
        encoding="utf-8"
    )

    immune_hierarchy = hierarchy[hierarchy["is_immune_level2"]].copy()
    level1_n = immune_hierarchy["level1"].nunique()
    level2_n = immune_hierarchy["level2"].nunique()
    marker_record_n = len(library)
    marker_gene_n = library["marker_gene"].nunique()
    core_n = library.loc[
        library["marker_class"].eq("A. Lineage-core"), "marker_gene"
    ].nunique()
    extended_n = library.loc[
        library["signature_eligible"]
        & ~library["marker_class"].eq("A. Lineage-core"),
        "marker_gene",
    ].nunique()
    available_n = detect_patient.loc[detect_patient["gene_in_matrix"], "gene"].nunique()
    detectable_n = detect_patient.loc[
        detect_patient["expression_status"].eq("detected"), "gene"
    ].nunique()

    all21 = coverage[coverage["scope"].eq("all_21_patients")].copy()
    coverage_wide = all21.pivot(
        index="lineage", columns="panel", values=["panel_gene_n", "detectable_gene_n"]
    )
    coverage_wide.columns = [f"{a}: {b}" for a, b in coverage_wide.columns]
    coverage_wide = coverage_wide.reset_index()

    previews = render_pdf_previews(figures, assets)
    qmd_path = scripts / "03_extended_immune_marker_report.qmd"
    rel_assets = [
        (title, Path(os.path.relpath(path, scripts)).as_posix())
        for title, path in previews
    ]

    body = [
        "---",
        'title: "Extended Immune Marker Library and Stereo-seq RNA Detectability Audit"',
        'author: "Stereo-seq analysis team"',
        "date: today",
        "format:",
        "  html:",
        "    toc: true",
        "    toc-float: true",
        "    code-fold: true",
        "    theme: bootstrap",
        "    embed-resources: true",
        "execute:",
        "  enabled: false",
        "---",
        "",
        "## Scope and approved analysis boundary",
        "",
        "This report summarizes the approved atlas marker-library construction and "
        "whole-tissue Bin50 raw-count detectability audit. It does not rerun RCTD, "
        "clustering, neighborhood scoring, disease testing, DEG/GSEA, CellChat, or "
        "frozen aggregate analyses.",
        "",
        "## Executive summary",
        "",
        f"- Immune atlas hierarchy: **{level1_n} Level 1** and **{level2_n} Level 2** categories.",
        f"- Immune marker library: **{marker_record_n:,} records** and **{marker_gene_n:,} unique genes**.",
        f"- Candidate panels: **{core_n} strict-core genes** and **{extended_n} eligible non-core genes**.",
        f"- Stereo-seq: **{available_n} genes available** and **{detectable_n} detected at least once**.",
        "- Increased panel coverage is evidence of improved RNA detectability only; it "
        "does not establish improved immune-cell identification accuracy.",
        "",
        "## Inputs and provenance",
        "",
        markdown_table(
            input_manifest[[
                column for column in input_manifest.columns
                if column in {
                    "absolute_path",
                    "file_name",
                    "size_bytes",
                    "modified_time",
                    "sha256",
                    "analysis_purpose",
                }
            ]]
        ),
        "",
        "## Core versus extended coverage",
        "",
        markdown_table(coverage_wide),
        "",
        "## Signature detectability summary",
        "",
        markdown_table(signature.head(60)),
        "",
        "## English figures",
        "",
    ]
    for title, rel_path in rel_assets:
        body.extend([f"### {title}", "", f"![]({rel_path})", ""])

    body.extend(
        [
            "## Chinese key findings",
            "",
            findings,
            "",
            "## Final audit",
            "",
            markdown_table(final_audit),
            "",
            "## Automated output validation",
            "",
            markdown_table(validation),
            "",
            "## Session information",
            "",
            f"- Python: `{platform.python_version()}`",
            f"- Platform: `{platform.platform()}`",
            f"- pandas: `{package_version('pandas')}`",
            f"- anndata: `{package_version('anndata')}`",
            f"- scipy: `{package_version('scipy')}`",
            f"- openpyxl: `{package_version('openpyxl')}`",
            f"- matplotlib: `{package_version('matplotlib')}`",
            f"- pypdfium2: `{package_version('pypdfium2')}`",
            "",
            "<!-- Exact render wrapper: powershell.exe -NoProfile -ExecutionPolicy Bypass "
            "-File scripts/04_render_extended_immune_report.ps1 -->",
        ]
    )
    qmd_path.write_text("\n".join(body) + "\n", encoding="utf-8")
    print(f"REPORT_SOURCE={qmd_path}")
    print(f"PDF_PREVIEWS={len(previews)}")
    print(
        "SESSION_INFO="
        f"Python {platform.python_version()}; pandas {package_version('pandas')}; "
        f"pypdfium2 {package_version('pypdfium2')}"
    )


if __name__ == "__main__":
    main()
