#!/usr/bin/env python3
"""Independent output validation and PDF rasterization for visual QA."""

from __future__ import annotations

import argparse
import hashlib
from pathlib import Path

import pandas as pd
import pypdfium2 as pdfium
from pypdf import PdfReader


EXPECTED_TABLES = [
    "01_ATLAS_MARKER_INPUT_AUDIT.tsv",
    "02_ATLAS_CELLTYPE_HIERARCHY.tsv",
    "03_EXTENDED_IMMUNE_MARKER_LIBRARY.tsv",
    "04_MARKER_SPECIFICITY_AND_OVERLAP_AUDIT.tsv",
    "05_ATLAS_TO_REFERENCE18_MAPPING.tsv",
    "06_STEREOSEQ_MARKER_DETECTABILITY_BY_PATIENT.tsv",
    "07_STEREOSEQ_MARKER_DETECTABILITY_BY_CHIP.tsv",
    "08_CORE_VS_EXTENDED_MARKER_COVERAGE.tsv",
    "09_EXTENDED_IMMUNE_SIGNATURE_CANDIDATES.tsv",
    "10_SIGNATURE_DETECTABILITY_SUMMARY.tsv",
]

EXPECTED_PDFS = [
    "01_ATLAS_IMMUNE_CELLTYPE_MARKER_COVERAGE_EN.pdf",
    "02_MARKER_DETECTION_BY_CHIP_EN.pdf",
    "03_IMMUNE_MARKER_PATIENT_COVERAGE_EN.pdf",
    "04_CORE_VS_EXTENDED_MARKER_COVERAGE_EN.pdf",
    "05_IMMUNE_SIGNATURE_OVERLAP_EN.pdf",
]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output-dir", required=True, type=Path)
    return parser.parse_args()


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        while chunk := handle.read(1024 * 1024):
            digest.update(chunk)
    return digest.hexdigest()


def check(rows: list[dict[str, object]], name: str, condition: bool, detail: str) -> None:
    rows.append({"check": name, "status": "PASS" if condition else "FAIL", "detail": detail})
    if not condition:
        raise AssertionError(f"{name}: {detail}")


def main() -> None:
    args = parse_args()
    output = args.output_dir.resolve()
    tables = output / "tables"
    figures = output / "figures"
    qc = output / "qc"
    render_dir = output / "tmp" / "pdfs"
    render_dir.mkdir(parents=True, exist_ok=True)
    rows: list[dict[str, object]] = []

    check(rows, "all_expected_tables_exist", all((tables / name).stat().st_size > 0 for name in EXPECTED_TABLES), ";".join(EXPECTED_TABLES))
    check(rows, "all_expected_pdfs_exist", all((figures / name).stat().st_size > 0 for name in EXPECTED_PDFS), ";".join(EXPECTED_PDFS))

    hierarchy = pd.read_csv(tables / EXPECTED_TABLES[1], sep="\t")
    library = pd.read_csv(tables / EXPECTED_TABLES[2], sep="\t")
    specificity = pd.read_csv(tables / EXPECTED_TABLES[3], sep="\t")
    mapping = pd.read_csv(tables / EXPECTED_TABLES[4], sep="\t")
    patient = pd.read_csv(tables / EXPECTED_TABLES[5], sep="\t")
    chip = pd.read_csv(tables / EXPECTED_TABLES[6], sep="\t")
    coverage = pd.read_csv(tables / EXPECTED_TABLES[7], sep="\t")
    signatures = pd.read_csv(tables / EXPECTED_TABLES[8], sep="\t")

    gene_n = patient["gene"].nunique()
    check(rows, "atlas_level2_count", len(hierarchy) == 59 and hierarchy["level2"].nunique() == 59, f"rows={len(hierarchy)}")
    check(rows, "immune_library_records", len(library) == 1450 and library["source_level2"].nunique() == 29, f"rows={len(library)}, subtypes={library['source_level2'].nunique()}")
    check(rows, "immune_library_unique_genes", library["marker_gene"].nunique() == len(specificity), f"library={library['marker_gene'].nunique()}, audit={len(specificity)}")
    check(rows, "marker_class_complete", library["marker_class"].notna().all(), str(library["marker_class"].value_counts().to_dict()))
    check(rows, "reference18_mapping_states_nonempty", mapping["reference18_state"].fillna("").str.len().gt(0).all(), f"rows={len(mapping)}")
    check(rows, "patient_table_cartesian", len(patient) == gene_n * 21, f"rows={len(patient)}, genes={gene_n}")
    check(rows, "chip_table_cartesian", len(chip) == gene_n * 3, f"rows={len(chip)}, genes={gene_n}")
    check(rows, "patient_count", patient["patient"].nunique() == 21, str(patient["patient"].nunique()))
    check(rows, "chip_set", set(patient["chip"]) == {"L3", "J2", "Y40102K8"}, ";".join(sorted(patient["chip"].unique())))
    check(rows, "all_patient_gene_keys_unique", not patient.duplicated(["gene", "chip", "patient"]).any(), "gene x chip x patient")
    check(rows, "all_chip_gene_keys_unique", not chip.duplicated(["gene", "chip"]).any(), "gene x chip")
    check(rows, "fractions_bounded", patient["detected_bin50_fraction"].between(0, 1).all() and chip["detected_bin50_fraction"].between(0, 1).all(), "0 <= fraction <= 1")
    unavailable = patient[~patient["gene_in_matrix"]]
    check(rows, "unavailable_counts_zero", (unavailable["total_raw_counts"] == 0).all() and (unavailable["detected_bin50_n"] == 0).all(), f"rows={len(unavailable)}")
    available = patient[patient["gene_in_matrix"]]
    check(rows, "available_status_not_unavailable", (available["expression_status"] != "unavailable").all(), f"rows={len(available)}")
    check(rows, "raw_counts_integer_nonnegative", (patient["total_raw_counts"] >= 0).all() and (patient["total_raw_counts"] % 1 == 0).all(), "patient totals")
    check(rows, "coverage_scope", set(coverage["scope"]) == {"patient", "all_21_patients"}, str(coverage["scope"].value_counts().to_dict()))
    check(rows, "three_signature_tiers", signatures["signature_type"].nunique() == 3, ";".join(sorted(signatures["signature_type"].unique())))

    pdf_rows = []
    for name in EXPECTED_PDFS:
        path = figures / name
        reader = PdfReader(str(path))
        page_count = len(reader.pages)
        check(rows, f"pdf_page_count_{name}", page_count == 1, f"pages={page_count}")
        document = pdfium.PdfDocument(str(path))
        for page_index in range(len(document)):
            page = document[page_index]
            bitmap = page.render(scale=1.6)
            image = bitmap.to_pil()
            png = render_dir / f"{path.stem}_page_{page_index + 1}.png"
            image.save(png)
            width, height = image.size
            pdf_rows.append({
                "pdf": name,
                "page": page_index + 1,
                "rendered_png": str(png.relative_to(output)).replace("\\", "/"),
                "width_px": width,
                "height_px": height,
                "png_nonempty": png.stat().st_size > 0,
                "pdf_size_bytes": path.stat().st_size,
            })
        document.close()
    pdf_qa = pd.DataFrame(pdf_rows)
    check(rows, "all_pdf_renders_nonempty", pdf_qa["png_nonempty"].all(), f"renders={len(pdf_qa)}")

    validation = pd.DataFrame(rows)
    validation.to_csv(qc / "00_OUTPUT_VALIDATION.tsv", sep="\t", index=False)
    pdf_qa.to_csv(qc / "00_PDF_RENDER_QA.tsv", sep="\t", index=False)

    manifest_rows = []
    for path in sorted(output.rglob("*")):
        relative_parts = path.relative_to(output).parts
        if (
            not path.is_file()
            or path.name == "output_manifest.tsv"
            or "tmp" in relative_parts
            or "__pycache__" in relative_parts
        ):
            continue
        manifest_rows.append({
            "relative_path": str(path.relative_to(output)).replace("\\", "/"),
            "purpose": "Final requested output, documentation, audit, reproducibility script, or run log",
            "status": "complete" if path.stat().st_size > 0 else "empty",
            "nonempty": path.stat().st_size > 0,
            "size_bytes": path.stat().st_size,
            "sha256": sha256(path),
        })
    pd.DataFrame(manifest_rows).to_csv(output / "output_manifest.tsv", sep="\t", index=False)
    print(f"Validation PASS: {len(validation)} checks; {len(pdf_qa)} PDF pages rendered for visual QA")


if __name__ == "__main__":
    main()
