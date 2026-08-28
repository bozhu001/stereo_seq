#!/usr/bin/env python3
"""Summarize final L3 7/7 BCR/CDR3 support and compare with frozen J2.

This post-finalizer step reads only the frozen 453 classification output and
the completed J2 reclassification table. It does not open original FASTQs,
generate candidates, or run TRUST4. Counts distinguish corrected independent
CID-UMI molecules from distinct raw UMI sequences so sequence repetition alone
is never described as clonal expansion.

Run:
    python3 scripts/454_summarize_l3_all7_clonal_support.py \
      results/reanalysis/l3_all7_bcr_cdr3_completion_260828_004934
"""

from __future__ import annotations

import argparse
import csv
import html
import platform
from collections import Counter
from datetime import datetime
from pathlib import Path
from typing import Any


ROOT = Path(__file__).resolve().parents[1]
J2_FINAL = ROOT / (
    "results/reanalysis/l3_j2_final_umi_audit_260827_160202/"
    "final_reparsed_all_chunks_260827_234500/"
    "CLONOTYPE_RECLASSIFICATION_AFTER_UMI_QC.tsv"
)
PATIENT_ORDER = [
    "HC/NL-66",
    "HC/NL-72",
    "IPF/FO23-1-06168",
    "IPF/FO23-1-06170",
    "SSC/05957/17B",
    "SSC/07998/15A",
    "SSC/15491/14",
]


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t", quoting=csv.QUOTE_NONE))


def write_tsv(path: Path, rows: list[dict[str, Any]]) -> None:
    if not rows:
        raise RuntimeError(f"Refusing to write empty table: {path}")
    with path.open("x", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=list(rows[0]),
            delimiter="\t",
            lineterminator="\n",
            quoting=csv.QUOTE_NONE,
            escapechar="\\",
        )
        writer.writeheader()
        writer.writerows(rows)


def integer(row: dict[str, str], key: str) -> int:
    return int(row[key])


def patient_summary(rows: list[dict[str, str]]) -> list[dict[str, Any]]:
    output: list[dict[str, Any]] = []
    for patient in PATIENT_ORDER:
        subset = [row for row in rows if row["patient"] == patient]
        if not subset:
            raise RuntimeError(f"No final L3 clonotypes for {patient}")
        molecules = [integer(row, "corrected_molecules") for row in subset]
        repeated = [row for row in subset if integer(row, "corrected_molecules") >= 2]
        distinct_umi = [row for row in subset if integer(row, "unique_raw_umi") >= 2]
        multi_space = [
            row
            for row in repeated
            if integer(row, "unique_cid") >= 2 and integer(row, "unique_bin") >= 2
        ]
        reliable = [
            row
            for row in subset
            if row["final_classification"] == "reliable expanded clonotype"
        ]
        maximum = max(molecules)
        max_rows = [
            row for row in subset if integer(row, "corrected_molecules") == maximum
        ]
        if reliable:
            evidence = "formal_reliable_rule_hit_but_provisional_not_unqualified_expansion"
        else:
            evidence = "no_reliable_clonal_expansion"
        output.append(
            {
                "patient": patient,
                "disease_group": subset[0]["disease_group"],
                "quality_corrected_productive_bcr_molecules": sum(molecules),
                "unique_productive_bcr_clonotypes": len(subset),
                "clonotypes_ge2_corrected_cid_umi_molecules": len(repeated),
                "clonotypes_ge2_distinct_raw_umi_sequences": len(distinct_umi),
                "maximum_corrected_molecules_per_clonotype": maximum,
                "maximum_clonotype_ids": "|".join(
                    row["clonotype_id"] for row in max_rows
                ),
                "maximum_clonotype_chains": "|".join(
                    row["chain"] for row in max_rows
                ),
                "maximum_clonotype_unique_cid": max(
                    integer(row, "unique_cid") for row in max_rows
                ),
                "maximum_clonotype_unique_bin": max(
                    integer(row, "unique_bin") for row in max_rows
                ),
                "repeated_clonotypes_with_ge2_cid_and_ge2_bin": len(multi_space),
                "formal_reliable_rule_clonotypes": len(reliable),
                "formal_reliable_rule_molecules": sum(
                    integer(row, "corrected_molecules") for row in reliable
                ),
                "reliable_expansion_interpretation": evidence,
            }
        )
    return output


def chip_metrics(
    chip: str,
    rows: list[dict[str, str]],
    molecule_key: str,
    raw_umi_key: str,
    patient_key: str,
) -> dict[str, Any]:
    molecules = [integer(row, molecule_key) for row in rows]
    repeated = [row for row in rows if integer(row, molecule_key) >= 2]
    distinct_umi = [row for row in rows if integer(row, raw_umi_key) >= 2]
    multi_space = [
        row
        for row in repeated
        if integer(row, "unique_cid") >= 2 and integer(row, "unique_bin") >= 2
    ]
    reliable = [
        row for row in rows if row["final_classification"] == "reliable expanded clonotype"
    ]
    return {
        "chip": chip,
        "patients": len({row[patient_key] for row in rows}),
        "quality_corrected_productive_bcr_molecules": sum(molecules),
        "unique_productive_bcr_clonotypes": len(rows),
        "clonotypes_ge2_corrected_cid_umi_molecules": len(repeated),
        "clonotypes_ge2_distinct_raw_umi_sequences": len(distinct_umi),
        "maximum_corrected_molecules_per_clonotype": max(molecules),
        "repeated_clonotypes_with_ge2_cid_and_ge2_bin": len(multi_space),
        "singleton": sum(row["final_classification"] == "singleton" for row in rows),
        "technical_candidate": sum(
            row["final_classification"] == "technical candidate" for row in rows
        ),
        "moderate_candidate": sum(
            row["final_classification"] == "moderate candidate" for row in rows
        ),
        "formal_reliable_rule_clonotypes": len(reliable),
        "interpretation": (
            "one_formal_hit_reported_as_spatially_supported_provisional_expanded_IGH_candidate"
            if chip == "L3"
            else "no_reliable_expanded_clonotype"
        ),
    }


def markdown_table(rows: list[dict[str, Any]], fields: list[str]) -> str:
    header = "| " + " | ".join(fields) + " |"
    rule = "| " + " | ".join("---" for _ in fields) + " |"
    body = [
        "| " + " | ".join(str(row[field]) for field in fields) + " |"
        for row in rows
    ]
    return "\n".join([header, rule, *body])


def html_table(rows: list[dict[str, Any]]) -> str:
    fields = list(rows[0])
    head = "".join(f"<th>{html.escape(field)}</th>" for field in fields)
    body = "".join(
        "<tr>"
        + "".join(f"<td>{html.escape(str(row[field]))}</td>" for field in fields)
        + "</tr>"
        for row in rows
    )
    return f"<div class='scroll'><table><thead><tr>{head}</tr></thead><tbody>{body}</tbody></table></div>"


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("run_dir")
    args = parser.parse_args()
    run_dir = (ROOT / args.run_dir).resolve()
    if not (run_dir / "L3_ALL7_FINAL_COMPLETE.ok").is_file():
        raise RuntimeError("453 final completion marker is absent")

    l3_path = run_dir / "L3_ALL7_CLONOTYPE_CLASSIFICATION.tsv"
    l3_rows = read_tsv(l3_path)
    j2_rows = [row for row in read_tsv(J2_FINAL) if row["chip_id"] == "Y40105J2"]
    if len({row["patient"] for row in l3_rows}) != 7:
        raise RuntimeError("L3 classification does not contain 7 patients")
    if len({row["sample_id"] for row in j2_rows}) != 7:
        raise RuntimeError("J2 classification does not contain 7 patients")

    patient_rows = patient_summary(l3_rows)
    chip_rows = [
        chip_metrics(
            "J2",
            j2_rows,
            "corrected_independent_cid_umi_molecules",
            "unique_raw_umi_sequences",
            "sample_id",
        ),
        chip_metrics(
            "L3",
            l3_rows,
            "corrected_molecules",
            "unique_raw_umi",
            "patient",
        ),
    ]

    patient_path = run_dir / "L3_ALL7_PATIENT_CLONAL_SUPPORT_SUMMARY.tsv"
    comparison_path = run_dir / "L3_J2_IDENTICAL_THRESHOLD_METRICS.tsv"
    write_tsv(patient_path, patient_rows)
    write_tsv(comparison_path, chip_rows)

    stamp = datetime.now().strftime("%y%m%d%H%M%S")
    report_stem = f"L3_ALL7_FINAL_METRICS_REPORT_{stamp}"
    report_md = run_dir / f"{report_stem}.md"
    report_html = run_dir / f"{report_stem}.html"
    key_fields = [
        "patient",
        "quality_corrected_productive_bcr_molecules",
        "unique_productive_bcr_clonotypes",
        "clonotypes_ge2_corrected_cid_umi_molecules",
        "clonotypes_ge2_distinct_raw_umi_sequences",
        "maximum_corrected_molecules_per_clonotype",
        "maximum_clonotype_unique_cid",
        "maximum_clonotype_unique_bin",
        "formal_reliable_rule_clonotypes",
        "reliable_expansion_interpretation",
    ]
    class_counts = Counter(row["final_classification"] for row in l3_rows)
    report = f"""# L3 7/7 final BCR/CDR3 support metrics

Generated: {datetime.now().astimezone().isoformat()}

## Definitions

- A quality-corrected molecule is one unique CID plus corrected UMI within patient, chain and exact CDR3 nt.
- Sequence repetition, a single UMI, or a single CID is not clonal expansion.
- The uniform reliable rule requires IGH, at least 2 corrected molecules, at least 2 CID, at least 2 Bin, and at least 2 positional and strict R2 fragment signatures.
- Distinct raw UMI sequence counts are shown separately because identical UMI sequences in distinct CID can still be distinct molecules but are not distinct UMI sequences.

## L3 patient metrics

{markdown_table(patient_rows, key_fields)}

## J2 versus L3 under the frozen identical classification threshold

{markdown_table(chip_rows, list(chip_rows[0]))}

## Interpretation

L3 has {len(l3_rows)} productive BCR clonotypes: {class_counts['singleton']} singleton, {class_counts['technical candidate']} technical, {class_counts['moderate candidate']} moderate, and {class_counts['reliable expanded clonotype']} formal reliable-rule hit. The sole formal hit is SSC_15491_14_TRUST4_0003 (IGH; 3 corrected molecules; 3 CID; 3 Bin), retained as a spatially supported provisional expanded IGH candidate because all supporting reads carry the anomalously frequent CGCTTGGCCT motif and the exact CDR3 has a cross-patient warning. It is not reported as unqualified reliable clonal expansion.

Python: {platform.python_version()} ({platform.platform()})
"""
    report_md.write_text(report, encoding="utf-8")
    report_html.write_text(
        "<!doctype html><html><head><meta charset='utf-8'><title>L3 7/7 final metrics</title>"
        "<style>body{font:14px Arial,sans-serif;max-width:1500px;margin:24px auto;line-height:1.5}"
        "table{border-collapse:collapse;font-size:11px}th,td{border:1px solid #bbb;padding:5px;vertical-align:top}"
        "th{background:#eaf2f8}.scroll{overflow:auto}.callout{border-left:5px solid #2b6cb0;background:#eef6ff;padding:12px}</style>"
        "</head><body><h1>L3 7/7 final BCR/CDR3 support metrics</h1>"
        "<div class='callout'>Single UMI, single CID, or sequence repetition is never interpreted as clonal expansion.</div>"
        "<h2>L3 patient metrics</h2>"
        + html_table(patient_rows)
        + "<h2>J2 versus L3: identical frozen threshold</h2>"
        + html_table(chip_rows)
        + "<h2>Full interpretation and provenance</h2><pre>"
        + html.escape(report)
        + "</pre></body></html>",
        encoding="utf-8",
    )
    print(f"PATIENT_SUMMARY={patient_path}")
    print(f"J2_L3_COMPARISON={comparison_path}")
    print(f"REPORT_MD={report_md}")
    print(f"REPORT_HTML={report_html}")


if __name__ == "__main__":
    main()
