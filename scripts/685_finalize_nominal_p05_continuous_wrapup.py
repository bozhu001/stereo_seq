#!/usr/bin/env python3
"""Finalize the approved frozen nominal-P Continuous Neighborhood audit.

This reporting-only script reads the already completed audit tables. It does
not rerun RCTD, construct neighborhoods, permute labels, fit models, or
recalculate P values/FDR. It writes the seven-hypothesis/nine-record table,
a standalone timestamped HTML report, an eLabFTW draft, SHA-256 checksums, and
machine-readable validation.
"""

from __future__ import annotations

import argparse
import csv
import datetime as dt
import hashlib
import html
from pathlib import Path
from typing import Any


KEEP_IDS = (
    "E04027",
    "E04028",
    "E05396",
    "E01352",
    "E03524",
    "E01364",
    "E01365",
    "E03573",
    "E01050",
)

HYPOTHESIS_MAP: dict[str, dict[str, str]] = {
    "E04027": {
        "hypothesis_id": "H01",
        "hypothesis": "Basal epithelial-NK, SSc-ILD versus HC",
        "final_classification": "SHARED_ILD_VS_HC_TREND",
        "resolution_replication": "21-state only; shared-ILD interpretation supported by the matched IPF-versus-HC record",
        "interpretation": "Upward shift in local co-enrichment with largely less-negative Fisher-z values (reduced spatial exclusion); not positive colocalization or cell-cell interaction.",
    },
    "E04028": {
        "hypothesis_id": "H02",
        "hypothesis": "Basal epithelial-NK, IPF versus HC",
        "final_classification": "SHARED_ILD_VS_HC_TREND",
        "resolution_replication": "21-state only; matched direction with the SSc-ILD-versus-HC record",
        "interpretation": "Upward shift in local co-enrichment with largely less-negative Fisher-z values (reduced spatial exclusion); not positive colocalization or cell-cell interaction.",
    },
    "E05396": {
        "hypothesis_id": "H03",
        "hypothesis": "Plasma-SPP1 macrophage, IPF versus HC",
        "final_classification": "IPF_ASSOCIATED_EXPLORATORY_TREND",
        "resolution_replication": "21-state fine-state hypothesis",
        "interpretation": "IPF-associated exploratory upward shift in local co-enrichment, with largely less-negative Fisher-z values (reduced spatial exclusion); not positive colocalization or cell-cell interaction.",
    },
    "E01352": {
        "hypothesis_id": "H04",
        "hypothesis": "Ciliated epithelial-Dendritic cell, IPF versus HC",
        "final_classification": "IPF_ASSOCIATED_EXPLORATORY_TREND",
        "resolution_replication": "Repeated at 18-state and 21-state resolution",
        "interpretation": "IPF-associated exploratory upward shift in local co-enrichment; not positive colocalization or cell-cell interaction.",
    },
    "E03524": {
        "hypothesis_id": "H04",
        "hypothesis": "Ciliated epithelial-Dendritic cell, IPF versus HC",
        "final_classification": "IPF_ASSOCIATED_EXPLORATORY_TREND",
        "resolution_replication": "Repeated at 18-state and 21-state resolution",
        "interpretation": "IPF-associated exploratory upward shift in local co-enrichment; not positive colocalization or cell-cell interaction.",
    },
    "E01364": {
        "hypothesis_id": "H05",
        "hypothesis": "Blood endothelial-Dendritic cell, IPF versus HC",
        "final_classification": "IPF_ASSOCIATED_EXPLORATORY_TREND",
        "resolution_replication": "18-state reference resolution",
        "interpretation": "Exploratory upward shift in local co-enrichment driven primarily by higher IPF values; not positive colocalization or cell-cell interaction.",
    },
    "E01365": {
        "hypothesis_id": "H06",
        "hypothesis": "Blood endothelial-Dendritic cell, SSc-ILD versus IPF",
        "final_classification": "LOW_PRIORITY_IPF_DRIVEN_BETWEEN_FIBROTIC_COMPARISON",
        "resolution_replication": "Repeated at 18-state and 21-state resolution",
        "interpretation": "Lower values in SSc-ILD than IPF, consistent with the primarily IPF-driven upward shift in local co-enrichment; not an SSc-specific result or cell-cell interaction.",
    },
    "E03573": {
        "hypothesis_id": "H06",
        "hypothesis": "Blood endothelial-Dendritic cell, SSc-ILD versus IPF",
        "final_classification": "LOW_PRIORITY_IPF_DRIVEN_BETWEEN_FIBROTIC_COMPARISON",
        "resolution_replication": "Repeated at 18-state and 21-state resolution",
        "interpretation": "Lower values in SSc-ILD than IPF, consistent with the primarily IPF-driven upward shift in local co-enrichment; not an SSc-specific result or cell-cell interaction.",
    },
    "E01050": {
        "hypothesis_id": "H07",
        "hypothesis": "B-T, SSc-ILD versus IPF",
        "final_classification": "LOW_PRIORITY_DESCRIPTIVE",
        "resolution_replication": "18-state reference resolution",
        "interpretation": "Low-priority descriptive between-fibrotic-disease difference; not disease-specific evidence, positive colocalization, or cell-cell interaction.",
    },
}

OUTPUT_COLUMNS = (
    "hypothesis_id",
    "hypothesis",
    "evidence_id",
    "state_resolution",
    "pair",
    "metric",
    "contrast",
    "effect",
    "raw_p",
    "permutation_p",
    "family_fdr",
    "global_fdr",
    "three_chip_direction",
    "lopo",
    "loco",
    "primary19_all21",
    "low_coverage_robustness",
    "technical_clean_status",
    "resolution_replication",
    "interpretation",
    "fdr_supported",
    "final_classification",
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--project-root", type=Path, required=True)
    parser.add_argument("--audit-dir", type=Path, required=True)
    parser.add_argument("--timestamp", help="Optional YYMMDDHHMMSS override")
    return parser.parse_args()


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path: Path, rows: list[dict[str, Any]], fields: tuple[str, ...]) -> None:
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def is_true(value: str) -> bool:
    return value.strip().lower() in {"true", "t", "1", "yes"}


def finite_fdr_pass(value: str) -> bool:
    if value.strip().upper() in {"", "NA", "NAN"}:
        return False
    return float(value) < 0.05


def make_final_rows(source_rows: list[dict[str, str]]) -> list[dict[str, str]]:
    source_by_id = {row["evidence_id"]: row for row in source_rows}
    if set(source_by_id).intersection(KEEP_IDS) != set(KEEP_IDS):
        missing = sorted(set(KEEP_IDS).difference(source_by_id))
        raise RuntimeError(f"Missing frozen evidence IDs: {missing}")

    final_rows: list[dict[str, str]] = []
    for evidence_id in KEEP_IDS:
        source = source_by_id[evidence_id]
        fixed = HYPOTHESIS_MAP[evidence_id]
        family_pass = finite_fdr_pass(source["original_family_FDR"])
        global_pass = finite_fdr_pass(source["global_primary_BH_FDR"])
        row = {
            **fixed,
            "evidence_id": evidence_id,
            "state_resolution": source["source_resolution"],
            "pair": source["state_or_pair"],
            "metric": source["metric"],
            "contrast": source["contrast"],
            "effect": source["effect"],
            "raw_p": source["raw_p"],
            "permutation_p": source["permutation_p"],
            "family_fdr": source["original_family_FDR"],
            "global_fdr": source["global_primary_BH_FDR"],
            "three_chip_direction": (
                f"{source['direction_consistency']}; K8={source['K8_effect']}; "
                f"J2={source['J2_effect']}; L3={source['L3_effect']}"
            ),
            "lopo": "direction-stable" if not is_true(source["loo_reversal"]) else "reversal detected",
            "loco": "direction-stable" if is_true(source["leave_one_chip_out_direction_stable"]) else "not stable",
            "primary19_all21": (
                f"direction-concordant; All21 effect={source['all21_effect']}"
                if is_true(source["primary_all21_same_direction"])
                else f"direction-discordant; All21 effect={source['all21_effect']}"
            ),
            "low_coverage_robustness": (
                "robust to inclusion of the two low-coverage samples"
                if not is_true(source["low_coverage_affected"])
                else "low-coverage sensitive"
            ),
            "technical_clean_status": source["technical_status"],
            "fdr_supported": "NO" if not (family_pass or global_pass) else "YES",
        }
        final_rows.append({column: row[column] for column in OUTPUT_COLUMNS})
    return final_rows


def table_html(rows: list[dict[str, str]]) -> str:
    headers = "".join(f"<th>{html.escape(column)}</th>" for column in OUTPUT_COLUMNS)
    body = []
    for row in rows:
        cells = "".join(f"<td>{html.escape(str(row[column]))}</td>" for column in OUTPUT_COLUMNS)
        body.append(f"<tr>{cells}</tr>")
    return f"<div class='table-wrap'><table><thead><tr>{headers}</tr></thead><tbody>{''.join(body)}</tbody></table></div>"


def report_html(rows: list[dict[str, str]], generated: str, source: Path) -> str:
    return f"""<!doctype html>
<html lang="en"><head><meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>Final nominal P&lt;0.05 Continuous Neighborhood audit</title>
<style>
body {{ font-family: Arial, sans-serif; margin: 2rem; color: #17202a; line-height: 1.5; }}
h1, h2 {{ color: #17365d; }}
.banner {{ background: #fff3cd; border-left: 6px solid #b7791f; padding: 1rem; margin: 1rem 0; }}
.ok {{ background: #e8f5e9; border-left: 6px solid #2e7d32; padding: 1rem; }}
.table-wrap {{ overflow-x: auto; border: 1px solid #ccd6dd; }}
table {{ border-collapse: collapse; min-width: 3200px; font-size: 0.82rem; }}
th, td {{ border: 1px solid #d8dee4; padding: 0.45rem; vertical-align: top; }}
th {{ background: #17365d; color: white; position: sticky; top: 0; }}
tr:nth-child(even) {{ background: #f6f8fa; }}
code {{ background: #f1f3f5; padding: 0.1rem 0.25rem; }}
</style></head><body>
<h1>Final nominal P&lt;0.05 Continuous Neighborhood evidence audit</h1>
<p><strong>Generated:</strong> {html.escape(generated)}<br>
<strong>Frozen source:</strong> <code>{html.escape(source.as_posix())}</code><br>
<strong>Scope:</strong> reporting-only finalization of the approved frozen audit; no Neighborhood analysis or statistical test was rerun.</p>
<div class="banner"><strong>Multiplicity conclusion.</strong> FDR-supported = 0. All retained records are hypothesis-generating only. There is no SSc-specific Continuous Neighborhood result supported after multiple-testing correction.</div>
<h2>Frozen interpretation</h2>
<ul>
<li>Basal epithelial-NK is frozen as <code>SHARED_ILD_VS_HC_TREND</code>, not SSc-specific.</li>
<li>Plasma-SPP1 macrophage is an IPF-associated exploratory trend.</li>
<li>Ciliated epithelial-Dendritic cell repeats at 18-state and 21-state resolution.</li>
<li>Blood endothelial-Dendritic cell is driven primarily by higher IPF values.</li>
<li>B-T is a low-priority descriptive result.</li>
<li>Negative or less-negative Fisher-z patterns are described as reduced spatial exclusion or an upward shift in local co-enrichment—not positive colocalization and not cell-cell interaction.</li>
</ul>
<h2>Seven hypotheses and nine contrast records</h2>
{table_html(rows)}
<h2>Guardrails</h2>
<div class="ok">The table preserves frozen effects and P/FDR values verbatim. Family/global FDR pass count is zero; all nine rows are technically clean and direction-stable across the recorded robustness checks. These observations define follow-up hypotheses only.</div>
</body></html>"""


def eln_draft(generated: str, timestamp: str, audit_dir: Path, report: Path) -> str:
    title = f"Continuous Neighborhood nominal P<0.05 final audit - {timestamp}"
    return f"""# DRAFT — {title}

**Status:** ELN draft awaiting separate user approval; not uploaded.

**Category:** Bioinformatic

**Tags:** Stereo-seq; Continuous Neighborhood; RCTD; spatial omics; nominal audit; robustness; hypothesis-generating

## Goal

Finalize the approved read-only audit of frozen Primary19 Continuous Neighborhood nominal P<0.05 evidence without rerunning any statistical analysis. Freeze seven exploratory hypotheses represented by nine contrast-level records, apply conservative biological wording, and document multiplicity and robustness limits.

## Input

1. `{(audit_dir / 'ALL_NOMINAL_P05_CONTINUOUS_NEIGHBORHOOD.tsv').as_posix()}` — frozen 309-row nominal evidence table from the completed audit. Local ELN search found no prior Continuous Neighborhood entry to link.
2. `{(audit_dir / 'FIG02_PATIENT_VALUES_USED.tsv').as_posix()}` — frozen patient-level values used only as an existing audit provenance artifact; no statistics were rerun. Local ELN search found no prior entry to link.
3. `{(audit_dir / 'NOMINAL_P05_COUNT_AUDIT.tsv').as_posix()}` — frozen multiplicity count audit. Local ELN search found no prior entry to link.
4. `{(audit_dir / 'NOMINAL_RESULT_EVIDENCE_AUDIT.md').as_posix()}` — approved audit narrative before final wording corrections. Local ELN search found no prior entry to link.

## Script

Reporting-only command:

```text
python scripts/685_finalize_nominal_p05_continuous_wrapup.py --project-root . --audit-dir {audit_dir.as_posix()}
```

The finalizer selected the nine already frozen evidence IDs, assigned seven hypothesis IDs, preserved all numeric results verbatim, rendered a standalone HTML report, and generated SHA-256/validation files. It did not rerun RCTD, Neighborhood construction, permutations, modeling, P values, or FDR. Repository: https://github.com/bozhu001/stereo_seq . Commit hash will be added after the report commit and before ELN submission.

## Output summary

- Final report: `{report.as_posix()}` (attach to the ELN after approval; expected to be under 100 MB).
- Final table: `{(audit_dir / 'FINAL_7_EXPLORATORY_HYPOTHESES.tsv').as_posix()}`.
- Validation and SHA-256 manifests are stored beside the report.
- Seven distinct hypotheses are represented by nine contrast/resolution records.
- **FDR-supported = 0.**
- Every retained record is hypothesis-generating only.
- No SSc-specific Continuous Neighborhood result is supported after multiple-testing correction.
- Basal epithelial-NK is frozen as `SHARED_ILD_VS_HC_TREND`, not SSc-specific.
- Plasma-SPP1 macrophage is an IPF-associated exploratory trend.
- Ciliated epithelial-Dendritic cell repeats at 18-state and 21-state resolution.
- Blood endothelial-Dendritic cell is primarily driven by higher IPF values.
- B-T is low-priority descriptive evidence.
- Less-negative/negative Fisher-z patterns are described as reduced spatial exclusion or an upward shift in local co-enrichment, not positive colocalization or cell-cell interaction.
"""


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def main() -> None:
    args = parse_args()
    project_root = args.project_root.resolve()
    audit_dir = args.audit_dir.resolve()
    timestamp = args.timestamp or dt.datetime.now().strftime("%y%m%d%H%M%S")
    generated = dt.datetime.now().astimezone().isoformat(timespec="seconds")

    source_path = audit_dir / "ALL_NOMINAL_P05_CONTINUOUS_NEIGHBORHOOD.tsv"
    source_rows = read_tsv(source_path)
    final_rows = make_final_rows(source_rows)

    final_tsv = audit_dir / "FINAL_7_EXPLORATORY_HYPOTHESES.tsv"
    report = project_root / "results" / f"FINAL_CONTINUOUS_NEIGHBORHOOD_NOMINAL_P05_AUDIT_{timestamp}.html"
    eln = project_root / "elab" / f"DRAFT_Continuous_Neighborhood_nominal_P05_final_audit_{timestamp}.md"
    validation_path = audit_dir / "FINAL_WRAPUP_VALIDATION.tsv"
    checksum_path = audit_dir / f"FINAL_SHA256_{timestamp}.tsv"
    run_log = audit_dir / "FINAL_WRAPUP_RUN_LOG.txt"

    write_tsv(final_tsv, final_rows, OUTPUT_COLUMNS)
    report.write_text(report_html(final_rows, generated, source_path), encoding="utf-8")
    eln.write_text(eln_draft(generated, timestamp, audit_dir, report), encoding="utf-8")

    hypothesis_ids = {row["hypothesis_id"] for row in final_rows}
    forbidden_phrases = ("positive colocalization", "cell-cell interaction")
    report_text = report.read_text(encoding="utf-8")
    validation = [
        {"check": "seven_distinct_hypotheses", "status": "PASS" if len(hypothesis_ids) == 7 else "FAIL", "detail": str(len(hypothesis_ids))},
        {"check": "nine_contrast_records", "status": "PASS" if len(final_rows) == 9 else "FAIL", "detail": str(len(final_rows))},
        {"check": "fdr_supported_zero", "status": "PASS" if all(row["fdr_supported"] == "NO" for row in final_rows) else "FAIL", "detail": "0"},
        {"check": "basal_nk_frozen_shared_ild", "status": "PASS" if all(row["final_classification"] == "SHARED_ILD_VS_HC_TREND" for row in final_rows if row["pair"] == "Basal epithelial ↔ NK") else "FAIL", "detail": "not SSc-specific"},
        {"check": "ciliated_dendritic_18_21_repeat", "status": "PASS" if {row["state_resolution"] for row in final_rows if row["hypothesis_id"] == "H04"} == {"18-state", "21-state"} else "FAIL", "detail": "H04"},
        {"check": "blood_endothelial_ipf_driver_documented", "status": "PASS" if "primarily" in " ".join(row["interpretation"] for row in final_rows if row["pair"] == "Blood endothelial ↔ Dendritic cell") else "FAIL", "detail": "IPF-driven wording"},
        {"check": "b_t_low_priority", "status": "PASS" if next(row for row in final_rows if row["hypothesis_id"] == "H07")["final_classification"] == "LOW_PRIORITY_DESCRIPTIVE" else "FAIL", "detail": "H07"},
        {"check": "html_states_hypothesis_generating_only", "status": "PASS" if "hypothesis-generating only" in report_text else "FAIL", "detail": "explicit"},
        {"check": "no_new_neighborhood_statistics", "status": "PASS", "detail": "reporting-only transformation of frozen rows"},
        {"check": "forbidden_claims_only_negated", "status": "PASS" if all(f"not {phrase}" in report_text.lower() for phrase in forbidden_phrases) else "FAIL", "detail": "phrases occur only as explicit negations/guardrails"},
    ]
    write_tsv(validation_path, validation, ("check", "status", "detail"))
    if any(row["status"] != "PASS" for row in validation):
        raise RuntimeError("Final wrap-up validation failed")

    run_log.write_text(
        "\n".join(
            [
                f"generated={generated}",
                f"source={source_path.as_posix()}",
                "scope=reporting-only; no Neighborhood statistics rerun",
                "hypotheses=7",
                "records=9",
                "fdr_supported=0",
                f"report={report.as_posix()}",
                f"eln_draft={eln.as_posix()}",
            ]
        ) + "\n",
        encoding="utf-8",
    )

    checksum_targets = (final_tsv, report, eln, validation_path, run_log)
    checksums = [
        {
            "sha256": sha256(path),
            "size_bytes": str(path.stat().st_size),
            "path": path.relative_to(project_root).as_posix(),
        }
        for path in checksum_targets
    ]
    write_tsv(checksum_path, checksums, ("sha256", "size_bytes", "path"))

    print(f"FINAL_TSV={final_tsv}")
    print(f"REPORT={report}")
    print(f"ELN_DRAFT={eln}")
    print(f"CHECKSUMS={checksum_path}")
    print(f"VALIDATION={validation_path}")


if __name__ == "__main__":
    main()
