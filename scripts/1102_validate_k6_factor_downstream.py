#!/usr/bin/env python3
"""Validate the frozen K6 factor downstream-analysis deliverables.

This script performs structural and numerical checks only. It does not refit
RCTD or NMF, change factor labels, or perform disease comparisons.
"""

from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd


def add_check(rows: list[dict[str, str]], check: str, passed: bool, detail: str) -> None:
    rows.append(
        {
            "check": check,
            "status": "PASS" if passed else "FAIL",
            "detail": detail,
        }
    )


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--run-dir", required=True, type=Path)
    args = parser.parse_args()
    run_dir = args.run_dir.resolve()

    factor_path = run_dir / "01_FACTOR_DEFINITIONS.tsv"
    metrics_path = run_dir / "02_PATIENT_FACTOR_METRICS.tsv"
    expression_path = run_dir / "03_PATIENT_FACTOR_RAW_EXPRESSION_VALIDATION.tsv"
    summary_path = run_dir / "04_FACTOR_RAW_EXPRESSION_SUMMARY.tsv"
    marker_path = run_dir / "05_MARKER_PANEL_AUDIT.tsv"
    final_path = run_dir / "06_FINAL_FACTOR_ANALYSIS_TABLE.tsv"

    factors = pd.read_csv(factor_path, sep="\t")
    metrics = pd.read_csv(metrics_path, sep="\t")
    expression = pd.read_csv(expression_path, sep="\t")
    summary = pd.read_csv(summary_path, sep="\t")
    markers = pd.read_csv(marker_path, sep="\t")
    final = pd.read_csv(final_path, sep="\t")

    checks: list[dict[str, str]] = []
    expected_pairs = {(chip, f"F{i:02d}") for chip in ("K8", "J2", "L3") for i in range(1, 7)}

    for name, table in (("factor definitions", factors), ("expression summary", summary), ("final table", final)):
        pairs = set(zip(table["chip"], table["factor"]))
        add_check(checks, f"{name}: 18 unique chip-factor rows", pairs == expected_pairs and len(table) == 18, f"rows={len(table)}; unique_pairs={len(pairs)}")

    for name, table in (("patient metrics", metrics), ("patient expression", expression)):
        pairs = set(zip(table["chip"], table["factor"]))
        add_check(checks, f"{name}: all factors represented", pairs == expected_pairs, f"unique_pairs={len(pairs)}")
        add_check(checks, f"{name}: 18 factors x 7 patients", len(table) == 126, f"rows={len(table)}")
        per_factor = table.groupby(["chip", "factor"])["patient_id"].nunique()
        add_check(checks, f"{name}: seven patients per factor", bool((per_factor == 7).all()), f"range={per_factor.min()}-{per_factor.max()}")

    add_check(checks, "all 21 patients represented", metrics["patient_id"].nunique() == 21, f"n={metrics['patient_id'].nunique()}")
    add_check(checks, "three chips represented", set(metrics["chip"]) == {"K8", "J2", "L3"}, ",".join(sorted(metrics["chip"].unique())))
    add_check(checks, "marker audit has 18 Reference18 states", len(markers) == 18 and markers["reference18_state"].nunique() == 18, f"rows={len(markers)}")
    add_check(checks, "every marker panel has raw-count genes", bool((markers["covered_marker_n"] > 0).all()), f"minimum={markers['covered_marker_n'].min()}")

    numeric_metrics = [
        "mean_usage",
        "high_usage_fraction",
        "usage_umi_spearman",
        "usage_n_genes_spearman",
        "usage_neighbor_n_spearman",
        "technical_residual_usage_moran_i",
    ]
    missing_metrics = int(metrics[numeric_metrics].isna().sum().sum())
    add_check(checks, "critical patient metrics are finite", missing_metrics == 0, f"missing_values={missing_metrics}")

    numeric_expression = [
        "neighbor_top3_marker_residual_spearman",
        "top3_composite_marker_residual_spearman",
        "top3_composite_marker_log2fc",
    ]
    missing_expression = int(expression[numeric_expression].isna().sum().sum())
    add_check(checks, "critical raw-expression metrics are finite", missing_expression == 0, f"missing_values={missing_expression}")

    usage_range_ok = bool(metrics["mean_usage"].between(0, 1).all())
    fractions_ok = bool(metrics["high_usage_fraction"].between(0, 1).all())
    add_check(checks, "usage metrics have valid ranges", usage_range_ok and fractions_ok, "expected range 0-1")

    figures = sorted((run_dir / "figures").glob("*.png"))
    add_check(checks, "four non-empty English audit figures", len(figures) == 4 and all(p.stat().st_size > 0 for p in figures), f"n={len(figures)}")
    add_check(checks, "run completion marker exists", (run_dir / "COMPLETE.ok").exists(), str(run_dir / "COMPLETE.ok"))
    add_check(checks, "input manifest exists", (run_dir / "00_INPUT_MANIFEST.tsv").exists(), str(run_dir / "00_INPUT_MANIFEST.tsv"))

    forbidden = [p.name for p in run_dir.iterdir() if "disease" in p.name.lower() or "aggregate" in p.name.lower()]
    add_check(checks, "no disease-comparison or aggregate output created", not forbidden, ",".join(forbidden) if forbidden else "none")

    validation = pd.DataFrame(checks)
    output = run_dir / "08_VALIDATION.tsv"
    validation.to_csv(output, sep="\t", index=False)
    failures = validation.loc[validation["status"] == "FAIL"]
    print(validation.to_string(index=False))
    print(f"\nValidation output: {output}")
    return 1 if len(failures) else 0


if __name__ == "__main__":
    raise SystemExit(main())
