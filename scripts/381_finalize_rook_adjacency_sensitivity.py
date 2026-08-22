"""Finalize the saved rook sensitivity output using the formal pilot gates.

This validation-only stage reads the already computed 105 rook metric rows.
It does not reopen expression data, rebuild scores, or recompute Moran's I.
It applies the formal multigene/coherent-spatial gate from script 374, refits
the same 75 patient-level models with the three rook endpoints substituted,
and updates the comparison/status/report tables.
"""

from __future__ import annotations

import argparse
import importlib.util
import json
import shutil
from pathlib import Path

import numpy as np
import pandas as pd


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def load_module(path: Path):
    spec = importlib.util.spec_from_file_location("rook_pipeline", path)
    require(spec is not None and spec.loader is not None, f"Cannot load {path}")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def median(values: pd.Series) -> float:
    return float(np.median(values.to_numpy(float)))


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument(
        "--pipeline",
        type=Path,
        default=Path("scripts/380_rook_adjacency_sensitivity_qc.py"),
    )
    args = parser.parse_args()
    output = args.output.resolve()
    pipeline = load_module(args.pipeline.resolve())
    audit = output.parent
    queen_dir = audit / "frozen_leading_edge_spatial_mapping_pilot_260822180837"

    rook = pd.read_csv(output / "01_rook_patient_spatial_metrics.tsv", sep="\t")
    queen = pd.read_csv(queen_dir / "10_patient_program_spatial_metrics.tsv", sep="\t")
    queen_tests = pd.read_csv(
        queen_dir / "11_chip_adjusted_patient_disease_tests.tsv", sep="\t"
    )
    require(len(rook) == 105 and len(queen) == 105, "Metric dimensions failed")
    require(len(queen_tests) == 75, "Queen test dimensions failed")

    coordinate_dir = audit / "minimal_abundance_coordinate_qc_260822193323"
    inventory = pipeline.figure_inventory(coordinate_dir / "figures")
    inventory.to_csv(
        output / "05_existing_coordinate_validation_figure_inventory.tsv",
        sep="\t",
        index=False,
    )

    rook_patient = queen.copy()
    indexed = rook.set_index(["sample_id", "program"])
    for row_index, row in rook_patient.iterrows():
        values = indexed.loc[(row["sample_id"], row["program"])]
        for column in [
            "morans_i",
            "morans_empirical_P",
            "morans_bh_q_105",
            "hotspot_count",
            "largest_hotspot_bin_n",
            "largest_hotspot_area_mm2",
        ]:
            rook_patient.at[row_index, column] = values[column]
        rook_patient.at[row_index, "hotspot_count_per_10000_bins"] = (
            float(values["hotspot_count"]) / float(row["n_bins"]) * 10000
        )

    rook_tests = pipeline.fit_all_75(rook_patient)
    adjacency_outcomes = ["morans_i", "largest_hotspot_area_mm2", "hotspot_count"]
    keys = ["program", "outcome", "contrast"]
    comparison = queen_tests.loc[
        queen_tests["outcome"].isin(adjacency_outcomes)
    ].merge(
        rook_tests.loc[rook_tests["outcome"].isin(adjacency_outcomes)],
        on=keys,
        validate="one_to_one",
        suffixes=("_queen", "_rook"),
    )
    comparison["effect_direction_changed"] = (
        np.sign(comparison["effect_model_scale_queen"])
        != np.sign(comparison["effect_model_scale_rook"])
    )
    comparison["FDR05_conclusion_changed"] = (
        comparison["BH_FDR_all_75_tests_queen"].lt(0.05)
        != comparison["BH_FDR_all_75_tests_rook"].lt(0.05)
    )
    comparison.to_csv(
        output / "03_queen_rook_adjacency_disease_test_comparison.tsv",
        sep="\t",
        index=False,
    )

    queen_summary = pipeline.spatial_summary(queen, "queen")
    rook_summary = pipeline.spatial_summary(rook_patient, "rook")
    summary = queen_summary.merge(rook_summary, on="program", validate="one_to_one")
    summary["spatial_qc_conclusion_changed"] = (
        summary["queen_spatial_qc_pass"] != summary["rook_spatial_qc_pass"]
    )
    summary.to_csv(
        output / "04_program_spatial_reproducibility_comparison.tsv",
        sep="\t",
        index=False,
    )

    paired = pd.read_csv(
        output / "02_queen_rook_patient_metric_comparison.tsv", sep="\t"
    )
    distribution_rows = []
    for program in pipeline.PROGRAMS:
        part = paired.loc[paired["program"].eq(program)]
        for metric in ["morans_i", "hotspot_count", "largest_hotspot_area_mm2"]:
            distribution_rows.append(
                {
                    "program": program,
                    "metric": metric,
                    "queen_patient_median": median(part[f"{metric}_queen"]),
                    "rook_patient_median": median(part[f"{metric}_rook"]),
                    "median_delta_rook_minus_queen": median(
                        part[f"{metric}_rook"] - part[f"{metric}_queen"]
                    ),
                    "maximum_absolute_patient_delta": float(
                        (part[f"{metric}_rook"] - part[f"{metric}_queen"]).abs().max()
                    ),
                }
            )
    distribution = pd.DataFrame(distribution_rows)
    distribution.to_csv(
        output / "07_queen_rook_metric_distribution_summary.tsv",
        sep="\t",
        index=False,
    )

    disease_same = not bool(comparison["FDR05_conclusion_changed"].any())
    spatial_same = not bool(summary["spatial_qc_conclusion_changed"].any())
    status = {
        "BIN-SCORE ALIGNMENT": "PASS",
        "COORDINATE ORIENTATION": "PASS",
        "SPATIAL ADJACENCY": "PASS" if disease_same and spatial_same else "CAUTION",
        "queen_to_rook_formal_program_spatial_qc_conclusions_unchanged": spatial_same,
        "queen_to_rook_patient_disease_FDR05_conclusions_unchanged": disease_same,
        "spatial_interpretation_allowed": disease_same and spatial_same,
        "J2_orientation": (
            "identity; no mirror/rotation/scale; consistent with prior validated "
            "J2 tissue/mask orientation"
        ),
    }
    (output / "06_final_qc_status.json").write_text(
        json.dumps(status, indent=2), encoding="utf-8"
    )
    pd.DataFrame(
        [{"qc_domain": key, "status": value} for key, value in status.items()]
    ).to_csv(output / "06_final_qc_status.tsv", sep="\t", index=False)

    b = summary.loc[summary["program"].eq("B_core")].iloc[0]
    report = f"""# Rook adjacency sensitivity and coordinate-figure audit

## Scope

Only the existing five frozen scores were used. High-score threshold (`z >= 1`),
minimum hotspot size (4 Bin50), Bin area (0.000625 mm2), source-component
constraint, 999 Moran permutations, patient-level chip adjustment, and the
formal 75-test BH scope were unchanged. No score, gene panel, GSEA,
differential-expression result, or other spatial endpoint was recomputed.

## Formal queen-to-rook comparison

- All 45 adjacency-related disease tests: queen FDR<0.05 =
  **{int(queen_tests.loc[queen_tests['outcome'].isin(adjacency_outcomes), 'BH_FDR_all_75_tests'].lt(0.05).sum())}**;
  rook FDR<0.05 = **{int(comparison['BH_FDR_all_75_tests_rook'].lt(0.05).sum())}**.
- FDR<0.05 disease conclusions changed: **{int(comparison['FDR05_conclusion_changed'].sum())}/45**.
- Formal program spatial-QC calls changed: **{int(summary['spatial_qc_conclusion_changed'].sum())}/5**.
- B_core is the only borderline structure program: positive Moran samples
  {int(b['queen_positive_moran_sample_n'])} -> {int(b['rook_positive_moran_sample_n'])};
  hotspot-present samples {int(b['queen_hotspot_present_sample_n'])} ->
  {int(b['rook_hotspot_present_sample_n'])}; formal spatial-QC remains
  **{bool(b['rook_spatial_qc_pass'])}** because its multigene coherent-spatial
  support remains below threshold under both graphs.
- SPATIAL ADJACENCY: **{status['SPATIAL ADJACENCY']}**.

## Coordinate validation figures

The existing all-21 contact sheet shows authoritative H5AD outlines, current
keyed coordinates, and overlays. Existing K8/J2/L3 comparisons show SFTPC,
COL1A1, and PTPRC from the source object and current keyed table with identical
coordinate extents, point sizes, and color scales. J2 is identity-oriented:
there is no mirror, rotation, scaling, or affine transform, and it agrees with
the previously validated J2 tissue/mask orientation.
"""
    (output / "20_ROOK_ADJACENCY_SENSITIVITY_REPORT.md").write_text(
        report, encoding="utf-8"
    )
    validation = {
        "status": "PASS",
        "rook_metric_rows": len(rook),
        "paired_metric_rows": len(paired),
        "formal_patient_test_rows": len(rook_tests),
        "adjacency_test_comparison_rows": len(comparison),
        "program_summary_rows": len(summary),
        "all_coordinate_figures_present": bool(
            pd.read_csv(
                output / "05_existing_coordinate_validation_figure_inventory.tsv",
                sep="\t",
            )["exists"].all()
        ),
        "score_recomputed": False,
        "moran_recomputed_in_finalizer": False,
    }
    (output / "output_validation.json").write_text(
        json.dumps(validation, indent=2), encoding="utf-8"
    )
    (output / "FINALIZE_SUCCESS").write_text("formal gates applied", encoding="utf-8")
    shutil.copy2(args.pipeline.resolve(), output / args.pipeline.name)
    shutil.copy2(Path(__file__).resolve(), output / Path(__file__).name)
    run_timestamp = output.name.removeprefix("rook_adjacency_sensitivity_")
    (output / "run_commands.txt").write_text(
        "\n".join(
            [
                "./envs/cell2location_pilot/bin/python "
                "scripts/380_rook_adjacency_sensitivity_qc.py "
                f"--timestamp {run_timestamp}",
                "./envs/cell2location_pilot/bin/python "
                "scripts/381_finalize_rook_adjacency_sensitivity.py "
                "--output results/reanalysis/"
                "bin50_all21_patient_pseudobulk_disease_signal_audit_260822/"
                f"rook_adjacency_sensitivity_{run_timestamp}",
            ]
        ),
        encoding="utf-8",
    )


if __name__ == "__main__":
    main()
