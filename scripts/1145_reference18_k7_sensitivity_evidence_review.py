#!/usr/bin/env python3
"""Review K=7 state-sensitivity evidence without one-to-one matching artifacts."""

from __future__ import annotations

import argparse
import hashlib
import importlib.util
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
ANALYSIS_SCRIPT = ROOT / "scripts/1143_reference18_neighborhood_nmf_k7_state_sensitivity.py"
BASE_SCRIPT = ROOT / "scripts/1093_reference18_per_chip_dual_nmf_rank_audit.py"
PRIMARY_RUN = ROOT / "results/reanalysis/REFERENCE18_NEIGHBORHOOD_NMF_K7_K12_EXTENSION_20261015_attempt5"
SENSITIVITY_RUN = ROOT / "results/reanalysis/REFERENCE18_NEIGHBORHOOD_NMF_K7_STATE_SENSITIVITY_20261016_210302_attempt5"
CONVERGENCE_RUN = ROOT / "results/reanalysis/REFERENCE18_NEIGHBORHOOD_NMF_K8_K7_CONVERGENCE_EXTENSION_20261016_212100"
CHIPS = ("K8", "J2", "L3")
SCENARIOS = ("no_basal_lec16", "no_b17")
DELETED_FAMILY = {
    "no_basal_lec16": "LEC-associated",
    "no_b17": "B-enriched/mixed-immune",
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output-dir", type=Path, required=True)
    return parser.parse_args()


def load_module(name: str, path: Path):
    spec = importlib.util.spec_from_file_location(name, path)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"Could not load {path}")
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def load_model(chip: str, scenario: str) -> dict[str, object]:
    if scenario == "main18":
        stem = f"neighborhood50_{chip}_K07"
        root = PRIMARY_RUN / "models"
    else:
        stem = f"{scenario}_{chip}_K07"
        root = SENSITIVITY_RUN / "models"
    loading_path = root / f"{stem}_LOADINGS.tsv"
    model_path = root / f"{stem}_MODEL_AND_USAGE.npz"
    frame = pd.read_csv(loading_path, sep="\t", index_col=0)
    with np.load(model_path, allow_pickle=True) as archive:
        usage = archive["relative_usage"].astype(float)
        entity = archive["entity_key"].astype(str)
    return {
        "states": frame.columns.tolist(),
        "loadings": frame.to_numpy(float),
        "usage": usage,
        "entity": entity,
        "paths": [loading_path, model_path],
    }


def main() -> None:
    args = parse_args()
    output = args.output_dir.resolve()
    output.mkdir(parents=True, exist_ok=False)
    analysis = load_module("state_sensitivity", ANALYSIS_SCRIPT)
    base = load_module("dual_nmf_base_review", BASE_SCRIPT)
    datasets, _, _ = base.load_datasets()
    datasets = {
        item.chip: item for item in datasets if item.route == "neighborhood50"
    }
    models = {}
    paths = [
        SENSITIVITY_RUN / "02_MODEL_ROBUSTNESS.tsv",
        SENSITIVITY_RUN / "07_CROSS_CHIP_FAMILY_MAP.tsv",
        SENSITIVITY_RUN / "09_B_FACTOR_ALL_CANDIDATES.tsv",
        CONVERGENCE_RUN / "01_CONVERGENCE_EXTENSION.tsv",
    ]
    for chip in CHIPS:
        for scenario in ("main18", *SCENARIOS):
            models[(scenario, chip)] = load_model(chip, scenario)
            paths.extend(models[(scenario, chip)]["paths"])

    rows = []
    for chip in CHIPS:
        patients = datasets[chip].metadata["patient_id"].astype(str).to_numpy()
        primary = models[("main18", chip)]
        for scenario in SCENARIOS:
            sensitivity = models[(scenario, chip)]
            common = sensitivity["states"]
            left = analysis.aligned_loadings(
                primary["loadings"], primary["states"], common
            )
            right = analysis.aligned_loadings(
                sensitivity["loadings"], sensitivity["states"], common
            )
            similarity = analysis.cosine_matrix(left, right)
            best_right = similarity.argmax(axis=1)
            best_left = similarity.argmax(axis=0)
            for primary_index, sensitivity_index in enumerate(best_right):
                global_rho, patient_rho, jaccard = analysis.patient_spatial_metrics(
                    primary["usage"][:, primary_index],
                    sensitivity["usage"][:, sensitivity_index],
                    patients,
                )
                family = analysis.factor_label(
                    primary["loadings"][primary_index], primary["states"]
                )
                rows.append(
                    {
                        "chip": chip,
                        "scenario": scenario,
                        "primary_factor": f"F{primary_index + 1:02d}",
                        "primary_family": family,
                        "expected_to_persist": family != DELETED_FAMILY[scenario],
                        "best_sensitivity_factor": f"F{sensitivity_index + 1:02d}",
                        "best_loading_cosine": float(
                            similarity[primary_index, sensitivity_index]
                        ),
                        "reciprocal_best_loading_match": bool(
                            best_left[sensitivity_index] == primary_index
                        ),
                        "same_bin_global_usage_spearman": global_rho,
                        "median_patient_usage_spearman": patient_rho,
                        "median_patient_q90_hotspot_jaccard": jaccard,
                    }
                )
    independent = pd.DataFrame(rows)
    retained = independent.loc[independent["expected_to_persist"]].copy()
    summaries = []
    for (scenario, chip), frame in retained.groupby(["scenario", "chip"]):
        summaries.append(
            {
                "scenario": scenario,
                "chip": chip,
                "retained_family_n": len(frame),
                "minimum_best_loading_cosine": float(frame["best_loading_cosine"].min()),
                "median_best_loading_cosine": float(frame["best_loading_cosine"].median()),
                "minimum_patient_usage_spearman": float(frame["median_patient_usage_spearman"].min()),
                "median_patient_usage_spearman": float(frame["median_patient_usage_spearman"].median()),
                "minimum_q90_hotspot_jaccard": float(frame["median_patient_q90_hotspot_jaccard"].min()),
                "median_q90_hotspot_jaccard": float(frame["median_patient_q90_hotspot_jaccard"].median()),
                "reciprocal_best_fraction": float(frame["reciprocal_best_loading_match"].mean()),
            }
        )
    summary = pd.DataFrame(summaries)
    cross = pd.read_csv(SENSITIVITY_RUN / "07_CROSS_CHIP_FAMILY_MAP.tsv", sep="\t")
    cross["minimum_pairwise_cosine"] = cross[
        ["K8_J2_cosine", "K8_L3_cosine", "J2_L3_cosine"]
    ].min(axis=1)
    main_cross = cross.loc[cross["scenario"].eq("main18")]
    convergence = pd.read_csv(
        CONVERGENCE_RUN / "01_CONVERGENCE_EXTENSION.tsv", sep="\t"
    )
    b_candidates = pd.read_csv(
        SENSITIVITY_RUN / "09_B_FACTOR_ALL_CANDIDATES.tsv", sep="\t"
    )
    b_max = b_candidates.groupby("chip", as_index=False).agg(
        maximum_candidate_patient_usage_spearman=("median_patient_usage_spearman", "max"),
        maximum_candidate_q90_hotspot_jaccard=("median_patient_q90_hotspot_jaccard", "max"),
        maximum_candidate_loading_cosine=("common_state_loading_cosine", "max"),
    )

    decisions = pd.DataFrame(
        [
            (
                "common_K7_cross_chip_composition",
                bool((main_cross["minimum_pairwise_cosine"] >= 0.80).all()),
                f"minimum pairwise cosine={main_cross['minimum_pairwise_cosine'].min():.6f}",
            ),
            (
                "retained_families_loading_persist",
                bool((summary["minimum_best_loading_cosine"] >= 0.80).all()),
                f"minimum independent best cosine={summary['minimum_best_loading_cosine'].min():.6f}",
            ),
            (
                "retained_families_spatial_scores_persist",
                bool((summary["minimum_patient_usage_spearman"] >= 0.70).all()),
                f"minimum median patient rho={summary['minimum_patient_usage_spearman'].min():.6f}",
            ),
            (
                "K8_5k_solution_confirmed_at_20k",
                bool((convergence["matched_5k_to_20k_min_cosine"] >= 0.95).all()),
                f"minimum matched cosine={convergence['matched_5k_to_20k_min_cosine'].min():.6f}",
            ),
            (
                "B_factor_has_no_coherent_single_noB_spatial_successor",
                bool(
                    (b_max["maximum_candidate_patient_usage_spearman"] < 0.30).all()
                    and (b_max["maximum_candidate_q90_hotspot_jaccard"] < 0.10).all()
                ),
                "descriptive gate: every chip max median-patient rho <0.30 and Q90 Jaccard <0.10",
            ),
            ("disease_blind", True, "No disease field supplied to NMF or matching"),
        ],
        columns=["decision", "pass", "detail"],
    )
    independent.to_csv(output / "01_INDEPENDENT_BEST_MATCHES.tsv", sep="\t", index=False)
    summary.to_csv(output / "02_RETAINED_FAMILY_SUMMARY.tsv", sep="\t", index=False)
    b_max.to_csv(output / "03_B_SUCCESSOR_SUMMARY.tsv", sep="\t", index=False)
    decisions.to_csv(output / "04_DECISIONS.tsv", sep="\t", index=False)
    pd.DataFrame(
        [{"path": str(path.resolve()), "sha256": sha256(path)} for path in paths]
    ).to_csv(output / "00_INPUT_MANIFEST.tsv", sep="\t", index=False)
    if not decisions["pass"].all():
        raise RuntimeError(
            "Evidence gate failed: "
            + ";".join(decisions.loc[~decisions["pass"], "decision"])
        )
    status = {
        "status": "ANALYSIS_COMPLETE_PENDING_USER_APPROVAL",
        "all_evidence_gates_pass": True,
    }
    (output / "FINAL_STATUS.json").write_text(
        json.dumps(status, indent=2), encoding="utf-8"
    )
    (output / "COMPLETE.ok").write_text("PASS\n", encoding="utf-8")
    (output / "RUN_LOG.txt").write_text(
        f"Script: scripts/{Path(__file__).name}\n"
        f"Command: {sys.executable} scripts/{Path(__file__).name} --output-dir {output}\n"
        "Scope: read-only evidence review of persisted K=7 models\n"
        "Disease used: no\nStatus: PASS\n",
        encoding="utf-8",
    )
    print(json.dumps(status, indent=2))


if __name__ == "__main__":
    main()
