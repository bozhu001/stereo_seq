#!/usr/bin/env python3
"""Extend K8 K=7 state-removal fits to 20,000 iterations.

This targeted disease-blind numerical audit reuses the exact seeds and balanced
reference from script 1143. It tests whether the 5,000-iteration medoid factors
persist when every K8 restart is allowed up to 20,000 coordinate-descent steps.
"""

from __future__ import annotations

import argparse
import hashlib
import importlib.util
import json
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd
from joblib import Parallel, delayed
from scipy.optimize import linear_sum_assignment


ROOT = Path(__file__).resolve().parents[1]
BASE_SCRIPT = ROOT / "scripts/1093_reference18_per_chip_dual_nmf_rank_audit.py"
SOURCE_RUN = ROOT / (
    "results/reanalysis/"
    "REFERENCE18_NEIGHBORHOOD_NMF_K7_STATE_SENSITIVITY_20261016_210302_attempt5"
)
SCENARIOS = {
    "no_basal_lec16": ("Basal-like epithelial", "Lymphatic endothelial", 1),
    "no_b17": ("B", 2),
}
RANK = 7
CHIP = "K8"
ROUTE = "neighborhood50"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--n-jobs", type=int, default=4)
    parser.add_argument("--restarts", type=int, default=50)
    parser.add_argument("--max-iter", type=int, default=20_000)
    return parser.parse_args()


def load_base():
    spec = importlib.util.spec_from_file_location("dual_nmf_base", BASE_SCRIPT)
    if spec is None or spec.loader is None:
        raise RuntimeError("Could not load script 1093")
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


def subset(matrix, states, removed):
    keep = [index for index, state in enumerate(states) if state not in removed]
    result = matrix[:, keep].copy()
    result /= result.sum(axis=1, keepdims=True)
    return result, [states[index] for index in keep]


def cosine_matrix(first, second):
    first = first / np.maximum(np.linalg.norm(first, axis=1, keepdims=True), 1e-12)
    second = second / np.maximum(np.linalg.norm(second, axis=1, keepdims=True), 1e-12)
    return first @ second.T


def main() -> None:
    args = parse_args()
    output = args.output_dir.resolve()
    output.mkdir(parents=True, exist_ok=False)
    (output / "models").mkdir()
    base = load_base()
    datasets, _, input_paths = base.load_datasets()
    dataset = next(
        item for item in datasets if item.route == ROUTE and item.chip == CHIP
    )
    patients = dataset.metadata["patient_id"].astype(str).to_numpy()
    reference = base.balanced_indices(
        patients, base.ROUTES.index(ROUTE), base.CHIPS.index(CHIP)
    )
    summaries = []
    restart_tables = []
    source_paths = list(input_paths)
    started = time.time()

    for scenario, specification in SCENARIOS.items():
        removed = specification[:-1]
        scenario_index = specification[-1]
        matrix, states = subset(dataset.matrix, dataset.state_names, removed)
        seed_base = (
            base.SEED
            + 80_000_000
            + scenario_index * 10_000_000
            + base.CHIPS.index(CHIP) * 1_000_000
            + RANK * 10_000
        )
        print(f"[{time.strftime('%H:%M:%S')}] {scenario}: fitting", flush=True)
        fits = Parallel(n_jobs=args.n_jobs, prefer="threads")(
            delayed(base.fit_nmf)(
                matrix[reference], RANK, seed_base + restart, args.max_iter
            )
            for restart in range(args.restarts)
        )
        medoid_index, restart = base.select_medoid(fits)
        components = fits[medoid_index].components
        _, loadings, projection_iter = base.project_fixed(
            matrix, components, args.max_iter
        )
        source_path = SOURCE_RUN / "models" / (
            f"{scenario}_{CHIP}_K07_LOADINGS.tsv"
        )
        source = pd.read_csv(source_path, sep="\t", index_col=0)
        source_paths.append(source_path)
        if source.columns.tolist() != states:
            raise RuntimeError(f"State order mismatch: {scenario}")
        similarity = cosine_matrix(source.to_numpy(float), loadings)
        rows, columns = linear_sum_assignment(-similarity)
        matched = similarity[rows, columns]
        summaries.append(
            {
                "scenario": scenario,
                "chip": CHIP,
                "rank": RANK,
                "max_iter": args.max_iter,
                "restarts": args.restarts,
                "convergence_fraction": float(
                    np.mean([fit.n_iter < args.max_iter for fit in fits])
                ),
                "median_n_iter": float(np.median([fit.n_iter for fit in fits])),
                "maximum_n_iter": int(max(fit.n_iter for fit in fits)),
                "restart_fraction_all_factors_ge_0_80": float(
                    (restart["min_factor_cosine_to_medoid"] >= 0.80).mean()
                ),
                "restart_min_factor_cosine": float(
                    restart["min_factor_cosine_to_medoid"].min()
                ),
                "matched_5k_to_20k_mean_cosine": float(matched.mean()),
                "matched_5k_to_20k_min_cosine": float(matched.min()),
                "projection_n_iter": projection_iter,
                "selected_seed": fits[medoid_index].seed,
            }
        )
        restart.insert(0, "scenario", scenario)
        restart_tables.append(restart)
        np.savez_compressed(
            output / "models" / f"{scenario}_{CHIP}_K07_MAXITER20000.npz",
            components=components.astype(np.float32),
            normalized_loadings=loadings.astype(np.float32),
            state_names=np.asarray(states, dtype=str),
        )
        pd.DataFrame(
            loadings,
            index=[f"F{index + 1:02d}" for index in range(RANK)],
            columns=states,
        ).to_csv(
            output / "models" / f"{scenario}_{CHIP}_K07_MAXITER20000_LOADINGS.tsv",
            sep="\t",
            index_label="factor",
        )
        print(f"[{time.strftime('%H:%M:%S')}] {scenario}: complete", flush=True)

    summary = pd.DataFrame(summaries)
    summary.to_csv(output / "01_CONVERGENCE_EXTENSION.tsv", sep="\t", index=False)
    pd.concat(restart_tables, ignore_index=True).to_csv(
        output / "02_RESTARTS.tsv.gz", sep="\t", index=False, compression="gzip"
    )
    pd.DataFrame(
        [{"path": str(path.resolve()), "sha256": sha256(path)} for path in source_paths]
    ).to_csv(output / "00_INPUT_MANIFEST.tsv", sep="\t", index=False)
    checks = pd.DataFrame(
        [
            ("two_scenarios", len(summary) == 2, str(len(summary))),
            ("all_50_restarts", all(len(frame) == 50 for frame in restart_tables), str([len(frame) for frame in restart_tables])),
            ("all_factors_reproduced_ge_0_80", bool((summary["restart_fraction_all_factors_ge_0_80"] == 1).all()), str(summary["restart_fraction_all_factors_ge_0_80"].tolist())),
            ("5k_20k_medoid_min_cosine_ge_0_95", bool((summary["matched_5k_to_20k_min_cosine"] >= 0.95).all()), str(summary["matched_5k_to_20k_min_cosine"].tolist())),
            ("disease_blind", True, "No disease field supplied to NMF"),
        ],
        columns=["check", "pass", "detail"],
    )
    checks.to_csv(output / "03_VALIDATION.tsv", sep="\t", index=False)
    if not checks["pass"].all():
        raise RuntimeError(
            "Validation failed: "
            + ";".join(checks.loc[~checks["pass"], "check"])
        )
    status = {
        "status": "ANALYSIS_COMPLETE_PENDING_EVIDENCE_REVIEW",
        "elapsed_seconds": time.time() - started,
    }
    (output / "FINAL_STATUS.json").write_text(
        json.dumps(status, indent=2), encoding="utf-8"
    )
    (output / "COMPLETE.ok").write_text("PASS\n", encoding="utf-8")
    (output / "RUN_LOG.txt").write_text(
        f"Script: scripts/{Path(__file__).name}\n"
        f"Command: {sys.executable} scripts/{Path(__file__).name} "
        f"--output-dir {output} --n-jobs {args.n_jobs} "
        f"--restarts {args.restarts} --max-iter {args.max_iter}\n"
        "Disease used: no\nStatus: PASS\n",
        encoding="utf-8",
    )
    print(json.dumps(status, indent=2), flush=True)


if __name__ == "__main__":
    main()
