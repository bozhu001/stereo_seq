"""Rook-adjacency sensitivity audit for five frozen Bin50 programs.

Purpose
-------
Recompute only Moran's I, hotspot count, and largest hotspot area from the
already frozen program-score table, replacing the historical queen graph with
a 50-native-unit rook graph. Compare every patient/program result and the
corresponding chip-adjusted patient tests with the completed queen analysis.
No gene panel, score, GSEA, differential-expression result, or other spatial
endpoint is created or modified.

Inputs
------
* keyed frozen score table from the completed coordinate audit
* completed queen patient metrics and 75 disease tests
* frozen topology HDF5 (for source-component labels)
* completed coordinate-QC status and existing validation figures

Outputs
-------
A timestamped directory containing paired queen/rook tables, final adjacency
status, a concise markdown report, provenance, software versions, and logs.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import platform
import shutil
import sys
import time
from collections import Counter
from datetime import datetime
from pathlib import Path

import h5py
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from statsmodels.formula.api import ols
from statsmodels.stats.multitest import multipletests


SEED = 20260822
PERMUTATIONS = 999
STEP = 50
BIN_AREA_MM2 = 0.000625
HIGH_Z = 1.0
MIN_HOTSPOT_BINS = 4
PROGRAMS = [
    "Fibroblast_ECM",
    "B_core",
    "Cilia_Axoneme",
    "AT2",
    "TNF_NFkB_Oxidative_Stress",
]
ROOK_HALF_OFFSETS = [(0, STEP), (STEP, 0)]


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def write_json(path: Path, value: object) -> None:
    path.write_text(
        json.dumps(value, ensure_ascii=False, indent=2, default=str),
        encoding="utf-8",
    )


def morans_i(values: np.ndarray, u: np.ndarray, v: np.ndarray) -> float:
    centered = np.asarray(values, dtype=float) - float(np.mean(values))
    denominator = float(np.dot(centered, centered))
    if denominator <= 0 or len(u) == 0:
        return float("nan")
    numerator = float(2.0 * np.sum(centered[u] * centered[v]))
    return float(len(centered) / (2.0 * len(u)) * numerator / denominator)


def empirical_p(observed: float, null: list[float]) -> float:
    finite = np.asarray([value for value in null if np.isfinite(value)])
    require(len(finite) > 0, "No finite Moran null values")
    return float((1 + np.sum(finite >= observed)) / (1 + len(finite)))


def permute_within_components(
    values: np.ndarray,
    components: np.ndarray,
    rng: np.random.Generator,
) -> np.ndarray:
    result = np.array(values, copy=True)
    for component in np.unique(components):
        indices = np.flatnonzero(components == component)
        result[indices] = values[rng.permutation(indices)]
    return result


def component_sizes(mask: np.ndarray, u: np.ndarray, v: np.ndarray) -> list[int]:
    selected = np.flatnonzero(mask)
    if len(selected) == 0:
        return []
    parent = np.arange(len(mask), dtype=np.int64)

    def find(value: int) -> int:
        while parent[value] != value:
            parent[value] = parent[parent[value]]
            value = int(parent[value])
        return value

    for left, right in zip(u, v, strict=True):
        if not (mask[left] and mask[right]):
            continue
        root_left = find(int(left))
        root_right = find(int(right))
        if root_left != root_right:
            parent[root_right] = root_left
    counts = Counter(find(int(index)) for index in selected)
    return sorted(counts.values(), reverse=True)


def attach_components(group: h5py.Group, n_bins: int) -> np.ndarray:
    source_to_g15 = np.asarray(group["source_to_g15_local"][:])
    source_component = np.asarray(group["source_component"][:])
    labels = np.full(n_bins, -1, dtype=np.int64)
    for source_index, local_g15 in enumerate(source_to_g15):
        if local_g15 >= 0:
            labels[int(local_g15)] = int(source_component[source_index])
    require(np.all(labels >= 0), "Missing source-component label")
    return labels


def rook_edges(coords: np.ndarray, components: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    index = {tuple(map(int, xy)): i for i, xy in enumerate(coords)}
    require(len(index) == len(coords), "Duplicated coordinates")
    left: list[int] = []
    right: list[int] = []
    for i, (x, y) in enumerate(coords.astype(int)):
        for dx, dy in ROOK_HALF_OFFSETS:
            j = index.get((int(x + dx), int(y + dy)))
            if j is not None and components[i] == components[j]:
                left.append(i)
                right.append(j)
    return np.asarray(left, dtype=np.int64), np.asarray(right, dtype=np.int64)


def fit_all_75(patient: pd.DataFrame) -> pd.DataFrame:
    contrasts = {
        "IPF_vs_HC": ("IPF", "HC"),
        "SScILD_vs_HC": ("SSc-ILD", "HC"),
        "SScILD_vs_IPF": ("SSc-ILD", "IPF"),
    }
    outcomes = {
        "program_median_raw_rank": ("program_median_raw_rank", "identity", False),
        "high_score_area_fraction": ("high_score_area_fraction", "identity", False),
        "morans_i": ("morans_i", "identity", False),
        "largest_hotspot_area_mm2": ("largest_hotspot_area_mm2", "log1p", False),
        "hotspot_count": ("hotspot_count", "log1p", True),
    }
    rows: list[dict[str, object]] = []
    for program in PROGRAMS:
        data = patient.loc[patient["program"].eq(program)].copy()
        require(len(data) == 21, f"Unexpected patient count: {program}")
        data["disease_group"] = pd.Categorical(
            data["disease_group"], categories=["HC", "IPF", "SSc-ILD"]
        )
        data["chip_id"] = pd.Categorical(data["chip_id"])
        for outcome, (column, transform, adjust_area) in outcomes.items():
            data["model_value"] = data[column].astype(float)
            if transform == "log1p":
                data["model_value"] = np.log1p(data["model_value"])
            formula = (
                "model_value ~ C(chip_id) + "
                "C(disease_group, Treatment(reference='HC'))"
            )
            if adjust_area:
                data["log_n_bins"] = np.log(data["n_bins"].astype(float))
                formula += " + log_n_bins"
            fit = ols(formula, data=data).fit(cov_type="HC3")
            names = list(fit.params.index)
            ipf = "C(disease_group, Treatment(reference='HC'))[T.IPF]"
            ssc = "C(disease_group, Treatment(reference='HC'))[T.SSc-ILD]"
            for contrast, (left, right) in contrasts.items():
                vector = np.zeros(len(names))
                if left == "IPF":
                    vector[names.index(ipf)] += 1
                elif left == "SSc-ILD":
                    vector[names.index(ssc)] += 1
                if right == "IPF":
                    vector[names.index(ipf)] -= 1
                elif right == "SSc-ILD":
                    vector[names.index(ssc)] -= 1
                test = fit.t_test(vector)
                rows.append(
                    {
                        "program": program,
                        "outcome": outcome,
                        "contrast": contrast,
                        "effect_model_scale": float(np.asarray(test.effect).ravel()[0]),
                        "standard_error_HC3": float(np.asarray(test.sd).ravel()[0]),
                        "raw_P": float(np.asarray(test.pvalue).ravel()[0]),
                        "model": formula,
                        "transform": transform,
                        "n_patient": 21,
                    }
                )
    result = pd.DataFrame(rows)
    require(len(result) == 75, "Disease-test count drift")
    result["BH_FDR_all_75_tests"] = multipletests(result["raw_P"], method="fdr_bh")[1]
    return result


def spatial_summary(patient: pd.DataFrame, label: str) -> pd.DataFrame:
    rows = []
    for program in PROGRAMS:
        part = patient.loc[patient["program"].eq(program)].copy()
        positive = part["morans_i"].gt(0) & part["morans_bh_q_105"].lt(0.05)
        multi_gene_support = (
            part["detected_gene_n_ge_0_1pct"].ge(2)
            & part["high_bins_two_gene_fraction"].ge(0.50)
        )
        coherent = positive & part["hotspot_count"].gt(0) & multi_gene_support
        by_chip = (
            part.assign(coherent=coherent)
            .groupby("chip_id", observed=True)["coherent"]
            .sum()
            .to_dict()
        )
        positive_n = int(positive.sum())
        hotspot_n = int(part["hotspot_count"].gt(0).sum())
        coherent_n = int(coherent.sum())
        passed = bool(
            coherent_n >= 11
            and len(by_chip) == 3
            and all(value >= 2 for value in by_chip.values())
        )
        rows.append(
            {
                "program": program,
                f"{label}_positive_moran_sample_n": positive_n,
                f"{label}_hotspot_present_sample_n": hotspot_n,
                f"{label}_coherent_spatial_sample_n": coherent_n,
                f"{label}_coherent_K8_n": int(by_chip.get("Y40102K8", 0)),
                f"{label}_coherent_J2_n": int(by_chip.get("Y40105J2", 0)),
                f"{label}_coherent_L3_n": int(by_chip.get("Y40105L3", 0)),
                f"{label}_spatial_qc_pass": passed,
            }
        )
    return pd.DataFrame(rows)


def figure_inventory(figure_dir: Path) -> pd.DataFrame:
    def project_relative(path: Path) -> str:
        parts = path.parts
        require("results" in parts, f"Figure is outside project results: {path}")
        return str(Path(*parts[parts.index("results") :]))

    rows = [
        {
            "chip": "ALL21",
            "content": "source H5AD outline; current keyed coordinates; overlay",
            "sample_id": "all 21",
            "path": project_relative(figure_dir / "coordinate_alignment_contact_sheet_all21.png"),
        },
        {
            "chip": "J2",
            "content": "source coordinates; current displayed coordinates; prior validated J2 overlay",
            "sample_id": "all J2 samples",
            "path": project_relative(figure_dir / "J2_source_current_prior_orientation_audit.png"),
        },
    ]
    representatives = {
        "K8": "SSC/24-1-18170A2",
        "J2": "SSC/23-105334B1",
        "L3": "SSC/15491/14",
    }
    safe = {
        "K8": "SSC_24_1_18170A2",
        "J2": "SSC_23_105334B1",
        "L3": "SSC_15491_14",
    }
    for chip in ["K8", "J2", "L3"]:
        rows.append(
            {
                "chip": chip,
                "content": "source-vs-current SFTPC, COL1A1, PTPRC maps; identical extent/point size/color scale",
                "sample_id": representatives[chip],
                "path": project_relative(
                    figure_dir / f"gene_map_comparison_{chip}_{safe[chip]}.png"
                ),
            }
        )
    result = pd.DataFrame(rows)
    result["exists"] = result["path"].map(lambda value: Path(value).is_file())
    result["size_bytes"] = result["path"].map(
        lambda value: Path(value).stat().st_size if Path(value).is_file() else 0
    )
    require(result["exists"].all(), "A required existing coordinate figure is missing")
    return result


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--timestamp", default=datetime.now().strftime("%y%m%d%H%M%S"))
    parser.add_argument(
        "--audit-base",
        type=Path,
        default=Path("results/reanalysis/bin50_all21_patient_pseudobulk_disease_signal_audit_260822"),
    )
    args = parser.parse_args()
    audit = args.audit_base.resolve()
    queen_dir = audit / "frozen_leading_edge_spatial_mapping_pilot_260822180837"
    coord_dir = audit / "minimal_abundance_coordinate_qc_260822193323"
    topology_path = Path(
        "results/reanalysis/rctd_bcell_spatial_niche_analysis_260809035034/prepared_niche_inputs.h5"
    ).resolve()
    score_path = coord_dir / "04_all21_frozen_program_scores_with_bin_id.tsv.gz"
    output = audit / f"rook_adjacency_sensitivity_{args.timestamp}"
    output.mkdir(parents=True, exist_ok=False)
    (output / "figures").mkdir()

    started = time.time()
    score = pd.read_csv(score_path, sep="\t")
    queen = pd.read_csv(queen_dir / "10_patient_program_spatial_metrics.tsv", sep="\t")
    queen_tests = pd.read_csv(queen_dir / "11_chip_adjusted_patient_disease_tests.tsv", sep="\t")
    require(len(score) == 152640 and score["bin_id"].is_unique, "Frozen score key failed")
    require(len(queen) == 105 and len(queen_tests) == 75, "Queen input dimensions failed")

    topology_by_sample: dict[str, h5py.Group] = {}
    rook_rows: list[dict[str, object]] = []
    with h5py.File(topology_path, "r") as handle:
        stored_ids = np.asarray(handle["bin_id"].asstr()[:])
        global_id_to_index = {value: index for index, value in enumerate(stored_ids)}
        require(len(global_id_to_index) == len(stored_ids), "Topology Bin IDs are duplicated")
        for sample_i, sample in enumerate(queen["sample_id"].drop_duplicates()):
            candidates = [
                handle[f"topology/{key}"]
                for key in handle["topology"].keys()
                if str(handle[f"topology/{key}"].attrs["sample_id"]) == sample
            ]
            require(len(candidates) == 1, f"Topology group mismatch: {sample}")
            group = candidates[0]
            part = score.loc[score["sample_id"].eq(sample)].copy()
            expected_global = np.asarray(group["g15_global_indices"][:], dtype=int)
            observed_global = part["bin_id"].map(global_id_to_index).to_numpy(int)
            require(
                np.array_equal(observed_global, expected_global),
                f"Frozen score/topology order mismatch: {sample}",
            )
            coords = part[["x", "y"]].to_numpy(float)
            components = attach_components(group, len(part))
            u, v = rook_edges(coords, components)
            require(len(u) > 0, f"No rook edges: {sample}")
            info = queen.loc[queen["sample_id"].eq(sample)].iloc[0]
            for program_i, program in enumerate(PROGRAMS):
                z = part[f"z__{program}"].to_numpy(float)
                high = z >= HIGH_Z
                sizes = component_sizes(high, u, v)
                hotspots = [size for size in sizes if size >= MIN_HOTSPOT_BINS]
                observed = morans_i(z, u, v)
                rng = np.random.default_rng(SEED + sample_i * 100 + program_i)
                null = [
                    morans_i(permute_within_components(z, components, rng), u, v)
                    for _ in range(PERMUTATIONS)
                ]
                rook_rows.append(
                    {
                        "sample_id": sample,
                        "patient_id": info["patient_id"],
                        "chip_id": info["chip_id"],
                        "disease_group": info["disease_group"],
                        "program": program,
                        "n_bins": len(part),
                        "rook_edge_n": len(u),
                        "rook_step_native": STEP,
                        "high_z_threshold_unchanged": HIGH_Z,
                        "minimum_hotspot_bins_unchanged": MIN_HOTSPOT_BINS,
                        "morans_i": observed,
                        "morans_empirical_P": empirical_p(observed, null),
                        "hotspot_count": len(hotspots),
                        "largest_hotspot_bin_n": int(hotspots[0]) if hotspots else 0,
                        "largest_hotspot_area_mm2": (
                            float(hotspots[0] * BIN_AREA_MM2) if hotspots else 0.0
                        ),
                    }
                )
                done = sample_i * len(PROGRAMS) + program_i + 1
                print(f"PROGRESS {done}/105 sample={sample} program={program}", flush=True)

    rook = pd.DataFrame(rook_rows)
    rook["morans_bh_q_105"] = multipletests(rook["morans_empirical_P"], method="fdr_bh")[1]
    rook.to_csv(output / "01_rook_patient_spatial_metrics.tsv", sep="\t", index=False)

    keys = ["sample_id", "patient_id", "chip_id", "disease_group", "program", "n_bins"]
    metric_names = [
        "morans_i",
        "morans_empirical_P",
        "morans_bh_q_105",
        "hotspot_count",
        "largest_hotspot_bin_n",
        "largest_hotspot_area_mm2",
    ]
    paired = queen[keys + metric_names].merge(
        rook[keys + metric_names + ["rook_edge_n"]],
        on=keys,
        how="outer",
        validate="one_to_one",
        suffixes=("_queen", "_rook"),
        indicator=True,
    )
    require(paired["_merge"].eq("both").all(), "Queen/rook patient keys differ")
    paired = paired.drop(columns="_merge")
    for metric in ["morans_i", "hotspot_count", "largest_hotspot_area_mm2"]:
        paired[f"delta_rook_minus_queen__{metric}"] = (
            paired[f"{metric}_rook"] - paired[f"{metric}_queen"]
        )
    paired["moran_significance_changed"] = (
        paired["morans_bh_q_105_queen"].lt(0.05)
        != paired["morans_bh_q_105_rook"].lt(0.05)
    )
    paired["hotspot_presence_changed"] = (
        paired["hotspot_count_queen"].gt(0) != paired["hotspot_count_rook"].gt(0)
    )
    paired.to_csv(output / "02_queen_rook_patient_metric_comparison.tsv", sep="\t", index=False)

    rook_patient = queen.copy()
    replacement = rook.set_index(["sample_id", "program"])
    for row_index, row in rook_patient.iterrows():
        values = replacement.loc[(row["sample_id"], row["program"])]
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

    rook_tests = fit_all_75(rook_patient)
    adjacency_outcomes = ["morans_i", "largest_hotspot_area_mm2", "hotspot_count"]
    queen_adj = queen_tests.loc[queen_tests["outcome"].isin(adjacency_outcomes)].copy()
    rook_adj = rook_tests.loc[rook_tests["outcome"].isin(adjacency_outcomes)].copy()
    test_keys = ["program", "outcome", "contrast"]
    test_compare = queen_adj.merge(
        rook_adj,
        on=test_keys,
        validate="one_to_one",
        suffixes=("_queen", "_rook"),
    )
    test_compare["effect_direction_changed"] = (
        np.sign(test_compare["effect_model_scale_queen"])
        != np.sign(test_compare["effect_model_scale_rook"])
    )
    test_compare["FDR05_conclusion_changed"] = (
        test_compare["BH_FDR_all_75_tests_queen"].lt(0.05)
        != test_compare["BH_FDR_all_75_tests_rook"].lt(0.05)
    )
    test_compare.to_csv(
        output / "03_queen_rook_adjacency_disease_test_comparison.tsv",
        sep="\t",
        index=False,
    )

    queen_summary = spatial_summary(queen, "queen")
    rook_summary = spatial_summary(rook_patient, "rook")
    summary = queen_summary.merge(rook_summary, on="program", validate="one_to_one")
    summary["spatial_qc_conclusion_changed"] = (
        summary["queen_spatial_qc_pass"] != summary["rook_spatial_qc_pass"]
    )
    summary.to_csv(output / "04_program_spatial_reproducibility_comparison.tsv", sep="\t", index=False)

    inventory = figure_inventory(coord_dir / "figures")
    inventory.to_csv(output / "05_existing_coordinate_validation_figure_inventory.tsv", sep="\t", index=False)

    coord_status = pd.read_csv(coord_dir / "11_final_qc_status.tsv", sep="\t")
    bin_pass = bool(
        coord_status.loc[coord_status["qc_domain"].eq("BIN-SCORE ALIGNMENT"), "status"].iloc[0]
        == "PASS"
    )
    orientation_pass = bool(
        coord_status.loc[coord_status["qc_domain"].eq("COORDINATE ORIENTATION"), "status"].iloc[0]
        == "PASS"
    )
    disease_conclusion_same = not bool(test_compare["FDR05_conclusion_changed"].any())
    reproducibility_same = not bool(summary["spatial_qc_conclusion_changed"].any())
    adjacency_status = (
        "PASS" if bin_pass and orientation_pass and disease_conclusion_same and reproducibility_same else "CAUTION"
    )
    status = {
        "BIN-SCORE ALIGNMENT": "PASS" if bin_pass else "FAIL",
        "COORDINATE ORIENTATION": "PASS" if orientation_pass else "FAIL",
        "SPATIAL ADJACENCY": adjacency_status,
        "queen_to_rook_program_spatial_qc_conclusions_unchanged": reproducibility_same,
        "queen_to_rook_patient_disease_FDR05_conclusions_unchanged": disease_conclusion_same,
        "spatial_interpretation_allowed": adjacency_status == "PASS",
        "J2_orientation": "identity; no mirror/rotation/scale; consistent with prior validated J2 tissue/mask orientation",
    }
    write_json(output / "06_final_qc_status.json", status)
    pd.DataFrame(
        [{"qc_domain": key, "status": value} for key, value in status.items()]
    ).to_csv(output / "06_final_qc_status.tsv", sep="\t", index=False)

    fig, axes = plt.subplots(1, 3, figsize=(15, 5))
    for axis, metric in zip(
        axes,
        ["morans_i", "hotspot_count", "largest_hotspot_area_mm2"],
        strict=True,
    ):
        x = paired[f"{metric}_queen"].to_numpy(float)
        y = paired[f"{metric}_rook"].to_numpy(float)
        low = min(float(x.min()), float(y.min()))
        high = max(float(x.max()), float(y.max()))
        axis.scatter(x, y, s=16, alpha=0.75)
        axis.plot([low, high], [low, high], linestyle="--", color="black", linewidth=1)
        axis.set_xlabel("queen")
        axis.set_ylabel("rook")
        axis.set_title(metric)
    fig.suptitle("Paired patient/program adjacency sensitivity (n=105)")
    fig.tight_layout()
    fig.savefig(output / "figures" / "queen_rook_patient_metric_comparison.png", dpi=220)
    fig.savefig(output / "figures" / "queen_rook_patient_metric_comparison.svg")
    plt.close(fig)

    max_deltas = {
        metric: float(paired[f"delta_rook_minus_queen__{metric}"].abs().max())
        for metric in ["morans_i", "hotspot_count", "largest_hotspot_area_mm2"]
    }
    report = f"""# Rook adjacency sensitivity and coordinate-figure audit

## Scope

Only the existing five frozen scores were used. The high-score threshold
(`z >= {HIGH_Z}`), minimum hotspot size ({MIN_HOTSPOT_BINS} Bin50), Bin area
({BIN_AREA_MM2} mm2), source-component constraint, 999 Moran permutations,
patient-level chip-adjusted models, and 75-test BH scope were unchanged. No
score, panel, GSEA, differential-expression result, or other spatial endpoint
was recomputed.

## Queen-to-rook conclusion

- Program spatial-QC conclusions unchanged: **{reproducibility_same}**
- Patient disease FDR<0.05 conclusions unchanged: **{disease_conclusion_same}**
- Maximum absolute paired changes: `{json.dumps(max_deltas)}`
- SPATIAL ADJACENCY: **{adjacency_status}**

## Coordinate validation figures

The existing all-21 contact sheet contains authoritative H5AD outlines,
current keyed coordinates, and overlays. K8/J2/L3 source-vs-current maps use
SFTPC, COL1A1, and PTPRC with identical coordinate limits, point sizes, and
color scales. J2 uses the identity orientation (no mirror, rotation, scaling,
or affine transform) and agrees with the previously validated J2 tissue/mask
orientation.
"""
    (output / "20_ROOK_ADJACENCY_SENSITIVITY_REPORT.md").write_text(report, encoding="utf-8")

    inputs = [score_path, queen_dir / "10_patient_program_spatial_metrics.tsv", queen_dir / "11_chip_adjusted_patient_disease_tests.tsv", topology_path, coord_dir / "11_final_qc_status.tsv"]
    pd.DataFrame(
        [
            {"path": str(path), "size_bytes": path.stat().st_size, "sha256": sha256_file(path)}
            for path in inputs
        ]
    ).to_csv(output / "00_input_provenance.tsv", sep="\t", index=False)
    write_json(
        output / "parameters_manifest.json",
        {
            "seed": SEED,
            "permutations": PERMUTATIONS,
            "rook_half_offsets": ROOK_HALF_OFFSETS,
            "native_step": STEP,
            "high_z_threshold_unchanged": HIGH_Z,
            "minimum_hotspot_bins_unchanged": MIN_HOTSPOT_BINS,
            "bin_area_mm2_unchanged": BIN_AREA_MM2,
            "score_recomputed": False,
            "gene_panel_modified": False,
            "gsea_or_de_run": False,
            "other_spatial_analysis_run": False,
            "elapsed_seconds": time.time() - started,
        },
    )
    write_json(
        output / "software_versions.json",
        {
            "python": sys.version,
            "platform": platform.platform(),
            "numpy": np.__version__,
            "pandas": pd.__version__,
        },
    )
    shutil.copy2(Path(__file__), output / Path(__file__).name)
    (output / "SUCCESS").write_text(datetime.now().isoformat(), encoding="utf-8")
    print(f"OUTPUT={output}", flush=True)


if __name__ == "__main__":
    main()
