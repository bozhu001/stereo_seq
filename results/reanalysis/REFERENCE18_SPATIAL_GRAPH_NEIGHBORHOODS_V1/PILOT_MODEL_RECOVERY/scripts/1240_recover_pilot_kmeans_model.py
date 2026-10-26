"""Deterministically recover and verify the frozen PN1-PN8 pilot KMeans model.

The original final script is imported as the authoritative implementation for
input loading, representative selection, graph construction, and sampling.
Only K=8 is fitted; K selection is not rerun. No 21-patient projection or
upstream analysis is performed.
"""

from __future__ import annotations

import hashlib
import importlib.util
import json
import logging
import platform
import sys
from datetime import datetime
from pathlib import Path

import joblib
import numpy as np
import pandas as pd
import scipy
import sklearn
from scipy.sparse import coo_matrix
from sklearn.cluster import KMeans
from sklearn.metrics import adjusted_rand_score, normalized_mutual_info_score


SEED = 20261025
K = 8
PER_PATIENT_TRAINING_N = 2500


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def logger_for(path: Path) -> logging.Logger:
    logger = logging.getLogger("pilot_model_recovery")
    logger.setLevel(logging.INFO)
    logger.handlers.clear()
    formatter = logging.Formatter("%(asctime)s\t%(levelname)s\t%(message)s")
    for handler in (logging.FileHandler(path, encoding="utf-8"), logging.StreamHandler()):
        handler.setFormatter(formatter)
        logger.addHandler(handler)
    return logger


def import_original(path: Path):
    spec = importlib.util.spec_from_file_location("frozen_pilot_1233", path)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"Cannot import original pilot script: {path}")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def compute_features_and_parameters(
    pilot: pd.DataFrame,
    weight_columns: list[str],
    pairs: np.ndarray,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    weights = pilot[weight_columns].to_numpy(float)
    adjacency = coo_matrix(
        (
            np.ones(2 * len(pairs)),
            (
                np.r_[pairs[:, 0], pairs[:, 1]],
                np.r_[pairs[:, 1], pairs[:, 0]],
            ),
        ),
        shape=(len(pilot), len(pilot)),
    ).tocsr()
    degree = np.asarray(adjacency.sum(axis=1)).ravel().astype(int)
    neighbor_sum = adjacency @ weights
    ego_mean = (neighbor_sum + weights) / (degree[:, None] + 1)
    hellinger = np.sqrt(ego_mean)
    patient_centers = np.zeros_like(hellinger)
    centered = hellinger.copy()
    for _, group in pilot.groupby("patient_id"):
        indices = group.index.to_numpy(int)
        center = hellinger[indices].mean(axis=0)
        centered[indices] -= center
        patient_centers[indices] = center
    raw_scale = centered.std(axis=0, ddof=1)
    divisor = np.where(raw_scale > 1e-8, raw_scale, 1.0)
    standardized = centered / divisor
    return standardized, degree, raw_scale, patient_centers


def provenance_row(item: str, expected: str, observed: str, status: str, evidence: str) -> dict[str, str]:
    return {
        "provenance_item": item,
        "expected": expected,
        "observed": observed,
        "status": status,
        "evidence": evidence,
    }


def main() -> None:
    script_path = Path(__file__).resolve()
    output = script_path.parents[1]
    project = script_path.parents[5]
    pilot_root = project / "results/reanalysis/REFERENCE18_SPATIAL_GRAPH_NEIGHBORHOODS_V1"
    original_script = project / "scripts/1233_reference18_spatial_graph_neighborhoods_pilot.py"
    archived_script = pilot_root / "scripts/1233_reference18_spatial_graph_neighborhoods_pilot.py"
    old_assignment_path = pilot_root / "tables/07_PILOT_NEIGHBORHOOD_ASSIGNMENTS.tsv.gz"
    old_edge_path = pilot_root / "tables/05_PILOT_DELAUNAY_EDGES.tsv.gz"
    representative_path = pilot_root / "tables/02_PILOT_REPRESENTATIVE_PATIENTS.tsv"
    k_selection_path = pilot_root / "tables/06_PILOT_K_SELECTION.tsv"
    models_dir = output / "models"
    models_dir.mkdir(parents=True, exist_ok=True)
    log_path = output / "logs/1240_recover_pilot_kmeans_model.log"
    logger = logger_for(log_path)
    logger.info("Starting deterministic K=8 pilot model recovery")

    original_hash = sha256(original_script)
    archived_hash = sha256(archived_script)
    if original_hash != archived_hash:
        raise RuntimeError("Working and archived final pilot scripts differ")
    module = import_original(original_script)
    if module.SEED != SEED:
        raise RuntimeError(f"Unexpected original seed: {module.SEED}")

    data, weight_columns, _ = module.load_frozen_inputs(project)
    input_checks = module.validate_inputs(data, weight_columns)
    if any(row["status"] != "PASS" for row in input_checks):
        raise RuntimeError("Frozen input validation failed")
    representatives = module.select_representatives(data)
    old_representatives = pd.read_csv(representative_path, sep="\t")
    representative_columns = ["chip_label", "disease", "patient_id", "bin_n"]
    representatives_exact = representatives[representative_columns].reset_index(drop=True).equals(
        old_representatives[representative_columns].reset_index(drop=True)
    )
    if not representatives_exact:
        raise RuntimeError("Representative patient selection is not identical")

    pilot, _, edges = module.build_patient_graphs(data, representatives)
    generated_pairs = edges["delaunay"]
    old_assignment = pd.read_csv(old_assignment_path, sep="\t")
    entity_order_exact = pilot.entity_key.astype(str).tolist() == old_assignment.entity_key.astype(str).tolist()
    if not entity_order_exact:
        raise RuntimeError("Pilot Bin order differs from the final assignment artifact")

    generated_edge_keys = {
        tuple(sorted(pair))
        for pair in zip(
            pilot.entity_key.to_numpy(str)[generated_pairs[:, 0]],
            pilot.entity_key.to_numpy(str)[generated_pairs[:, 1]],
        )
    }
    old_edges = pd.read_csv(old_edge_path, sep="\t")
    old_edge_keys = {
        tuple(sorted(pair))
        for pair in old_edges[["left_entity_key", "right_entity_key"]].itertuples(index=False, name=None)
    }
    edge_set_exact = generated_edge_keys == old_edge_keys
    if not edge_set_exact:
        raise RuntimeError("Reconstructed Delaunay edge set differs from the final artifact")

    feature, degree, raw_scale, patient_centers_per_bin = compute_features_and_parameters(
        pilot, weight_columns, generated_pairs
    )
    original_feature, _, _, original_degree = module.local_features(
        pilot, weight_columns, generated_pairs
    )
    feature_exact = np.array_equal(feature, original_feature)
    feature_max_abs_difference = float(np.max(np.abs(feature - original_feature)))
    degree_exact = np.array_equal(degree, original_degree) and np.array_equal(
        degree, old_assignment.delaunay_degree.to_numpy(int)
    )
    if not feature_exact or not degree_exact:
        raise RuntimeError("Recovered feature matrix or graph degree does not exactly match original code")

    eligible = np.flatnonzero(degree >= 2)
    rng = np.random.default_rng(SEED)
    training = module.balanced_sample_indices(
        pilot, eligible, PER_PATIENT_TRAINING_N, rng
    )
    # The original script drew a validation subset here. It did not affect the
    # KMeans fit because each model had its own explicit random_state.
    _ = rng.choice(eligible, min(6000, len(eligible)), replace=False)

    logger.info("Fitting K=8 only on %s balanced Pilot bins", len(training))
    model = KMeans(
        n_clusters=K,
        random_state=SEED,
        n_init=10,
        max_iter=500,
        algorithm="lloyd",
    ).fit(feature[training])
    recovered_zero_based = model.predict(feature[eligible])
    recovered_labels = np.array([f"PN{label + 1}" for label in recovered_zero_based])
    old_labels = old_assignment.pilot_neighborhood.to_numpy(str)[eligible]

    exact_matches = recovered_labels == old_labels
    agreement = float(exact_matches.mean())
    ari = float(adjusted_rand_score(old_labels, recovered_labels))
    nmi = float(normalized_mutual_info_score(old_labels, recovered_labels))
    exact_recovered = bool(np.all(exact_matches))
    decision = "MODEL_RECOVERED_EXACT" if exact_recovered else "MODEL_NOT_EXACTLY_RECOVERED"
    logger.info("Label agreement=%.12f ARI=%.12f NMI=%.12f", agreement, ari, nmi)

    states = [column.removeprefix("weight__") for column in weight_columns]
    selected_k = pd.read_csv(k_selection_path, sep="\t").query("pilot_selected == True")
    formal_serialized = [
        path for path in pilot_root.rglob("*")
        if path.is_file() and path.suffix.lower() in {".pkl", ".pickle", ".joblib", ".npz", ".npy", ".sav", ".rds", ".qs"}
    ]

    provenance = [
        provenance_row("original_script_sha256", original_hash, archived_hash, "PASS" if original_hash == archived_hash else "FAIL", str(original_script.relative_to(project))),
        provenance_row("formal_preexisting_serialized_model", "recoverable original fitted object", f"candidate_n={len(formal_serialized)}", "NOT_FOUND" if not formal_serialized else "FOUND", "Search restricted to the authoritative analysis tree; no original model was found before refitting."),
        provenance_row("frozen_input_rows", "126449", str(len(data)), "PASS" if len(data) == 126449 else "FAIL", "Current K8/J2/L3 frozen Reference18 inputs."),
        provenance_row("pilot_patient_n", "9", str(representatives.patient_id.nunique()), "PASS" if representatives.patient_id.nunique() == 9 else "FAIL", "Representative patient selection rerun from the final script."),
        provenance_row("pilot_bin_order", f"{len(old_assignment)} rows in exact entity_key order", f"rows={len(pilot)}; exact={entity_order_exact}", "PASS" if entity_order_exact else "FAIL", "Compared with 07_PILOT_NEIGHBORHOOD_ASSIGNMENTS.tsv.gz."),
        provenance_row("delaunay_edge_set", f"{len(old_edge_keys)} unique stored edges", f"generated={len(generated_edge_keys)}; exact={edge_set_exact}", "PASS" if edge_set_exact else "FAIL", "Rebuilt with final source, within patient/ROI/tissue component and maximum 25*sqrt(2) um."),
        provenance_row("reference18_state_order", "18 frozen columns in stored order", "|".join(states), "PASS" if len(states) == 18 else "FAIL", "Column order read directly from frozen inputs."),
        provenance_row("ego_and_hellinger", "ego=(self+neighbor sum)/(degree+1); sqrt", "exact final source implementation", "PASS", "Feature calculation duplicated and compared exactly to original local_features()."),
        provenance_row("patient_centering", "subtract unlabeled within-patient Hellinger mean", "9 patient-specific 18-vectors recovered", "PASS", "No PN labels used for centering."),
        provenance_row("feature_standardization", "pooled Pilot ddof=1 SD; zero guard 1e-8", f"18 scales; exact_feature_matrix={feature_exact}; max_abs_diff={feature_max_abs_difference:.17g}", "PASS" if feature_exact else "FAIL", "Compared elementwise with original local_features()."),
        provenance_row("eligible_bin_rule", "Delaunay degree >= 2", f"eligible={len(eligible)}; not_evaluable={len(pilot)-len(eligible)}; degree_exact={degree_exact}", "PASS" if len(eligible) == 50972 and degree_exact else "FAIL", "Compared with final assignment graph degrees."),
        provenance_row("balanced_training_bins", "up to 2500 per patient using default_rng(20261025)", f"training_n={len(training)}", "PASS", "Indices regenerated before the original validation draw."),
        provenance_row("k", "8", str(K), "PASS" if len(selected_k) == 1 and int(selected_k.k.iloc[0]) == K else "FAIL", "K selection was not rerun."),
        provenance_row("random_seed", "20261025", str(SEED), "PASS", "Used for balanced sample and KMeans random_state."),
        provenance_row("kmeans_configuration", "n_init=10; max_iter=500; algorithm=lloyd", f"n_init={model.n_init}; max_iter={model.max_iter}; algorithm={model.algorithm}", "PASS", f"scikit-learn={sklearn.__version__}"),
        provenance_row("final_recovery_decision", "100% exact agreement required", decision, "PASS" if exact_recovered else "FAIL", "No label permutation or post-hoc remapping was used."),
    ]
    pd.DataFrame(provenance).to_csv(output / "01_MODEL_PROVENANCE.tsv", sep="\t", index=False)

    parameter_rows: list[dict[str, object]] = []
    for state_index, state in enumerate(states):
        parameter_rows.append({"parameter_type": "feature_raw_sd_ddof1", "pn": "", "patient_id": "", "reference18_state": state, "value": raw_scale[state_index]})
        parameter_rows.append({"parameter_type": "feature_scale_divisor", "pn": "", "patient_id": "", "reference18_state": state, "value": raw_scale[state_index] if raw_scale[state_index] > 1e-8 else 1.0})
    patient_centers = pd.DataFrame(patient_centers_per_bin, columns=states).assign(patient_id=pilot.patient_id).groupby("patient_id", sort=True)[states].first()
    for patient_id, row in patient_centers.iterrows():
        for state in states:
            parameter_rows.append({"parameter_type": "pilot_patient_hellinger_center", "pn": "", "patient_id": patient_id, "reference18_state": state, "value": row[state]})
    for cluster_index, center in enumerate(model.cluster_centers_):
        for state_index, state in enumerate(states):
            parameter_rows.append({"parameter_type": "standardized_cluster_center", "pn": f"PN{cluster_index + 1}", "patient_id": "", "reference18_state": state, "value": center[state_index]})
    for name, value in {
        "k": K,
        "seed": SEED,
        "n_init": model.n_init,
        "max_iter": model.max_iter,
        "n_iter_observed": model.n_iter_,
        "inertia": model.inertia_,
        "training_bin_n": len(training),
        "eligible_bin_n": len(eligible),
    }.items():
        parameter_rows.append({"parameter_type": f"model_{name}", "pn": "", "patient_id": "", "reference18_state": "", "value": value})
    pd.DataFrame(parameter_rows).to_csv(output / "02_RECOVERED_MODEL_PARAMETERS.tsv", sep="\t", index=False)

    reproducibility_rows = [
        {"metric": "exact_label_agreement", "scope": "all_evaluable", "value": agreement, "numerator": int(exact_matches.sum()), "denominator": len(exact_matches)},
        {"metric": "mismatched_bin_n", "scope": "all_evaluable", "value": int((~exact_matches).sum()), "numerator": int((~exact_matches).sum()), "denominator": len(exact_matches)},
        {"metric": "adjusted_rand_index", "scope": "all_evaluable", "value": ari, "numerator": "", "denominator": len(exact_matches)},
        {"metric": "normalized_mutual_information", "scope": "all_evaluable", "value": nmi, "numerator": "", "denominator": len(exact_matches)},
    ]
    for pn in [f"PN{index}" for index in range(1, K + 1)]:
        mask = old_labels == pn
        reproducibility_rows.append(
            {
                "metric": "per_old_pn_exact_agreement",
                "scope": pn,
                "value": float((recovered_labels[mask] == pn).mean()),
                "numerator": int((recovered_labels[mask] == pn).sum()),
                "denominator": int(mask.sum()),
            }
        )
    pd.DataFrame(reproducibility_rows).to_csv(output / "03_PILOT_LABEL_REPRODUCIBILITY.tsv", sep="\t", index=False)

    labels_order = [f"PN{index}" for index in range(1, K + 1)]
    confusion = pd.crosstab(
        pd.Categorical(old_labels, categories=labels_order, ordered=True),
        pd.Categorical(recovered_labels, categories=labels_order, ordered=True),
        dropna=False,
    )
    confusion.index.name = "old_pilot_label"
    confusion.columns = [f"recovered_{column}" for column in confusion.columns]
    confusion.reset_index().to_csv(output / "04_PN_CONFUSION_MATRIX.tsv", sep="\t", index=False)

    training_table = pilot.loc[training, ["entity_key", "patient_id", "disease", "chip_label", "roi_id", "tissue_component_id", "x_um", "y_um"]].copy()
    training_table.insert(0, "training_rank", np.arange(len(training)))
    training_table.insert(1, "pilot_row", training)
    training_table.to_csv(models_dir / "TRAINING_BIN_INDICES.tsv.gz", sep="\t", index=False, compression="gzip")
    pd.DataFrame({"cluster_index_zero_based": np.arange(K), "pn_label": labels_order}).to_csv(models_dir / "CLUSTER_TO_PN_MAPPING.tsv", sep="\t", index=False)
    pd.DataFrame({"feature_index_zero_based": np.arange(len(states)), "reference18_state": states, "weight_column": weight_columns}).to_csv(models_dir / "REFERENCE18_FEATURE_ORDER.tsv", sep="\t", index=False)
    recovered_assignment = old_assignment.copy()
    recovered_assignment["recovered_pn"] = "NOT_EVALUABLE"
    recovered_assignment.loc[eligible, "recovered_pn"] = recovered_labels
    recovered_assignment["exact_match_to_original"] = recovered_assignment.recovered_pn.eq(recovered_assignment.pilot_neighborhood)
    recovered_assignment.to_csv(models_dir / "RECOVERED_PILOT_LABELS.tsv.gz", sep="\t", index=False, compression="gzip")

    bundle = {
        "model": model,
        "cluster_centers": model.cluster_centers_.copy(),
        "feature_raw_sd_ddof1": raw_scale.copy(),
        "feature_scale_divisor": np.where(raw_scale > 1e-8, raw_scale, 1.0),
        "reference18_states": states,
        "weight_columns": weight_columns,
        "cluster_to_pn": {index: f"PN{index + 1}" for index in range(K)},
        "seed": SEED,
        "k": K,
        "eligible_degree_min": 2,
        "delaunay_max_distance_um": 25 * np.sqrt(2),
        "patient_centering": "subtract new patient's unlabeled Hellinger feature mean",
        "original_script_sha256": original_hash,
        "training_entity_keys": pilot.entity_key.to_numpy(str)[training],
        "decision": decision,
    }
    joblib.dump(bundle, models_dir / "PN1_PN8_RECOVERED_KMEANS_BUNDLE.joblib", compress=3)

    decision_report = f"""# Pilot K-means模型恢复结论

## 结论

**{decision}**

在最终Pilot脚本、冻结Reference18输入、原9位患者、原Delaunay构图规则、原特征转换、K=8及随机种子20261025下，仅重新拟合K=8一次。没有重选K，没有使用21位患者，也没有运行RCTD、NMF、edgeR、GSEA、LR或Marker分析。

## 可重复性结果

- 原Pilot可评估Bin50：{len(eligible):,}
- 新旧标签完全一致：{int(exact_matches.sum()):,}/{len(exact_matches):,}（{agreement:.12%}）
- 不一致Bin：{int((~exact_matches).sum()):,}
- ARI：{ari:.12f}
- NMI：{nmi:.12f}
- Delaunay边集合完全一致：{edge_set_exact}
- 标准化特征矩阵与原函数逐元素完全一致：{feature_exact}；最大绝对差={feature_max_abs_difference:.17g}
- 未进行标签置换或事后PN重命名；cluster 0–7直接映射PN1–PN8。

## 对21例External Projection的判断

{"现在可以严格支持下一阶段21例External Projection。恢复包已保存训练中心、scale、训练Bin索引、Reference18顺序和PN映射；下一阶段仍须明确区分9位in-sample与12位out-of-sample。" if exact_recovered else "当前仍不能严格支持21例External Projection。新旧标签未达到100%一致，恢复模型不得冒充原模型。"}
"""
    (output / "05_MODEL_RECOVERY_DECISION_ZH.md").write_text(decision_report, encoding="utf-8")

    metadata = {
        "decision": decision,
        "timestamp": datetime.now().astimezone().isoformat(),
        "exact_label_agreement": agreement,
        "matching_bin_n": int(exact_matches.sum()),
        "mismatched_bin_n": int((~exact_matches).sum()),
        "evaluable_bin_n": len(exact_matches),
        "ari": ari,
        "nmi": nmi,
        "training_bin_n": len(training),
        "original_script_sha256": original_hash,
        "old_assignment_sha256": sha256(old_assignment_path),
        "old_delaunay_edges_sha256": sha256(old_edge_path),
        "sklearn": sklearn.__version__,
        "numpy": np.__version__,
        "pandas": pd.__version__,
        "scipy": scipy.__version__,
        "python": sys.version,
        "platform": platform.platform(),
    }
    (models_dir / "MODEL_METADATA.json").write_text(json.dumps(metadata, indent=2, ensure_ascii=False), encoding="utf-8")

    logger.info("Recovery decision: %s", decision)
    for handler in logger.handlers:
        handler.flush()
    manifest_rows = []
    for path in sorted(output.rglob("*")):
        if path.is_file() and path != output / "06_OUTPUT_MANIFEST.tsv":
            manifest_rows.append({"relative_path": str(path.relative_to(output)), "bytes": path.stat().st_size, "sha256": sha256(path)})
    pd.DataFrame(manifest_rows).to_csv(output / "06_OUTPUT_MANIFEST.tsv", sep="\t", index=False)


if __name__ == "__main__":
    main()
