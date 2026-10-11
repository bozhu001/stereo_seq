#!/usr/bin/env python3
"""Aggregate frozen K3 neighborhood composition by patient and niche.

The source HDF5 contains the same 18-state neighborhood-composition vectors
used to derive the disease-blind shared K3 atlas. This script only averages
those frozen vectors within each patient/niche; it does not rerun RCTD or K3.
"""

from __future__ import annotations

import hashlib
from pathlib import Path

import h5py
import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
ATLAS = (
    ROOT
    / "results/reanalysis/ALL21_SHARED_UNBIASED_DISCRETE_NICHE_ATLAS_20261011_191652_attempt4"
)
OUT = (
    ROOT
    / "results/reanalysis/ALL21_K3_NICHE_PSEUDOBULK_DEG_20261011_201008"
)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def decode(values: np.ndarray) -> list[str]:
    return [value.decode() if isinstance(value, bytes) else str(value) for value in values]


def main() -> None:
    assignment_path = ATLAS / "02_ALL21_SHARED_NICHE_ASSIGNMENTS.tsv.gz"
    assay_path = ATLAS / "12_SHARED_NICHE_ASSAY.h5"
    sample_path = OUT / "01_PSEUDOBULK_SAMPLE_METADATA.tsv"
    assignments = pd.read_csv(assignment_path, sep="\t")
    samples = pd.read_csv(sample_path, sep="\t")
    with h5py.File(assay_path, "r") as handle:
        composition = handle["composition"][:]
        state_names = decode(handle["state_names"][:])
        entity_keys = decode(handle["entity_key"][:])

    if entity_keys != assignments["entity_key"].astype(str).tolist():
        raise RuntimeError("Frozen assay and assignment order differ")
    keep = assignments["niche"].isin(["CN1", "CN2", "CN3"]).to_numpy()
    if not np.allclose(composition[keep].sum(axis=1), 1.0, atol=1e-5):
        raise RuntimeError("Evaluable composition rows do not sum to one")
    meta = assignments.loc[
        keep, ["patient_id", "chip_label", "disease", "niche"]
    ].reset_index(drop=True)
    values = pd.DataFrame(composition[keep], columns=state_names)
    frame = pd.concat([meta, values], axis=1)
    grouped = (
        frame.groupby(
            ["patient_id", "chip_label", "disease", "niche"], observed=True
        )[state_names]
        .mean()
        .reset_index()
    )
    counts = (
        frame.groupby(
            ["patient_id", "chip_label", "disease", "niche"], observed=True
        )
        .size()
        .rename("bin_n")
        .reset_index()
    )
    grouped = grouped.merge(
        counts,
        on=["patient_id", "chip_label", "disease", "niche"],
        validate="one_to_one",
    )
    grouped.insert(
        0,
        "pseudobulk_id",
        grouped["patient_id"].astype(str) + "::" + grouped["niche"],
    )
    grouped["airway_epithelial_sum"] = grouped[
        [
            "Ciliated epithelial",
            "Secretory/other airway epithelial",
            "Basal-like epithelial",
        ]
    ].sum(axis=1)
    grouped["ciliated_of_airway"] = grouped["Ciliated epithelial"] / grouped[
        "airway_epithelial_sum"
    ].replace(0, np.nan)
    grouped = samples[["pseudobulk_id", "low_bin_coverage_lt20"]].merge(
        grouped, on="pseudobulk_id", validate="one_to_one"
    )
    if len(grouped) != 63 or grouped["patient_id"].nunique() != 21:
        raise RuntimeError("Expected 63 patient/niche rows and 21 patients")
    if not (grouped[state_names].sum(axis=1).sub(1).abs() < 1e-5).all():
        raise RuntimeError("Aggregated composition rows do not sum to one")

    grouped.to_csv(
        OUT / "19_PATIENT_NICHE_RCTD_COMPOSITION.tsv", sep="\t", index=False
    )
    manifest = pd.DataFrame(
        [
            {
                "role": "frozen_k3_assignments",
                "path": str(assignment_path.relative_to(ROOT)),
                "sha256": sha256(assignment_path),
            },
            {
                "role": "frozen_k3_neighborhood_composition",
                "path": str(assay_path.relative_to(ROOT)),
                "sha256": sha256(assay_path),
            },
            {
                "role": "patient_niche_pseudobulk_metadata",
                "path": str(sample_path.relative_to(ROOT)),
                "sha256": sha256(sample_path),
            },
        ]
    )
    manifest.to_csv(OUT / "20_COMPOSITION_INPUT_MANIFEST.tsv", sep="\t", index=False)
    validation = pd.DataFrame(
        [
            ("patient_niche_rows", len(grouped), 63, "PASS"),
            ("patients", grouped["patient_id"].nunique(), 21, "PASS"),
            ("niches", grouped["niche"].nunique(), 3, "PASS"),
            ("reference18_states", len(state_names), 18, "PASS"),
            ("evaluable_bins", int(grouped["bin_n"].sum()), 126_199, "PASS"),
            (
                "composition_sum_max_abs_error",
                grouped[state_names].sum(axis=1).sub(1).abs().max(),
                "<1e-5",
                "PASS",
            ),
        ],
        columns=["check", "observed", "expected", "status"],
    )
    validation.to_csv(OUT / "21_COMPOSITION_PREPARATION_VALIDATION.tsv", sep="\t", index=False)
    print(OUT)


if __name__ == "__main__":
    main()
