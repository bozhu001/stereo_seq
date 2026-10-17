#!/usr/bin/env python3
"""Map frozen harmonized K=7 Neighborhood-NMF families onto ssDNA.

This is a visualization-only workflow. It reuses the frozen disease-blind K=7
models, their cross-chip family matching, the frozen Bin50 metadata, and the
previously validated chip-specific ssDNA coordinate transforms. It does not
refit RCTD or NMF and does not perform a disease test.

Typical command:
    envs/stereo_seq_py/Scripts/python.exe \
        scripts/1157_render_reference18_k7_ssdna_maps.py
"""

from __future__ import annotations

import argparse
import datetime as dt
import hashlib
import json
import platform
import sys
from contextlib import ExitStack
from pathlib import Path

import matplotlib
import numpy as np
import pandas as pd
from PIL import Image

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.cm import ScalarMappable
from matplotlib.colors import Normalize


ROOT = Path(__file__).resolve().parents[1]
MODEL_ROOT = ROOT / (
    "results/reanalysis/"
    "REFERENCE18_NEIGHBORHOOD_NMF_K7_K12_EXTENSION_20261015_attempt5/models"
)
MAP_PATH = ROOT / (
    "results/reanalysis/"
    "REFERENCE18_NEIGHBORHOOD_NMF_K7_STATE_SENSITIVITY_20261016_210302_attempt5/"
    "07_CROSS_CHIP_FAMILY_MAP.tsv"
)
SPATIAL_ROOT = ROOT / (
    "results/reanalysis/REFERENCE18_UMI50_VOLTRON_PER_CHIP_NICHE_20261005_153036"
)
COORDINATE_AUDIT = ROOT / "results/ppt_qc/ssDNA_coordinate_audit.json"

CHIPS = {
    "K8": {
        "chip_id": "Y40102K8",
        "tif": Path(
            r"D:/Lung ST/chip/chip 1/QC/Y40102K8_20260428_204657/"
            r"Y40102K8_SC_20260428_204657_4.2.2/ssDNA/"
            r"Experiment-64-ApoTome-01-Stitching-02_ORG.tif"
        ),
    },
    "J2": {
        "chip_id": "Y40105J2",
        "tif": Path(
            r"D:/Lung ST/chip/chip 2/QC/Y40105J2_20260520_160413/"
            r"Y40105J2_SC_20260520_160413_4.2.2/ssDNA/Y40105J2.tif"
        ),
    },
    "L3": {
        "chip_id": "Y40105L3",
        "tif": Path(
            r"D:/Lung ST/chip/chip 3/QC/Y40105L3_20260522_205054/"
            r"Y40105L3_SC_20260522_205054_4.2.2/ssDNA/YY40105L3.tif"
        ),
    },
}

CHIP_ORDER = ("K8", "J2", "L3")
DISEASE_ORDER = ("HC", "IPF", "SSc-ILD")
FAMILY_ORDER = (
    "Lymphatic-endothelial-associated",
    "Blood-endothelial-associated",
    "AT1-associated",
    "Perivascular-stromal-associated",
    "B-enriched/mixed-immune",
    "Ciliated-associated",
    "Fibroblast-associated",
)
SOURCE_FAMILY = {
    "LEC-associated": "Lymphatic-endothelial-associated",
    "Blood-endothelial": "Blood-endothelial-associated",
    "AT1-associated": "AT1-associated",
    "Perivascular-stromal": "Perivascular-stromal-associated",
    "B-enriched/mixed-immune": "B-enriched/mixed-immune",
    "Ciliated-associated": "Ciliated-associated",
    "Fibroblast-associated": "Fibroblast-associated",
}
FAMILY_DISPLAY = {
    "Lymphatic-endothelial-associated": "Lymphatic endothelial",
    "Blood-endothelial-associated": "Blood endothelial",
    "AT1-associated": "AT1",
    "Perivascular-stromal-associated": "Perivascular stromal",
    "B-enriched/mixed-immune": "B-enriched / mixed immune",
    "Ciliated-associated": "Ciliated",
    "Fibroblast-associated": "Fibroblast",
}
FAMILY_COLORS = {
    "Lymphatic-endothelial-associated": "#2A9D8F",
    "Blood-endothelial-associated": "#4575B4",
    "AT1-associated": "#F4A261",
    "Perivascular-stromal-associated": "#8E6C8A",
    "B-enriched/mixed-immune": "#D73027",
    "Ciliated-associated": "#00A6CA",
    "Fibroblast-associated": "#E6AB02",
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--output-dir",
        type=Path,
        help="New output directory. Defaults to a dynamically timestamped path.",
    )
    return parser.parse_args()


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def output_path(argument: Path | None) -> Path:
    if argument is not None:
        return argument.resolve()
    stamp = dt.datetime.now().strftime("%y%m%d%H%M%S")
    return ROOT / "results/reanalysis" / f"REFERENCE18_K7_SSDNA_MAPS_{stamp}"


def load_coordinate_transforms() -> dict[str, np.ndarray]:
    audit = json.loads(COORDINATE_AUDIT.read_text(encoding="utf-8"))
    transforms = {}
    for chip, configuration in CHIPS.items():
        raw_to_registered = np.asarray(
            audit["chips"][configuration["chip_id"]]["ipr_transform"][
                "raw_to_registered_bin1_matrix"
            ],
            dtype=float,
        )
        transforms[chip] = np.linalg.inv(raw_to_registered)
    return transforms


def load_spatial_data() -> tuple[pd.DataFrame, dict[str, Path]]:
    family_map = pd.read_csv(MAP_PATH, sep="\t")
    family_map = family_map.loc[family_map["scenario"].eq("main18")].copy()
    family_map["family"] = family_map["family_label"].map(SOURCE_FAMILY)
    if set(family_map["family"]) != set(FAMILY_ORDER):
        raise ValueError("Frozen cross-chip family map is not the expected K=7 map")

    transforms = load_coordinate_transforms()
    frames = []
    inputs: dict[str, Path] = {
        "cross_chip_family_map": MAP_PATH,
        "coordinate_audit": COORDINATE_AUDIT,
    }
    for chip in CHIP_ORDER:
        model_path = MODEL_ROOT / f"neighborhood50_{chip}_K07_MODEL_AND_USAGE.npz"
        metadata_path = SPATIAL_ROOT / chip / "01_FROZEN_REFERENCE18_INPUT.tsv.gz"
        inputs[f"model_{chip}"] = model_path
        inputs[f"metadata_{chip}"] = metadata_path
        inputs[f"ssdna_{chip}"] = CHIPS[chip]["tif"]

        archive = np.load(model_path, allow_pickle=True)
        entity_key = archive["entity_key"].astype(str)
        relative_usage = archive["relative_usage"].astype(float)
        dominant = archive["dominant_factor"].astype(int)
        if relative_usage.shape != (len(entity_key), 7):
            raise ValueError(f"Unexpected K=7 usage shape for {chip}: {relative_usage.shape}")
        if dominant.min() != 1 or dominant.max() != 7:
            raise ValueError(f"Dominant factors are not one-based F01-F07 for {chip}")

        chip_map = family_map.set_index("family")
        factor_to_family = {
            chip_map.loc[family, f"{chip}_factor"]: family for family in FAMILY_ORDER
        }
        factor_names = np.asarray(
            [f"F{index:02d}" for index in range(1, 8)], dtype=object
        )
        factor_to_column = {factor: index for index, factor in enumerate(factor_names)}

        metadata = pd.read_csv(
            metadata_path,
            sep="\t",
            usecols=["entity_key", "patient_id", "disease", "x_raw", "y_raw"],
        ).drop_duplicates("entity_key")
        metadata["disease"] = metadata["disease"].replace(
            {"SSC": "SSc-ILD", "SSc": "SSc-ILD"}
        )
        model = pd.DataFrame(
            {
                "entity_key": entity_key,
                "dominant_factor": factor_names[dominant - 1],
            }
        )
        model["family"] = model["dominant_factor"].map(factor_to_family)
        b_factor = chip_map.loc["B-enriched/mixed-immune", f"{chip}_factor"]
        model["b_mixed_usage"] = relative_usage[:, factor_to_column[b_factor]]
        model = model.merge(metadata, on="entity_key", how="left", validate="one_to_one")
        if model[["family", "patient_id", "disease", "x_raw", "y_raw"]].isna().any().any():
            raise ValueError(f"K=7 model-to-spatial metadata join failed for {chip}")

        homogeneous = np.column_stack(
            [model["x_raw"].to_numpy(float), model["y_raw"].to_numpy(float), np.ones(len(model))]
        )
        raw = (transforms[chip] @ homogeneous.T).T
        model["ssdna_x"] = raw[:, 0]
        model["ssdna_y"] = raw[:, 1]
        model["chip"] = chip
        frames.append(model)

    return pd.concat(frames, ignore_index=True), inputs


def panel_entries(data: pd.DataFrame) -> dict[str, list[dict[str, object]]]:
    entries = {disease: [] for disease in DISEASE_ORDER}
    for disease in DISEASE_ORDER:
        disease_data = data.loc[data["disease"].eq(disease)]
        pairs = disease_data[["chip", "patient_id"]].drop_duplicates()
        pairs["chip"] = pd.Categorical(pairs["chip"], CHIP_ORDER, ordered=True)
        pairs = pairs.sort_values(["chip", "patient_id"])
        for row in pairs.itertuples(index=False):
            subset = disease_data.loc[
                disease_data["chip"].eq(str(row.chip))
                & disease_data["patient_id"].eq(row.patient_id)
            ].copy()
            entries[disease].append(
                {"chip": str(row.chip), "patient": row.patient_id, "data": subset}
            )
    return entries


def crop_bounds(frame: pd.DataFrame, image: Image.Image) -> tuple[int, int, int, int]:
    padding = 160
    xmin = max(0, int(np.floor(frame["ssdna_x"].min() - padding)))
    xmax = min(image.width, int(np.ceil(frame["ssdna_x"].max() + padding)))
    ymin = max(0, int(np.floor(frame["ssdna_y"].min() - padding)))
    ymax = min(image.height, int(np.ceil(frame["ssdna_y"].max() + padding)))
    if xmax <= xmin or ymax <= ymin:
        raise ValueError("Invalid ssDNA crop bounds")
    return xmin, ymin, xmax, ymax


def configure_axis(axis: plt.Axes, entry: dict[str, object]) -> None:
    axis.set_axis_off()
    patient = str(entry["patient"])
    short_patient = patient.split("/", 1)[-1]
    axis.set_title(f"{entry['chip']} | {short_patient}", fontsize=8.2, pad=3)


def render_dominant_map(
    entries: dict[str, list[dict[str, object]]],
    images: dict[str, Image.Image],
    output: Path,
) -> None:
    figure, axes = plt.subplots(3, 9, figsize=(30, 11.2), squeeze=False)
    for row_index, disease in enumerate(DISEASE_ORDER):
        for column_index, axis in enumerate(axes[row_index]):
            axis.set_axis_off()
            if column_index >= len(entries[disease]):
                continue
            entry = entries[disease][column_index]
            frame = entry["data"]
            image = images[str(entry["chip"])]
            xmin, ymin, xmax, ymax = crop_bounds(frame, image)
            crop = np.asarray(image.crop((xmin, ymin, xmax, ymax)).convert("L"))
            axis.imshow(crop, cmap="gray", vmin=20, vmax=225, alpha=0.72, origin="upper")
            for family in FAMILY_ORDER:
                subset = frame.loc[frame["family"].eq(family)]
                if subset.empty:
                    continue
                axis.scatter(
                    subset["ssdna_x"] - xmin,
                    subset["ssdna_y"] - ymin,
                    s=5.2,
                    c=FAMILY_COLORS[family],
                    alpha=0.88,
                    linewidths=0,
                    rasterized=True,
                )
            configure_axis(axis, entry)
        axes[row_index, 0].text(
            -0.18,
            0.5,
            disease,
            transform=axes[row_index, 0].transAxes,
            rotation=90,
            ha="center",
            va="center",
            fontsize=13,
            fontweight="bold",
        )

    legend = [
        Line2D(
            [0],
            [0],
            marker="o",
            linestyle="",
            markerfacecolor=FAMILY_COLORS[family],
            markeredgecolor="none",
            markersize=7,
            label=FAMILY_DISPLAY[family],
        )
        for family in FAMILY_ORDER
    ]
    figure.legend(
        handles=legend,
        loc="lower center",
        ncol=4,
        frameon=False,
        bbox_to_anchor=(0.5, 0.014),
        fontsize=9,
    )
    figure.suptitle(
        "Harmonized K=7 neighborhood-factor map across all 21 lung samples",
        fontsize=17,
        fontweight="bold",
    )
    figure.text(
        0.5,
        0.935,
        "Each Bin50 is colored by its dominant frozen K=7 factor; grayscale background = ssDNA.",
        ha="center",
        fontsize=9.5,
        color="#4D5963",
    )
    figure.subplots_adjust(
        left=0.043, right=0.995, top=0.90, bottom=0.115, wspace=0.035, hspace=0.24
    )
    for extension in ("png", "pdf"):
        figure.savefig(
            output / f"FIG01_ALL21_K7_DOMINANT_FACTORS_ON_SSDNA.{extension}",
            dpi=260,
            bbox_inches="tight",
        )
    plt.close(figure)


def render_b_mixed_map(
    entries: dict[str, list[dict[str, object]]],
    images: dict[str, Image.Image],
    output: Path,
) -> None:
    figure, axes = plt.subplots(3, 9, figsize=(30, 11.2), squeeze=False)
    color_map = plt.get_cmap("magma")
    normalization = Normalize(vmin=0, vmax=1)
    for row_index, disease in enumerate(DISEASE_ORDER):
        for column_index, axis in enumerate(axes[row_index]):
            axis.set_axis_off()
            if column_index >= len(entries[disease]):
                continue
            entry = entries[disease][column_index]
            frame = entry["data"]
            image = images[str(entry["chip"])]
            xmin, ymin, xmax, ymax = crop_bounds(frame, image)
            crop = np.asarray(image.crop((xmin, ymin, xmax, ymax)).convert("L"))
            axis.imshow(crop, cmap="gray", vmin=20, vmax=225, alpha=0.82, origin="upper")
            axis.scatter(
                frame["ssdna_x"] - xmin,
                frame["ssdna_y"] - ymin,
                s=5.2,
                c=frame["b_mixed_usage"],
                cmap=color_map,
                norm=normalization,
                alpha=0.90,
                linewidths=0,
                rasterized=True,
            )
            configure_axis(axis, entry)
        axes[row_index, 0].text(
            -0.18,
            0.5,
            disease,
            transform=axes[row_index, 0].transAxes,
            rotation=90,
            ha="center",
            va="center",
            fontsize=13,
            fontweight="bold",
        )

    color_axis = figure.add_axes((0.36, 0.045, 0.28, 0.018))
    color_bar = figure.colorbar(
        ScalarMappable(norm=normalization, cmap=color_map),
        cax=color_axis,
        orientation="horizontal",
    )
    color_bar.set_label("B-enriched / mixed-immune continuous usage", fontsize=9)
    figure.suptitle(
        "B-enriched / mixed-immune K=7 neighborhood usage on ssDNA",
        fontsize=17,
        fontweight="bold",
    )
    figure.text(
        0.5,
        0.935,
        "Fixed 0–1 scale across all patients and chips; this is neighborhood exposure, not absolute B-cell abundance.",
        ha="center",
        fontsize=9.5,
        color="#4D5963",
    )
    figure.subplots_adjust(
        left=0.043, right=0.995, top=0.90, bottom=0.115, wspace=0.035, hspace=0.24
    )
    for extension in ("png", "pdf"):
        figure.savefig(
            output / f"FIG02_ALL21_B_MIXED_USAGE_ON_SSDNA.{extension}",
            dpi=260,
            bbox_inches="tight",
        )
    plt.close(figure)


def main() -> None:
    args = parse_args()
    output = output_path(args.output_dir)
    output.mkdir(parents=True, exist_ok=False)
    Image.MAX_IMAGE_PIXELS = None

    data, inputs = load_spatial_data()
    entries = panel_entries(data)
    expected_counts = {"HC": 6, "IPF": 6, "SSc-ILD": 9}
    observed_counts = {disease: len(entries[disease]) for disease in DISEASE_ORDER}
    if observed_counts != expected_counts:
        raise ValueError(f"Unexpected disease patient counts: {observed_counts}")

    inside_rows = []
    with ExitStack() as stack:
        images = {
            chip: stack.enter_context(Image.open(CHIPS[chip]["tif"])) for chip in CHIP_ORDER
        }
        for chip in CHIP_ORDER:
            frame = data.loc[data["chip"].eq(chip)]
            image = images[chip]
            inside = (
                frame["ssdna_x"].between(0, image.width - 1)
                & frame["ssdna_y"].between(0, image.height - 1)
            )
            inside_rows.append(
                {
                    "chip": chip,
                    "n_bins": len(frame),
                    "inside_ssdna_fraction": float(inside.mean()),
                    "ssdna_width": image.width,
                    "ssdna_height": image.height,
                }
            )
        if min(row["inside_ssdna_fraction"] for row in inside_rows) < 0.999:
            raise ValueError("One or more transformed Bin50 centers fall outside ssDNA")
        render_dominant_map(entries, images, output)
        render_b_mixed_map(entries, images, output)

    manifest = (
        data.groupby(["chip", "disease", "patient_id"], observed=True)
        .size()
        .rename("mapped_bin50_n")
        .reset_index()
        .sort_values(["disease", "chip", "patient_id"])
    )
    manifest.to_csv(output / "01_SAMPLE_MAP_MANIFEST.tsv", sep="\t", index=False)
    pd.DataFrame(inside_rows).to_csv(
        output / "02_COORDINATE_VALIDATION.tsv", sep="\t", index=False
    )
    input_manifest = pd.DataFrame(
        [
            {
                "input": name,
                "path": str(path),
                "exists": path.exists(),
                "sha256": sha256(path) if path.exists() and path.suffix.lower() != ".tif" else "not_computed_large_tiff",
            }
            for name, path in inputs.items()
        ]
    )
    input_manifest.to_csv(output / "00_INPUT_MANIFEST.tsv", sep="\t", index=False)

    figure_paths = [
        output / "FIG01_ALL21_K7_DOMINANT_FACTORS_ON_SSDNA.png",
        output / "FIG01_ALL21_K7_DOMINANT_FACTORS_ON_SSDNA.pdf",
        output / "FIG02_ALL21_B_MIXED_USAGE_ON_SSDNA.png",
        output / "FIG02_ALL21_B_MIXED_USAGE_ON_SSDNA.pdf",
    ]
    validation = {
        "patient_n": int(manifest["patient_id"].nunique()),
        "disease_patient_counts": observed_counts,
        "mapped_bin50_n": int(manifest["mapped_bin50_n"].sum()),
        "family_n": int(data["family"].nunique()),
        "families": list(FAMILY_ORDER),
        "minimum_inside_ssdna_fraction": min(
            row["inside_ssdna_fraction"] for row in inside_rows
        ),
        "four_figures_exist_and_nonempty": all(
            path.exists() and path.stat().st_size > 0 for path in figure_paths
        ),
        "nmf_or_rctd_refit": False,
        "disease_test": False,
    }
    (output / "03_VALIDATION.json").write_text(
        json.dumps(validation, indent=2), encoding="utf-8"
    )
    (output / "RUN_LOG.txt").write_text(
        "\n".join(
            [
                "Workflow: frozen harmonized K=7 maps on ssDNA",
                f"Date: {dt.datetime.now().isoformat(timespec='seconds')}",
                "Command: envs/stereo_seq_py/Scripts/python.exe scripts/1157_render_reference18_k7_ssdna_maps.py",
                f"Output: {output}",
                "NMF/RCTD refit: no",
                "Disease test: no",
                "Status: PASS",
            ]
        )
        + "\n",
        encoding="utf-8",
    )
    (output / "SESSION_INFO.txt").write_text(
        "\n".join(
            [
                f"Python: {sys.version}",
                f"Platform: {platform.platform()}",
                f"numpy: {np.__version__}",
                f"pandas: {pd.__version__}",
                f"matplotlib: {matplotlib.__version__}",
            ]
        )
        + "\n",
        encoding="utf-8",
    )
    (output / "COMPLETE.ok").write_text("COMPLETE\n", encoding="utf-8")
    print(output)


if __name__ == "__main__":
    main()
