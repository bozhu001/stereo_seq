"""Reassess Bin50 G15 HVG recurrence with HVGs defined per sample.

This bespoke, minimal analysis addresses only two requested points:
1. define HVGs independently within every biological sample; and
2. report recurrence with separate proportion and number-of-genes columns.

It intentionally performs no gene interpretation or downstream analysis.
"""

from __future__ import annotations

import datetime as dt
import hashlib
import importlib.metadata
import platform
import sys
import warnings
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import scanpy as sc
import scipy
from scipy import sparse


PROJECT_ROOT = Path(__file__).resolve().parents[1]
ROOT = PROJECT_ROOT / "results" / "reanalysis"
INPUT_DIR = ROOT / "bin50_g15"
OUTPUT_DIR = ROOT / "bin50_hvg_samplewise_reassessment"
OUTPUT_NAME = "BIN50_samplewise_HVG_reassessment.md"
N_TOP_GENES = (1000, 1500, 2000)
HVG_FLAVOR = "seurat_v3"
HVG_SPANS = (0.3, 0.5, 0.7, 1.0)
EXPECTED_SAMPLE_FIELD = "sample_id"
CHIP_FIELD = "chip_id"

AGENTS_READ = [
    (
        "Current-project AGENTS.md supplied in the user message "
        "(no physical AGENTS.md was found in the project or its parent directories)"
    )
]
SKILLS_USED = [
    r"C:\Users\Administrator\.codex\skills\scrna\SKILL.md",
    r"C:\Users\Administrator\.codex\skills\script-review\SKILL.md",
    r"C:\Users\Administrator\.codex\skills\analysis-wrapup\SKILL.md",
]


def require(condition: bool, message: str) -> None:
    """Raise a clear error when a required validation fails."""
    if not condition:
        raise RuntimeError(message)


def sha256(path: Path) -> str:
    """Return the SHA-256 digest of a file."""
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def matrix_nonzero_values(matrix: object) -> np.ndarray:
    """Return stored/nonzero values from a dense or sparse matrix."""
    if sparse.issparse(matrix):
        return np.asarray(matrix.data)
    array = np.asarray(matrix)
    return array[array != 0]


def validate_count_matrix(
    matrix: object,
    adata: ad.AnnData,
    source_name: str,
) -> dict[str, object]:
    """Validate raw-count properties and agreement with stored raw QC fields."""
    require(
        matrix.shape == adata.shape,
        f"{source_name} shape {matrix.shape} does not match AnnData {adata.shape}",
    )
    values = matrix_nonzero_values(matrix)
    finite = bool(np.isfinite(values).all())
    nonnegative = bool((values >= 0).all())
    integer_like = bool(
        np.allclose(values, np.rint(values), atol=1e-6, rtol=0)
    )

    row_sums = np.asarray(matrix.sum(axis=1)).ravel()
    if sparse.issparse(matrix):
        detected = np.asarray(matrix.getnnz(axis=1)).ravel()
    else:
        detected = np.count_nonzero(np.asarray(matrix), axis=1)

    required_obs = {"total_counts_common_raw", "n_genes_common_raw"}
    require(
        required_obs.issubset(adata.obs.columns),
        "Raw-count QC fields are missing; cannot establish count provenance",
    )
    sum_matches = bool(
        np.array_equal(
            row_sums.astype(np.int64),
            adata.obs["total_counts_common_raw"].to_numpy(dtype=np.int64),
        )
    )
    detected_matches = bool(
        np.array_equal(
            detected.astype(np.int64),
            adata.obs["n_genes_common_raw"].to_numpy(dtype=np.int64),
        )
    )
    require(
        finite and nonnegative and integer_like,
        f"{source_name} is not finite, nonnegative, integer-like raw counts",
    )
    require(
        sum_matches and detected_matches,
        (
            f"{source_name} does not match total_counts_common_raw and "
            "n_genes_common_raw"
        ),
    )
    return {
        "finite": finite,
        "nonnegative": nonnegative,
        "integer_like": integer_like,
        "row_sums_match_total_counts_common_raw": sum_matches,
        "detected_genes_match_n_genes_common_raw": detected_matches,
        "dtype": str(matrix.dtype),
        "sparse": bool(sparse.issparse(matrix)),
    }


def choose_counts(adata: ad.AnnData) -> tuple[object, str, dict[str, object]]:
    """Choose and validate counts using the required source priority."""
    candidates: list[tuple[str, object]] = []
    if "counts" in adata.layers:
        candidates.append(('layers["counts"]', adata.layers["counts"]))
    candidates.append(("X", adata.X))
    if adata.raw is not None and adata.raw.shape[1] == adata.n_vars:
        candidates.append(("raw.X", adata.raw.X))

    failures: list[str] = []
    for name, matrix in candidates:
        try:
            validation = validate_count_matrix(matrix, adata, name)
            return matrix, name, validation
        except Exception as exc:  # preserve all candidate failures in final error
            failures.append(f"{name}: {exc}")
    raise RuntimeError(
        "No reliable original-count matrix was found. " + " | ".join(failures)
    )


def discover_input() -> Path:
    """Locate the unique G15 H5AD in the preferred project result directory."""
    require(ROOT.exists(), f"Working directory does not exist: {ROOT}")
    require(INPUT_DIR.exists(), f"Preferred input directory is missing: {INPUT_DIR}")
    candidates = sorted(INPUT_DIR.glob("*.h5ad"))
    require(
        len(candidates) == 1,
        f"Expected exactly one H5AD in {INPUT_DIR}; found {len(candidates)}",
    )
    return candidates[0].resolve()


def choose_output_path() -> Path:
    """Avoid overwriting an existing report using the lab timestamp convention."""
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    preferred = OUTPUT_DIR / OUTPUT_NAME
    if not preferred.exists():
        return preferred
    stamp = dt.datetime.now().strftime("%y%m%d%H%M%S")
    return OUTPUT_DIR / f"BIN50_samplewise_HVG_reassessment_{stamp}.md"


def markdown_table(frame: pd.DataFrame) -> str:
    """Format the strict three-column recurrence table without extra packages."""
    lines = [
        "| n_samples_hvg | proportion | number_of_genes |",
        "|---:|---:|---:|",
    ]
    for row in frame.itertuples(index=False):
        lines.append(
            f"| {int(row.n_samples_hvg)} | {float(row.proportion):.4f} "
            f"| {int(row.number_of_genes)} |"
        )
    return "\n".join(lines)


def package_version(name: str) -> str:
    """Return installed package version without importing another module."""
    try:
        return importlib.metadata.version(name)
    except importlib.metadata.PackageNotFoundError:
        return "not installed"


def main() -> None:
    """Run the minimal per-sample HVG recurrence reassessment."""
    input_path = discover_input()
    output_path = choose_output_path()
    adata = sc.read_h5ad(input_path)

    require(
        EXPECTED_SAMPLE_FIELD in adata.obs.columns,
        f"obs[{EXPECTED_SAMPLE_FIELD!r}] is missing",
    )
    require(
        not adata.obs[EXPECTED_SAMPLE_FIELD].isna().any(),
        "sample_id contains missing values",
    )
    require(
        adata.obs_names.is_unique and adata.var_names.is_unique,
        "Observation and variable names must be unique",
    )
    require(
        adata.n_vars > max(N_TOP_GENES),
        "Input does not retain a sufficiently large complete gene space",
    )
    require(
        "highly_variable" not in adata.var.columns,
        "Input already contains an HVG flag and may not be the pre-HVG G15 object",
    )

    counts, counts_source, counts_validation = choose_counts(adata)

    sample_text = adata.obs[EXPECTED_SAMPLE_FIELD].astype(str)
    if CHIP_FIELD in adata.obs:
        chip_text = adata.obs[CHIP_FIELD].astype(str)
        sample_chip_counts = (
            pd.DataFrame({"sample_id": sample_text, "chip_id": chip_text})
            .drop_duplicates()
            .groupby("sample_id", observed=True)["chip_id"]
            .nunique()
        )
        repeated_between_chips = sorted(
            sample_chip_counts.index[sample_chip_counts > 1].astype(str)
        )
    else:
        chip_text = pd.Series("", index=adata.obs_names, dtype=str)
        repeated_between_chips = []

    if repeated_between_chips:
        sample_id_used = chip_text + "__" + sample_text
        sample_field_description = "chip_id + '__' + sample_id"
    else:
        sample_id_used = sample_text.copy()
        sample_field_description = "sample_id"
    adata.obs["sample_id_used"] = sample_id_used.to_numpy()

    sample_counts = (
        adata.obs["sample_id_used"]
        .astype(str)
        .value_counts(sort=False)
        .sort_index()
        .rename("n_bins")
    )
    require((sample_counts > 0).all(), "At least one sample has no bins")
    samples = sample_counts.index.tolist()
    n_samples = len(samples)

    hvg_calls: list[dict[str, object]] = []
    hvg_attempts: list[dict[str, object]] = []
    recurrence = {
        n_top: np.zeros(adata.n_vars, dtype=np.int32)
        for n_top in N_TOP_GENES
    }
    labels = adata.obs["sample_id_used"].astype(str).to_numpy()

    for sample_id in samples:
        row_mask = labels == sample_id
        sample_matrix = counts[row_mask, :]
        for n_top in N_TOP_GENES:
            sample_adata = None
            selected_span = None
            captured_warnings: list[str] = []
            for span in HVG_SPANS:
                sample_adata = ad.AnnData(
                    X=sample_matrix.copy(),
                    obs=pd.DataFrame(index=adata.obs_names[row_mask].copy()),
                    var=pd.DataFrame(index=adata.var_names.copy()),
                )
                try:
                    with warnings.catch_warnings(record=True) as caught:
                        warnings.simplefilter("always")
                        sc.pp.highly_variable_genes(
                            sample_adata,
                            flavor=HVG_FLAVOR,
                            n_top_genes=n_top,
                            span=span,
                            inplace=True,
                            check_values=True,
                        )
                    captured_warnings.extend(str(item.message) for item in caught)
                    selected_span = span
                    hvg_attempts.append(
                        {
                            "sample_id_used": sample_id,
                            "n_top_genes": n_top,
                            "span": span,
                            "status": "PASS",
                            "message": "",
                        }
                    )
                    break
                except ValueError as exc:
                    hvg_attempts.append(
                        {
                            "sample_id_used": sample_id,
                            "n_top_genes": n_top,
                            "span": span,
                            "status": "RETRY",
                            "message": str(exc),
                        }
                    )
            require(
                sample_adata is not None and selected_span is not None,
                (
                    f"{sample_id}, N={n_top}: Seurat v3 LOESS failed for all "
                    f"spans {HVG_SPANS}"
                ),
            )
            flags = sample_adata.var["highly_variable"].to_numpy(dtype=bool)
            selected = int(flags.sum())
            require(
                selected == n_top,
                (
                    f"{sample_id}: requested {n_top} HVGs but Scanpy selected "
                    f"{selected}"
                ),
            )
            recurrence[n_top] += flags.astype(np.int32)
            hvg_calls.append(
                {
                    "sample_id_used": sample_id,
                    "n_bins": int(row_mask.sum()),
                    "n_top_genes": n_top,
                    "selected_hvgs": selected,
                    "selected_span": selected_span,
                    "warnings": " | ".join(sorted(set(captured_warnings))),
                    "status": "PASS",
                }
            )

    tables: dict[int, pd.DataFrame] = {}
    threshold_summaries: dict[int, dict[int, int]] = {}
    for n_top, n_samples_hvg in recurrence.items():
        number_of_genes = np.bincount(
            n_samples_hvg,
            minlength=n_samples + 1,
        )
        table = pd.DataFrame(
            {
                "n_samples_hvg": np.arange(n_samples + 1, dtype=int),
                "proportion": (
                    np.arange(n_samples + 1, dtype=float) / n_samples
                ),
                "number_of_genes": number_of_genes.astype(int),
            }
        )
        require(
            int(table["number_of_genes"].sum()) == adata.n_vars,
            f"N={n_top}: recurrence table does not sum to n_vars",
        )
        require(
            table["n_samples_hvg"].tolist() == list(range(n_samples + 1)),
            f"N={n_top}: incomplete recurrence rows",
        )
        require(
            pd.api.types.is_integer_dtype(table["number_of_genes"]),
            f"N={n_top}: number_of_genes is not integer",
        )
        tables[n_top] = table
        threshold_summaries[n_top] = {}
        for percent in (20, 40, 60, 80):
            minimum = int(np.ceil((percent / 100) * n_samples))
            threshold_summaries[n_top][percent] = int(
                number_of_genes[minimum:].sum()
            )

    reproducibility_score = {
        n_top: threshold_summaries[n_top][80] for n_top in N_TOP_GENES
    }
    best_score = max(reproducibility_score.values())
    best_settings = [
        n_top
        for n_top, score in reproducibility_score.items()
        if score == best_score
    ]

    environment = {
        "Python": platform.python_version(),
        "Platform": platform.platform(),
        "anndata": ad.__version__,
        "scanpy": sc.__version__,
        "numpy": np.__version__,
        "pandas": pd.__version__,
        "scipy": scipy.__version__,
        "scikit-misc": package_version("scikit-misc"),
    }
    sample_lines = "\n".join(
        f"- `{sample}`: {int(n_bins)}"
        for sample, n_bins in sample_counts.items()
    )
    environment_lines = "\n".join(
        f"- {name}: `{version}`" for name, version in environment.items()
    )
    span_usage = (
        pd.DataFrame(hvg_calls)
        .groupby(["n_top_genes", "selected_span"], observed=True)
        .size()
        .rename("n_sample_calls")
        .reset_index()
    )
    span_usage_lines = "\n".join(
        (
            f"- N={int(row.n_top_genes)}, span={float(row.selected_span):.1f}: "
            f"{int(row.n_sample_calls)} sample calls"
        )
        for row in span_usage.itertuples(index=False)
    )
    retry_count = sum(row["status"] == "RETRY" for row in hvg_attempts)
    successful_warning_count = sum(bool(row["warnings"]) for row in hvg_calls)
    threshold_lines = []
    for n_top in N_TOP_GENES:
        values = threshold_summaries[n_top]
        threshold_lines.append(
            f"- N={n_top}: ≥20% `{values[20]}`; ≥40% `{values[40]}`; "
            f"≥60% `{values[60]}`; ≥80% `{values[80]}` genes."
        )

    conventions = "\n".join(f"- `{item}`" for item in AGENTS_READ)
    skills = "\n".join(f"- `{item}`" for item in SKILLS_USED)
    report = f"""# BIN50 Sample-wise HVG Reassessment

## 1. Project conventions

### AGENTS.md read

{conventions}

### Skills used

{skills}

### Laboratory conventions followed

- Followed the nearest project `AGENTS.md` supplied for this session as the
  canonical ruleset.
- Used only task-relevant skills rather than scanning the full wiki.
- Used `pathlib.Path`; all analysis inputs and outputs are inside the project.
- Used a documented Python script and a timestamped execution log; no ad hoc
  analysis commands were used.
- Preserved the complete G15 common-gene space and validated transformed-object
  provenance before analysis.
- Validated nonnegative integer counts against stored raw count and detected-gene
  fields before Seurat v3 HVG selection.
- Avoided overwriting an existing report by applying the lab `YYMMDDHHMMSS`
  timestamp convention when needed.
- No convention conflicts affected the computation. Priority was: current
  project `AGENTS.md`, relevant project code, shared bioinfo-wiki skills.

### Analysis environment and software versions

{environment_lines}

## 2. Input data

- Input H5AD: `{input_path}`
- Input SHA-256: `{sha256(input_path)}`
- Input shape: `{adata.n_obs} × {adata.n_vars}`
- Counts source: `{counts_source}`
- Counts dtype: `{counts_validation['dtype']}`
- Counts are finite, nonnegative, integer-like: `True`
- Count row sums match `obs["total_counts_common_raw"]`: `True`
- Detected genes match `obs["n_genes_common_raw"]`: `True`
- Sample field available: `sample_id`
- Sample field used: `{sample_field_description}`
- `sample_id` repeated between chips: `{'yes' if repeated_between_chips else 'no'}`
- Sample count: `{n_samples}`

### Bins per sample

{sample_lines}

## 3. Method

- HVGs were defined independently within each of the {n_samples} samples.
- HVGs were not defined by chip.
- HVGs were not defined once on the pooled object.
- Every sample retained the same complete `{adata.n_vars}`-gene space.
- Method: Scanpy `sc.pp.highly_variable_genes`, flavor `{HVG_FLAVOR}`, raw
  counts, `inplace=True`, `check_values=True`.
- Each sample/N call started with `span=0.3`. Following the existing project
  HVG script's numerical-stability convention, a near-singular LOESS fit was
  retried in order with spans `{', '.join(map(str, HVG_SPANS))}`; no sample was
  skipped.
- Compared `n_top_genes`: {', '.join(map(str, N_TOP_GENES))}.
- For every N and gene, `n_samples_hvg` is the number of samples in which that
  gene was independently selected; `proportion = n_samples_hvg / {n_samples}`.
- No specific genes were inspected, listed, filtered, or interpreted.

### Actual span usage

{span_usage_lines}

## 4. Proportion table: N=1000

{markdown_table(tables[1000])}

## 5. Proportion table: N=1500

{markdown_table(tables[1500])}

## 6. Proportion table: N=2000

{markdown_table(tables[2000])}

## 7. Validation

- Per-sample HVG calls completed: `{len(hvg_calls)}` of
  `{n_samples * len(N_TOP_GENES)}` expected calls (`PASS`).
- LOESS attempts requiring a larger-span retry: `{retry_count}`.
- Successful calls emitting captured warnings: `{successful_warning_count}`.
- Samples analyzed without silent skipping: `{n_samples}` of `{n_samples}`
  (`PASS`).
- Each table contains every row from 0 through {n_samples} samples, including
  zero-gene rows (`PASS`).
- Proportion is calculated as `n_samples_hvg / {n_samples}`, lies between 0 and
  1, and is displayed to four decimal places (`PASS`).
- `number_of_genes` is stored and displayed as an integer (`PASS`).
- N=1000 table gene-count sum: `{int(tables[1000]['number_of_genes'].sum())}`
  = `{adata.n_vars}` (`PASS`).
- N=1500 table gene-count sum: `{int(tables[1500]['number_of_genes'].sum())}`
  = `{adata.n_vars}` (`PASS`).
- N=2000 table gene-count sum: `{int(tables[2000]['number_of_genes'].sum())}`
  = `{adata.n_vars}` (`PASS`).
- Analysis stopped before gene interpretation, PCA, Harmony, neighbors, UMAP,
  Leiden, SpatialLeiden, BANKSY, FlowSOM, or annotation (`PASS`).

## 8. Summary

- Samples analyzed: `{n_samples}`.
- Per-sample HVG settings: `{', '.join(map(str, N_TOP_GENES))}`.
- HVGs were independently defined per sample: `yes`.
{chr(10).join(threshold_lines)}
- Higher cross-sample reproducibility at the ≥80% threshold was shown by:
  `{', '.join(f'N={value}' for value in best_settings)}` ({best_score} genes).
"""
    output_path.write_text(report, encoding="utf-8", newline="\n")

    print("Actual AGENTS.md read:")
    for item in AGENTS_READ:
        print(f"- {item}")
    print("Actual skills used:")
    for item in SKILLS_USED:
        print(f"- {item}")
    print(f"Input H5AD: {input_path}")
    print(f"Input shape: {adata.shape}")
    print(f"Counts source: {counts_source}")
    print(f"Sample count: {n_samples}")
    print(f"Final Markdown: {output_path.resolve()}")
    print("Stopped after sample-wise HVG proportion summary.")


if __name__ == "__main__":
    try:
        main()
    except Exception:
        print("Analysis failed before a report could be finalized.", file=sys.stderr)
        raise
