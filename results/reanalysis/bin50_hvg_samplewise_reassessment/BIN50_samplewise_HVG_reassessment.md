# BIN50 Sample-wise HVG Reassessment

## 1. Project conventions

### AGENTS.md read

- `Current-project AGENTS.md supplied in the user message (no physical AGENTS.md was found in the project or its parent directories)`

### Skills used

- `C:\Users\Administrator\.codex\skills\scrna\SKILL.md`
- `C:\Users\Administrator\.codex\skills\script-review\SKILL.md`
- `C:\Users\Administrator\.codex\skills\analysis-wrapup\SKILL.md`

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

- Python: `3.12.13`
- Platform: `Windows-11-10.0.26200-SP0`
- anndata: `0.13.2`
- scanpy: `1.12.3`
- numpy: `2.4.6`
- pandas: `3.0.5`
- scipy: `1.18.0`
- scikit-misc: `0.5.2`

## 2. Input data

- Input H5AD: `D:\bb\博士课题组资料\新建文件夹\自己课题相关\Lung ST\codex agent\stereo_seq\stereo_seq\results\reanalysis\bin50_g15\BIN50_joint_reanalysis_G15_raw_counts.h5ad`
- Input SHA-256: `99696a573577f5ef452ec073346cd6e611badc71d6fc96aa89f8565a42ad4fc7`
- Input shape: `152640 × 27618`
- Counts source: `X`
- Counts dtype: `uint32`
- Counts are finite, nonnegative, integer-like: `True`
- Count row sums match `obs["total_counts_common_raw"]`: `True`
- Detected genes match `obs["n_genes_common_raw"]`: `True`
- Sample field available: `sample_id`
- Sample field used: `sample_id`
- `sample_id` repeated between chips: `no`
- Sample count: `21`

### Bins per sample

- `HC/NL-50`: 5708
- `HC/NL-53`: 3249
- `HC/NL-55`: 10636
- `HC/NL-63`: 5774
- `HC/NL-66`: 7135
- `HC/NL-72`: 9919
- `IPF/FO22-1-09404`: 6439
- `IPF/FO23-1-03100`: 9275
- `IPF/FO23-1-06168`: 5358
- `IPF/FO23-1-06170`: 7156
- `IPF/FO23-1-06474`: 1095
- `IPF/FO23-1-09473`: 2831
- `SSC/05957/17B`: 7325
- `SSC/07998/15A`: 6477
- `SSC/15275/13`: 1287
- `SSC/15491/14`: 10213
- `SSC/23-105334B1`: 26898
- `SSC/24-1-18170A2`: 12837
- `SSC/3342/13`: 701
- `SSC/FO20-I-08293B1`: 4450
- `SSC/FO21-1-01756B1`: 7877

## 3. Method

- HVGs were defined independently within each of the 21 samples.
- HVGs were not defined by chip.
- HVGs were not defined once on the pooled object.
- Every sample retained the same complete `27618`-gene space.
- Method: Scanpy `sc.pp.highly_variable_genes`, flavor `seurat_v3`, raw
  counts, `inplace=True`, `check_values=True`.
- Each sample/N call started with `span=0.3`. Following the existing project
  HVG script's numerical-stability convention, a near-singular LOESS fit was
  retried in order with spans `0.3, 0.5, 0.7, 1.0`; no sample was
  skipped.
- Compared `n_top_genes`: 1000, 1500, 2000.
- For every N and gene, `n_samples_hvg` is the number of samples in which that
  gene was independently selected; `proportion = n_samples_hvg / 21`.
- No specific genes were inspected, listed, filtered, or interpreted.

### Actual span usage

- N=1000, span=0.3: 15 sample calls
- N=1000, span=0.5: 5 sample calls
- N=1000, span=0.7: 1 sample calls
- N=1500, span=0.3: 15 sample calls
- N=1500, span=0.5: 5 sample calls
- N=1500, span=0.7: 1 sample calls
- N=2000, span=0.3: 15 sample calls
- N=2000, span=0.5: 5 sample calls
- N=2000, span=0.7: 1 sample calls

## 4. Proportion table: N=1000

| n_samples_hvg | proportion | number_of_genes |
|---:|---:|---:|
| 0 | 0.0000 | 14262 |
| 1 | 0.0476 | 8170 |
| 2 | 0.0952 | 3440 |
| 3 | 0.1429 | 1233 |
| 4 | 0.1905 | 369 |
| 5 | 0.2381 | 107 |
| 6 | 0.2857 | 28 |
| 7 | 0.3333 | 6 |
| 8 | 0.3810 | 2 |
| 9 | 0.4286 | 0 |
| 10 | 0.4762 | 0 |
| 11 | 0.5238 | 0 |
| 12 | 0.5714 | 0 |
| 13 | 0.6190 | 0 |
| 14 | 0.6667 | 1 |
| 15 | 0.7143 | 0 |
| 16 | 0.7619 | 0 |
| 17 | 0.8095 | 0 |
| 18 | 0.8571 | 0 |
| 19 | 0.9048 | 0 |
| 20 | 0.9524 | 0 |
| 21 | 1.0000 | 0 |

## 5. Proportion table: N=1500

| n_samples_hvg | proportion | number_of_genes |
|---:|---:|---:|
| 0 | 0.0000 | 10515 |
| 1 | 0.0476 | 8589 |
| 2 | 0.0952 | 4809 |
| 3 | 0.1429 | 2264 |
| 4 | 0.1905 | 939 |
| 5 | 0.2381 | 333 |
| 6 | 0.2857 | 121 |
| 7 | 0.3333 | 40 |
| 8 | 0.3810 | 5 |
| 9 | 0.4286 | 1 |
| 10 | 0.4762 | 1 |
| 11 | 0.5238 | 0 |
| 12 | 0.5714 | 0 |
| 13 | 0.6190 | 0 |
| 14 | 0.6667 | 0 |
| 15 | 0.7143 | 1 |
| 16 | 0.7619 | 0 |
| 17 | 0.8095 | 0 |
| 18 | 0.8571 | 0 |
| 19 | 0.9048 | 0 |
| 20 | 0.9524 | 0 |
| 21 | 1.0000 | 0 |

## 6. Proportion table: N=2000

| n_samples_hvg | proportion | number_of_genes |
|---:|---:|---:|
| 0 | 0.0000 | 7932 |
| 1 | 0.0476 | 8006 |
| 2 | 0.0952 | 5682 |
| 3 | 0.1429 | 3176 |
| 4 | 0.1905 | 1681 |
| 5 | 0.2381 | 703 |
| 6 | 0.2857 | 281 |
| 7 | 0.3333 | 103 |
| 8 | 0.3810 | 41 |
| 9 | 0.4286 | 11 |
| 10 | 0.4762 | 1 |
| 11 | 0.5238 | 0 |
| 12 | 0.5714 | 0 |
| 13 | 0.6190 | 0 |
| 14 | 0.6667 | 0 |
| 15 | 0.7143 | 0 |
| 16 | 0.7619 | 0 |
| 17 | 0.8095 | 0 |
| 18 | 0.8571 | 0 |
| 19 | 0.9048 | 1 |
| 20 | 0.9524 | 0 |
| 21 | 1.0000 | 0 |

## 7. Validation

- Per-sample HVG calls completed: `63` of
  `63` expected calls (`PASS`).
- LOESS attempts requiring a larger-span retry: `21`.
- Successful calls emitting captured warnings: `0`.
- Samples analyzed without silent skipping: `21` of `21`
  (`PASS`).
- Each table contains every row from 0 through 21 samples, including
  zero-gene rows (`PASS`).
- Proportion is calculated as `n_samples_hvg / 21`, lies between 0 and
  1, and is displayed to four decimal places (`PASS`).
- `number_of_genes` is stored and displayed as an integer (`PASS`).
- N=1000 table gene-count sum: `27618`
  = `27618` (`PASS`).
- N=1500 table gene-count sum: `27618`
  = `27618` (`PASS`).
- N=2000 table gene-count sum: `27618`
  = `27618` (`PASS`).
- Analysis stopped before gene interpretation, PCA, Harmony, neighbors, UMAP,
  Leiden, SpatialLeiden, BANKSY, FlowSOM, or annotation (`PASS`).

## 8. Summary

- Samples analyzed: `21`.
- Per-sample HVG settings: `1000, 1500, 2000`.
- HVGs were independently defined per sample: `yes`.
- N=1000: ≥20% `144`; ≥40% `1`; ≥60% `1`; ≥80% `0` genes.
- N=1500: ≥20% `502`; ≥40% `3`; ≥60% `1`; ≥80% `0` genes.
- N=2000: ≥20% `1141`; ≥40% `13`; ≥60% `1`; ≥80% `1` genes.
- Higher cross-sample reproducibility at the ≥80% threshold was shown by:
  `N=2000` (1 genes).
