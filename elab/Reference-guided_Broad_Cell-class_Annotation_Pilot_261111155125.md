# Reference-guided Broad Cell-class Annotation Pilot — 261111155125

Status: approved by the user on 2026-11-11. Local final copy. A remote
eLabFTW entry could not be created because `ELABFTW_API_KEY` was unavailable
in the process, Windows user, and Windows machine environments. No entry ID
or URL was invented.

Category: Bioinformatic

Tags: `Stereo-seq`, `CellBin`, `Reference18`, `SingleR`,
`broad-class annotation`, `pilot`, `abstention`, `FFPE`

## Goal

Assess whether label transfer from the frozen Reference18 with SingleR was
suitable for exploratory broad cell-class annotation of the G10/U20 CellBins
from all 21 Stereo-seq FFPE samples. The prespecified design required a useful
Unknown/abstention output and allowed expansion to all 359,329 assigned
G10/U20 CellBins only after a small pilot passed input, reference-validation,
and uncertainty gates.

## Input

Each input was checked against local eLabFTW copies. No prior logged entry was
found that exactly represented this combination of frozen Reference18 inputs,
Joint CellBin V1, and CellBin QC metrics; therefore no provenance link was
invented.

1. Frozen Reference18 source counts:
   `results/reanalysis/RCTD_SIX_DATASET_LEVEL1_ALL21_20260916_223802/prepared_inputs/six_dataset_level1_rctd_inputs.h5`
   - Size: 427,977,861 bytes
   - SHA-256: `0ca7d516402b8c1e41eaabeaa3b02811bfffce13451ea5cab7c17638a47d84f9`
2. Selected reference-cell metadata:
   `results/reanalysis/RCTD_SIX_DATASET_LEVEL1_ALL21_20260916_223802/prepared_inputs/selected_reference_cells.tsv.gz`
   - Size: 517,681 bytes
   - SHA-256: `6a81bc6563eb67367e5380091f0e29bde783b5566cf929a62e0c3cf33c781707`
3. Frozen Reference18 cell-to-state mapping:
   `results/reanalysis/RCTD_LEVEL2_REFERENCE18_ALL21_20260919_150957/prepared_inputs/REFERENCE_18_CELL_MAPPING.tsv`
   - Size: 4,244,517 bytes
   - SHA-256: `c34939305996aa66018e699b2798c2730b74636e1a47b344a5ed705db868588b`
4. Frozen Reference18 18-state definitions:
   `results/reanalysis/RCTD_LEVEL2_REFERENCE18_ALL21_20260919_150957/02_REFERENCE18_STATE_DEFINITION.tsv`
   - Size: 1,435 bytes
   - SHA-256: `2bb380b15ae179293318d0422c191b2e0ce98f50dadf7e321681af0cf1173c1e`
5. Joint CellBin V1:
   `results/cell bin/CellBin_joint_raw_21samples/Joint_CellBin_raw_with_sampleID.h5ad`
   - Size: 183,428,028 bytes
   - SHA-256: `c9fde3462856f5df70a36435554fe9dea7bf3dbdc82313e870f98cb0237f593b`
6. CellBin QC metrics:
   `results/cell bin/CellBin_QC_diagnostics_21samples/CellBin_QC_metrics.parquet`
   - Size: 28,785,865 bytes
   - SHA-256: `2a28d56b123c3d78bb1ef00850000af8eea34b1df098437ff957b9052d6a5108`

The complete machine-readable provenance table is
`results/cell bin/CellBin_BroadClass_pilot_21samples/Input_provenance.tsv`.

## Script

Pipeline commands, run in sequence from the repository root:

```text
python scripts/1248_cellbin_broadclass_reference_preflight.py
python scripts/1249_prepare_cellbin_broadclass_singler_pilot.py
Rscript scripts/1250_run_cellbin_broadclass_singler_pilot.R
python scripts/1251_evaluate_cellbin_broadclass_singler_pilot.py
powershell -ExecutionPolicy Bypass -File scripts/render_cellbin_broadclass_pilot_report.ps1
```

Key parameters and safeguards:

- The original 18 Reference18 states were unchanged and were mapped
  traceably to four broad classes for this pilot only.
- Reference cells were aggregated into 889 donor × Reference18-state profiles
  from 124 donors.
- The query pilot contained 500 assigned G10/U20 CellBins from each of all 21
  samples (10,500 total).
- SingleR 2.12.0 default pruning was retained as Unknown; no spatial or
  nearest-neighbor filling was allowed.
- Reference self-validation used donor-grouped five-fold cross-validation.
- Canonical markers were orthogonal validation evidence, not the classifier
  or the sole gate.

Repository: `https://github.com/bozhu001/stereo_seq`

Analysis/report commit:
`b36f44a98dd756cf3ccca27c409f8f97861e02b9`

## Output summary

Output directory:
`results/cell bin/CellBin_BroadClass_pilot_21samples/`

Rendered standalone report:
`results/cellbin_broadclass_pilot_21samples_261111155125.html`

- Size: 2,034,390 bytes
- SHA-256: `B242ED2A2B8748272205508EFFAE975DFFF85C87124790FAD87D6AFDE3731744`
- Render status: PASS

The report is under 100 MB and should be attached if this approved local entry
is later recreated in the eLabFTW web interface.

Reference18 donor-grouped broad-class self-validation accuracy was 99.8875%,
but default SingleR pruning produced no Unknown query CellBins. The query
median maximum score was 0.0591 versus 0.3627 for the fifth percentile among
correctly classified held-out reference profiles, giving a query/reference
score-domain ratio of 0.163. Only 32.96% of forced predictions had at least
two independent canonical markers consistent with the predicted broad class.
Maximum score had Spearman correlations of 0.564 and 0.560 with nGenes and
nUMI, respectively.

The pilot therefore stopped at the prespecified uncertainty gate. It was not
expanded to all 359,329 CellBins, and its predictions are not formal or
provisional cohort cell-type labels. The workflow is retained as an
experimental, project-reusable workflow and is not frozen as the formal
annotation method.

Depth-matched Reference18 validation and Route A unsupervised clustering were
not started in this analysis. Route A is explicitly reserved for a separate
future task.
