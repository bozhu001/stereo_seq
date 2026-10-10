# B-lineage coverage mapping to immune subniches — 261010225750

- Category: `Bioinformatic`
- Tags: `Stereo-seq`, `Bin50`, `B-lineage`, `immune-subniche`, `K-resolution`, `spatial-QC`
- Approval: user approved the ELN draft on 2026-10-10
- Remote eLabFTW ID/URL: unavailable; `ELABFTW_API_KEY` was not provisioned in the agent environment
- Report attachment for manual upload: `results/B_LINEAGE_SUBNICHE_MAPPING_REPORT_261010225750.html`

## Goal

Map the frozen 2,669-Bin50 B-lineage RNA coverage mask onto the frozen
chip-specific K=3/4/5 immune-subniche assignments. Determine whether
B-lineage support follows reproducible cross-resolution branches while
preserving patient-level replication and patient-imbalance warnings.

## Input

1. Final calibrated B-lineage coverage mask:
   `results/reanalysis/B_LINEAGE_CALIBRATED_COVERAGE_MASK_GE5_ONLY_20261010_220617/01_FINAL_CALIBRATED_B_COVERAGE_MASK.tsv.gz`.
   Prior local ELN:
   `elab/Bin50_calibrated_B-lineage_RNA_coverage_mask_261010215502.md`.
   No remote eLabFTW ID exists, so an API provenance link could not be
   created.
2. Frozen K=3/4/5 immune-subniche assignments:
   `results/reanalysis/EXTENDED_IMMUNE_SIGNATURE_VALIDATION_20261010_192044/work/01_FROZEN_K345_ASSIGNMENTS.tsv.gz`.
   No standalone remote eLabFTW entry was found for this input.

## Script

Analysis command:

```powershell
powershell.exe -ExecutionPolicy Bypass `
  -File scripts/1055_run_b_lineage_subniche_mapping.ps1 `
  -OutputDir results/reanalysis/B_LINEAGE_SUBNICHE_MAPPING_20261010_225151
```

Report rendering:

```bash
tmux new-session -d -s b_subniche_report \
  'bash scripts/1056_render_b_lineage_subniche_mapping_report.sh'
```

Analysis/report commit:
`6e94986`
(`https://github.com/bozhu001/stereo_seq/commit/6e94986`).

No RCTD, neighborhood construction, clustering, K selection, marker
calibration, or disease comparison was rerun. Patient, rather than Bin50, was
treated as the biological replication unit.

## Output summary

- Reproduced 29,645 immune-candidate Bin50 and the final 2,669-bin B-lineage
  mask.
- L3 B-enriched lineage: K3 I1 → K4 I4 → K5 I3; 6/7 patients showed
  enrichment. K3 captured 74.5% of the L3 B mask; median patient enrichment
  ratio 3.08.
- J2 B-enriched lineage: K3 I3 → K4 I1 → K5 I3; 6/7 patients showed
  enrichment. K3 captured 45.7% of the J2 B mask; median patient enrichment
  ratio 1.91.
- Y40102K8 B-enriched lineage: K3 I1 → K4 I4 → K5 I3; 5/7 patients showed
  enrichment. K3 captured 99.5% of the chip B mask, but the result was
  strongly patient-imbalanced.
- All recurrent mappings retained a patient-imbalance warning.
- Mapping labels are not final immune-subniche names. Final naming awaits
  lineage-specific T, NK, Plasma, Myeloid, and Mast/Basophil coverage audits.
- All automated final checks passed.

Final output directory:
`results/reanalysis/B_LINEAGE_SUBNICHE_MAPPING_20261010_225151/`.

Rendered HTML:
`results/B_LINEAGE_SUBNICHE_MAPPING_REPORT_261010225750.html`;
980,312 bytes; SHA256
`F94F68D53E8D85095306B32083C694F96B7C291EDE5C24300BC106232D3558D1`.
