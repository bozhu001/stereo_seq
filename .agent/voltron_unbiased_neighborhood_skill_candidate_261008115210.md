# Skill candidate: reusable per-chip VoltRon unbiased neighborhood analysis

Proposed issue title: `Create reusable per-chip VoltRon unbiased neighborhood niche workflow`

## Why this should be a skill

The same independently fitted workflow has been used for L3 and J2 and is
approved for Y40102K8. It has stable scope, repeated inputs/outputs, and a
substantial audit burden that should not be reimplemented chip by chip.

## Proposed skill name

`voltron-unbiased-neighborhood`

## Purpose

Run an independent per-chip Bin50-level spatial neighborhood discovery from
frozen Reference18 continuous weights, including graph audit, native CLR,
K-means resolution scanning, stability diagnostics, conservative composition
annotation, cross-patient consistency, and UMI sensitivity.

## Required inputs

- Chip ID.
- Frozen Bin50 table with patient, ROI, physical coordinates, UMI, and 18
  continuous Reference18 weights.
- UMI threshold, radius, K range, random seeds, and output root.
- Reference18 state order and patient metadata.

## Core workflow

1. Audit inputs and eligible/evaluable Bin counts.
2. Build or validate a within-patient/ROI radius neighborhood without turning
   empty neighborhoods into zero profiles.
3. Compute neighborhood composition and VoltRon-native CLR features.
4. Fit K-means over a declared K range with repeated fixed seeds.
5. Audit convergence, ARI/NMI, niche size, patient coverage, UMI association,
   spatial distribution, and Reference18 interpretability.
6. Freeze a primary K plus lower/higher comparators without using disease
   separation or aggregate labels.
7. Generate primary-K composition, within-patient enrichment, cross-patient
   direction consistency, UMI sensitivity, conservative annotations, final
   audit, figures, logs, and report-ready tables.

## Generalized parameters

- `chip_id`, `umi_threshold`, `radius_um`, `center_bin_included`.
- `k_min`, `k_max`, `seeds_per_k`, seed schedule, convergence controls.
- Minimum reporting thresholds for tiny/single-patient niches.
- Primary/lower/higher resolution recommendation rules.
- Output naming, color keys, and report paths.

## Guardrails

- Fit each chip independently.
- Do not use disease labels or frozen aggregate labels in neighborhood
  construction, clustering, K selection, or annotation.
- Do not map niche numbers across chips.
- Do not silently replace official VoltRon behavior with a custom transform.
- Keep not-evaluable Bin separate from zero-valued features.

## Source implementations

- L3 resolution workflow around `scripts/997_l3_unbiased_niche_resolution.R`
  and its follow-up scripts.
- J2 implementation: `scripts/1013_j2_unbiased_niche_resolution.R` through
  `scripts/1017_j2_unbiased_niche_resolution_report.qmd`.

## Expected outputs

Input/neighborhood audits, seed-level and K-level diagnostics, frozen
assignments, primary-K Reference18 tables, patient-specific enrichment,
cross-patient consistency, UMI sensitivity, conservative annotations, English
figures, Chinese findings, final audit, logs, and rendered HTML.

