# Session summary — CellBin SingleR broad-class pilot wrap-up

Date: 2026-11-11

## Scope

Completed the gated Reference18-guided SingleR broad cell-class annotation
pilot for the 21-sample Stereo-seq CellBin project and archived the user-
approved stop decision. No full-cohort prediction, new annotation method, or
Route A clustering was started during wrap-up.

## Inputs protected

- Joint CellBin V1 remained unchanged.
- Frozen Reference18 counts, state definitions, and cell mappings remained
  unchanged.
- Existing G10/U20 QC metrics and decisions remained unchanged.
- Six input paths, sizes, and SHA-256 values are recorded in
  `Input_provenance.tsv` and the approved local ELN copy.

## Pilot design

- SingleR 2.12.0 with 889 donor × unchanged-Reference18-state profiles from
  124 reference donors.
- Four traceable broad classes; no changes to the 18 frozen states.
- 10,500 query CellBins: 500 assigned G10/U20 CellBins from each of 21
  samples.
- Donor-grouped five-fold reference cross-validation.
- Default SingleR pruning retained as Unknown; no spatial filling.

## Decision evidence

- Reference donor-grouped accuracy: 99.8875%.
- Query Unknown fraction under default pruning: 0%.
- Query median maximum score: 0.0591.
- Correct held-out reference maximum-score P5: 0.3627.
- Query/reference score-domain ratio: 0.163.
- Predicted-class support from at least two independent markers: 32.96%.
- Maximum-score Spearman correlation with nGenes/nUMI: 0.564/0.560.

The uncertainty gate failed. The workflow stopped before expansion to all
359,329 CellBins. Pilot labels are diagnostic only and are not formal cell
types. The workflow is retained as an experimental project-reusable workflow,
not a frozen annotation method.

## Archive

- HTML report:
  `results/cellbin_broadclass_pilot_21samples_261111155125.html`
- Report SHA-256:
  `B242ED2A2B8748272205508EFFAE975DFFF85C87124790FAD87D6AFDE3731744`
- Analysis/report commit:
  `b36f44a98dd756cf3ccca27c409f8f97861e02b9`
- Approved local ELN copy:
  `elab/Reference-guided_Broad_Cell-class_Annotation_Pilot_261111155125.md`
- Remote eLabFTW ID/URL: unavailable because `ELABFTW_API_KEY` was absent from
  process, Windows user, and Windows machine environments; no ID was invented.

## Deferred work

- Do not start Reference18 depth-matched validation in this session.
- Route A (conventional unsupervised G10/U20 CellBin clustering pilot using one
  representative sample from K8, J2, and L3) is the next priority, but must be
  executed as an independent task with separate analysis and ELN records.
