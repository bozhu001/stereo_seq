# Session summary — Reference18 harmonized K=7 robustness

Date: 2026-10-17  
Project: Stereo-seq lung spatial transcriptomics  
Analysis/report commit: `f2ce760d6b2c0ce85905ce0b88815969cc272864`

## Objective and decisions

The session resolved whether K8, J2 and L3 could use a common K=7
Reference18 Neighborhood NMF representation for cross-chip analysis. J2/L3
K=12 was explicitly deferred; it will be revisited only if the harmonized K=7
model fails a downstream requirement.

The user approved K=7 as the primary common coarse resolution. The seven
frozen cross-chip families are LEC-associated, Blood-endothelial,
AT1-associated, Perivascular-stromal, B-enriched/mixed-immune,
Ciliated-associated and Fibroblast-associated.

The user further confirmed that B-enriched/mixed-immune matches the project's
immune focus. It is retained as the primary immune-focused neighborhood axis,
while its score must not be reported as pure B-cell abundance.

## Analyses completed

1. `scripts/1143_reference18_neighborhood_nmf_k7_state_sensitivity.py`
   refitted disease-blind K=7 models after removing Basal+LEC or B, using the
   same patient-balanced reference, unpenalized objective, 50 restarts and
   three LOPO restarts.
2. `scripts/1144_k8_k7_state_sensitivity_convergence_extension.py` repeated
   the two K8 sensitivity models with 20,000 allowed iterations and the same
   50 seeds.
3. `scripts/1145_reference18_k7_sensitivity_evidence_review.py` used
   independent best factor matches to avoid forced Hungarian mismatches after
   one axis disappeared.
4. `scripts/1146_plot_reference18_k7_robustness_summary.py` froze the
   interpretation and created the final English six-panel PNG/PDF/SVG.
5. `scripts/1147_reference18_k7_robustness_report.ipynb` and
   `scripts/1148_render_reference18_k7_robustness_report.sh` generated a
   timestamped standalone HTML report in WSL tmux. Code input was hidden in
   the final HTML; the notebook preserves complete executable code.

## Key results

- The seven main K=7 families had a minimum pairwise cross-chip loading cosine
  of 0.912.
- Across both state-removal models, the six expected retained families had a
  minimum independent loading cosine of 0.926, minimum median-patient Bin50
  score Spearman of 0.852, and minimum median-patient Q90 hotspot Jaccard of
  0.526.
- Both K8 20,000-iteration models reached 100% restart convergence. The
  minimum 5,000-to-20,000-iteration factor cosine was 0.9986.
- The original B-enriched factor had no coherent single successor after B
  removal. Maximum candidate median-patient score correlations were 0.257,
  0.056 and 0.125 in K8, J2 and L3, respectively; maximum Q90 Jaccards were
  below 0.076.
- The B factor is interpreted as real B-associated signal plus co-localized
  immune composition and possible immune-state deconvolution cross-talk. It is
  not interpreted as pure B abundance or as merely a mislabeled T/NK factor.
- LEC is retained but sensitivity-qualified; absence of separate T/NK/Mast/
  Basal/AT2 factors at K=7 is not absence of those signals.

## Final artifacts

- Figure directory:
  `results/reanalysis/REFERENCE18_K7_ROBUSTNESS_SUMMARY_20261017_165039/`.
- HTML report:
  `results/reanalysis/REFERENCE18_K7_ROBUSTNESS_REPORT_20261017_165550/REFERENCE18_K7_ROBUSTNESS_REPORT_261017165550.html`.
- Report SHA-256:
  `25875d2c41ca164542338b644d8b6a1a069251e8a4d94a8f0feae2f57f931cca`.
- Browser screenshot visual QA passed; report content, tables, embedded figure
  and session information rendered without truncation or execution errors.

## eLabFTW and script-review status

The user separately approved the eLabFTW draft. Approved local copy:
`elab/Harmonized_K7_Reference18_Neighborhood_NMF_robustness_261017165550.md`.
No remote entry ID/URL exists because `ELABFTW_API_KEY` is unavailable. The
HTML is under 100 MB and should be attached through the web UI.

Script-review classification remains pending user choice between recurring
workflow and one-off audit; no skill-candidate issue was filed.

## Operational notes

The first Windows parallel run failed because the sandbox denied thread-pool
creation; the formal run used the established WSL `stereo` Python environment.
Two early HTML render attempts were superseded: the first used an incorrect
working-directory heuristic and the next exposed code inputs. The final
timestamped report at `261017165550` corrected both issues. A verified stale,
zero-byte `.git/index.lock` from a timed-out read-only Git status command was
removed after confirming no Git process was active; no analysis data were
deleted.
