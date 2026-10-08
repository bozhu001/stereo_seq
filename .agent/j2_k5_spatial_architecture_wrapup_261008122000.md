# J2 K=5 spatial architecture wrap-up — 2026-10-08

The approved Y40105J2 K=5 spatial architecture audit was finalized locally.
It reused 58,772 evaluable Bin50 and the frozen K=5 assignment, constructed
218,614 unique within-patient/ROI Queen edges, and analyzed all five niches and
all ten heterotypic pairs symmetrically. The 999-permutation patient/ROI null
was completed for all seven patients; a downstream ordering error was fixed by
reusing the complete permutation checkpoint rather than rerunning it.

N2, N4, and N1 had the highest median largest-component fractions, but all
five label fields were fragmented. Nine of ten heterotypic pairs were
direction-consistent in at least six of seven patients, all negatively relative
to the exchangeable-label null; no recurrent positive pair was found. Results
are spatial association evidence, not independent validation or communication.

The standalone report is
`results/J2_K5_SPATIAL_ARCHITECTURE_REPORT_261008121443.html`; analysis/report
commit is `793085cee696e2432ce055212a286229d4429232`. The approved local ELN is
`elab/J2_K5_spatial_niche_architecture_261008121443.md`. No remote eLabFTW
entry exists because `ELABFTW_API_KEY` was unavailable.

J2 discovery and spatial results were frozen using
`results/reanalysis/J2_UNBIASED_NICHE_RESOLUTION_20261008_110522/J2_NICHE_RESULTS_FREEZE_20261008.tsv`
and `J2_FROZEN.ok`. The user directed that no further J2 analysis be performed.
