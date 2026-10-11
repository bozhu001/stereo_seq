# Reference18 UMI50 overlapping composition-neighborhood analysis

## Scope

This analysis reused the frozen 50-um, center-included Reference18 UMI>=50 assay. It did not rerun RCTD, rebuild physical neighborhoods, or choose a discrete K.

## Construction

- Patient-balanced reference: 10,500 bins (500 per patient).
- Initial/refined unique anchors: 1,050/679.
- Primary overlapping composition-neighborhood size: 500 nearest bins; sensitivity: 300 and 800.
- Primary covered bins: 64,178/126,449 (50.8%); median memberships among covered bins: 3.

## Patient-level inference

- Primary tested composition neighborhoods: 262.
- Tests: 786 across three disease contrasts.
- Global HC3 FDR <0.05: 0.
- Global within-chip permutation FDR <0.05: 0.
- Fully supported findings: 0.

Disease labels were never used to construct anchors or memberships. Interpretation of the Reference18 composition was attached only after testing.
