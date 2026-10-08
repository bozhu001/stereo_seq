# J2 K=5 spatial architecture audit

This directory reuses the frozen J2 K=5 assignment and evaluates all five niches symmetrically.

- Queen adjacency: same patient and ROI; |dx| <= 25 um and |dy| <= 25 um; self-pairs excluded.
- Each physical edge is counted once.
- The patient/ROI null uses 999 label permutations and preserves niche counts.
- Frozen aggregate inputs are not read.
- Disease is metadata only; no disease-significance test is performed.
- Adjacency is spatial association, not cell-cell communication.
