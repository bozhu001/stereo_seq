# L3 K=8 spatial architecture

本目录复用冻结K=8 assignment和x_um/y_um物理坐标。

- Queen adjacency: same patient and ROI; |dx| <= 25 um and |dy| <= 25 um; not the same Bin.
- Only the four unique offsets (25,0), (0,25), (25,25), (25,-25) are constructed, so each physical edge is counted once.
- Fragmentation index = 1 - largest component fraction.
- Spatial coverage = niche Bin count / all evaluable Bin count in that patient.
- Bin50 area is reported as 625 um2 because the frozen physical conversion defines each Bin50 as 25 x 25 um.
- Row-normalized adjacency uses incident edge ends: a homotypic edge contributes two ends to its source niche; a heterotypic edge contributes one end to each source niche.
- Null: 999 label permutations within each patient/ROI; seed 2026100711; log2 offset 0.5.
- Empirical P is two-sided from upper/lower permutation tails with +1 correction; BH is performed separately within each patient for 28 heterotypic and 8 homotypic pairs.
- The exchangeable-label null does not preserve niche-label spatial autocorrelation or overlap among adjacent 50 um input neighborhoods; enrichment therefore describes departure from label exchangeability, not independent validation of a biological interface.
- Frozen aggregate inputs were not read. Disease was not used in graph construction or testing.
