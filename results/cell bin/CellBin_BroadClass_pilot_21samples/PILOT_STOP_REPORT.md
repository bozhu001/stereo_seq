# Broad-class pilot stop report

## Decision

**STOP before full 21-sample expansion.** The primary SingleR pilot passed ID and reference self-validation checks but failed the predeclared uncertainty gate. No formal full-cohort prediction Parquet or 21 spatial prediction maps were created.

## Evidence

- Pilot: 10,500 G10/U20 CellBins (500 from each of all 21 samples; all three chips).
- Frozen Reference18: 57,492 cells aggregated into 889 donor × unchanged-state profiles from 124 donors.
- Donor-grouped Reference18 broad-class self-validation accuracy: 99.8875%.
- Default SingleR pruning marked **0.0000%** of query CellBins Unknown (0 rows), so all sparse queries were forced into one of four classes.
- Median query maximum score: 0.0591; fifth percentile among correctly classified held-out reference profiles: 0.3627; ratio: 0.163. These scores are method-specific, not probabilities, but the separation demonstrates a large reference/query domain shift.
- Only 32.96% of forced predictions had at least two independent canonical markers for the predicted broad class. This is orthogonal supporting evidence, not the sole decision rule.
- Spearman correlation of maximum score with nGenes/nUMI was 0.564/0.560. Depth association was moderate rather than the sole failure mechanism.

## Interpretation

Reference18 broad classes are internally separable, but default SingleR confidence behavior does not transfer to these very sparse FFPE CellBins. Introducing an unvalidated fixed score cutoff or using space/density to fill labels would violate the task constraints. The pilot labels therefore remain diagnostic only and are not promoted as provisional cohort annotations.

## Possible recovery path requiring a new reviewed method decision

1. Build target-depth synthetic reference profiles by multinomially downsampling donor/state pseudobulks to observed CellBin UMI strata.
2. Estimate class- and depth-stratified error/abstention behavior using donor-held-out simulation.
3. Predefine an abstention rule from those held-out errors, then repeat a small pilot.
4. Expand only if the calibrated pilot yields a meaningful Unknown set, acceptable marker direction, and spatial plausibility without chip-depth domination.

This recovery is a new calibration analysis and was not run automatically after the primary-method gate failed.
