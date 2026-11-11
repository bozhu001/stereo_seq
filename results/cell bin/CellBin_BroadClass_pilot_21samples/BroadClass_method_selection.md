# Broad-class method selection (pre-execution)

## Selected primary method

**SingleR 2.12.0 with donor/state-aggregated frozen Reference18 profiles**, followed by SingleR pruning. Reference18 labels remain unchanged; only the documented 18-state-to-four-broad-class mapping is used for the pilot.

## Why this method

- The CellBin query is FFPE and very sparse. SingleR uses rank/correlation evidence across genes rather than requiring a dense full matrix or de novo clustering.
- Aggregating reference cells within donor × Reference18 state reduces reference noise and memory while preserving independent donor profiles when the frozen metadata supports them.
- `pruned.labels`, `delta.next`, and per-class scores provide explicit uncertainty outputs. Pruned predictions are retained as `Unknown`; no spatial or nearest-neighbor filling is allowed.
- Scores/deltas are method-specific evidence, **not calibrated probabilities or accuracy estimates**.
- Query raw counts remain unchanged. Library-size log normalization is used only as a transient input representation required for comparison.
- The main risks are FFPE sparsity, ambient RNA, reference/query technology mismatch, and chip-depth effects. These are evaluated in the pilot and full validation.

## Alternatives reviewed

- **scANVI:** not selected because `scvi-tools` is not installed in the project environment, joint training would be substantially more resource intensive, and it would add a learned latent model whose depth/batch behavior requires a separate validation program.
- **RCTD/cell2location:** not selected for the primary method because their continuous mixture outputs should not be converted directly into single-CellBin identities; this task requests provisional broad labels with Unknown.
- **Canonical-marker scoring alone:** retained only as orthogonal validation, not as the classifier.

## Frozen inputs and audit results

- Reference18: 57,492 cells, 18 unchanged states, 18,966 raw-count genes.
- Query: 359,329 assigned G10/U20 CellBins from all 21 samples.
- Common unique gene symbols: 18,966.
- Reference donor/sample candidate columns: donor_id.

## Pilot gate

The small pilot will include balanced cells from all three chips and the six prespecified caution samples where possible. Expansion is allowed only if:

1. all IDs and score rows match the pilot input exactly;
2. reference self-validation is materially above chance at broad-class level;
3. pruning produces a nonzero Unknown set rather than forcing every label;
4. predictions are not nearly mono-class within a chip/sample without marker support;
5. score/delta behavior is not overwhelmingly explained by nGenes/nUMI;
6. sparse processing completes without constructing the 359,329 × 32,089 matrix densely.

## Planned outputs

Each provisional query row will retain the unique Joint V1 ID, original cellbin ID, chip/sample/disease fields, x/y coordinates, predicted broad class, all broad-class scores, delta/uncertainty fields, `unknown_flag`, and the explicit label `Exploratory / Provisional annotation`.
