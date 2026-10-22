# P0-B figure sources and limits

Generated: 2026-10-22T16:49:56.461087+00:00

These are saved corrected results, not new searches or field solves. Raw result files and the original ledgers are read-only.

## Current outputs

- 06_p0b_calibration_and_rhs: GENERATED_NOT_YET_VISUALLY_REVIEWED
- 07_p0b_offline_material_comparison: GENERATED_NOT_YET_VISUALLY_REVIEWED
- 08_p0b_fixed_probe8: GENERATED_NOT_YET_VISUALLY_REVIEWED

Only outputs listed as generated in the current manifest belong to this snapshot; earlier files may remain after a skipped/failed rebuild.

## Interpretation and completeness

- Images use a 256×192 display grid; original theta coordinates are exported separately. Continuous known truth is image-only and never selects a candidate.
- Standalone RHS uses method_rhs+80 for baselines and +160 for TASK_AWARE. These are offline attributions, not sums of actual shared ledger work.
- Fixed probe8 covers 1/24 tasks and uses the same target in both worlds; the recorded centroid differences are not a joint-compensation experiment.

- [excluded run] Earlier results/p0b/search (p0b-search-001) is invalid for method comparison and excluded: its corrector continued from a rejected proposal instead of the accepted best/no-op anchor. Its paid costs remain in the original ledger; these plots do not reassign them.
- [NOT_RUN] Seed 2026102202 / narrow: Insufficient remaining search quota for all four full method allowances after retained repair costs; remaining_rhs=480. No missing method results are filled in.
- [07_p0b_offline_material_comparison] Exported 18432 original theta entries separately. The coarse display does not recompute the original offline material scores.
- [08 raw arrays] Exported 2880 existing frequency/source/receiver rows. Plotted distances remain exactly the saved SCORED values, not recomputed values.
