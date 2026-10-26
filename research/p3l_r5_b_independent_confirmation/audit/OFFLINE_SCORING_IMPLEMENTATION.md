# Offline scoring implementation — prepared, not executed

`src/score_predictions.py` is prepared for the parent orchestrator. It has not been executed, and no fresh truth file has been opened by this task. The frozen statistical protocol and evaluator are unchanged.

Before any offline truth access the scorer requires:

- `raw/predictions/ALL_PREDICTIONS_SEALED.json`, both group `ONLINE_COMPLETE.json` files, and `truth_read=false` in all three seals.
- Global `groups` exactly A/B, pointing to the corresponding group seal paths. Paths may be absolute or relative to the confirmation root.
- Each group has exactly the 13 canonical prediction archives (physics seed 2026102601, four learned methods × three seeds); optional `det_frozen` must have all three seeds in both groups.
- Global `artifacts` is the exact union of group artifact entries. Each entry has `path`, integer `bytes`, and full lowercase BLAKE2b-512 `blake2b` (128 hexadecimal characters, Python `hashlib.blake2b(content).hexdigest()`). All bytes are read, checked, and loaded from the same in-memory snapshot.
- Each group contains exactly the frozen 256 parent IDs, with no duplication or A/B overlap. Arrays must have finite `[256,576]` raw/clipped predictions and Boolean `[256]` fallback vectors.
- All anchor q0/chi0 are finite with frozen shapes; the sealed physics raw prediction equals the saved chi0 exactly.

Group seal `parents` accepts either count 256 or the explicit full parent ID list. The output directory must not already exist; interrupted or previous outputs are not implicitly overwritten. A `PRE_GT_VALIDATION.json` is written only after all online artifacts/anchors validate and immediately before opening the first truth file. Only then are `chi` and `q` loaded from `raw/offline_only/{parent}.npz`.

Expected outputs under `raw/offline_scores/`:

| File | Purpose |
|---|---|
| `PARENT_LOSSES.npz` | Frozen statistics schema: 512 parents, three seeds, five canonical methods, raw/clipped full-grid complex MSE, separate ensemble losses |
| `PARENT_METRICS.csv` | Per parent/method seed, seed-mean, and ensemble rows; raw/clipped full/ROI complex MSE/MAE, fallback flag and anchor phase cell counts |
| `PARENT_METRICS.npz` | Compact full/ROI metrics, phase counts, ROI proportions, phase squared-error sums, phase conditional per-cell MSE and full-grid contributions, fallback flags |
| `OFFLINE_SCOREPACK.npz` | Explicit offline-only full GT chi/q, anchor chi/q, raw/clipped predictions, ensemble predictions, anchor field/phase labels and ROI mask for plots |
| `SCORING_SUMMARY.json` | Descriptive metrics and phase decomposition; no confidence or scientific decision |

Phase labels depend only on the frozen anchor and original R4 MaterialPrior static config/points: `field=A@q0[:64]/halfwidth`. The original `field_mean=0` is checked. Outside ROI is its own label; within ROI, host plateau is `field<=-1`, transition is `-1<field<1`, inclusion plateau is `field>=1`. The four full-grid MSE contributions sum to the full-grid MSE. Empty phase conditional MSEs use zero storage plus `phase_cell_mse_valid=false` and cannot be interpreted as observed zero error. Summary phase MSE aggregates squared errors over phase cell counts. These are anchor-defined saturation diagnostics, not evidence of a Maxwell-specific transition mechanism.

`det_frozen`, when present, enters descriptive CSV/scorepack/summary outputs only. It never enters the five-method `PARENT_LOSSES.npz` or the frozen four-comparison family. Physics has one prediction repeated over three seed slots. No solver, training or model inference modules are imported.

After explicit parent authorization and all seals exist:

```sh
../../.venv_nn/bin/python research/p3l_r5_b_independent_confirmation/src/score_predictions.py
../../.venv_nn/bin/python research/p3l_r5_b_independent_confirmation/src/statistics_confirm.py \
  research/p3l_r5_b_independent_confirmation/raw/offline_scores/PARENT_LOSSES.npz \
  --output research/p3l_r5_b_independent_confirmation/raw/offline_scores/STATISTICAL_RESULTS.json
```

Validation status for this added scorer: `NOT_RUN_PARENT_REQUESTED_NO_EXECUTION_YET`. Existing statistics implementation tests are separate and remain as previously recorded.
