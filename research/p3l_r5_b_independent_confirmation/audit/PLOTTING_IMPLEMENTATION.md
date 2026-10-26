# Offline plotting — prepared, not executed

Execution update (2026-10-26 16:51 UTC): parent authorized plotting after all512 predictions were sealed and offline scoring/statistics completed. All eight figures were rendered and visually inspected; see `FIGURE_VISUAL_QA.json`. The preparation notes below preserve the pre-execution interface and constraints.

`src/plot_confirmation.py` is prepared for execution only after parent-authorized offline scoring. It has not been run, loaded the scorepack, or opened fresh GT. No model, physics, scoring or statistical code was changed.

Inputs are the completed `raw/offline_scores/SCORING_SUMMARY.json`, `OFFLINE_ACCESS_STATUS.json`, `OFFLINE_SCOREPACK.npz`, `PARENT_METRICS.npz`, and root `RESULTS_AND_STATISTICS.json`. The latter accepts the frozen statistics output either directly or under `statistics`. It must report `CONFIRMATORY_ANALYSIS` with 512 parents. Canonical methods/parent order are checked across scorepack/metrics. No bootstrap or loss is recomputed for scientific decisions; plots read frozen statistical intervals and descriptive scoring products.

Outputs default to the new, previously nonexistent `figures/confirmation/` directory:

| Figure | Selection / meaning |
|---|---|
| `FIXED_PARENT_MATERIAL_A/B` | First two manifest-ordered parents in each batch; GT and five raw prediction methods at seed2026102601; separate real/imag scales shared across both batches and all shown methods |
| `RISK_AND_PAIRED_BATCH_EFFECTS` | Pooled raw seed-mean risks; pooled/A/B paired effects and frozen98.75% intervals, with batches explicitly descriptive |
| `ANCHOR_PHASE_RISK_AND_GAIN` | A/B ROI phase shares, per-phase cell-weighted raw seed-mean MSE, full-grid gain contributions versus StrongCNN |
| `RETENTION_AND_PAIRED_CONTRASTS` | Four frozen contrasts, raw/clipped separated, retention one-sided lower bound explicitly marked, ratio only when denominator qualification permits |
| `RAW_CLIPPED_SEEDMEAN_AND_ENSEMBLE` | Separate raw/clipped panels and separate mean-seed-loss / ensemble-prediction-loss values |
| `OFFLINE_ERROR_SELECTED_DIAGNOSTICS` | Two worst distilled-minus-Strong raw seed-mean parent losses, explicitly selected after GT for explanatory diagnosis; ensemble images labeled separately |
| `CPU_PROVIDER_AND_PREDICTION_COST` | Only if completed `audit/timing/CONFIRMATION_TIMING.json` exists; provider/full prediction/load-plus-first timing categories and batch1/64 separate |

Every figure is written as PNG/PDF, with compact source CSV or NPZ beside it. Captions are saved in `FIGURE_CAPTIONS.md`; fixed IDs and cost availability are saved in `FIGURE_RECEIPT.json`. Optional `det_frozen` does not enter primary five-method figures. No generated figure has yet been visually reviewed.

The existing timing schema excludes physical acquisition and anchor construction, so these values are explicitly inference costs. Anchor-inclusive end-to-end cost is not fabricated or obtained by subtraction; the receipt records its absence. A new measured end-to-end table, if later supplied, needs a separately explicit plotting schema.

```sh
../../.venv_nn/bin/python research/p3l_r5_b_independent_confirmation/src/plot_confirmation.py
```

Execution/visual QA: `NOT_RUN_PARENT_REQUESTED_NO_EXECUTION_BEFORE_OFFLINE_AUTHORIZATION`.
