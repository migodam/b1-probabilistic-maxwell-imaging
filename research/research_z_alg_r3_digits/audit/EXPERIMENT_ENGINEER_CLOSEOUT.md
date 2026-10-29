# Experiment engineer closeout

Core numerical pipeline completed, with all failures preserved. Scientific advantage claims are **UNDER AUDIT** after the user challenged the custom Native/GN baseline. The frozen label “Full-J native” means this implementation's native material coordinates, not a reproduction of Route A or Route B.

## Main nine targets

- 24×24 real material pixels; all 576 unknown; 12 sources/48 receivers; k=2π; spacing .04; pulse radius .018. Full nonlinear scalar-TM disk VIE, noiseless primary observations.
- Local font generation only, no MNIST/downloaded data/template input to solver. Source/config freeze before holdout generation. Main solver inputs contain G/H/U0 and opaque Y IDs only. Private generator entropy remains excluded from seedless public truth.
- Original 27 predictions froze at 20:04:41 UTC. Remaining-budget supplement 27 froze at 20:08:25. Matched ordinary data-adjoint control 9 froze at 20:13:14. All 63 snapshots and source equality passed before first main GT scoring.
- Primary means: custom GN(m) RMSE 4.45928, reciprocal GN 4.45988, Gram dual + Full-J 3.95593. IoU .30143/.30082/.38014; SSIM .30565/.30477/.40911. These are poor on many digits. Native 4/6 and Gram 4/6/9 are visually legible; this is visual review, not a classification metric or all-nine success.
- Remaining-budget means: 4.46860/4.47003/3.95798. Ordinary data-adjoint control RMSE 4.86834. Its final576 initialization matches the original common coarse initialization exactly for all nine. The finite paired result does not establish superiority over mature Native implementations.
- Born-vs-full disagreement 3.54–8.32 confirms strong multiple scattering in this discrete study. Independent total-field endpoint checks have residuals ≤2.46e-15; condition estimates and opaque ZGECON work are recorded, not called free. No continuous-PDE or cross-discretization certification.
- Retained archives are optimizer iterates, not certified branches. All primary RMS≥1 greedy distinct counts are 1. Oracle best@K is retrospective and separately labeled.

Key files: `results/holdout_scores.json`, `budget_filled_scores.json`, `data_adjoint_scores.json`, `ordinary_attribution_comparison.json`, `holdout_coverage.json`; `audit/EXPERIMENT_ALL_PREDICTIONS_FROZEN_GT_RELEASE.json`.

## Requested image filenames

- `figures/digits_1to9_reconstruction.png`: exact copy of primary original 9×4.
- `figures/cost_vs_image_quality.png`: exact copy of primary cost/error scatter.
- `figures/holdout_with_ordinary_control_supplement.png`: clearly labeled 9×5 supplement/control view.
- `figures/contrast_regime_sensitivity.png`: post-GT physical regime diagnostic only.

GT-self-IoU at material threshold5 versus stored alpha≥.5 averages .89215. It is an auxiliary reference, not a global ceiling or a revised acceptance gate.

## Paired noise extension: partial, no aggregate score

Noise source/config froze at 20:12:52 UTC, before main GT. Parent GO 20:15:19. All same nine materials received new private circular random directions normalized to exact relative norm1e-3; no new PDE was needed because exact clean public Y was reused. This is not eighteen independent materials.

Parent reprioritized the Native audit. Queued launches paused at 20:25:04; both in-flight jobs completed. Fourteen result JSON files and all model counters were persisted; thirteen queued jobs were explicitly NOT_RUN. After verifying zero running numerical children, the stopped launcher was terminated at 20:31:46. No partial noise performance metric or `digits_noisy.png` was produced. Outer timing includes pause/cancellation and is not a normal completed batch runtime.

Definitive paid partial counters: `audit/EXPERIMENT_NOISE_PARTIAL_PAUSED.json` (14 runs, 1107.897175 CPU seconds, 1,610,748 physical complex RHS). Status: `audit/EXPERIMENT_NOISE_CANCELLATION.json`. Do not add partial totals again if reading individual root result accounting.

## Post-GT controlled contrast sensitivity

Two public failure-selected morphologies, original IDs a63/42, were transformed on all pixels by `m_low=1+1.3*(mGT-min)/(max-min)`. This changes the physical regime to [1,2.3]; it does not replace the original 8–18 main test. New generator, solver and configuration were snapshotted before the new anonymous Y. Generation paid 2 LU / 24 RHS / .0186423 seconds.

Custom GN, same generic prior and 576 unknowns, on digits8/9 gave residuals 2.1250e-4 / 1.8830e-4, RMSE .231967/.231460, normalized by 1.3 range .178436/.178046. Actual walls 8.59315/.55777 seconds under 45-second caps. Prediction snapshots precede diagnostic truth scoring. Both shapes are visible with amplitude bias; this does not prove the high-contrast implementation bug-free. Main IoU threshold5 is inapplicable to materials below2.3; diagnostic threshold1.65 is separately labeled and never substituted into main results.

`results/contrast_control_scores.json`; `audit/EXPERIMENT_PREDICTIONS_FREEZE_CONTRAST_CONTROL.json`. Independent validation adds 2 LU/24 explicit RHS + 2 opaque ZGECON calls/.0298687 seconds.

## Preservation

`audit/EXPERIMENT_FINAL_SOURCE_EQUALITY.json`: all 45 source/config comparison entries across primary, noise, contrast, budget and ordinary control match frozen bytes. No hashes were computed. No old algorithm files were edited after their freezes. Detailed costs are owned by the separate cost worker; this document does not replace the COST_LEDGER.
