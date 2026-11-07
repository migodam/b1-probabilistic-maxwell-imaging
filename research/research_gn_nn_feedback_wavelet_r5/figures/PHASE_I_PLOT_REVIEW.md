# Phase I rendered figure review

All outputs are **Phase I, historical strict-GN inputs, DEVELOPMENT VAL used for
checkpoint/model selection**. Image/mother exclusion does not imply unseen
classes/fonts. No fresh TEST is represented and Phase II noise-aware COMMON
inputs are excluded.

The plotter ran with the existing Gaussian `.venv_nn/bin/python` and BLAS thread
limits of one. Only plotting and metadata were performed; no model fitting,
inference or physical corrector was rerun. TRAIN truth and TEST data were not
opened. Completed ordered predictions were available before VAL reference truth
was read.

## Rendered artifacts

- `01_actual_architecture.png` / `.pdf`: actual imported classes, all seven
  wavelet bands, zero inverse-layer LU, separate full nonlinear loss/adjoint,
  fixed-prior and free-control distinctions.
- `02_train_val_curves.png`: separate TRAIN and VAL panels, actual completed
  epochs. Four arms have 12 history entries; `haar_closed_no_fullfield` has 7.
- `03_prespecified_val_reconstructions.png`: all five arms plus GT/GN for the
  first reference-ordered A=.2 digit2, digit8, islands and primitive_composites
  examples, each clean and noisy. Indices are 56,57,72,73,88,89,80,81; no missing
  prescribed group and no success-based replacement. Fixed display [.02,.55],
  saturation counts and exact exposure flags are shown.
- `gallery_01.png` through `gallery_26.png`: all 156 ordered raw observations,
  with all five arms and failures retained.
- `04_absolute_haar_error_energies.png` and
  `HAAR_ERROR_BAND_ENERGIES.csv`: absolute contributions divided by 576, with all
  coefficients retained. Maximum absolute band-sum minus voxel-MSE difference:
  **5.960464477539063e-8**.
- `05_quality_residual_actual_cost.png`: separate clean/noisy quality, actual
  physical residual/qualification, full-LU TRAIN+VAL receipts and reported
  prediction-forward wall time, with distinct cost scopes.
- `07_raw_vs_corrected_matched_subset.png`: exactly matched **72 observations
  per arm / 360 total**, using the recorded chosen online corrected result.
  Best-K/oracle scores are not used. Raw full cohort has **0/780 qualified**;
  matched raw subset has **0/360**; the separately charged corrected subset has
  **360/360**. Qualification is not conflated with material quality.
- `06_depth_diagnostics.png`, `_clean.png`, `_noisy.png` and
  `DEPTH_PLOT_SCOPE.json`: declared six-observation subset per arm (three clean,
  three noisy), depths 1,2,4,8,16. This is not all156 VAL. Only the depth plotting
  command was rerun after adding the explicit subset labels and noise split.

## Visual checks and limits

Viewed the architecture, main raw reconstruction panel, Haar energies,
TRAIN/VAL curves, raw quality/cost, corrected comparison, all three depth
summaries, and gallery page01. Fixed actual overlap between the reconstruction
header and first row, and between full-LU annotations and subplot title.
The remaining 25 gallery pages were generated but not individually viewed.

Raw NN qualification is never inferred from corrected qualification. The depth
plots distinguish fixed-energy/constitutive/state-change diagnostics from
actual physical residual and material error, and explicitly withhold a shared
fixed-energy guarantee from the free cascade. Scientific claims and acceptance
remain with the parent review.

`FIGURES.json` records scope, all outputs, actual history lengths, input byte
sizes/modification times, matched corrected IDs and the band-energy discrepancy.
Use `plot_development.py --depth-only` for a later depth CSV update; it does not
read material truth, model checkpoints, raw predictions or run physical solves.
