# Frozen R4 shapes, data and independent scoring

`experiments/generate_diversity.py` constructs nine new glyphs and three analytic shapes without downloads. Odd digits 1/3/5/7/9 use local Arial Bold; even digits 2/4/6/8 use local Georgia Bold. Font names are in `FONTS` in that script. A glyph is cropped to its actual text bounds, fit with preserved aspect ratio inside width 0.64 and height 0.70 of a centered 384² canvas, then pixel-area averaged to 24². Ring, two disjoint ellipses and a thick circular arc use the same canvas and area averaging. This produces a fixed soft alpha and a fixed mask `alpha>=0.5` per shape.

For each shape a single background is sampled uniformly from `[0.08,0.15]`; its alpha and background remain fixed across contrast amplitudes `A=0.4,2.5,12`. Material is `m=background+A*alpha`. All methods share bounds `[.02,20]` and initial material `.1`. Shapes are deterministic from source/font; private random entropy controls backgrounds, complex noise directions and public record order. Noise uses one independent normalized complex Gaussian direction per shape, shared across amplitudes, scaled to `1e-3*||Y_clean||` at each amplitude. Each clean/noisy pair shares one exact nonlinear scalar-TM disk-pulse VIE forward result. The 36 generation solves use 36 full 576-dimensional LU factorizations and 432 source RHS, with 12 source / 48 receiver frozen R3 operators at `k=2pi`. The inversion and data generation do not use a Born approximation. For each already-computed clean state, the generator also forms `Y_Born=S@(m*E0)` by matrix multiplication and records `||Y_Born-Y_full||/||Y_full||` only in the private mapping. This adds no LU/forward solve and allows measured nonlinear stratification; low amplitude is not assumed to be Born-valid.

The 72 public NPZ files contain only `Y` and opaque `observation_id`. Shared operators contain no truth. The public manifest contains only anonymous observation IDs and GPU eligibility. Eligibility selects digit 2, digit 7 and ring, all three amplitudes and both noise levels: 18 paired observations. Public order is randomized and no name encodes a class, font, contrast, background or noise label. Do not upload `private_truth/`, which contains material/alpha/mask/cleanY and the private mapping including entropy. The remote upload must be limited to solver source, runner, public data and frozen run configuration.

Generation has not been run by preparation. After the parent's solver/config freeze and explicit generation GO:

```sh
python experiments/generate_diversity.py --go
```

Outputs are `public_data/{operators.npz,manifest.json,obs_*.npz}`, `private_truth/{mapping_PRIVATE.json,obs_*.npz}` and the label-free `audit/GENERATION_RECEIPT.json`. Existing observations are preserved. The generator prints costs and state residuals, not truth or private labels. No truth is sent to the solver.

After the parent closes execution, every planned task is recorded as complete, failed or NOT_RUN before offline truth scoring. The seal command records existing receipt paths and marks missing planned tasks NOT_RUN; it performs no copies, retries or extra physics. The parent preserves the closed prediction files. The independent scorer imports no solver and performs no new physical solves:

```sh
python score_results.py seal --manifest hardware/run_manifest.json --results remote_results --output audit/ALL_PREDICTIONS_SEALED.json
python score_results.py score --seal audit/ALL_PREDICTIONS_SEALED.json --truth private_truth --output results/scores.json
python make_figures.py --scores results/scores.json --truth private_truth --out figures
```

The scorer reports independent pixel RMSE, NRMSE=`RMSE/A`, MAE, and reconstruction-only Otsu IoU. Otsu uses 256 equal-width bins between the reconstruction's own minimum and maximum; constant reconstruction predicts an empty foreground. The fixed GT mask is never used to select this threshold. The offline known-background/known-amplitude threshold `background+A/2` gives a separately labeled midpoint IoU. These measures distinguish material accuracy and shape segmentation. SSIM is omitted. The primary endpoint is the solver's frozen public-physical-score top-1; GT never changes this selection or the solver. A separate offline oracle-best-K RMSE/NRMSE uses finite correction-phase optimizer candidates and explicitly does not count these as independent branches. Data success requires relative residual <=.001 clean or <=.003 noisy; material success requires NRMSE<=.15 and Otsu IoU>=.7; joint success requires both. Thresholds are read from configs/CONTRACT.json, and missing unrun tasks remain NOT_RUN.

Results include finite/failed counts, separately aggregated CPU and CUDA metrics, LU/RHS/wall/CPU costs including known failed-task costs, clean/noisy paired differences, and within-backend pairwise NRMSE wins/ties. The tie criterion is `abs(delta)<=max(.02,.10*min(NRMSE_left,NRMSE_right))`; every continuous paired delta is also retained. Missing/failed pairs remain explicit. No pooled CPU/GPU algorithm ranking or hardware speedup inference is made. Plates cover all CPU digits at every contrast and noise level, with non-digits in separate plates; CUDA plates show its fixed subset. Material scales are shared within a plate, with noise pairing and win/cost figures separately provided.
