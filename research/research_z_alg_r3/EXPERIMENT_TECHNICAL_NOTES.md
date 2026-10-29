# Frozen experiment and replay notes

The final finite result is one additional fresh V7 data fit: GuardedGlobalDual fits 3/4 NEWCONFIRM12 targets at the predeclared all-data relative residual threshold of 1e-6, versus 2/4 for GuardedReciprocal. Both meet the stricter optimizer stop goal of 1e-10 on 2/4. No material-RMS success gate was predeclared. Report the actual RMS values; do not replace the primary residual-selected top1 with offline oracle best@K. Independent-paper novelty/support remains HOLD.

All physical work stopped at the science freeze, 2026-10-29 18:19:17 UTC, 33.4 minutes into the 90-minute sprint. Later operations are documentation, plotting, disclosed-truth scoring and publication packaging. No N24, scale study, new neural network, GPU or Full66 experiment was run.

## Evidence entry points

- [RESULTS_SUMMARY.json](RESULTS_SUMMARY.json) aggregates V6/V7, DEV costs, definitions and scientific limits.
- [V6 summary](results/HOLDOUT_V6_SUMMARY.md), [V6 per-case CSV](results/HOLDOUT_V6_METRICS.csv), [V7 summary](results/CONFIRM_V7_SUMMARY.md) and [V7 per-case CSV](results/CONFIRM_V7_METRICS.csv) preserve every image outcome.
- [Cost JSON](COST_LEDGER.json), [cost CSV](COST_LEDGER.csv) and [cost overview](COST_LEDGER.md) include every original paid run, target generation, independent validation, phase work and known missing timers.
- [V6 protocol](BLIND_EXPERIMENT.md), [V6 frozen contract](experiments/HOLDOUT12_CONTRACT.json) and [fresh V7 contract](experiments/CONFIRM_V7_CONTRACT.json) establish what was fixed before targets were generated.
- [V6 image comparison](figures/holdout_materials.png) uses the fixed reciprocal FullJ control for every row and retains the GlobalDual easy-case failure. The earlier native-control presentation is archived in `figures/v6_native_control/`.
- [V7 image comparison](figures/confirm_v7_materials.png) uses the fixed GuardedReciprocal control, the actual GuardedGlobalDual top1 and a clearly labeled GT-selected oracle-best column. Both V6 and V7 have separate cost/coverage figures and vector PDF exports.

## What is blind, and what is now disclosed

The geometry, public material box [0.5,30], acquisition and spatial-blob prior are known. Every one of the 12 material values is unknown to each solver. There is no hidden support supplied to a decoder, no mask-specific image fit counted as an all-data fit, and no target material or target-generation seed input to any solver.

The two DEV targets are exploratory. V6 used four fresh targets after its source/config freeze. Once scored, their outcomes were development evidence for the separately frozen V7 guarded revision. V7 used four new targets generated only after the V7 freeze; V6 targets are never counted as fresh V7 evidence. All 20 V6 predictions and all eight V7 predictions froze before their independent validation and material scoring. Original predictions and scores remain intact.

After V7 scoring, ten material-only NPZ files were exported under [truth_after_freeze/](truth_after_freeze/README.md). Each contains exactly the array `m`. Random target-generation seeds remain in the ignored local `private_truth/` tree; the export code never reads the seed values. Public deterministic optimizer seeds in configs are separate. A reader who has this disclosed package is replaying exposed targets, not conducting a new blind trial.

## Replay the recorded scores without a physics rerun

Run from this R3 directory using Python with NumPy installed. The exact local versions are recorded in [audit/EXPERIMENT_RUNTIME.json](audit/EXPERIMENT_RUNTIME.json); [experiments/requirements.txt](experiments/requirements.txt) records NumPy, SciPy, threadpoolctl and Matplotlib versions used in the sprint. No package installation was needed for this work.

```sh
python experiments/score_predictions.py \
  --split confirm \
  --freeze-manifest audit/CONFIRM_V7_FINAL_FREEZE.json \
  --truth-root truth_after_freeze \
  --disclosed-truth-replay \
  --output /tmp/r3_v7_replay_scores.json
```

The analogous V6 command is:

```sh
python experiments/score_predictions.py \
  --split holdout \
  --freeze-manifest audit/HOLDOUT_V6_FINAL_FREEZE.json \
  --truth-root truth_after_freeze \
  --disclosed-truth-replay \
  --output /tmp/r3_v6_replay_scores.json
```

The output path must be new. The replay flag permits changed modification times after a Git clone, checks recorded prediction byte lengths, and explicitly marks the result as disclosed replay. It cannot establish fresh blindness. No SHA checks are computed. The V7 disclosed replay has been executed: every scored row is exactly identical to the original frozen score, with zero physics RHS. Its receipt is [results/confirm_v7_public_replay_scores.json](results/confirm_v7_public_replay_scores.json).

Re-running a solver on these now-disclosed observations is a new numerical run on exposed data. The public CLI and frozen config are retained for a reader who chooses to do so; this sprint did not rerun physics after the science freeze. For example, `run_portfolio.py` accepts `--public`, `--method guarded_globaldual` or `guarded_reciprocal`, `--config`, and a new `--output`. Its source never imports a scoring module or opens either truth directory.

## Independent validation and costs

`experiments/validate_predictions.py` does not import the solver, `Model`, reduced kernel or a retained factorization. For every archived material candidate it independently forms `A=I-diag(m)G`, solves `A j=diag(m)E0` with NumPy, and evaluates `S j`. It checks the public box, original all-data residual, full resolvent singular gap, state equation and constitutive residual. The physical operator construction is shared through the public cached `G`, `E0` and `S`; this is an independent full solve, not an independent discretization or continuous-PDE validation.

All 3,500 V6 and 968 V7 archived candidates pass these checks. Exact raw dual endpoints that had merged with nearby archive entries needed 49 additional independent solves for mechanistic analysis. The independent validation cost is 13,551 full complex RHS in total. A valid forward evaluation does not make an inaccurate image successful.

Every full RHS count is an actual complex column for a 12-dimensional full system. Reduced algebraic RHS and GlobalDual inner least-squares RHS are separate. GlobalDual inner QR solves are 12-by-9 complex least-squares systems with three receiver RHS per call. Their cost is explicit; no conversion to full-PDE-equivalent RHS is asserted. Full-J correction work is included, not hidden behind proposal generation. [The clarified V6 mechanism analysis](results/holdout_globaldual_mechanism_v2.json) distinguishes raw endpoints, archived near-duplicate representatives and exact paid-trial improvements.

Across all DEV/V6/V7 solver work, target generation and independent validation, the paid total is 1,120,095 full complex RHS, plus 23,316 reduced RHS and 289,302 algebraic complex least-squares RHS. These categories are not interchangeable. The all-solver numerical/outer wall total is 63.008937 s. Original exploratory runs, failed attempts and both public/private score operations remain charged. The first target-generation attempt failed before generation and consumed zero RHS; its elapsed time was not retained. Legacy DEV scoring elapsed times were also not retained.

The V7 total outer budget is 3 s, including public loading, the common guard, repeated setup, continuation and wrapper bookkeeping. Continuation receives only measured remaining time. Actual terminal/bookkeeping overshoot is retained, with maximum 0.0128075 s. Full subprocess wall includes interpreter/imports and result file IO and is separately recorded; it contains the numerical wall and must not be added to it. Generation Green-kernel setup and receiver SVD are charged once in generation and shared by the methods. Stage subtotals and the all-solver total are overlapping views and must not be added together.

The shared V7 guard produces identical materials, residuals and operation counts in each paired case. All four fresh cases invoke continuation; none stops at the guard. Every V7 root accounting ledger exactly equals the sum of its phase ledgers. This confirms the paired comparison and cost boundary, not a general performance claim.

## Interpretation limits

This is a noiseless finite scalar-TM disk-pulse VIE with 12 unknown real cell values and nine complex observations. It retains nonlinear multiple scattering. Four fresh V7 cases, one additional data fit and mixed V6 outcomes do not establish global enumeration, generic material recovery, a global optimum, performance at scale, a continuous Maxwell result or an independent-paper novelty claim. Oracle coverage is a diagnostic, while the displayed top1 images remain fixed by original-data residual ranking.
