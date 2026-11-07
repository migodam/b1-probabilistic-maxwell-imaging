# Both-block matrix-feedback exposed DEV control

The isolated R11 copy replaces only the current residual H/T maps and matched adjoint, while the material block retains the full moving R10 matrix readout. Data metrics, native Fourier/APCS, Green FFT, original PR directions, exact quadratic line, material degrees of freedom, stage betas, online public-Phi selection and corrector remain inherited. Current beta=1 and material beta=1/.2/.2/.2 are deliberately different: no common alternating objective is claimed. No frozen six-arm or R8 source file was edited.

The actual production helper passed its one necessary wiring check: adjoint relative error 2.87e-15 and current-objective directional derivative relative error 1.52e-7. The check uses no Maxwell forward, adjoint or GT. Every one of the 3000 saved current steps in each case decreased its fixed-m quadratic objective; this is a saved-history audit, not an extra optimization or physical evaluation.

All three local DEV runs completed and sealed predictions before GT scoring. The first two used two lanes; the third entered the ready lane. BLAS threads=1. The 480 CPU cap did not trigger, no run failed and no retry occurred.

|Case|Amplitude NRMSE|Fixed IoU|Phi|Material / CN99|Exact child CPU s|Corrector stop|
|---|---:|---:|---:|---|---:|---|
|r10_twod_0001|0.105166|0.713287|5.2520686e-05|PASS / PASS|101.890481|STOP: TOTAL NO. OF ITERATIONS REACHED LIMIT (3000)|
|r10_twod_0004|0.515137|0.086957|0.17393422|FAIL / FAIL|60.735483|CONVERGENCE: NORM OF PROJECTED GRADIENT <= PGTOL (536)|
|r10_twod_0022|0.100878|0.706349|5.5032405e-05|PASS / PASS|62.382625|CONVERGENCE: RELATIVE REDUCTION OF F <= FACTR*EPSMCH (611)|

The bridge succeeds in this exposed three-object DEV control after both lagged-only changes failed it. The shell remains at a high data-residual stopping point. The annulus remains successful. These observations support only this scoped comparison; they are not universal recovery, global basin, novelty or convergence claims. The annulus corrector reached its 3000-iteration limit. The bridge stopped on relative-F reduction, which is not by itself a projected-gradient certificate.

Exact external wait4 CPU for the three children totals 225.008589 seconds, including import/setup/inversion/serialization/exit. Per-run near-end process_time and 60/120/240/480-second paid-incumbent checkpoints are also retained. Their clocks serve different purposes and are not added. Controller and wiring CPU are separately recorded. Physical operations total {'cached_operator_imports': 3, 'FULL_FORWARD': 4314, 'actual_complex_RHS': 103380, 'LU': 4314, 'source_RHS': 51768, 'full_factorization_dimension': 2484864, 'DATA_ADJOINT': 4301, 'adjoint_RHS': 51612}. No physical operation occurs during scoring.

Current H construction reuses one retained 19x19 Woodbury LU within each current step. The actual post-material current/material history losses use their appropriate moving matrix denominator; equal betas share one factor. 3000 current feedback factors and 5500 history small factors per run are included in costs. No new dense Maxwell factorization was added for logging.

Artifacts: configs/twod_dev_bothblocks_tasks.json; src/cie_feedback_current.py; experiments/diagnose_bothblocks_run.py; RESULTS/exposed_dev/*__fullGreenBothBlocks__480CPU/; BOTHBLOCKS_ALL3_SEALED.json; BOTHBLOCKS_WIRING_CHECK.json; BOTHBLOCKS_CURRENT_STEP_AUDIT.json; BOTHBLOCKS_QUEUE.json; BOTHBLOCKS_DEV_SCORES.json/.csv; BOTHBLOCKS_DEV_CHECKPOINT_SCORES.json; BOTHBLOCKS_DEV_COST_LEDGER.json. Startup source/config snapshots are preserved in bothblocks_run_logs/.
