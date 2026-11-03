# R5 trajectory diagnostic implementation and execution handoff

Status: complete for the frozen 8 replay runs and 6 reused historical states. Truth access: NONE. No corrector, formula identity tests, parameter scan, GPU use or publication. Scientific interpretation remains with the parent Codex thread.

- 168 replay checkpoints (21 per run), including 24 iteration500 pre-readout states with complete current physical Jacobian spectra and exact actual-step directional HVP.
- 6 reused historical fixed states: self_absorption and onehop_phase_pi, stages0/1/2. Their material, auxiliary J and Q are available; actual previous material, material step and PR trajectory are NOT_AVAILABLE. All actual-step relationships and directional HVP fields are null. A stage-to-stage difference is never substituted for a single step.
- 30 total endpoint/state spectra and 120 counterfactual material readouts (96 replay +24 historical), each evaluated with unchanged original physical G. Ordinary uses its original mature scalar CIE readout; self/plus/minus use the original onehop readout with Goff multiplier0/+1/-1. They share the public maximum60 inner steps and original compatible seed policies; they are not claimed to solve one identical nonlinear subproblem. Self is a formula-equivalent onehop implementation, not a byte replay of the historical blockself path.
- Frozen material LS uses all sources jointly with real material columns and fixed source metric diag(1/sum_cells |Einc|²). Physical data whitening is separately1/||Y||. No per-source material solve or material-truth score is used.
- Ordinary/self frozen metric transport and its driving RHS are analytically zero. Their canonical arrays/norms are zero and direction cosines/dots are null; raw subtraction/projection roundoff arrays remain for audit. Plus/minus transport vectors retain measured directions. No arbitrary cutoff was added.
- Endpoint exact full-Phi curvature, prior curvature, exact data curvature, data GN curvature and residual feedback are separate. No full physical Hessian was constructed. The inherited HVP reuses the retained Maxwell LU.
- All168 replay signed cross sums are positive in this sample (range0.000298917650938734 to0.10793615734972448). This is an observation, not a scientific gate or material-recovery statement.

## Artifacts

`SUMMARY.json`, `COST_TOTALS.json`, `RUNS.csv`, `ALL_CHECKPOINTS.csv`, `LEGACY_STATES.csv`, `COUNTERFACTUALS.csv`, and one JSON/NPZ pair per state. Per-run `CHECKPOINTS.csv` retains the wide checkpoint table. `ZERO_PHYSICS_CORRECTIONS_RECEIPT.json` records interpretation/label corrections using stored arrays, and `OUTPUT_AUDIT.json` confirms174 state summaries and missing-data labels.

## Cost and failed attempt

All successful replay diagnostics:408 LU /8928 complex RHS; all6 reused fixed states:30 LU /720 complex RHS. A failed first devglobal attempt caused by a local NumPy indexing error added1 LU /24 complex RHS, retained in that run's `prior_attempts` and `failed_attempts`. Its auxiliary counter was not returned before the exception, so that portion is explicitly incomplete rather than reconstructed. Total charged new physical diagnostics:439 LU /9672 complex RHS. Measured run/attempt wall time:9.685399s; CPU time:9.346890s. Process import startup is outside those timers and unavailable. Money:null. Extra Green, small4 Gram, source/block RHS and inner readout counts remain separate in per-state JSON and aggregate receipts.

The callback replay itself is charged separately. No historical source/result file was edited.
