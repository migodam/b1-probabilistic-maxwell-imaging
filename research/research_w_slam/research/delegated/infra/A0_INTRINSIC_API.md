# A0 intrinsic nuisance stress: code ready, execution NOT_RUN

Frozen root definition: `A0_INTRINSIC_SUPPLEMENT.json`, created without loading any main A0 result or offline performance file. Runner: `tools/a0_intrinsic.py --root ROOT [--output PATH] [--ledger PATH]`, default output `results/a0/intrinsic_stress`. It must run in the parent-serialized A0 wrapper. Syntax/JSON parsed only by worker; no imports, tests, or physics executed.

Panel is fixed parents0,3,6,9, both noise repeats: eight conditions. Initialization examines saved selections from **all** main arms and recomputes the exact hierarchical online objective from saved base prediction, observed base data, unchanged proper support/grid-level prior and gain/phase prior. No GT or holdout informs selection. Intrinsic starts contrast3, radius.12; fixed loss.02. Main A0 files are read-only.

Eight physical/statistical coordinates: epsinf,ell,px,py,gain_offset,phase,contrast,radius. Contrast uniform bounds[2,5], radius[.08,.18]; other bounds/gain priors match main A0. Uniform normalizers are common and cancel in comparisons; no numerical damping is a prior. Support remains the saved main support. Public scaling is explicit. No jitter or clipping.

FullJ TRF has max_nfev6. Each newly requested physical point pays one full shared J (2 forward RHS +16 transpose RHS, one LU); its data/J are cached algebraically for gain/phase changes. Physical J columns `[eps,ell,px,py,contrast,loss,radius]` become residual columns `[gamma*first4, exp(i phase)*F, i*gamma*F, gamma*contrastJ, gamma*radiusJ]`; fixed loss column is omitted. All paid points, residual attempts, residual Jacobians and failure markers persist. The selected online candidate is the lowest unchanged-prior data objective among paid residual trials and the original saved-main no-op. Nonconverged optimizer status is retained even if a usable no-op or trial remains.

Before any fitting, headroom must be at least900RHS and60FE in the same A0 ledger, with positive CPU remainder. Otherwise write NOT_RUN; no borrowing. Frozen maximum: training864RHS48FE, eight independent outer-source holdouts16RHS8FE, finer parent0 both repeats4RHS2FE, total884RHS58FE. Finer audit uses n384/grid3, qualified primary n256/grid2, all h.17/order3/f1.

All eight online results and `ONLINE_COMPLETE.json` are terminal and durable **before** first offline label or holdout access. Independent full2TX16RX holdout scoring uses the original frozen64-real noise envelope. Separate `OFFLINE_SCORES.json` includes pose/material errors, geometry, support IoU/area, intrinsic errors, paired train/holdout compatibility and fixed finer prediction differences.

Outputs: root `RESULTS.json`, `ONLINE_COMPLETE.json`, `OFFLINE_SCORES.json`; each `parent_XX_repeat_R/` has `ONLINE_INITIALIZATION.json`, `PAID_XX.npz`, `ALL_TRIALS.npz`, `CURRENT_ONLINE_SELECTION.json`, `ONLINE_RESULT.json`, `PRIMARY_HOLDOUT.npz`, and fixed parent0 finer files. Pending physical markers remain on interrupted/failed actions. Existing output is refused; paid work is not automatically retried.

This is a **target-intrinsic nuisance stress** with fixed-intrinsic truth generator. It is **not unknown-m prior SBC**, a calibrated posterior claim, or scientific acceptance. Parent review remains required.
