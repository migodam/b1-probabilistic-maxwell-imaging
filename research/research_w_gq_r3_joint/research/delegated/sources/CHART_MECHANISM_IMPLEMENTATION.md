# Public 2×2 chart mechanism audit — prepared, not launched

Own implementation: `tools/chart_mechanism.py`.
Frozen before any PDE: `research/delegated/sources/CHART_MECHANISM_FROZEN.json`.
Status: **READY / PDE NOT_RUN**. Syntax compilation only; no physics smoke test.

## Predeclared comparison

Four matched public states: uniform/changed wall × absent/strong target. Strong geometry is identical in all four states; only target material changes from `(contrast,loss)=(0,0)` to `(6.8,.075)`. The absent state is a forward diagnostic/reference, not a legal reconstructed target under a positive-contrast inversion prior. Public changed wall and target numbers are frozen in JSON. Three shared frequencies and current common-mesh JointModel.

Rank16 chart is built solely at uniform+absent. The reference approximate Jacobian is reused from the chart's construction audit. Native full reference Jacobian uses 16 receiver adjoints per frequency, no additional FE factorization. Other states use a shared full response, full receiver adjoints, and fixed chart propagation. Expected **12 FE / 432 RHS**, hard local caps **18 FE / 600 RHS**; actual pre-action reservations remain authoritative.

State files save full/fixed complex J, identical exact responses, scaled real-white J, state coordinates, parameter scales, scaled prior precision, data Schur, step Schur, step Hessians/gradients, unconstrained steps and F29 residual consequence. Reference files save public U/spectrum/mass diagonal/audit arrays. Measured development y is copied into the run evidence; no private data paths are used.

Data Schur uses the public independent wall prior in scaled coordinates, zero target precision, no numerical damping. Step Schur/steps add **1e-3 I on all19 scaled coordinates**. These steps are unconstrained algebra diagnostics, not feasible/accepted reconstructions and not posterior precisions. Data-Schur eigenvalues are retained without clipping. F27 bound is computed from the full measured audit matrices; it is an offline bound, not a cheap deployable certificate.

The final four factorial contrasts concern relative error metrics; they do not assert additive independent scattering physics. Parent interprets mechanisms and validity.

## Launch only after parent queue grant

From R3 root, use existing guarded single-thread/RSS wrapper (the wrapper owns polled4GiB receipt; model methods own shared `.fem.lock`):

```text
python tools/run_suite.py chart_mechanism_public PYTHON tools/guarded_online.py --audit-log research/delegated/sources/CHART_MECHANISM_FILE_ACCESS.jsonl tools/chart_mechanism.py --run
```

Use the project Python environment as needed to invoke run_suite. Its child token `PYTHON` resolves the existing A25 environment. No new dependency installation.

Output: `research/delegated/sources/CHART_MECHANISM_RUN/`. Failure leaves partial arrays, per-state summaries and local receipt. Wrapper adds its independent receipt under `results/jobs/chart_mechanism_public/` and JOB_LEDGER; physics reservations also enter the existing PHYSICS_LEDGER. No additional rank or parameter sweep, Monte Carlo, or tiny-identity test.
