# R2 CIE global-source implementation and execution audit

## Reused source and exact change

`src/cie_global.py::run_cie(model,Y,m0,config)` is a narrow extension of
`research/research_fbg_r3_experiments/src/cie_tsom_baseline.py::run_cie`.
The existing `sourcewise_cg=False` already supplies global scalar PR+ and
scalar exact line minimization. This work does not invent or rename that
algorithm. It supplies the missing full fixed source-metric contractions and
a source-rotation invariant coefficient stopping diagnostic. Non-equivariance
of a sourcewise update is not evidence of worse inversion.

The extension directly imports the mature `fourier_apcs`, `_regularizer` and
`_readout`. Receiver DPCS, projected four-corner IDFT APCS, native possibly
redundant/nonorthogonal Fourier coefficients, exact coefficient zero padding
between masks, PR+ and non-descent resets, physical real material box, dual-beta
transport and physical material prior are preserved. It is explicit dense
Fourier work, not a claim of FFT complexity. FFT calls are zero.

New `data_metric` and `state_metric` config entries are fixed Hermitian SPD
source matrices. Legacy `data_weights` and `state_weights` remain accepted.
Defaults are the same per-source energy normalizations, evaluated once from
original Y/E. For residual matrices D/R the objective is
`0.5 Re tr(D* D Omega_d)+0.5 rho Re tr(R* R Omega_s)` and the current gradient is
`SU* D Omega_d + rho H* R Omega_s`. Full-metric material readout reuses mature
real/prior readout after source Cholesky whitening: for `Omega_s=L L*`, pass
`J L,Q L` with unit diagonal weights. This realizes the full weighted shared
real least squares, not sourcewise material fits or a prior on transformed r.

For `E'=EU,Y'=YU` with unitary U, caller supplies
`Omega_d'=U*Omega_d U` and similarly for state metric. It must not renormalize
mixed source columns. Global PR and global `||Delta A||F^2/||A||F^2` diagnostic
are covariant under these unitary coordinate changes. The helper
`transport_source_metric` computes the general residual-metric congruence
`B^-1 Omega B^-H`; that metric identity does **not** assert optimizer covariance
for nonunitary B under an untransported Frobenius coefficient metric.

The full-metric extension rejects `sourcewise_cg=True` with off-diagonal
metrics; it does not pretend an independent-column line search solves a
coupled objective. Main native comparison instead invokes the unmodified
mature source itself. A temporary forwarding-only `_readout` observer copies
J/Q/m at steps 500/1000/1500 and restores the original function afterward;
it does not change returned numbers. The extension returns copied
`stage_states`, including J/Q/m/coefficients, and supports `stage_callback`.

An optional `current_solver='exact'` is a **tiny diagnostic only**: whitened
stacked Kronecker least squares solved with rank-aware LAPACK GELSD. It retains
native Fourier coordinates and guards coefficient unknowns at 256 by default.
Its SVD, RHS, matrix entries and stationarity are recorded. No such solve is
silently added to the main PR arms, and this code is not a scalable exact
current solver.

## Frozen main protocol and runner

`experiments/run_blind.py` reads only the public manifest and an explicit
allowlist of payload arrays. It never imports the generator or opens truth,
clean observations, shape parameters, mapping or labels. It preserves any
existing terminal disposition rather than retrying or overwriting it.

The parent froze three Fourier stages: mf=3/beta=1; mf=6 and10 with
beta_current=1,beta_material=.2; 500 iterations each, miniter=50,
dpcs_rtol=1e-3, material inner maxiter=60, coefficient_tol=0 for both native and
global. Thus both run all 1500 alternating iterations and coefficient change
is diagnostic only. There is no wall cap. Native/global capacity and iteration
count match; wall and actual operation counts can differ.

Four sequential CPU1 arms use the same public uniform host initialization,
physical real box [.02,20], original scalar-TM nonlinear Maxwell model and
historical smooth/TV material prior:

- Native mature sourcewise dual-beta CIE.
- Global PR flag with the same original diagonal objective metrics and full
  metric capability; no rotated data are substituted in the main inversion.
- Direct original physical L-BFGS-B, up to1000 iterations.
- Known full-current variable-projection control supplied in
  `src/profile_current.py`, fixed beta=1, up to250 L-BFGS-B iterations. Its
  uniform total-energy source metric differs from native/global and is
  explicitly a separate control, not an attribution solely to representation.

Each arm selects its pre-corrector incumbent by the original physical objective
including the unchanged prior, **always including the common initial material**.
Each then gets the same original physical L-BFGS-B corrector allowance of1000
iterations/4000 function evaluations. Direct therefore has up to1000+1000,
with both passes paid. Stage endpoints/selected predictions are stored; the
runner does no GT selection, scoring or image acceptance.

## Numerical check and accounting

The one requested new profile-gradient wiring check uses existing genuine R1
4x4 public operators, one disclosed constant root .55, initial .1, and one
direction at step1e-5. Its analytic directional derivative is
2.79450336446818e-4 versus finite difference2.79450335788356e-4 (absolute
difference6.58e-13). `audit/PROFILE_GRADIENT_WIRING.json` retains 3 exact-current
Cholesky factorizations/18 current RHS and the paid disclosed-root physical
forward (1 LU/6 source RHS). It is a new objective-wiring check, not a repeated
historical identity suite. This does not validate optimization superiority.

Every main case writes its own `results/blind/<opaque_id>/<arm>/RECEIPT.json`,
auxiliary/corrector histories, stage arrays and sealed selected prediction.
Receipts include process wall/CPU, original physical objective and data
residual, the supplied99% whitened noise radius, physical model LU/RHS,
physical operation wall/CPU, auxiliary current/inner costs, and failures.
No material error is computed before independent scoring.

The global extension records full source-metric products, metric Cholesky,
current-objective and gradient evaluations, PR/restarts, real-prior inner
evaluations, IDFT/SVD/basis columns and original Maxwell model deltas. The
legacy `reduced_operator_columns` is a historical work counter rather than an
exhaustive count of additional diagnostic contractions; added objective
logging costs remain in full wall/CPU. Physical RHS counts are actual retained
LU/adjoint solves and are not conflated with auxiliary LS RHS. Failed execution
receipts attach to the original exception and can be written/callbacked.

Main timers include public load, setup, auxiliary and corrector work, physical
candidate reranking, stage/prediction serialization and intermediate JSON.
They exclude imports and their own final receipt write. Geometry/operator
assembly is frozen historical work, not rerun; fresh generator cost is reported
separately by the parent. No hashes, GPU, installation or historical edits are
performed.

First native case elapsed10.58 wall/10.50 CPU seconds and sealed a prediction;
its original whitened data residual24.3966 is below supplied radius25.1654.
That is a data-only observation, not material reconstruction acceptance. Main
outcomes and final scientific interpretation remain the parent's review.

All12 main dispositions completed as SEALED_PREDICTION, with no main failures.
Readback confirms every native/global trajectory performed exactly1500 steps
and three MAXITER stage endings; this is not convergence. Summed sequential
case receipts are90.779540584 wall seconds/91.659207 process CPU seconds,
9279 physical LU and222048 actual complex physical RHS (111348 source plus
110700 adjoint). Generator, gradient-wiring and historical assembly costs are
separate; source/auxiliary work must not be hidden behind these physical counts.
The sum excludes the Python-process imports/startup and final receipt writes.
Parent owns the global cost closure including any separate diagnostics.
