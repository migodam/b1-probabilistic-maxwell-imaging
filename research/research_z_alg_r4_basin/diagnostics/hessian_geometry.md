# Exact Hessian and box KKT mechanism experiment

Status: offline public-Y-only diagnostic. No GT is opened by this program. Exact data Hessian includes the residual-Hessian term; the regularizer Hessian differentiates smoothed TV exactly. The eigenproblem is restricted to coordinates more than 1e-6 from both box bounds. No exact-rank, second-order global optimality, or formula accuracy claim is made.

| Case | Free/lower/upper | Projected gradient infinity norm | Smallest algebraic eigenvalues | Eigensolver status | Best ray objective decrease |
|---|---|---:|---|---|---:|
| V1_holdout7_data_lbfgs | 31/537/8 | 9.537e-09 |  | MAX_250_MATVEC_PAID | 0 |

## V1_holdout7_data_lbfgs

Original objective=0.0532222738, relative residual=0.326231075. KKT lower/upper violations at gradient threshold1e-8: 0/0. The reported KKT residual is a raw metric rather than a newly imposed acceptance gate.

| Eigenpair | Exact curvature | GN data + exact regularizer | Residual-Hessian contribution |
|---|---:|---:|---:|

Paid cost: 1 LU, 6024 complex RHS; 250 HVP. Wall=1.3750s; CPU=1.3159s. Eigensolver matvec cap=250, tol=1e-3. Failed/partial work remains paid.
| old_holdout9_FullJ | 200/366/10 | 7.808e-05 |  | MAX_250_MATVEC_PAID | 0 |

## old_holdout9_FullJ

Original objective=0.0894633054, relative residual=0.42297396. KKT lower/upper violations at gradient threshold1e-8: 0/0. The reported KKT residual is a raw metric rather than a newly imposed acceptance gate.

| Eigenpair | Exact curvature | GN data + exact regularizer | Residual-Hessian contribution |
|---|---:|---:|---:|

Paid cost: 1 LU, 6024 complex RHS; 250 HVP. Wall=1.0339s; CPU=0.9756s. Eigensolver matvec cap=250, tol=1e-3. Failed/partial work remains paid.
| old_holdout8_FullJ | 43/533/0 | 0.0002001 | -5.6588e-06, -2.4991e-07, 1.2941e-07 | CONVERGED_TOL_1E_MINUS3 | 0.00028002 |

## old_holdout8_FullJ

Original objective=0.405040531, relative residual=0.900042514. KKT lower/upper violations at gradient threshold1e-8: 0/0. The reported KKT residual is a raw metric rather than a newly imposed acceptance gate.

| Eigenpair | Exact curvature | GN data + exact regularizer | Residual-Hessian contribution |
|---|---:|---:|---:|
| 0 | -5.65881e-06 | 8.643652e-05 | -9.209533e-05 |
| 1 | -2.499134e-07 | 1.210991e-06 | -1.460904e-06 |
| 2 | 1.29408e-07 | 3.766548e-07 | -2.472468e-07 |

A negative free-coordinate eigenvalue is a feasible local negative-curvature direction, up to the eigensolver tolerance. A strict saddle interpretation additionally requires sufficient first-order stationarity; the KKT metrics are retained so that this is not assumed. Finite rays evaluate the actual original regularized objective and use the scalar box-feasible interval, never coordinate clipping.

| requested signed RMS | actual signed alpha | actual objective decrease | relative residual |
|---|---:|---:|---:|
| -0.10 | -0.1 | -0.0004617285 | 0.9005554 |
| -0.50 | -0.33225 | -0.00143741 | 0.9016381 |
| -1.00 | -0.33225 | -0.00143741 | 0.9016381 |
| -2.00 | -0.33225 | -0.00143741 | 0.9016381 |
| +0.10 | +0.057502 | +0.0002800224 | 0.8997313 |
| +0.50 | +0.057502 | +0.0002800224 | 0.8997313 |
| +1.00 | +0.057502 | +0.0002800224 | 0.8997313 |
| +2.00 | +0.057502 | +0.0002800224 | 0.8997313 |

Paid cost: 9 LU, 3564 complex RHS; 142 HVP. Wall=0.7323s; CPU=0.7168s. Eigensolver matvec cap=250, tol=1e-3. Failed/partial work remains paid.

## Interpretation boundary

- V1 holdout7 and old holdout9: `INCONCLUSIVE_MAX_250_MATVEC`. The capped eigensolver returned no converged eigenpair; no nonnegative or positive curvature conclusion is drawn.
- Old holdout8: negative free-coordinate curvature is observed while projected gradient remains 2.00e-4. It is not a stationary-saddle certificate. The feasible positive direction ends at RMS .05750, and its finite objective improvement has a substantial first-order contribution. Nonlocal escape is not established.
- The eight requested signed-RMS probes collapse to three unique feasible alpha values. All eight actual forward/LU costs remain paid in the ledger; repeated boundary-limited values are not independent evidence.
- No GT was opened and no post-hoc GT scoring was performed.

## Authorized dense fixed-face representation followup

V1_holdout7_data_lbfgs: explicit 31×31 fixed-face Hessian from coordinate HVP, followed by one symmetric dense eigh. Smallest three eigenvalues: 1.760459765e-08, 2.411887526e-08, 3.546555008e-08.

Added paid work: 31 HVP, 1 base LU, 768 complex RHS; wall 0.52570s, CPU 0.47376s. The previous process had exited; this includes its necessary new base forward/adjoint, with LU reused by every HVP.

The original 250-matvec ARPACK failure remains in the ledger. This is one representation change, not a tolerance sweep, derivative accuracy comparison, or formula test. A nonnegative fixed-free-face Hessian does not certify full box-constrained local or global optimality. No further ray or GT test is performed.

old_holdout9_FullJ: NOT_RUN, free dimension 200 exceeds80. No added physical work.
