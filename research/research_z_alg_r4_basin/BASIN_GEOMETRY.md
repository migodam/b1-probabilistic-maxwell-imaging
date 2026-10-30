# Basin geometry: DEV7 versus holdout7, with 9/8 representatives

Status: `OFFLINE_DIAGNOSTIC`; truth was already exposed in frozen history. These runs describe four finite 576-cell scalar-TM problems. They are not fresh blind evidence, rank certificates, exact alternative-solution certificates, or global basin proofs.

The strongest contrast is between DEV7 and holdout7: DEV7 reaches a small-residual endpoint while all three historical holdout7 methods stop near residual 0.335. The GT line for holdout7 has a larger residual rise and a later adverse-gradient region. On holdout9 the historical GlobalDual continuation reaches a small-residual endpoint while Full-J fails, while holdout8 remains a failure for all three methods.

## Frozen actual endpoints

| Case | Method | Relative residual | GT material RMSE | cos(-physical gradient, GT-m) | cos(-small-Gram gradient, GT-m) | resolvent norm estimate |
|---|---|---:|---:|---:|---:|---:|
| DEV7 | fullj_native | 0.000266983 | 0.8711 | 0.019 | 0.005 | 9.429 |
| DEV7 | fullj_reciprocal | 0.000266918 | 0.8719 | 0.020 | 0.005 | 9.429 |
| DEV7 | globaldual | 0.000225079 | 0.8136 | 0.017 | 0.004 | 9.462 |
| holdout7 | fullj_native | 0.334945 | 5.4023 | -0.346 | 0.028 | 6.293 |
| holdout7 | fullj_reciprocal | 0.334969 | 5.4017 | -0.346 | 0.028 | 6.293 |
| holdout7 | globaldual | 0.335544 | 5.3928 | -0.346 | 0.027 | 6.261 |
| holdout9 | fullj_native | 0.422974 | 4.7774 | -0.251 | 0.140 | 12.778 |
| holdout9 | fullj_reciprocal | 0.424073 | 4.7637 | -0.249 | 0.140 | 12.781 |
| holdout9 | globaldual | 0.000228314 | 0.8849 | 0.021 | 0.003 | 12.250 |
| holdout8 | fullj_native | 0.900043 | 6.6797 | -0.499 | -0.253 | 5.011 |
| holdout8 | fullj_reciprocal | 0.900022 | 6.6779 | -0.499 | -0.252 | 5.001 |
| holdout8 | globaldual | 0.801796 | 6.9250 | -0.330 | -0.242 | 8.626 |

The residuals are against the original anonymous Y. Every listed nonzero residual is an approximate fit, never another solution of the original Y. Small residual and material recovery remain separate: DEV7 and successful holdout9 still have RMSE about 0.81–0.89.

## Six-node GT-line mechanism table

At t=0,.2,.4,.6,.8,1, m(t)=1+t(GT-1). No optimizer follows this prescribed line. The gradient alignment is undefined at the exact GT endpoint and is recorded as null. The maximum principal angle by itself can be dominated by one direction, so the mean sine is also retained in the data.

| Case | Path residuals, increasing t | Physical descent/GT cosine, t<1 | Gram condition range | Resolvent estimate range | Max GT-line finite-step relative remainder |
|---|---|---|---:|---:|---:|
| DEV7 | 0.999, 1.003, 0.914, 0.706, 0.373, 0.000 | -0.155, +0.164, +0.463, +0.596, +0.621 | 630–787 | 4.13–9.84 | 0.131 |
| holdout7 | 0.858, 0.881, 0.968, 0.892, 0.624, 0.000 | +0.112, -0.291, +0.069, +0.182, +0.660 | 630–899 | 4.13–13.44 | 0.254 |
| holdout9 | 1.452, 1.589, 1.547, 1.359, 0.932, 0.000 | -0.385, -0.174, +0.374, +0.543, +0.679 | 567–931 | 4.13–13.26 | 0.254 |
| holdout8 | 1.216, 1.313, 1.290, 1.220, 0.808, 0.000 | -0.315, +0.241, -0.292, +0.702, +0.369 | 556–1720 | 4.13–44.18 | 0.630 |

1. **DEV/holdout digit identity does not fix the basin.** The source and receiver geometry and m=1 Jacobian are shared. Yet the targets have distinct backgrounds, foreground contrasts and shapes. DEV7 foreground=12.300, background=1.372; holdout7 foreground=17.016, background=1.377, with a different font and stroke. Hence this is a compound shape/contrast change, not an isolated contrast experiment.
2. **Local descent is not necessarily GT-directed.** At frozen failed Full-J endpoints, physical descent/GT cosines are −0.346 (holdout7), −0.251 (holdout9), −0.499 (holdout8). Their RMS=.1 negative-gradient probes still decrease the physical loss, with actual/predicted decrease near 1.00–1.01. This supports locally coherent descent away from the truth direction; it does not certify stationarity or a mathematical local minimum. This probe uses the unregularized physical loss; it does not determine descent of the original regularized solver objective or explain its stopping criterion by itself.
3. **Gram collapse does not explain this contrast.** Source Gram conditions range from about 5.6e2 to 1.72e3 across the 24 path nodes, with no QR fallback. Full-J leading-32 squared singular-value energy stays at 97.4%–99.5% along paths. Singular tails reach numerical precision, but no exact nullspace or rank claim follows.
4. **The observable directions rotate substantially.** Top-32 spaces relative to m=1 reach near-right maximum angles for successful and failed cases, so rotation alone does not separate success. Full spectra and all 32 principal angles are retained rather than reduced to a rank assertion.
5. **Finite-step linear predictions can get the sign wrong.** With h=.05(GT−1), the predicted decrease is negative but the actual decrease positive at holdout7 t=0 and .4 and holdout8 t=.2. Their signed actual/predicted ratios are −1.298, −1.373 and −1.638. These are finite-step curvature effects, not a derivative accuracy test. Small RMS=.1 gradient probes remain locally predictable.
6. **Resolvent growth is case dependent, not sufficient.** Estimated ||A^-1||1 increases from 4.13 at m=1 to 9.84 (DEV7 GT), 13.1 (holdout7 GT), 13.3 (holdout9 GT), and 44.2 (holdout8 GT). The difficult 8 has the largest GT resolvent and remainder, but low-resolvent failed endpoints also exist. These four cases do not prove a universal threshold.

## Actual historical traces versus GT paths

`diagnostics/frozen_actual_traces.json` preserves the original scalar best-objective trace and stage events for every case/method. The traces have seconds, residuals and costs but no material vector at each point; therefore no unrecorded Jacobian or GT-distance trajectory is invented. Historical endpoint vectors are available and were diagnosed directly. The last panel of the figure is the actual frozen Full-J best trace; the other panels are the prescribed offline GT line.

## Definitions and numerical scope

The physical loss is .5||F(m)−Y||²/||Y||². Real analytic J comes directly from the frozen physics module and shared public G/E0/S. Spectrum uses a real 1152×576 dense SVD at each of 36 nodes; the first 32 material-space right singular vectors define the observational subspace summary. Gram eigenvalues refer to E^H E (12×12), not J^T J. The small-Gram objective and analytic gradient reuse the frozen implementation.

A=I−diag(m)G. LAPACK zgecon estimates its inverse 1-norm from the already retained LU and ||A||1; no explicit resolvent inverse or extra inverse columns are formed. Its internal RHS work is opaque and explicitly not counted as zero. The true remainder is formed directly as F(m+h)−F(m)−J(m)h. There are two probes: case-fixed h=.05(GT−1), and negative physical gradient scaled to RMS=.1; both respect [.5,20] by clipping. No finite differences, standard formula accuracy tests, hashes, source readback comparisons or optimization reruns were added.

## Paid cost

New run: 108 LU, 108 full forward calls, 36 analytic Full-J builds, 36 dense real-J SVDs, 36 source-Gram factorizations, 36 source-adjoint gradients, 36 zgecon estimates. Exposed complex RHS=3456: source=1296, receiver-transpose=1728, source-adjoint=432. Small-Gram RHS=1728 are separately counted small-system RHS, not full physical RHS. Wall=2.6203 s; process CPU=2.5146 s. Environment requests single-thread BLAS via five standard variables and threadpool_limits(1); threadpoolctl reported no discoverable pools on this runtime, so its list is empty rather than an independent thread-count measurement. Frozen historical solver costs are retained separately in the endpoint rows.

## Files

- `src/diagnostics.py`: reusable executable offline diagnostic.
- `diagnostics/mechanism_all.json`: all raw spectra, step probes, angles, costs and definitions.
- `diagnostics/mechanism_table.csv`: readable scalar mechanism table.
- `diagnostics/frozen_actual_traces.json`: original actual trace records.
- `diagnostics/COST_LEDGER.json`: paid operation ledger.
- `figures/basin_geometry_dev_vs_holdout.{png,pdf,svg}`: six-panel mechanism and historical trace figure.
