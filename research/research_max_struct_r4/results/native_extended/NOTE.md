# Extended ordinary native CIE control

Status: `SEALED_DEVELOPMENT_PREDICTION`; truth access `NONE`. This is one authorized capacity/budget control on the disclosed development observation, not a parameter sweep or a new method. The frozen mature sourcewise CIE implementation is unchanged. The runner supplies the identical historical prior value/gradient and omits constructing its unused dense IRLS matrix.

The original three 500-update stages (mf3, mf6, mf10; beta-current 1, beta-material 1/.2/.2) were followed by 6500 mf12 updates with beta-current 1 and beta-material .2. Coefficient tolerance was zero. All 8000 updates completed. mf12 has 576 native Fourier coefficient columns and rank559 APCS, with the same rank17 locked receiver DPCS. Exact comparison of the first three material/current snapshots against frozen R2 native snapshots found zero difference in every entry; see `HISTORICAL_PREFIX_COMPARISON.json`.

| Candidate | Updates | Original Phi | Physical relative residual |
|---|---:|---:|---:|
| Host .1 | 0 | .523523272161 | 1.02325292295 |
| mf3 stage | 500 | .267734534836 | .731751747771 |
| mf6 stage | 1000 | .175174650928 | .591890476712 |
| mf10 stage | 1500 | .0864275665506 | .415735546726 |
| mf12 stage | 8000 | .0177079178671 | .188076732096 |
| Selected mf12 + common original corrector | 8000 + 1000 | .00479822370453 | .0977065589361 |

Only original Phi with the unchanged prior selected among the host and the four actual stage endpoints. The selected mf12 material received exactly one call to the same `cooperative_paths.refine` used for A: maxiter1000, objective scale1e6, ftol1e-12, scaled gtol.01. It used 1011 objective evaluations and stopped at the iteration limit. Final projected-gradient infinity norm is 5.4923885e-7 in the unscaled original objective. Whitened physical residual2344.96335 exceeds the original noise radius25.1654222. This runner makes no material-recovery assertion and has no GT access.

Full execution cost is 73.5941804 s wall and 68.743973 s CPU, including startup/imports, auxiliary updates, common correction, original physical evaluations and saved artifacts. Auxiliary execution was 59.8068902 s wall / 55.847282 s CPU. The corrector was 12.9347008 s wall. Total physical cost is 1015 original LU and 24324 complex RHS, including five auxiliary candidate validations (5 LU/60 RHS), common correction (1010 LU/24252 RHS) and one final diagnostic adjoint (12 RHS). All solves use frozen public operators, CPU1 and complex128; no assembly or GPU setup is charged because none was run.

Native counters record 8000 material readouts, 40000 material-inner evaluations, 8000 material-inner iterations, 8000 PR-CG steps, 576000 reduced-operator columns, four APCS SVDs, one receiver SVD and 1156 explicitly generated Fourier columns. The transparent mode-work proxy sum(ns × native Fourier columns × updates) is 48,408,000, compared with 3,480,000 for the original 1500-stage schedule: 13.9103 times the mode-work proxy and 5.3333 times the outer updates. The per-update nine dense current matrix products correspond to an estimated 148,709,376,000 complex multiply-adds; this estimate excludes setup, elementwise work, SVDs, material readouts and the corrector, and is not a FLOP or elapsed-time equivalence claim.

For context, A's separate receipt reports 59.4237267 s total wall, 56.562409 s CPU, 1033 physical LU/24708 physical RHS, and 61440 material-inner evaluations, in addition to its explicitly charged 18,342,720 small 4x4 factorizations / 261,872,640 small complex RHS. This extended native control uses more ordinary modes and outer iterations and takes longer wall time in these local runs; it does not exactly match A's work, material-inner effort, or operation mix. Neither run certifies convergence, and no general speed comparison follows.

`AUXILIARY.json` preserves every auxiliary iteration and stage score. `checkpoints/iteration_00500.npz` through `iteration_08000.npz` retain every 500th material/current/field without extra physics. `CHECKPOINTS.json` records contemporaneous wall/CPU and native costs; its native counter snapshot occurs immediately before the mature routine increments the current iteration's PR/footer counters. These interim checkpoints were not ranked as additional physical candidates. `RECEIPT.json`, `CORRECTOR.json`, `ALL_PHYSICAL_CANDIDATES.npz` and `prediction.npz` retain the final disposition. No failed attempt occurred. The numerical run cost excludes the later read-only snapshot comparison and this note.
