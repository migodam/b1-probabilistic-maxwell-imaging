# Fixed block-loop absorption: implementation and recorded result

One configured 24×24 disclosed curvedbar development case was run, without reading truth. The new prediction and the global common-corrector control are sealed. No parameter tuning, CIE rerun for the control, formula identity suite, commit, or push was performed. Scientific conclusions and material scoring remain with the parent.

## Implementation and exact comparison boundary

`src/cie_block_absorption.py` retains the R2 global PR loop, fixed source metrics, native Fourier coefficient coordinates, DPCS offset, stage schedule, bounds, and original material prior. It imports the mature Fourier/DPCS helper and shared-real scalar CIE readout. Disjoint 2×2 spatial-cell groups give 144 actual 4×4 Green submatrices including self terms. These groups and their alignment are fixed.

The parent-supplied replacement is B=beta I-Gin, K=I+diag(m)B, R=beta K^-1(J-diag(m)Q). Current H/T use this same fixed-material transform. Material evaluation uses a=K^-H(R Omega), z=B R+beta Q, gradient=-rho Re sum(conj(a)*z)+prior gradient. The scalar CIE readout supplies a seed only; both previous and seed are scored with the exact block loss plus prior before L-BFGS-B(maxiter60), with fallback retained. This adopts known block absorption and CIE algebra and makes no originality assertion.

All three phases ran 500 PR updates: mf3/beta1, mf6/current1/material0.2, mf10/current1/material0.2. DPCS rtol=1e-3, host=0.1, box=[0.02,20]. The prior uses the exact original `model.reg.evaluate` values/gradients; the unused dense Hessian built by the historical `model.regularizer` wrapper is omitted. Thus wall comparisons also include this implementation difference, and no speed claim is made.

The control loads the minimum-original-Phi R2 global auxiliary candidate (stage 2, Phi 0.011269442424001625), applies the identical R4 `cooperative_paths.refine(maximum=1000)`, and retains the published R2 global prediction as an incumbent. Both new paths use original Phi including prior for selection. Historical native and R3 scaled-native references have different cost histories and are context, not equal-compute arms.

## Public physical scores

| Arm | Original Phi including prior | Relative physical residual | Physical LU | Physical complex RHS |
|---|---:|---:|---:|---:|
| block_CIE_plus_common_corrector | 6.50627883468e-06 | 0.000967014434494 | 1033 | 24708 |
| R2_global_plus_common_corrector | 0.00583248209351 | 0.107820545567 | 1013 | 24288 |
| published_R2_global | 0.00583348384709 | 0.107830606924 | 689 | 16452 |
| published_R2_native | 0.0050962193293 | 0.100734862445 | 1020 | 24396 |
| R3_scaled_native_restart | 0.00509610269221 | 0.100733686717 | 748 | 17964 |

Block auxiliary Phi was 0.2645905236382217,0.030312610837575287,0.0007066885666717203 across stages; matched R2 global auxiliary Phi was 0.2767745966664839,0.15300752495451386,0.011269442424001625. Both new physical correctors reached the 1000-iteration cap. All three CIE outer stages reached 500. Of 1500 material inner solves, 601 reached 60 iterations, 734 stopped on projected gradient, and 165 on relative function reduction. No endpoint is asserted to be a certified local minimum or basin result.

## Validation, costs and provenance

The sole new gradient wiring check passed at beta0.2 with one fixed public-operator direction: relative error=2.60446043399e-08 at h=1e-6. Its three exact value/gradient calls cost 864 small4×4 factorizations/solves and 10368 complex RHS, with zero physical LU. There was no extra identity or error-floor suite.

Full block auxiliary work used 61440 material value/gradient evaluations and 53467 L-BFGS inner iterations. Including wiring: 127386 batched calls, 18343584 small4×4 factorizations, 18343584 small4×4 solves, and 261883008 small-system complex RHS. NumPy's dense library solve refactors per call; this cost is explicitly separate from Maxwell LU/RHS. There were no inner fallbacks. Existing iteration receipts show decreasing fixed-current block material loss and decreasing fixed-material current objective at every recorded step; these are recorded run facts, not additional experiment tests.

Block process wall/CPU=59.423727/56.562409s; common-control new wall/CPU=10.607888/10.154329s. These include startup/imports and the artifacts produced inside each runner, with the summary postprocessing timed separately. Historical operator assembly was not rerun. New runs total 2046 physical LU/48996 complex RHS. Both startup records imported successfully; failed_attempts=[] in both receipts. No costs for hypothetical unexecuted work are reported as measured.

`STARTUP.json`, `GRADIENT_WIRING.json`, `AUXILIARY_EXECUTION.json`, `AUXILIARY.json`, `CORRECTOR.json`, `RECEIPT.json`, `stage_0.npz` through `stage_2.npz`, `ALL_PHYSICAL_CANDIDATES.npz`, and `prediction.npz` preserve the evidence. `results/global_common_corrector/` holds its independent receipt and prediction. Executed source snapshots are saved; the current module differs from its executed copy only by explanatory docstring corrections. Frozen R2 CIE, mature baseline, and physical model source files were not modified.

Reproduce in a fresh output directory using the supplied `.venv_nn/bin/python`, CPU thread count 1, and `experiments/block_absorption.py`, then `experiments/global_common_corrector.py`; both runners refuse to overwrite an existing disposition. `experiments/summarize_block_absorption.py` reads receipts only. No GT-bearing score files or truth payloads are read by any of these runners.
