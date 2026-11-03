# DCTMC-style control: fixed 30-step disposition

Status: **SEALED_DEVELOPMENT_PREDICTION**. No GT arrays or shape labels were read. One public observation, constant material .1 initialization, exactly 30 iterations, one selected-candidate corrector, no parameter tuning. Full parameters were saved in `FROZEN_CONFIG.json` before numerical work. This is a known-method control, not a new algorithm or a reproduction of the paper's numerical experiment.

The source implements the T/V maps and singular-coordinate overwrite described in [Levinson–Markel, DCTMC I](https://arxiv.org/html/1401.3319v3), equations 10, 11, 23, 24 and 26. Pure diagonal extraction is a permitted special case of its locality operator. Intermediate T/V need not be symmetric. Our declared variations are public-host initialization, real diagonal projection clipped to [.02,20], relative product-singular cutoff .001, physical-score selection with the historical prior, and one original-physics corrector. No symmetry is forced; no weighted off-diagonal aggregation or convergence claim is used.

The completed nonlocal T and legal physical material have very different residuals:

| Stage | Original physical relative data residual | Original Phi | Completed nonlocal T data residual | V offdiag/full |
|---|---:|---:|---:|---:|
| Public host | 1.023252923 | .5235232722 | .0008317563 | not evaluated |
| Iteration 1 | 1.023035364 | .5233006788 | .0008319086 | .944984479 |
| Iteration 4, selected by physical Phi | 1.022327184 | .5225764422 | .0008324988 | .930988543 |
| Iteration 30 | 1.039494744 | .5402747589 | .0008397992 | .732086182 |
| Selected iteration 4 plus common corrector | **.09586391932** | **.004614443850** | not evaluated | not evaluated |

The cutoff retained 185 of 576 measured complex singular coordinates, leaving 331591 of 331776 T coordinates unconstrained by the overwrite. Discarded observed data has relative norm .0008343813. The retained-coordinate residual is approximately 3.1e-15 after overwrite, while the omitted coordinates remain free. The .001 cutoff is a parent-declared DEV noise-level truncation, not a certified rank/noise-optimal threshold or full-data compatibility certificate.

At iteration30, nonlocal V has offdiag/diag=1.074681824, imaginary-diagonal relative norm=.0121460714, legal-local-projection relative change=.7323854301, and transpose-asymmetry ratio=.405875912. Projection clipped 264 real diagonal entries to the lower material bound. Decreasing offdiagonal fraction did not make the material data residual decrease. These are model-locality/data diagnostics only; material recovery was not scored.

The one corrector used the existing `experiments/cooperative_paths.py::refine`, original objective scaled1e6, maxiter1000, original-equivalent gtol1e-8. It performed 1000 iterations/1051 evaluations and stopped at the iteration limit. Final original projected-gradient infinity norm is 6.35426116e-6, so this is not a first-order convergence certificate. No extension was run.

| Cost scope | LU | Complex RHS | Details |
|---|---:|---:|---|
| Host plus30 lifting cycles | 61 | 35508 | 31 original physical LU;30 nonlocal locality LU;31×576 dense T RHS;30×576 locality RHS;31×12 source RHS |
| One common corrector | 1051 | 25224 | 12612 source and12612 adjoint RHS |
| Final exact gradient report | 0 | 12 | retained final physical LU |
| **Total** | **1112** | **60744** | **1082 original physical LU plus30 nonlocal locality LU** |

Measured run wall=15.04138483s, CPU=13.84945s. Lifting/setup/31-candidate persistence through pre-corrector checkpoint took3.09613633s. Detailed SVD, matrix products, dense solves, source/adjoint counts and operation timers are in `RECEIPT.json`. The shell's first Python lacked numpy and failed before scientific work; `STARTUP_FAILURE.json` retains zero physical LU/RHS and unknown startup wall/CPU. The successful run used the existing `Gaussian/.venv_nn/bin/python` and CPU1; no dependencies were installed.

All31 physical materials and predictions are saved separately under `candidates/`; `ALL_PHYSICAL_CANDIDATES.npz` also contains their material arrays and original scores. `FINAL_LIFTING_STATE.npz` contains the last completed/physical T and the nonlocal V immediately before the last projection. `SELECTED_BEFORE_CORRECTOR.npz`, `CORRECTOR.json`, `prediction.npz`, and `ITERATIONS.json` preserve the selected seed, complete correction history, sealed endpoint, and per-step defects. No publication, commit, GT scoring, or further iteration was performed by this worker.
