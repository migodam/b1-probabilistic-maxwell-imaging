# Self-only absorption control

The single frozen self-only run is complete and its prediction is sealed. No truth payload, material score, parameter tuning, additional derivative/identity suite, commit, or push was used.

The scientific change is `gin_mode='self'`: Gin contains only diag(G), with all intercell entries removed. The operator applies P_i=beta/[1+(beta-G_ii)m_i] by direct complex division. Its adjoint is also diagonal division. These are charged as actual scalar divisions, with zero small4×4 factorizations or solves. Default `gin_mode='block4'` retains the previous block branch. The original block execution source, frozen configuration, receipts and prediction remain untouched.

The generic exact material value/gradient and mature scalar CIE seed plus L-BFGS readout functions are AST-identical to the original executed block implementation. The proposed exact rational self readout was not implemented. All matched scientific fields listed in `COMPARISON.json` are unchanged: source metrics, prior, bounds, initialization, Fourier/DPCS, three 500-update stages, max 60 material iterations, and common scaled original-Phi corrector(maxiter1000). Candidate selection uses original Phi including prior and retains the original host incumbent. The one prior block-gradient wiring receipt is referenced; no new numerical check was run.

| Arm | Original Phi including prior | Relative physical residual |
|---|---:|---:|
| self_only | 0.00583250389712 | 0.107820651645 |
| block4 | 6.50627883468e-06 | 0.000967014434494 |
| R2_global_common_corrector | 0.00583248209351 | 0.107820545567 |

Self auxiliary stage Phi values were 0.26587463551255636,0.06756084411575954,0.010801894717538736. All three stages reached 500 updates. The common corrector reached 1000 iterations with 1048 evaluations. Among 1500 material solves, 676 hit the 60-iteration cap, 711 stopped on projected gradient and 113 on relative function reduction; zero fallbacks occurred. These caps are preserved, not convergence claims.

This run did not reproduce the block4 physical-fit improvement using only the self term. It does not establish that intercell blocks are universally necessary or prove a basin or material-recovery result. Parent scientific scoring and conclusions are separate.

Measured process wall/CPU: 37.285989/34.663160s. Physical Maxwell cost: 1055 LU and 25236 complex RHS. Self diagonal work: 144602 application calls and 1166529024 actual complex scalar divisions (682336512 forward, 484192512 adjoint). Material evaluations: 70051; inner L-BFGS iterations: 61172. Source-metric Cholesky and receiver/APCS SVD counts remain separately reported. Startup imported successfully and failed_attempts=[]; operator assembly was not rerun.

Files include `FROZEN_CONFIG.json`, `STARTUP.json`, `AUXILIARY_EXECUTION.json`, `AUXILIARY.json`, `CORRECTOR.json`, `RECEIPT.json`, three stage-current NPZs, `SELECTED_BEFORE_CORRECTOR.npz`, `ALL_PHYSICAL_CANDIDATES.npz`, `prediction.npz`, `RUN.log` and `executed_source_snapshot/`. The original block note predates this extension; its executed snapshot is authoritative for the old run, while this run has its own complete source snapshot.

Command: the scientific `.venv_nn/bin/python` with CPU thread count 1 runs `experiments/block_absorption.py --variant self_absorption`. An existing result receipt prevents overwrite. `experiments/summarize_self_absorption.py` reads saved receipts and source only, and its postprocessing timing is separate.
