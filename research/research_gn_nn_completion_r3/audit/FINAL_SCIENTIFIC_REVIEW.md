# Final scientific review — GN–NN R3

Review by the primary Codex investigator after prediction closure and TEST scoring. The earlier independent source reviews are preserved in PRETRAIN_REVIEW.md and EVALUATION_REVIEW.md. **This final interpretation is not a new independent external replication.** No extra training or physics sweep was commissioned after seeing TEST.

## Disposition

| Claim | Decision | Evidence and boundary |
|---|---|---|
| End-to-end full-VIE NN completion was implemented and run | PASS | 1024TRAIN mothers,64VAL,64TEST;192TEST observations;17models;3264model–observation evaluations,384controls and128additional filter evaluations |
| Learned material priors can improve these GN images | SUPPORTED WITH SCOPE | Clean image/voxel gains; noisy J-full gains in all3seeds; candidates independently re-evaluated through inherited full VIE |
| Full receiver-null projection is generally superior to full current | NOT SUPPORTED | J-null is modestly better in clean NRMSE, substantially worse in noisy NRMSE across all3seeds |
| V-physics is the consistently best architecture | NOT SUPPORTED | seed7 good;17/27 poor already onVAL; no best-seed replacement |
| Two-pass voxel Cascade is promising | EXPLORATORY | seed7 improves noisy joint recovery12→23/64;17/27Cascade NOT_RUN |
| Current-target accuracy establishes material recovery | REFUTED AS A SUFFICIENT INTERPRETATION HERE | J-null-current better current-target accuracy but worse material readout/top1 than material-loss arms |
| Large nullspace energy establishes a useful material search representation | NOT SUFFICIENT | Clean98.09%/noisy91.85%current-error null energy coexist with divergent readout and physical qualification |
| Hard-null material space has only one missing direction in general | NOT PROVEN | Conditional annihilator derivation; numerical575rank at two public anchors only |
| Neural output discovers probability mass or new physical branches | NOT CLAIMED | Finite deterministic candidates; no posterior, root certification or connected-component analysis |
| Universal gain over strong regularized inverse solvers | HOLD | One weak-reg baseline; discrepancy stopping and regularization path absent;A2.5 favors weak-reg |
| IEEE TAP-level independent algorithm novelty | HOLD | CNN/current learning/nullspace projection and correction have direct prior art; controlled comparison and mechanism may justify further work |

## What was actually frozen

Prediction receipt closure includes all17planned models and two controls before first numerical TEST truth scoring. The secondary32-observation subset and its four filter runs also closed first. Numerical truth scoring occurred at approximately2026-11-07 01:59UTC. Release packaging copied raw truth bytes earlier without loading arrays or passing truth to inference. This is process isolation, not an operating-system security guarantee.

Training normalization uses TRAIN only. Validation chooses checkpoints and authorized the256→1024 expansion. TEST images, font labels, masks, amplitude and family do not enter inference, correction or online selection. Source inspection nevertheless found a generator-level digit/amplitude correlation before TEST scoring; this limits distributional claims and is prominently reported.

## Independent quantities remain separate

- Physical feasible: full nonlinear data residual meets the frozen gate.
- Joint success: physical feasible plus materialNRMSE and OtsuIoU gates.
- Online top1: fixed data-feasible nearest-NN rule, without truth.
- best@K: GT-scored finite candidate coverage, never called an online choice.
- Natural current oracle: L2 projection of true current onto affine receiver constraint, not the optimal attainable material readout.
- Auxiliary constitutive defect: a physical consistency diagnostic, not a substitute for full-material forward validation.

All192primaryTEST anchors meet the physical gate. Consequently a100%feasible selection rate alone provides little evidence; it must be read alongside fallback, displacement and material gains. The study does so. Most learned proposals are corrected into the gate and retain substantial material displacement; this is measured finite-step behavior, not a transverse-contraction theorem.

## Fairness and remaining limitations

The voxel and current heads differ in output DOF and fixed dense algebra; only J-full/J-null directly isolate projection. The same backbone/epochs does not imply equal FLOPs. Data costs, training costs, reference setup, warm inference and correction are separated. There is no cold end-to-end GPU speed claim. Failure counters and noninterruptible wall overruns remain; strict-budget success counts are reported in addition to completed endpoints.

The weak-reg45 candidate archive includes common GN for offline comparison; its online selection does not use that candidate. Its offline best@K must not be interpreted as a free cold-start multistart result. Source and report explicitly disclose this distinction.

Three TRAIN numerical SVD stops remain usable finite accepted endpoints, not converged solves. One initial crash has unknownLU/RHS; exact all-attempt operation totals cannot be asserted. No precision/fidelity claim is manufactured from repeated error-floor tests.

The frozen Otsu vsalpha≥.5 definition is imperfect for continuous materials: GT itself passes the shape gate47/64. Posthoc calibration is reported without changing gates. This is not a mathematical47/64upperbound on reconstructions. Original joint rates, NRMSE, data residual and family results all remain visible.

The confirmed scope is low-contrast synthetic same-grid scalar-TM with fixed acquisition and A2.5 extrapolation. High contrastA12, physical mismatch, experimental data, general unseen geometries, fully noise-optimized conventional baselines and multi-seedCascade remain NOT_RUN.

## Final scientific judgment

The completed experiment supports learned material completion under full physical verification. It does **not** support making the exact receiver-null current space the sole NN search domain. The strongest next scientific question is how material readout and noise-consistent freedom in receiver-visible currents should be weighted, while keeping a stable voxel baseline and better conventional noise regularization. Further novelty depends on answering that question with new frozen data, not on adding more names to the existing solver.
