# R4 execution checkpoint — 2026-11-07 around06:12UTC

The six-hour job remains active. Deadline10:44:27UTC. New TEST truth arrays/metrics have NOT been read. All GT access so far was historical R3 or TRAIN/VAL generation; no new network has trained yet at this checkpoint.

- E0 complete:17 historical runs ×192observations =3264. Actual currents6532LU/79920RHS. New figures separate sequential Twofold, intersections, auxiliary heads and real physics.
- Data generated:504/78/108 TRAIN/VAL/TEST mothers; clean/noisy1380observations. Strict primary GN preparation ongoing. VAL156ready. Primary TRAIN still incomplete, so training must wait for actual `audit/DATA_READY.json` release.
- Actual processes after a local recovery: finite training waiter69065 (detached own process), M4 preparation controller69950/runner70433, fetch watcher70081, VAL public-evaluation controller65905. The i9 mainPython1292 continued throughout, so no remote restart was made. These PIDs are snapshots: inspect before acting.
- Old tool session handles were lost and old local preparation/fetch/waiter processes had exited. Root cause unknown. Recovered scripts reuse complete pairs, append logs and preserve old ledgers. Earlier17 misconfigured completed attempts and up-to4 interrupted attempts are recorded; recovery adds up-to2 interrupted local attempts with UNKNOWN_NOT_ZERO cost.
- Resource limits: only own i9Pcores4/6/12/14 withBLAS1; M4max2numericalCPUworkers; MPS for sequential NN training. RTX4060 not allocated because another session owns it. Never interrupt unrelated R11/other processes.

## Automatic work already scheduled

1. `delegated/data/watch_preparation.py` fetches complete pairs, publishes full TRAIN1008/VAL156 release, later starts remote independent VAL/TEST weakreg and discrepancy anchors after all strict anchors.
2. `delegated/networks/run_initial_five.py` waits for the full release, then trains from scratch the five nonsoft seed7 arms and immediately predicts public VAL per model. It must not be launched again while the real PID is live.
3. `delegated/data/watch_val_evaluation.py` waits until both M4 preparation workers have ended, then processes explicit runs in `audit/VAL_EVALUATION_REQUEST.json` using two BLAS1 partitions. It reads no truth. Parent owns adding runs/closing request.

## Parent's next operations

- After Jfull VAL predictions: prepare posthoc alphas0,.1,.3,.6 via `experiments/alpha_scan.py --source-run j_full_seed7_cascade2 --alphas 0 .1 .3 .6`; alpha1 reuses existing Jfull evaluation instead of repeating it. Add the four posthoc run names to the VAL request.
- Seal complete VAL evaluations using `experiments/score_predictions.py seal --split val --runs ...`, then scoreVAL. Apply `experiments/select_validation.py` under the predeclared `audit/VAL_SELECTION_CONTRACT.json`.
- Interior alpha: from-scratch soft model; alpha0/1: report endpoint degeneracy, do not fabricate sixth model. Choose constrained arm C* using VAL only.
- Train seeds17/27 forVoxel/Jfull/C*. Train const/internal .01 arms onJfull/C*, all fromscratch. Predict/evaluateVAL. Optional fullVIE matched-pair fine-tuning only if declared gate and time pass.
- New user steering: A17-inspired current→real shared material→fixedG physics residual→current RNN, avoiding whole L inverse. Theory derived in `theory/A17_RECURRENT_FEEDBACK.md`: exact finite readout gap; fixed-m descent; profile majorization; physical-normal Hessian and transverse local geometry. Do not confuse fixed-step derivative theorem with implemented direction-dependent exact-line steps.
- Conditional supplementary training under `audit/RECURRENT_SUPPLEMENT_CONTRACT.json`: only refreshed_q and feedback1, rho1, seed7; reuse existing Jfull as nofeedback control. `experiments/train_recurrent.py --authorized --rho 1 --modes refreshed_q feedback1 --predict-val`. No training has occurred yet. Feedback2 and small-A17-Schur training NOT_RUN. Preserve main six-space priority.
- Prepared public VAL16-iteration depth diagnostic: `experiments/recurrent_iterations.py`; frozen6VALobservations in `audit/RECURRENT_ITERATION_DIAGNOSTIC.json`. Out-of-trained-depth empirical stability only, no fixed-point proof.

## TEST and publication still pending

Freeze all method choices before new TEST GT. Public inference/correction must finish and be sealed before scoring. Keep planned216observations/108mothers, NOT_RUN denominators, top1 vs bestK, physical feasibility vs quality vs55s cost separate. Independent cold CIE reference is preselected13mothers×2noise unless explicitly expanded, never call that full108 coverage.

`experiments/test_energy.py` supports new true-current stages plus aux→physical realization gap and threshold sensitivity without new SVD. `experiments/plot_generalization.py` uses sealed tables only; `src/plots.py` makes predeclared13mother figure and full216observation gallery. Every image/class/font/rule exposure must be labelled.

Final source/data/weights/results publication remains NOT_DONE. Existing B1 repo/site only, no new repo. Local initial commit51aa92fed exists; origin contains unrelated R11 updates. Preserve frozen sources; only scoped new-study integration. Big arrays/models go to GitHubRelease, compact scientific report/figures/CSV to existing site. Verify real direct report links after publishing; do not claim the future R4 URL is live before verification.
