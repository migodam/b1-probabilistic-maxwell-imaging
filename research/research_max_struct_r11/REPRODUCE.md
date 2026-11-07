# R11 reproduction and accounting entry points

This is an initial, file-backed handoff. It distinguishes reading completed evidence from executing an inverse or opening confirmation GT. No command below was run to regenerate science while compiling this document. Existing output directories are immutable: a scientific repeat needs a separate manifest/output directory and the same frozen numerical contract. Do not restart the active Windows48 or vector queues.

## Existing environment and repository root

The actually used local interpreter is:

```
/Volumes/migodam's-external-brain/Research/Inv_SLAM/Gaussian/.venv_nn/bin/python
```

Keep that path without resolving its symlink. Example shell setup, from the actual R11 worktree repository root:

```sh
TASK_PY="/Volumes/migodam's-external-brain/Research/Inv_SLAM/Gaussian/.venv_nn/bin/python"
cd "/Volumes/migodam's-external-brain/Research/Inv_SLAM/research/delegated/max_struct_r1/r11_worktree"
```

Portable use substitutes the repository root and an existing compatible interpreter; preserve repository-relative public/config/output paths. CPU runs use numpy/scipy/threadpoolctl and BLAS/OpenMP threads=1. Vector GPU execution additionally uses the existing compatible Torch/CUDA environment selected in its execution config; this handoff does not install a new backend, change precision or authorize new GPU allocation. Historical transport scripts contain site-specific connection configuration: reuse the configured connection rather than copying credentials into reports. No SSH private-key contents are needed or provided here.

## Refresh costs without science or GT

These real scripts read local metadata/receipts. They do not retrieve missing remote arrays or interrupt a running process:

```sh
"$TASK_PY" research/research_max_struct_r11/experiments/twod_compile_costs.py
"$TASK_PY" research/research_max_struct_r11/experiments/twod_independent_costs.py
"$TASK_PY" research/research_max_struct_r11/experiments/summarize_completed_costs.py
```

The last script is stdlib-only and accepts `--study-root /absolute/path/to/research_max_struct_r11` and `--output-dir /absolute/report/directory`. It writes COST_LEDGER.json/MD. Rerun it after producers synchronize terminal receipts. It never reads NPZ/GT, executes an optimizer, opens SSH, or stops/restarts a queue. It prefers the authoritative upstream ledgers and replaces same-scope nested vector CPU with whole-child PHASE_COSTS when available.

Known additive CPU is a lower bound. Exact child CPU, own controller CPU, nested phase/readout timing, wall intervals, GPU stream intervals, memory peaks and physical counters are distinct quantities. 576-cell scalar-TM LU, 5184-dimensional vector-current LU and 17496-dimensional fine-grid iterative solves are separate. Small Woodbury LU is not a physical Maxwell LU. Counter aliases and portfolio overlap must not be summed twice.

## Completed scalar-TM DEV controls

The common runner is `experiments/twod_run.py`, with actual solver config `configs/twod_dev_solver.json`: uniform m=.1, [.02,20] real box, inherited prior/source metrics/noise99 DPCS; mf3/6/10/12 auxiliary 500/500/500/1500; corrector at most3000; whole-child CPU480 and paid prefixes60/120/240/480. Each task's PUBLIC_TASK/FROZEN_CONFIG/NUMERICAL_SCHEDULE and AUXILIARY_EXECUTION are the exact run receipts.

Existing task manifests:

- `configs/twod_dev_lagged_tasks.json`: five exposed freeze-A tasks; remote DEV ledger also retains one full replay.
- `configs/twod_dev_contracted_tasks.json`: three exposed contracted-T tasks, ordinary scalar contraction retained.
- `configs/twod_dev_bothblocks_tasks.json`: three fullGreenBothBlocks tasks; use its isolated wrapper.
- `configs/twod_dev_current_only_tasks.json`: three Green-current/ordinary-material factor controls; arm=global_dual retains the original prior and has_reg material branch.

Actual single-task command shapes, for use only with a new repeat destination:

```sh
"$TASK_PY" research/research_max_struct_r11/experiments/twod_run.py --manifest research/research_max_struct_r11/configs/twod_dev_contracted_tasks.json --config research/research_max_struct_r11/configs/twod_dev_solver.json --task-id r10_twod_0001__contractedLagged__480CPU
"$TASK_PY" research/research_max_struct_r11/experiments/diagnose_bothblocks_run.py --manifest research/research_max_struct_r11/configs/twod_dev_bothblocks_tasks.json --config research/research_max_struct_r11/configs/twod_dev_solver.json --task-id r10_twod_0022__fullGreenBothBlocks__480CPU
"$TASK_PY" research/research_max_struct_r11/experiments/diagnose_current_only_run.py --manifest research/research_max_struct_r11/configs/twod_dev_current_only_tasks.json --config research/research_max_struct_r11/configs/twod_dev_solver.json --task-id r10_twod_0022__GreenCurrentOrdinaryMaterial__480CPU
```

The independent BothBlocks wrapper patches dispatch only inside its child. Ordinary CIE/Dual/Global and material-only G19 must use the ordinary runner. Never patch a mixed-arm controller globally. The executed source snapshots, standalone H wiring check, full diff and portable interface are in `delegated/diagnostics/bothblocks_run_logs/`, `BOTHBLOCKS_WIRING_CHECK.json`, `BOTHBLOCKS_SOURCE_DIFF.patch`, and `BOTHBLOCKS_PORTABLE_INTERFACE.md`. Current-only source/diff/costs are separately preserved in `delegated/diagnostics/current_only_dev/`. Current beta=1 and material beta=1/.2/.2/.2 differ; these records do not claim a common alternating objective.

## Completed independent16 x five-arm confirmation

Frozen files: `configs/twod_bothblocks16_{generator,solver,scoring}.json`, the generated `twod_bothblocks16_all_tasks.json`, and `delegated/diagnostics/bothblocks_confirmation/FROZEN_CONTRACT.json`. Seed202611071116, eight material families x high/mixed, independently spawned geometry/noise streams, balanced high12/18, inherited 24x24/12source/48receiver acquisition. Public observations are under `data_public/bothblocks16/`; GT under `sealed_truth/bothblocks16/` was generator/scorer-only.

The actual finite generation/80-run entry is `experiments/diagnose_bothblocks16_queue.py --execute`. It performs generation with `diagnose_bothblocks16_generate.py`, preserves startup source/config snapshots, runs two local lanes in case-index cyclic arm order, uses ordinary dispatch for four arms and the isolated BothBlocks wrapper only for that arm, and captures exact external child wait4 CPU. Existing paths intentionally cause a repeat to refuse overwriting.

Actual postseal sequence, already completed after explicit root authorization:

```sh
"$TASK_PY" research/research_max_struct_r11/experiments/diagnose_bothblocks16_seal.py
"$TASK_PY" research/research_max_struct_r11/experiments/diagnose_bothblocks16_physics_prepare.py
"$TASK_PY" research/research_max_struct_r11/experiments/diagnose_bothblocks16_score.py
"$TASK_PY" research/research_max_struct_r11/experiments/diagnose_bothblocks16_frontier.py
"$TASK_PY" research/research_max_struct_r11/experiments/diagnose_bothblocks16_plot.py
```

The first two are GT-free selection/sealing. All80 terminal outcomes, 80 predictions, 64 predeclared paid-portfolio selections and 80 common-LU selections were sealed before scoring. The last three require authorized GT access and perform no new physical forward. The portfolio rule and actual costs are in `PORTFOLIO_RULE_BEFORE_DATA_SCORE.json`;240+240 is an allowed-budget view, not equal actual spend. Object-level high/mixed-stratified bootstrap is10000 replicates, seed202611071117. Source/seed/rules were not tuned after scores.

## Independent48 x six-arm campaign: still active

Use `delegated/twod/independent48/RUNNING_HANDOFF.md` and the existing transport/monitor state. Do not rerun launchers from this document. Original scientific files are `configs/twod_independent48_{generator,solver,scoring}.json`, `twod_independent48_all_tasks.json`, and public `data_public/independent48/manifest.json`. Original baseline-first dispatch completed five tasks before its explicit stop; the remaining283 use `twod_independent48_interleaved_remaining_tasks.json`. PLAN/FROZEN_MANIFEST/FROZEN_SOLVER and DISPATCH_REORDER_RECEIPT retain that chronology rather than rewriting history as initially interleaved.

`experiments/twod_independent_score.py` is the real primary scorer; it requires manifest/queue arguments, the frozen scoring contract, explicit postseal authorization/note, and the public secondary selection seal. Do not invoke it while any of288 attempts remains pending/running. Secondary `twod_secondary_score.py` only reindexes already scored primary results after its public seal. Its policy was registered after independent16 scoring and before any48GT read, not as the original primary pre-registration. The cost compiler reads only receipts, never48GT.

## Vector DEV and ten new observations: active finite phases

Real entry points and frozen contracts:

- `experiments/vector_source_basis.py`, `vector_basis_audit.py`, `vector_derivative_check.py`, `vector_z_derivative_check.py`: saved/public-background preparation and diagnostics; receipts state which were zero physical solve.
- `experiments/vector_generate_confirmation.py`: generation-only, hard-wired to `configs/vector_confirmation_contract.json` and `vector_confirmation_material_design.json`; existing private/public outputs refuse overwrite. Eight coarse5184 solves and two fine17496 iterative observation generations are distinct.
- `experiments/vector_run.py`: portable flags `--config`, `--queue`, `--output`, `--max-tasks`, optional `--offset`; these paths are interpreted relative to the R11 study directory.
- `experiments/vector_remaining_controller.py`: actual finite sequence of dev_denominator6 then confirmation50. It captures each scientific child's Windows GetProcessTimes into PHASE_COSTS; its own PROCESS_COST is only controller CPU. It is already active; do not start another copy.

The controller's actual phase command shapes are:

```sh
python research/research_max_struct_r11/experiments/vector_run.py --config configs/vector_dev_contract.json --queue configs/vector_dev_denominator_queue.json --max-tasks 6 --output delegated/vector/dev_denominator_process
python research/research_max_struct_r11/experiments/vector_run.py --config configs/vector_confirmation_execution.json --queue configs/vector_confirmation_queue.json --max-tasks 50 --output delegated/vector/confirmation_process
```

They are historical reproducibility records, not instructions to consume an active GPU allocation. Use a new approved destination and coordinated resource only for a new repeat. Ordinary six DEV arm receipts are nested inside `dev_ordinary_process/PROCESS_COST.json`; do not add their arm CPU sum. Gate-helper CPU is a separate process and counted once. The engineering six-LU library comparison is already inside the dev_denominator child. `CURRENT_RUN.json`, `IMPLEMENTATION_NOTES.md` and retained source overlays describe active provenance and measurement gaps. No vector confirmation GT has been opened for this cost compilation.

## Cached post-hoc diagnostics and figures

These actual scripts reuse saved public/exposed data and do not rerun an inverse:

- `diagnose_fixed_states.py` plus the one-time `diagnose_repair_ordinary.py`: final combined receipt contains50FW including the archived initial seed-only attempt and its correction. Re-running the correction is intentionally refused.
- `diagnose_metric_split.py`, `current_curvature_diagnostic.py`, `current_remainder_diagnostic.py`: diagnostic formulas/counts and scope in their receipts; their auxiliary small solves are not dense physical LU.
- `diagnose_plot_exposed.py`: exposure-only material and iteration plots under FIGURES/exposed_development, no16/48GT access.
- `diagnose_bothblocks16_noise_prior.py`: root-authorized cached Y_clean/Y_prediction finite-noise/prior audit, not online GT use; see its own REPRODUCE.md. Missing cached fields cause failure rather than a new forward.

Re-running a report is not required to inspect its saved JSON/CSV/PNG/PDF. Scientific gate decisions, GT exposure, failed starts, capped outcomes, setup, timing uncertainty and negative results remain separate from successful artifact rendering or publication.


## Post-hoc first-CN99-hit diagnostic

`experiments/diagnose_bothblocks16_first_hit.py --phase select` uses only stored public-Phi incumbents, public sigma/Y and the original prior; it seals the first-hit mapping/materials before `--phase score` reads authorized GT. Existing mappings refuse reselection. This is POSTHOC_OFFLINE_DIAGNOSTIC because confirmation GT was already historically exposed, not a new confirmation algorithm or change to the48 contract. CPU prefix timestamps are proxies, not measured rerun early-stop costs. All recorded outcomes and costs are under `delegated/diagnostics/bothblocks_confirmation/first_cn99_hit/`. No new forward, adjoint, threshold scan or optimization is performed.

Vector early DEV-only postseal export/scoring is another producer: `delegated/vector/dev_early_postseal/PROCESS_COST.json` is its own controller; `REMOTE_PACK_COST.json` is a separate remote archive CPU; its external scorer step replaces nested `dev_postseal_scores/SCORING_COST.json` CPU. These costs remain running/null until actual receipts arrive, and cannot authorize fresh vector-confirmation GT access.


## Post-hoc G19-reference CPU budget control

The executed fixed rule/32pair manifest is `delegated/diagnostics/bothblocks_confirmation/posthoc_budget_control/FROZEN_MANIFEST.json`. `diagnose_budget_control_freeze.py` reads only public original child CPU and stop reasons. `diagnose_budget_control_queue.py` executed only its7 triggered pairs through `diagnose_budget_control_run.py`, using unchanged r6 physical corrector plus R11 whole-child observer. Existing destinations refuse overwrite; L-BFGS history is reset.

`diagnose_budget_control_seal_score.py --phase seal` freezes all32 public-Phi selections before `--phase score` reads authorized GT. Original candidates and original80/48 results remain intact. Actual old+new CPU, unused allowance, charged overrun and eligible paid-state CPU are distinct fields. These are not naturally equal-spend executions; no idle-fill or new threshold was used. New7-child cost is additive, while old32/G19 deployment reference costs are not added to the original80 campaign again.
