# R11 2D execution wiring

`src/twod_budget.py`: no scientific mathematics; hooks genuine R8 PaidStateObserver, keeps all paid physical evaluation metadata and float64 new-incumbent trajectory, checks absolute whole-child CPU. CPU-budget control uses a dedicated BaseException so inherited Exception handlers do not misclassify it as a numerical failure. The inherited solver finally receipt remains available.

`experiments/twod_run.py`: fresh R11 wrapper imports unchanged actual R8 run_case/benchmark_solver. Noise99 selection uses unchanged R10 public-only adapter helpers. Runtime-only dispatch adds budget_check and supports frozen task `material_readout={module_file,callable,prepare}`; prepare(model,numerical) returns config, then material_readout is called by unchanged R8 CIE. Root candidate should use arm=a_plus. No R10 adapter output-path bypass or source modification.

Outputs: inherited snapshots and costs; R11_PAID_PHYSICAL_STATES.json; R11_INCUMBENT_TRAJECTORY.npz (float64 m, originalPhi, absoluteCPU, physicalLU/RHS); R11_CPU_CHECKPOINTS.json + material NPZ; R11_CAP_PREDICTION exact already-paid m/j/Q/Y; R11_ACTUAL_PAID_INCUMBENT and LAST_PAID_STATE preserve complete states including atomic overshoot. Final receipt explicitly separates cap eligibility from charged actual cost.

`experiments/twod_queue.py`: copied mature R10 finite queue with R11-owned paths and budget-stop success handling. Continues using one-process-per-task, BLAS thread limits, own-worker affinity/below-normal priority, STOP_DISPATCH, no retry, Windows GetProcessTimes. Do not launch without an authorized free CPU allocation. GETProcessTimes/controller/internal views are nested and never summed twice.

Verification: Python compile and deterministic fake-observer clock test, no physics/no GT. No data generation or inversion started. A real DEV run remains required once root freezes candidates and resources are free. No changes outside assigned R11 paths.

Known scope limits: R8 rejects cross-grid geometry; final internal CPU timestamp cannot include the write that records it, so whole-child external total is authoritative. Cooperative caps can overrun by one atomic call and final serialization. Checkpoint prediction uses strict completion times with no interpolation. Natural mature stopping retains actual spent work.

## Completed exposed DEV

All5 lagged DEV plus1 fullG19 replay sealed and synced locally under RESULTS/exposed_dev. No GT scoring by this worker. FullG19 replay m/j/Q/Y_prediction exactly equals R10 public prediction. TWOD_DEV_COSTS.json/csv/md preserve actual disjoint GetProcessTimes and known CPU lower bound; unknown control/setup overheads are explicit. All four inherited auxiliary stages completed in every task, with natural corrector stopping or3000 limit; no cap exhaustion.

For custom material readouts, HISTORY.state_term_material remains the legacy onehop diagnostic. Read material_final_block_loss and DENOMINATOR_SETUP.derivative_policy for the actual fixed-reference block objective.

## Current independent48 and supplementary execution

The earlier compile-only verification describes the initial wiring stage. Independent48 is now an authorized real finite288-task campaign, with all six frozen methods, whole-child480CPU allowance, natural mature stopping, and public60/120/240/480 checkpoints. Current counts come from `dev_transport/INDEPENDENT48_STATUS.json`, with complete public arrays synced in batches. The canonical active-attempt registry, not a raw count of all process attempts, defines the288 scientific outcomes.

Four original controller-interrupted attempts remain preserved with incomplete endCPU; three subsequent prep-only failures have exact billed CPU and zero physics. The latter use explicitly registered retry03 output/prepared directories. No finished original132 outcome is rerun. `independent48/COSTS.json` separates canonical scientific prefixes from all incurred attempts and unknown interrupted costs. Main48 keeps all objects; object0022 receives an additional labelled exclusion sensitivity and cannot support an exact end-to-end same-cost assertion.

Primary GT, secondary reindex and plots must use `twod_postprocess_run.py --kind primary|secondary|plots -- <unchanged authorized child arguments>`. Its external wait4 child CPU replaces nested self timers, and `POSTPROCESS_EXECUTIONS` is the authoritative disjoint postprocess ledger. The public selector children are already measured externally by the sync controller. These commands have been prepared, not authorized GT execution.

The supplementary64 controller uploads only the frozen public paired-grid8 archive. It waits for all canonical288 terminal outcomes including retry03, all main public seals/full synchronization, and the natural exit of the original own queue/workers. It then uses the same four physical cores and ordinary solver. A coarse inverse model always receives coarse operators/geometry; cross48to24 denotes independently generated fine-grid observations, not an unsupported cross-grid inverse adapter. No GT is transported or scored by this controller.
