# Cost ledger: completed recorded scientific work

Final receipt inventory: `audit/cost_ledger/COST_LEDGER_FULL_V2.json` and `.csv`. Rebuild with `python research/research_max_struct_r4/experiments/build_cost_ledger.py --full`. This reads saved JSON only and executes no physics or tests. Historical `COST_LEDGER.json/csv` and `COST_LEDGER_FULL_V1.json/csv` are preserved. V2 is the current accounting reference; older files retain their historical status.

The inventory contains 66 scientific/audit receipt rows. Known physical LU subtotal is 38757; model complex RHS excluding nonlocal locality solves is 980383. Counts are not a uniform-compute metric:

| Physical dimension | Physical LU | Model complex RHS excluding locality |
|---|---:|---:|
| 576 | 36698 | 931161 |
| 1024 | 2054 | 49116 |
| 16, independent port wiring | 5 | 106 |
| Known subtotal | 38757 | 980383 |

Nonlocal T-matrix locality solves are separate: 30 LU/17280 RHS. Its original top counter 1112 LU therefore splits into 1082 physical plus 30 locality LU. Full-operator retained solves, dense T reconstruction, rays and receiver resolvents remain recorded in their native counters; RHS totals do not imply identical matrix dimensions or operations.

Small4x4 totals (including the one recorded generic-block wiring direction) are 26710128 factorizations and 295349184 complex solve RHS; cached matrix applications are 252054144 RHS. Scalar complex divisions total 13611979008 and onehop block matrix applications 2313709056 RHS. Sparse partition-average neighbor applications remain in the separate auxiliary operation dictionaries. These arithmetic counters must not be equated to a 576-cell physical factorization.

## Included independent work

The preserved development inventory includes the failed initial T-matrix import, compute probe startup failure, mature/native/block/self/cached/onehop/phase arms, curvature/rays, port reference and envelope wiring, invariants and recorded postprocessing. Historical R2 reference 689 LU/16452 RHS is identified rather than recharged. Common-global is a new corrector on a saved historical auxiliary candidate, not a new end-to-end auxiliary run. Replays and failed attempts are chargeable even when scientific settings are unchanged.

The six-case confirmation includes 18 downloaded task `result.physical_model_cost` totals once each, plus the shared generator 6 LU/72 RHS once. Grid32 includes one generator 1 LU/12 RHS and the two 1024-dimensional solves 1027 and 1026 LU. Windows independent replay adds 1032 LU/24684 RHS. Post-confirmation partition-average adds two independent solves totaling 2062 LU/49320 RHS; endpoint stationarity adds 6 LU/144 RHS; onehop/native endpoint continuations add respectively 267 LU/6432 RHS and 3385 LU/81264 RHS. New scoring, plot and saved-data comparison/stop audits explicitly record zero new physics where their receipts say so. All these scientific executions are complete; V2 itself adds no science.

Each top receipt contributes once. Nested native `auxiliary_cost.model`, adapter stage costs, `completion_cost`, `corrector_cost`, operation timers and duplicate `CORRECTOR.json` are nonadditive views. Generic block/onehop embedded material wiring counts are added once to their auxiliary small-operation counts, not recharged through duplicate standalone files. Port reference and dimension16 envelope wiring are independent paid invocations. Box rays reuse already charged cooperative full geometry; the free-face curvature diagnostic is a distinct run.

## Timing, failures and unknowns

Confirmation actual CUDA worker count is zero. Its independently measured orchestrator elapsed is 1179.0318317911588 seconds, including remote startup, finite queue tasks and shutdown. Worker wall totals are not parallel elapsed. Windows solver wall is 376.1966163000325 seconds; its wrapper execution elapsed is separately recorded in `audit/confirmation/WINDOWS_VALIDATION_RUN.json` (377.40292970929295 seconds), with upload/archive/download overhead in transport records. Do not add confirmation and Windows intervals and call that total parallel time. No continuous timer covers the entire session's overlapping local/remote activity, so total session elapsed remains unknown/null.

Caller, solver, stage, operation and once-per-process import times retain their distinct scopes. Grid32 caller timers include more work than adapter solver timers. Original bootstrap imports and operator assembly are not invented as zero cost. Transport and ZIP normalization remain separate evidence; their elapsed is not added again to a containing orchestration timer.

The first Windows wrapper launch failed before science because the existing vendor path was missing: `WINDOWS_VALIDATION_RUN_STARTUP_FAILURE.json` retains paid 0.6626318749040365-second execution and both uploads. The later successful execution and its transport are independently retained. Windows ZIP backslash extraction failures and normalization are documented in the original audit receipts. Missing startup timers, absent postprocess counters and incompletely recorded repeats remain unknown. `missing_fields` lists unavailable values, while retained failure details preserve errors/refused attempts where recorded. The inventory cannot prove unrecorded retries absent.

Project monetary amount, currency and LLM fees are unknown/null. Git base is ordinary provenance; hashes were not verified. This is the known recorded-cost inventory, not a complete financial bill or a scientific judgment. No confirmation GT/private arrays were opened by the ledger builder.

Report/archive packaging is subsequent delivery work, with zero physical LU/RHS; its own wall/CPU receipt is `audit/EXPORT_RECEIPT.json`. Publication/readback is recorded separately in `audit/PUBLICATION_RECEIPT.json` after it occurs, rather than changing the sealed scientific cost denominator.
