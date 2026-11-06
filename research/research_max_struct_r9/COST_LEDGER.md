# R9 actual cost ledger

All113 inverse receipts (DEV47 +CONF66) are terminal. This ledger reads public cost receipts only: no GT arrays, scores, new solves or plot reruns. It reports a known CPU lower bound, not exact total research cost.

Known measured CPU including retained failed subsets and cost-compiler invocations: **4353.102408s**. Physical LU **179912**, complex RHS **4315260**. Currency and energy remain null.

| Category | Receipt count | Known CPU(s) | Known task-wall sum(s) | Physical LU | Complex RHS | Current Cholesky | FFT batch calls |
|---|---:|---:|---:|---:|---:|---:|---:|
| development_inversions | 47 | 1804.348707 | 1845.049701 | 76987 | 1845624 | 2829 | 1512420 |
| confirmation_inversions | 66 | 2534.792004 | 2597.139873 | 102900 | 2466660 | 1768 | 2592720 |
| curvature_diagnostics | 2 | 0.179027 | 0.189135 | 14 | 2760 | 0 | 0 |
| gradient_wiring | 1 | 0.083160 | 0.000000 | 5 | 144 | 3 | 0 |
| confirmation_generation | 1 | 0.476146 | 0.483195 | 6 | 72 | 0 | 0 |
| analysis_scoring_plots | 8 | 13.083116 | 14.349461 | 0 | 0 | 0 | 0 |
| measured_archived_import_failure_subset | 12 | 0.098573 | 0.101129 | 0 | 0 | 0 | 0 |
| cost_compiler_measured_subsets | 2 | 0.041675 | 0.070566 | 0 | 0 | 0 | 0 |

Nested setup, auxiliary/corrector timers, per-object generation timers and scoring cost ledgers are not added again. Inverse totals come directly from raw execution receipts, so curvature diagnostic cost is counted once in its own category rather than twice via attributed supplement summaries. Profile first-stage cost is already paid within full profile runs; weak standalone reuse and portfolio recipe sums are excluded from actual totals. Historical R8 endpoints are not new R9 jobs.

## Timing and logging boundaries

Inherited baselines use inclusive r9_total_cpu_seconds (imports, preparation, solver callbacks, paid snapshots and serialization); landscape uses process_CPU_seconds. A+ readout microsecond counters, physical-operation timers, inner auxiliary/corrector CPU and landscape startup are retained per receipt in ALL_COSTS.json. They are overlapping subsets. Isolated logging CPU is unknown; the inclusive baseline is not mislabeled pure arithmetic. Peak counters are maxima across receipts, not summed memory occupancy.

## Failures and missing metering

The28 initial wrong-interpreter imports used zero physics. Only12 landscape children recorded CPU; their measured subtotal is included, while16 baseline failure CPU entries remain null. The launcher-path typo and first scorer-schema error each retain a source artifact and unknown CPU/wall. The layout-repair plot attempt is included when present in the plot-attempt array. Literature/orchestration, seal checks, warm-initializer exports, plotting imports explicitly outside timers, and other unmetered tooling are not free: total cost remains a lower bound. New training count is0; monetary/energy cost is null, not0.

## Queue elapsed time and overlaps

Queue interval union=1192.657607s; exactly one queue active=1192.657607s; multiple queues active=0.000000s. First-to-last queue span=1848.434893s, including orchestration gaps. Sum of measured queue makespans=1192.613760s is separate from UTC intervals; the two clock boundaries are not asserted identical. The recorded queues are nonoverlapping in this run. Task-wall sums above are concurrent workload, not makespan or total research elapsed.

| Queue | Jobs | Measured makespan(s) | Exclusive nonoverlap(s) | Overlap with another queue(s) |
|---|---:|---:|---:|---:|
| confirmation_queue | 54 | 596.761000 | 596.776339 | 0.000000 |
| confirmation_weak_handoff_queue | 6 | 23.869012 | 23.868437 | 0.000000 |
| cost_controls_queue | 7 | 101.677308 | 101.676257 | 0.000000 |
| curvature_queue | 2 | 0.737080 | 0.737063 | 0.000000 |
| development_queue | 28 | 3.726989 | 3.726427 | 0.000000 |
| development_queue_02 | 36 | 371.844859 | 371.876526 | 0.000000 |
| development_supplements_queue | 3 | 16.240015 | 16.239583 | 0.000000 |
| final_controls_queue | 7 | 77.757497 | 77.756975 | 0.000000 |

UTC interval segments and actual timestamps are preserved in ALL_COSTS.json. Queue-active overlap is not a claim about SM/core utilization. Queue-controller wall or child CPU is not charged twice.

## Weak handoff standalone recipes

These are method costs only; actual new continuation cost is already in the inverse table. Initialization includes exact paid first-stage prefix plus measured startup, before unneeded poststage details factorization. Necessary current Cholesky/RHS are retained.

| Task | Actual new CPU(s) | Reused init CPU(s) | Standalone CPU(s) | Standalone LU | Standalone complex RHS | Standalone current Cholesky |
|---|---:|---:|---:|---:|---:|---:|
| obs_40ac5b79e8a74ac3__weak_profile_handoff | 16.352524 | 3.281207 | 19.633731 | 1372 | 32916 | 105 |
| r9_confirmation_0000__weak_profile_handoff | 6.462406 | 2.113123 | 8.575529 | 309 | 7416 | 75 |
| r9_confirmation_0001__weak_profile_handoff | 22.397152 | 2.626942 | 25.024094 | 1344 | 32244 | 102 |
| r9_confirmation_0002__weak_profile_handoff | 3.781059 | 3.001512 | 6.782571 | 181 | 4332 | 108 |
| r9_confirmation_0003__weak_profile_handoff | 22.534557 | 2.574194 | 25.108751 | 1355 | 32508 | 104 |
| r9_confirmation_0004__weak_profile_handoff | 1.241010 | 2.869254 | 4.110264 | 35 | 840 | 86 |
| r9_confirmation_0005__weak_profile_handoff | 13.879602 | 3.059822 | 16.939424 | 804 | 19284 | 106 |

Full additive physical/auxiliary counters, FFT transforms, extra rank-selection SVD counts, timing subsets, unknown costs, source receipt paths and source helper are retained in audit/ALL_COSTS.json. Rebuild with experiments/compile_cost_ledger.py; rebuilding only re-reads metadata, appends its measured compiler invocation, and refreshes these derived reports. No raw receipt is changed.

完整扩展、远端、三维和发布费用见[COST_LEDGER_EXTENSION.md](COST_LEDGER_EXTENSION.md)。本文件保留第一阶段冻结口径，不覆盖其原始结果。
