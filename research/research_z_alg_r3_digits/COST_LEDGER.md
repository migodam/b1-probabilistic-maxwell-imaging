# R3-DIGITS paid-work ledger

State: **CLOSED_RECORDED_SCOPE_NOISE_PARTIAL14_NOT_RUN13**. Total: **KNOWN LOWER BOUND, not complete total**.

Reproduce using `python aggregate_costs.py`; close only after all receipts using `python aggregate_costs.py --close`. No numerical solver, public observation payload, private truth or seed is read.

| Recorded quantity | Known subtotal |
|---|---:|
| physical_LU | 163699 |
| physical_complex_RHS | 6381144 |
| source_RHS | 1964388 |
| adjoint_RHS | 233652 |
| FullJ_receiver_RHS | 4183104 |
| FullJ_builds | 87148 |
| small_Gram_factorizations | 14608 |
| small_Gram_RHS | 701184 |
| normal_factorizations | 142744 |
| green_kernel_evaluations | 1099008 |
| solver_numerical_wall_seconds | 5002.701193999499 |
| solver_cpu_seconds | 4366.794058 |
| worker_process_wall_seconds | 3650.20861328952 |
| setup_wall_seconds | 0.9214688739739358 |
| generation_forward_wall_seconds | 0.24320608424022794 |
| validation_wall_seconds | 3.7896218781825155 |
| scoring_envelope_wall_seconds | 6.905566042754799 |
| audit_envelope_wall_seconds | 0.756445625 |
| prototype_process_wall_lower_bound_seconds | 1.002159583 |
| condition_estimator_calls | 167 |
| condition_estimator_wall_seconds | 0.3947177450172603 |
| kernel_benchmark_wall_seconds | 0.1448685829527676 |
| kernel_benchmark_cpu_seconds | 0.10388399999999998 |
| noise_generation_wall_seconds | 0.055320790968835354 |
| diagnostic_wall_seconds | 0.11961629311554134 |
| diagnostic_cpu_seconds | 0.11177500000000018 |
| validation_cpu_seconds | 0.028067000000000064 |

## Limits

- KNOWN_LOWER_BOUND: obsolete failed prototype physical work and CPU unrecovered. Actual complete scientific cost is unknown.
- Opaque LAPACK ZGECON estimator internal RHS are extra work beyond explicit PDE RHS and remain unknown, not zero.
- Solver numerical wall includes precompute, normalization reference and coarse search within each method budget; excludes public file load, imports/threadpool management and JSON write.
- Worker process wall, summed solver numerical wall, batch elapsed wall, scoring envelope and contained validation wall are distinct measures. Do not sum concurrent worker wall as experiment elapsed, or add validation again to scoring envelope.
- Noise closed PARTIAL:14completed charged,13queued NOT_RUN; no subset noise performance metrics, batch elapsed is pause-inclusive/cancelled and not normal runtime. POSTGT TRF/lowcontrast diagnostics are development calibration only.
- Paired noise repeats the same nine material targets with different observations; it is not 18 independent target materials.
- Budget-filled supplement is isolated after observation generation and before first GT scoring; only incremental accounting contributes, merged root costs are excluded.
- Physical LU/source/adjoint/receiver RHS count once per top-level attempt; nested stages/phases, batch receipt copies and frozen prediction copies are not additional attempts.
- CPU totals are recorded solver CPU only; generation, audits, validation and failed prototype CPU are incomplete.
- No model-priced token or monetary billing receipt available; monetary cost unknown.

## Batches

| Split | Workers | BLAS threads/worker | Batch elapsed s | Worker-process wall sum s |
|---|---:|---:|---:|---:|
| dev | 2 | 1 | 227.987606167 | 437.522363875 |
| dev_v2 | 2 | 1 | 227.377426542 | 451.008863958 |
| holdout | 2 | 1 | 959.661794084 | 1912.814902332 |
| budget_filled | 2 | 1 | 199.994073250 | 388.654575211 |
| data_adjoint | 2 | 1 | 231.407351875 | 460.207907914 |

## Main original plus paid continuation

| Method | Numerical wall sum s | CPU sum s | Physical RHS | Physical LU |
|---|---:|---:|---:|---:|
| fullj_native | 739.050413004 | 645.535452000 | 1054932 | 25787 |
| fullj_reciprocal | 739.043234208 | 649.717392000 | 1075716 | 26255 |
| globaldual | 761.293236790 | 673.943377000 | 1041180 | 28578 |

The table adds each original attempt to its continuation increment; it excludes the cumulative filled-root copy. Actual convergence or no accepted step can leave authorized time unused. These sums are paid work, not concurrent elapsed time.

Per-attempt evidence, scopes and unknowns are in `COST_LEDGER.csv` and `COST_LEDGER.json`. Main outputs can be present before the enclosing batch receipt is written; their process wall remains unknown until that receipt exists.
