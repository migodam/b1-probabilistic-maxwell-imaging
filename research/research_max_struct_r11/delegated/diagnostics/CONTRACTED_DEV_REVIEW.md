# Contracted-lagged exposed DEV execution

All three requested main DEV objects completed locally and sealed public predictions before scoring. The actual contracted-readout closure passed one direction derivative check including the inherited prior (relative error 5.25e-7, zero FW). No runtime failure, retry, ULP guard or budget extension occurred. The original 3000-step auxiliary schedule completed; the 480 CPU cap did not trigger. Two local lanes were used for annulus/shell, followed by bridge. Remote jobs and v1 code were not touched.

|Case|Amplitude NRMSE|Fixed IoU|Phi|CN99|CPU seconds|Corrector stop|
|---|---:|---:|---:|---|---:|---|
|r10_twod_0001|0.103394|0.718310|5.2568317e-05|PASS|68.49|STOP: TOTAL NO. OF ITERATIONS REACHED LIMIT (3000)|
|r10_twod_0004|0.529608|0.081871|0.2141742|FAIL|40.27|CONVERGENCE: NORM OF PROJECTED GRADIENT <= PGTOL (477)|
|r10_twod_0022|0.437098|0.128788|0.0103552|FAIL|61.19|STOP: TOTAL NO. OF ITERATIONS REACHED LIMIT (3000)|

Only the annulus passes the inherited material gate (amplitude NRMSE≤.15 and fixed IoU≥.7) and noise gate. The shell and bridge remain failures. This is an exposed three-case DEV control; it does not establish general performance or novelty. Annulus and bridge correctors reached the 3000-iteration limit. Shell stopped on the projected-gradient criterion with large data residual; do not describe all three as converged.

R11_CPU_CHECKPOINTS and CONTRACTED_DEV_CHECKPOINT_SCORES retain paid incumbents at 60/120/240/480 seconds. For runs ending earlier, later checkpoints retain the completed endpoint, not an interpolated state. All selections use original public Phi. No GT enters the online solver.

Whole-child CPU sampled near exit sums to 169.950853 seconds. Exact external exit rusage was not captured; final serialization/exit tail is unavailable and is not invented. Per-run receipts retain operation counts, setup, auxiliary/corrector timings and stop reasons. Aggregate physical costs: {'cached_operator_imports': 3, 'FULL_FORWARD': 6759, 'actual_complex_RHS': 162024, 'LU': 6759, 'source_RHS': 81108, 'full_factorization_dimension': 3893184, 'DATA_ADJOINT': 6743, 'adjoint_RHS': 80916}. Scoring performs no forward or adjoint.

Artifacts: CONTRACTED_ALL3_SEALED.json; CONTRACTED_WIRING_CHECK.json; CONTRACTED_DEV_SCORES.json/.csv; CONTRACTED_DEV_CHECKPOINT_SCORES.json; CONTRACTED_DEV_COST_LEDGER.json; configs/twod_dev_contracted_tasks.json; RESULTS/exposed_dev/*__contractedLagged__480CPU/. The post-run source snapshot and source timestamps are preserved in CONTRACTED_RUN_SOURCE_SNAPSHOT.py and CONTRACTED_SOURCE_PROVENANCE.json.
