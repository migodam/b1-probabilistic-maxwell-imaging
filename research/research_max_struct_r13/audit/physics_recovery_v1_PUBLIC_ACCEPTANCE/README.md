# Public prediction acceptance audit

Campaign `physics_recovery_v1`; 42/42 expected task rows retained; 0 rows have reported wiring/availability issues.

No truth, new Maxwell factorization, optimization or numerical-defect pass gate. Numerical field residuals and exact selector/index issues remain separate. See AUDIT.json for per-view units, relative norms, initial-state max differences, selector/prefix checks, missing counters and stage termination metadata.

Validation matrix actions/CPU are a separate cost, never added to solver physical LU/RHS or used to retroactively create candidates. Shared starts require exact public operators and material initializer; NN groups are separate. Missing tasks/costs are retained as missing/null, never as zero.

Reported issues:
- No exact wiring/index/availability issue detected in these inputs; numerical defects are recorded without a threshold.
