# Infrastructure API

Status: implementation complete; tests deliberately NOT_RUN for parent execution. No physics, providers, network, or old-ledger imports. Owned files: `src/wslam/{geometry,accounting,factors}.py`, `tools/run_job.py`, `tests/test_infra.py`.

## Geometry

- `expand_support(support)` accepts binary 2×2, 3×3, or 6×24 arrays; expands coarse regions into 6 y rows × 24 x columns. `GRID_ID='G_GRID_0125_V1'`, 0.125 cells, ROI x=[−1.5,1.5], y=[0.5,1.25].
- `validate_geometry(support, epsinf=2., ell=0., target_center=None, target_radius=None)` returns `legal,reasons,components,air_components,holes,area,gap,support,epsinf,ell,grid_id`. Empty support is legal. Nonempty support must be 4-connected. Every padded 2×2 block is checked for checkerboard diagonal contacts. Outside air is included when counting 4-connected air components; any number of holes is allowed. Boundary material contact is allowed.
- Scalars or matching coarse/grid material arrays apply only on occupied cells; returned air is exactly epsinf=1, ell=0. Occupied bounds epsinf=[1.5,6], ell=[0,0.25]. Area≤2.25.
- `disk_wall_gap(mask, center, radius)` computes exact point-to-occupied-rectangle distances minus radius, clipped at zero; empty wall gives infinity. Target bounds x=[−1,1], y=[1.6,2.6], r=[0.08,0.18], gap≥0.15. No true-mask input.
- `enumerate_supports(size=2, probability=.65)` returns `(coarse_binary_support, normalized_conditional_mass)` list for legal support. Size 2 or 3. Bernoulli product mass is normalized over legal support, including empty.

## Accounting

`BudgetBook(path=None)` / `Book` use a project-local atomic JSON ledger; default `results/budget_ledger.json`, overridden by `WSLAM_LEDGER`. Avoid `.jsonl` extension because content is a single JSON object. Schema/caps mismatch fails closed. Lock + atomic file replacement serializes each mutation; failed/abandoned physical reservations remain charged and cannot be refunded.

- `reserve(phase,label,rhs=0,factorizations=0,**metadata)`: checks caps before action, returns charged event. CPU exhaustion blocks even zero-count action. Counts must be nonnegative Python integers.
- `event(phase=None,label='',**metadata)`: append diagnostics or `process_cpu_seconds=...`; phase can come from `stage=` or `WSLAM_STAGE`. Events cannot change RHS/LU counts.
- `snapshot(phase=None)`: mapping canonical stage → `{cap,used,remaining}`, CPU key `process_cpu_seconds`.
- `totals(stage=None)`: `{rhs,factorizations,cpu_seconds}`.
- `BudgetExceeded` is raised before denied action. Canonical stages: E0E1 4000/1200/7200; A0, E2E3, E4 each 12000/3000/14400; E5 8000/2000/10800 (RHS/LU/process CPU seconds). E0/E1 and E2/E3 aliases share their canonical stage budgets.

`python tools/run_job.py STAGE LABEL COMMAND...` sets cwd to this new project; PYTHONPATH to src, vendor, sibling research-w-continuation/.deps; standard BLAS/OpenMP thread counts to 1; PYTHONDONTWRITEBYTECODE=1. A project-wide job lock serializes wrappers. Remaining stage CPU becomes child RLIMIT_CPU, and live process-group CPU samples provide an additional aggregate guard. Final CPU charge uses wrapper process CPU plus RUSAGE_CHILDREN CPU, bounded below by peak observed group CPU. Wall and child maximum RSS are separate diagnostics. Nonzero exits and interrupted/spawn-failed jobs are recorded. Descendants remaining in the job group are killed before releasing the lock. CPU guards have polling/kernel time-granularity overshoot; CPU count is measured, never silently clipped. RUSAGE_CHILDREN records descendants reaped within the child tree; detached or unreaped processes are not a supported job pattern. RSS is a wrapper-lifetime child high-water mark, not a sum or incremental allocation.

## Factor ownership

- `FactorRegistry.register(factor_id,kind,version=1,observation_id=None,payload=None)` supports prior/transition/likelihood. A likelihood requires one globally unique observation ID; only the same factor can replace it with a strictly increasing version and unchanged kind/observation. `get`, `factors`, and `len` inspect current registry.
- `refresh_logweights(logweights,old_loglikelihood,new_loglikelihood,normalize=True)` / `likelihood_refresh` computes `log w + newLL − oldLL` then stable normalization. Likelihood logs must be finite; weight logs can be −∞. Changed zero-support likelihoods require reconstruction from the base measure, so the helper fails closed there.

Tests cover boundary/empty/ring/passive materials, connected-but-pinched supports, exact disk-to-rectangle corners and target bounds, conditional Bernoulli mass, multiprocess reservation atomicity/caps/failure accounting, factor uniqueness/versioning, and exact grid posterior replacement. Parent should execute tests in its managed dependency environment.

## Integration corrections

Backend diagnostic `event(..., rhs=..., factorizations=...)` fields become `observed_rhs` and `observed_factorizations`; these never affect charged totals. Atomic ledger serialization handles NumPy scalar `.item()`, array `.tolist()`, and complex `[real, imag]` metadata, including complex arrays. A regression test checks informational counts do not double-charge and scientific diagnostics persist.

Job output now persists under `results/jobs/LABEL/OUTPUT.log` (append across reruns) and an atomically written `RECEIPT.json` (latest attempt). The wrapper prints one compact JSON pointer to both paths. Completion events also store both paths. Labels must be a single directory name. Tests remain NOT_RUN.

Infrastructure tests now use standard-library `unittest` and temporary directories, preserving all cases, spawned-process atomicity checks, and exception assertions. No pytest dependency or package installation. Tests remain NOT_RUN by this worker.

Compilation qualification correction: earlier source checks using `ast.parse` established parsing only. They did not catch finite_bank.py's repeated `outer=` keyword. Parent preserved the failed pre-PDE launch, fixed it to `outer_description`, and compiled all sources. Evidence export preserves that distinction and failed receipt; no prior parse-only statement should be read as successful compilation or execution.
