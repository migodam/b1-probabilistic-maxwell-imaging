# Data, blindness and cost review

This worker changes only the new `research/research_max_struct_r1` tree. It reads frozen R3/R4 reports and a small set of named receipts. No historical algorithm, original evidence, provider, precision, prior, acquisition or protocol is changed. No GPU, NN, SHA check, push or new blind imaging is run.

## Data and disclosure

The executed diagnostic loader opened **only `G`, `E0`, `S`, `grid_shape`** from frozen `research/research_fbg_r3_experiments/data/dev_calibration.npz`. It did not open that file's `Y` or any stored material field. Those four public arrays are now copied into `data_public/dynamics_operators_4x4.npz` for standalone packaging. The current reproducer uses that local file and included `frozen_r3_physics.py`/`frozen_diversity_physics.py`; only import/input routing changed after execution, with no new numerical run. It chooses a disclosed nonconstant 16-cell material, computes its physical root, and uses this development truth for local mechanics; a second disclosed material puts one cell on the common lower bound. Consequently D1/D3 are **OFFLINE_DIAGNOSTIC**, with no blind prediction/seal or online image ranking. Their final physical residual is inapplicable to an imaging gate; they are tests of iteration maps.

Existing genuine Hankel–Green matrices are reused unchanged. The R3 adapter directly depends on `research_z_alg_r4_diversity/src/physics.py`, absent from the isolated checkout. An unmodified file copy is kept at `research/delegated/audit/frozen_diversity_physics.py`, copied from the existing publication repository at `/Volumes/migodam's-external-brain/Research/Inv_SLAM/research/delegated/z_alg_r2_publication/repository/research/research_z_alg_r4_diversity/src/physics.py`. The executed diagnostic bound this frozen module before importing the tracked R3 adapter. For standalone packaging, `frozen_r3_physics.py` now preserves that adapter with only its import rerouted to the included module; the original Git source pointer and change are in `audit/D_INPUT_PROVENANCE.json`. Its retained LU, real material and common prior are reused; this is not a new Maxwell backend.

For historical E3, the scorer's post-closure truth-open record is preserved, but the objects themselves were already disclosed. They remain historical development evidence. The older E1 eager-truth access before complete closure remains a protocol failure; later lazy loading cannot repair the past. New confirmation sets and R4-BASIN fresh predictions have their own separate disclosure/closure contracts and cannot confer blindness on E3 or this diagnostic.

## Executed new work and costs

Machine/interpreter: existing Apple CPU scientific venv `Gaussian/.venv_nn/bin/python`, NumPy/SciPy, threadpool limit1. No package installation. Reproduce from the packaged study root:

```sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 /path/to/existing/scientific/python src/dynamics_audit.py
```

The executed successful numerical process records **.04976 s wall, .04155 s CPU**, including input/root setup, diagnostic execution and D1/D3 output serialization within `main`; interpreter imports/startup and final cost-JSON serialization are outside this timer. These very small values are bookkeeping, not a speed claim.

| Operation | Paid count |
|---|---:|
| Original physical LU / original physical complex RHS |2 / 12|
| Exact-current small Cholesky / complex RHS |450 / 450|
| Analytic-root-JVP small Cholesky / complex RHS |96 / 96|
| Exact-current blocks |75|
| PR current steps |105|
| Coupled sourcewise line solves (6×6 real least squares) |64|
| Material readouts |180|
| Common-prior material inner solves / evaluations / iterations |3 / 24 / 9|
| Explicitly counted full source-metric products |490|

The original physical counts include one shared root forward and one active-box root forward; all six sources are RHS of one retained LU per root. Cholesky counts are **not** physical LU or new forward solves. Dense `G/S`, Gram and regularizer products are paid in measured wall/CPU but not exhaustively enumerated as separate FLOP counts; absence from the operation dictionary does not imply zero cost. No FFT executes (explicit matrices). A16×16 spectral norm, one6×6 unitary QR, an eigenspectrum and64 small least-squares line solves are paid; their internal LAPACK factorization/SVD actions are not separately instrumented and are not reported as zero. Historical assembly is not recharged or called free: its wall is `null` with an explicit not-rerun reason.

Raw run receipts are in `research/delegated/audit/D1.json`, `D3.json`, and `dynamics_cost.json`; complete D3 material/current trajectories are in `D3_trajectories.npz`. Receipts retain code/config hashes as `null` under the user's no-hash instruction and identify source paths and parent Git provenance. Blind seal, candidate ranking and GPU fields are `null` with reasons. Physical residual gates are inapplicable: no final imaging task is attempted. Parent consolidation must retain its actual frozen base/commit path provenance.

Three startup failures precede the successful numerical run: an unclosed parenthesis (corrected before any import), system Python lacking NumPy, and the scientific venv import failing because the isolated checkout omitted the frozen diversity dependency. None reached root preparation or any numerical task; original physical LU/RHS are known zero for them. Their wall/CPU were not instrumented and remain **null/UNMEASURED**, not zero. Additional interpreter/path discovery probes used no physics. Read-only historical extraction and report composition were not separately timed and are **UNMEASURED**.

## Historical E3 cost audit: preserve actual consumed compute

Re-summing the 48 already-scored rows produces:

| Historical arm | Joint pass/12 | Worker wall sum/s | Physical LU | Actual complex RHS |
|---|---:|---:|---:|---:|
| Native TRF |0|296.8810|1192|71520|
| Fixed beta=.2 feedback |8|421.3715|20391|489024|
| CIE single + corrector |8|540.1542|15926|381660|
| CIE dual + corrector |10|540.2338|11825|283272|

These are referenced historical costs, **not new paid MAX-STRUCT-R1 runs**. Single CIE is wall-truncated in9/12; dual in11/12. Only3/12 single and1/12 dual finish their scheduled stages, and `SCHEDULE_COMPLETE` still permits MAXITER stages. The historical finite-set result is a45-second budget point, not a convergence ceiling. Two-lane shared-resource execution, order/warmup, setup sharing and differing stopping/corrector work prevent calling nominal wall envelopes matched actual compute. Worker wall sums are not parallel elapsed wall. Historical unmeasured failures and startup attempts remain separate from complete receipt sums.

## Review limits and required report language

Do not write that prior-free beta-invariant root derivatives imply identical global dynamics; D1 gives finite-step/prior/PR counter-controls. Do not assign a Fréchet spectrum to the one-PR directional map without differentiability. Do not call active clipping a measured beta-invariance failure when this experiment preserves first-order coincidence. Do not call full objective covariance alone complete sourcewise-iteration covariance. Do not attribute the sourcewise D3 difference only to PR history when its first step already differs. Do not reconstruct missing historical PR states or turn a wall-cut endpoint into a basin fixed point.

D1/D3 are adequate for these bounded diagnostic claims. They neither establish a new globally convergent method nor improve the historical 10/12 image score. Parent F/M results require independent output/receipt review before scientific consolidation.
