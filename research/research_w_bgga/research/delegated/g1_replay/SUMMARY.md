# G1 exposed finite24 replay implementation handoff

Status: **implementation written, NOT_RUN**. This worker did not execute Python,
tests, a replay, a physical action, or a new experiment. It did not open locked,
truth, inference, or existing POSTERIORS files. The parent owns execution through
the metered wrapper, mathematical interpretation, and scientific gate decisions.

## Files

- `src/bgga/replay.py`: pure online entry points, finite supplied-response
  methods, and separate offline scoring functions.
- `tools/g1_replay.py`: strict old-input whitelist, PLAN-first output, global
  online seal, independent reference, cost/memory receipts.
- `tests/test_replay.py`: supplied-array protocol/algebra tests, **NOT_RUN**.

## Frozen methods and numerical labels

All seven online method functions accept exactly `(prior, y, oracle)`. The prior
object includes fixed G0 numerical qualification metadata but no response table,
reference, GT, or file API. The oracle accepts one atom index and returns only
that atom's base TX0/1 by 16RX complex response. Every method has a separate
query session; queries never borrow another method's paid full table.

| Method | Implementation | Logical response queries |
| --- | --- | --- |
| M0 | Independent exact finite Bayes | 2950 |
| M1 | Own MAP wall-group, retain its target conditional, Dirac wall | 2950 |
| M2 | Legal nonempty wall maximizing minimum finite-bank material/pose secant retention | 2950 if qualified; 0 and NOT_RUN if qualification is absent |
| M3 | Legal wall medoid minimizing own-posterior matched-pose whitened squared response risk | 2950; all risk work included in its method timer |
| M4 | Unknown likelihood [0,1], descending prior times width per unit cost; no requery | Stop only when unknown-mass upper <=.02; otherwise reach 2950 |
| M5 | Stable first-fit complete-link wall clusters; representative likelihood and original prior conditional lift at each target pose | 2950; all cluster work included in its method timer |
| M6 | Exact own full-table cache/Schur logical control | 2950; measured Schur speed NOT_RUN |

M2 uses nearest legal finite-grid neighbors at each candidate wall's material
nodes and the **same own joint MAP pose** for all candidates. It saves actual
secant atom pairs, spans, scales, ranks, singular values and retention spectra.
The error-resolution rank rule is strictly `s > e`, where `e` is the maximum
over the three new-G0 fixtures of
`||packj(J-J_reference) diag(1,.1,.3,.2,1,.1,.04)||2`.
The fixture IDs are frozen as `strong_center`, `right_large`, and
`low_loss_high_contrast`. Only their `J` and `J_reference` NPZ arrays are loaded.
Absent/partial fixture qualification gives `NOT_RUN`; zero effective material
rank and empty walls remain `NOT_DEFINED`. No machine threshold is substituted.
The label is `EMPIRICAL_BANK_SECANT_RETENTION`, with the additional legacy label
`EMPIRICAL_BANK_DIFFERENCE_GEOMETRY`. It is neither a local-Jacobian claim nor a
uniform derivative/projector/retention certificate. Selection ties use higher
own wall posterior mass, then the stable wall-group ID.

M5 checks **every pair** in a cluster at all 25 prior poses: current potential
diameter <=atanh(.02) and future max-pose whitened response diameter <=.1.
The stable-ID-first legal group is the representative. Quotient cells retain
target pose as well as wall cluster; the lift uses the original within-cell
prior conditional, preserving a full 2950-atom joint output for scoring.
G0 radius/contrast stress does not extend this task's fixed-intrinsic 25-pose
scope. No chain merge or zero-joint-error claim is made.

## Data and phase separation

The CLI fixes source paths to the three entries audited in
`research/delegated/source_audit/CACHE_MAP.json`, checks their exact historical
directory/filenames and rejects locked/inference/POSTERIORS paths. It loads:

- GRID: prior weights, pose, wall group, support ID, four-cell occupancy,
  epsinf, ell. It validates an explicit Cartesian lookup rather than assuming a
  wall-major reshape; the old ordering is pose-major.
- BANK: only `data` (2950x4x16). The oracle slices TX0/1.
- OBSERVATIONS: only `data` (12x2x4x16). Cases are labeled by their positional
  `parent_index`/`repeat_index`; ID, truth, clean-mean, noise and posterior
  columns are not opened.

The first file in a new result directory is `PLAN.json`, before numerical data
loading or any method. Each case writes all seven method NPZ/JSON files and a
seven-row full-atom joint table, retaining an explicit NaN row plus status for
undefined methods. After **all 24 cases x7 methods** are sealed, the CLI writes
`ALL_ONLINE_OUTPUTS_SEALED.json`. Only then does the independent offline stage
compute exact finite references and read the sealed outputs for scoring. It
never calls an online method again after the barrier.

M4 stores the predeclared checkpoint schedule 16/64/256/1024/2048/2949/2950,
all reached checkpoint atom/log-likelihood bounds, full log-evidence bounds,
query traces and unknown-mass upper bounds. An early stop records unreached
checkpoints explicitly. A partial-query displayed estimate is labeled
`QUERIED_CONDITIONAL_POINT_ESTIMATE`; every unqueried original-prior atom remains
in the interval state. It is not relabeled as the true posterior.

## Scores and accounting

For each parent/repeat x7 the independent scorer saves coordinate-CDF sup,
the fixed events x>=0/y>=2.1, full joint TV, the predeclared joint partition
cell masses `(supportID, x<0 versus x>=0, y>=2.1)`, reference-important cells
with mass >=.05, support-only cells separately, and squared-position Bayes
regret `||mean_method-mean_ref||^2` against `.05*R_ref+.005^2`.
Zero nodes follow exact original-coordinate inequalities, with no numerical
zero tolerance. These are predeclared partition cells, not a proof that all
continuous posterior modes were discovered. GT loss remains
`NOT_SCORED_NO_GT_ACCESS`.

The offline diagnostic reports logZ_exact, p0min, and
`p0min/(Z_exact+p0min)` in the log domain. This is a non-tightening diagnostic
for the default [0,1] unknown-likelihood upper envelope; it is **not a lower
bound on true omitted posterior mass**. Bound comparisons save the original
signed atom/likelihood/unknown-mass violations; tiny positive values are not
silently clipped. Floating-point roundoff is explicitly uncertified.

Every artifact carries `REPLAY_NOT_NEW_PDE_SPEED_EVIDENCE`; actual new FE/RHS
are zero. The historical bank's 2950 FE, 11800 RHS and approximately 9288 CPU
seconds remain paid setup in PLAN/input/run receipts. Shared new-G0 setup is
separately read from this campaign's G0 paid reservations plus authoritative
wrapper completion CPU; it is not silently made free or allocated per M2 case.
Per-method timers include their own risk/clustering calculations but exclude
file writes; total phase/process receipts include output IO. No new scientific
gates are created and no physical or Schur speed saving is claimed.

The CLI uses NumPy small matrices and checks actual process high-water RSS
against 1,000,000,000 bytes. Memory estimates, syntax, tests, numerical outputs
and measured peak RSS are **NOT_RUN** until the parent executes the wrapper.

## Parent execution (not run by this worker)

Use the project's existing Python/runtime selection and metered wrapper. Because
this is zero-new-PDE replay, `PREFLIGHT` is the wrapper phase that does not demand
a new physics grant:

```text
python tools/run_job.py PREFLIGHT <unique-test-label> python -m unittest discover -s tests -p test_replay.py
python tools/run_job.py PREFLIGHT <unique-replay-label> python tools/g1_replay.py --output results/<fresh-replay-directory>
```

The replay reads new-G0 qualification files if all three are already available.
Do not alter its frozen rule or scales after observing the finite24 outcomes.
The parent should inspect `FAILURE.json` if execution stops, then the full
`SCORE_TABLE.json`, M4 checkpoint violations, undefined M2 rows, and all paid
setup/CPU receipts before interpreting any G1 result. Exact bound violations
require parent isolation/review under the existing GATES contract.
