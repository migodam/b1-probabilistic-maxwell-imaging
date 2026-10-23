# Wide-frequency continuation driver: build-only delivery

Implementation: `src/wr/frequency_continuation.py`. AST parsing PASS.
Physical execution, imports, and new RHS/factorizations: NOT_RUN / zero.
The campaign owner must start this driver in the existing serial physical queue.
Current running files and their data have not been changed.

## Frozen inputs and schedule

The driver reads only `inputs/mean_<seed>_wide.npz`, public original configuration,
and the shared measurement-only initialization APIs. It calls exactly the same
`born_initialize(y, full_wide_frequencies, grid)` and `selected_groups` as
`online_inference`. One `GroupPixels` support and parameter vector are reused
throughout, with no GT/oracle import and no active-support changes.

| Stage | Frequencies | TR iteration limit | Weak numerical prior |
|---|---|---:|---:|
| 1 | 0.8 | 8 | 0.002 |
| 2 | 0.8, 0.9 | 6 | 0.002 |
| 3 | 0.8, 0.9, 1.0 | 6 | 0.002 |
| 4 | 0.8, 0.9, 1.0, 1.1, 1.2 | 10 | 0 |

Bounds `[0,1]`, shared box trust-region solver, `l2_radius=False`, radius 0.15,
maximum radius 0.5. Each stage starts from the preceding stage's retained
solution. `iteration_scale` may reduce all stage limits, never enlarge them;
stages and frequencies remain unchanged. The original material prior and noise
model remain unchanged. Only the explicitly requested weak optimization penalty
is absent in the final stage. No residual-only convergence/floor stop is enabled.

Callable: `main(seed, iteration_scale=1., folder=None)`.
CLI: `python -m wr.frequency_continuation SEED [--iteration-scale SCALE]
[--folder RESULTS_RELATIVE_FOLDER]`, using the established campaign runtime and
Python paths. The CLI exits with status 2 if all four stages did not complete.
This report does not authorize a concurrent physical launch.

## Budget and independent receipts

Before constructing any `Forward`, a whole-sequence affordability gate checks
all four planned stages against remaining R RHS/factorization caps. Each stage
is then rechecked immediately before execution. CPU checks require remaining
time; elapsed runtime is still enforced by `Book` at individual reservations.

Each stage allows `max_calls=2*iterations+10`. A driver-side guard enforces at
most `1+iterations` full-J calls **before** calling the physical adapter; all
remaining allowed optimizer calls are costed as F-only. Every physical evaluator
exception is fatal for that stage (`fatal_exceptions=(Exception,)`), so failed
FJ refreshes cannot quietly initiate repeated derivative attempts. This includes
numerical failures, shape/type errors, and budget stops. The forecast is:

`frequency_count * [4*max_calls + 12*(1+iterations) + 4 final_forward]`.

The final potentially uncached F is included. This yields 216, 352, 528, and
1280 complex RHS (total **2376**) and 27, 46, 69, and 155 factorizations (total
297). Normal cache reuse may reduce actual cost; the call/derivative caps prevent
it from exceeding the forecast, aside from independently enforced CPU limits.

If the whole sequence is unaffordable, all stages are skipped before physics.
If a later stage becomes unaffordable, it and all subsequent stages are explicitly skipped in their
`RESULTS.json` with empty histories. No partial stage is quietly substituted,
no new budget is created, and the driver does not precharge/double-count its
forecast. Runtime `Book.reserve` remains authoritative. The serial queue is
required because a forecast does not lock out another campaign process.
CPU duration cannot be predicted; runtime budget rejection or physical failure
is recorded as a failed stage, retains anchors, skips all later stages, and
never initiates recovery physics. The existing launcher must bill process CPU
once; the driver reports measured CPU without charging it a second time.

Unique run/stage ledger phases isolate real physical receipts from other runs.
`RHS_RECEIPT.json` independently sums reserve IDs, charged complex RHS,
factorizations, confirmed `backend_solve` RHS, and confirmed factor events.
Reserved cost and confirmed completed backend work are kept separate, including
failed calls. Optimizer call counts are not used as physical cost receipts.

Each stage has a separate `Forward` endpoint folder, `FWI_HISTORY.json`,
`FINAL.npz` (or `RETAINED_ANCHOR.npz` on failure), `RESULTS.json`,
`BUDGET_FORECAST.json`, and `RHS_RECEIPT.json`. Top-level
`SEQUENCE_BUDGET_FORECAST.json`, `FINAL.npz`,
`STAGES.json`, `RESULTS.json`, and `RHS_RECEIPT.json` preserve the overall result.
If continuation stops early, final data are explicitly labelled with their
actual frequency subset; all-wide distance is null and status is incomplete.
Completing iteration allocations is not claimed as projected-gradient
convergence, geometry validation, model adequacy, or scientific success.

## Restricted-support limitation and separate proposed extension

A fixed measurement-selected mask cannot represent a legal feature outside its
selected tiles. Frequency continuation does not repair that representational
restriction. No support exchange is implemented or executed here.

If a subsequent, independently budgeted study is authorized, consider the
following bound-aware **outer** proposal mechanism; never modify active support
inside a corrector:

```text
finish bounded corrector on current support S; retain its accepted anchor
compute measurement-only residual and inactive-tile directional scores
rank inactive additions by feasible one-sided objective decrease at x_add=0
rank removable active tiles using feasible bound-aware loss/gradient evidence
propose a small area-neutral exchange S' (legal tiles, fixed area/feature cap)
map retained parameters to S', initialize additions at feasible zero
run a separate fixed-support bounded corrector under its own frozen budget
accept only a preregistered measured-objective improvement; otherwise retain S
record proposal, all trials, geometry checks, and actual RHS costs
```

This is a recommendation, not a claim that such exchanges improve the present
inverse problem or satisfy physical geometry constraints. Changing support
would change the restricted inverse class and needs its own reported protocol.
