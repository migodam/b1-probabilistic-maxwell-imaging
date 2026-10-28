# R7 independent physics adapter

Implementation: `research/p3l_r7_material_solution_fibers/src/physics.py`.
No closed R6 Book is imported and no historical ledger is written.

```python
from physics import Book, Target
with Book('E2', 'unique-job', parent='parent-id') as book:
    target = Target(book)
    p = target.evaluate(x)                    # x: finite Gaussianized (66,)
    Jx = target.jacobian(p)                    # (512,66), full whitened data
    DchiDx = target.material_tangent(x)        # complex128 (576,66), zero RHS
    next_current = target.solve_anchor(rhs, p, role='jet_1')
    target.close()                            # release sole LU owner
```

`p.x` is Gaussianized; `p.q` is native Gaussian64 × Uniform2.
`evaluate_native(q)` explicitly transforms an interior native point.
`jacobian_native(p)` returns whitened `Jq`; `Jx = Jq * probit_derivative(x)`
column by column. `phi_gradient` returns the likelihood gradient in x.
`objective = Phi + 0.5*||x||²`, including the Gaussianized property coordinates.
`white` uses all 32 receivers × 8 sources: source-major Re32 then Im32,
scaled by `sqrt(2)/sigma`. Original gain is 1; sigma is
`0.00119211798514519`. No acquisition mask or per-parent rescaling is admitted.

The backend is the original sibling R4 frozen P1 decoder/target and original
A30 saved weak background, SciPy CPU, complex128. New geometry/background
construction, alternate solver roots, and historical LU-cache imports are
rejected. The backend is imported lazily only when needed.

## Exact paid charges

| Operation | Complex RHS | States / LU |
|---|---:|---:|
| Fresh `evaluate` attempt | 8 | 1 / 1 |
| First `jacobian` with missing receiver T | 32 | 0 / 0 |
| Reuse J or reconstruct J from an existing T | 0 | 0 / 0 |
| `solve_anchor` with `(576,k)` RHS | k | 0 / 0 |
| `material_tangent`, decoder/whitening/matrix products | 0 | 0 / 0 |

The 32-RHS charge is the original receiver adjoint solve against `S.conj().T`,
yielding `T = S L^-1`. A jet with eight source columns costs 8 RHS for each
solve. Reservations precede factorization/solve; failure reservations remain.
The adapter supplies a physical wrapper for `solve_anchor`; direct naked
backend solver calls are rejected by the Book callback.

## Material derivative and jet convention

For raw complex material chi, `L = I - diag(chi) G`,
`L j = diag(chi) incident`, `e = incident + G j`, and `y = S j`.
Thus for a direction dchi, `L dj = diag(dchi) e` and
`dy = T @ (dchi[:,None] * e)`. `material_tangent(x)` is the decoder derivative
`Dchi/Dq * diag(dq/dx)`, with no extra volume/mass factor.
For a curve chi(t), the second current derivative solves
`L j'' = diag(chi'') e + 2 diag(chi') G j'`.
These identities describe the frozen convention; the fake tests establish
chain-rule/layout bookkeeping only, not Maxwell accuracy or scientific results.

## Memory and persistence

Before any fresh evaluation, the previous Point's state is replaced by endpoint
arrays without matrix/factor/native references, and receiver backend caches are
cleared. The Target retains at most two cached Points and one live LU owner.
Older Points with an already computed J/T remain usable for array-only work;
an older Point missing T cannot request a solver action without reevaluation.
Do not retain a separate reference to the original `p.state` across evaluations.
`threadpoolctl` plus environment settings enforces one BLAS thread; the only
installed dependency was `threadpoolctl==3.7.0` in the existing `.venv_nn`.

Every successful paid endpoint atomically persists x, native q, complex128
material/chi, prediction, current, field, and residual. Once computed, Jx and
Jnative are immediately written into that endpoint, with receiver T included.
Attempts and failures retain their own records. Persistence stores no LU.

## Authorization and recovery

Book requires nonempty valid JSON `configs/E2_FREEZE.json` and
`configs/AUTHORIZATION.json` inside the R7 research directory. Authorization
requires `explicit_user_authorization: true`, optional matching campaign ID
`P3LR7-20261028`, integer `total_complex_RHS <= 50000`, and stage caps
`E0=0, E1=0, E2<=50000`. E0/E1 cannot construct a physical Book.
Optional `controlled_CPU_wall_seconds` bounds job wall reservations.

One nonblocking flock covers a complete worker job. Journal append+fsync
precedes atomic snapshot replacement. Partial/malformed/foreign journals,
snapshot mismatch, missing snapshot after a paid journal, interrupted actions,
and interrupted jobs all block new work. Unique job names prevent replay.
`totals()` validates the journal and snapshot without mutating either. Opening
a Book performs this complete replay once. Live events then check the locked
latest journal row and snapshot against its last committed state, along with
file identity/size/modification time. `_check` uses the same live continuity
check, including forward-credit consumption callbacks. Outside changes fail closed. This avoids
replaying the whole journal at every physics action.

`reconcile_snapshot()` explicitly regenerates only the R7 snapshot from its
valid journal and preserves the displaced snapshot as an audit backup.
`reconcile_interrupted(job, reason)` requires a dead recorded PID and no live
worker; it appends a reconciliation event, retains every RHS/LU reservation,
and conservatively charges the interrupted job's full reserved wall. CPU time
for a dead worker is marked unavailable. It does not repair a damaged journal.
No recovery imports historical R4/R6 snapshots.

The wall deadline is checked before every new paid action; an individual
already-running numerical call is not forcibly interrupted. Signal/KILL may
leave an interrupted journal job, which fails closed on the next invocation.

## Nonphysical verification

Run from the repository root:
`../../.venv_nn/bin/python research/p3l_r7_material_solution_fibers/audit/physics_adapter/test_fake_adapter.py`.
Seven tests cover: chain rule/layout and persistence; one live LU/two Points;
8/32/8 RHS callback charges; cap/failure accounting; exclusive lock; zero E0/E1
and missing freeze; snapshot/journal corruption and explicit dead-job recovery;
live constant-sized continuity checks without full replay.
Close receipts include stdlib `resource` process peak RSS as `peak_RSS_bytes`
(macOS reports bytes; Linux kilobytes are converted to bytes).
The fake suite asserts that `a30.physics` is never imported. All fake charges
are in temporary directories; the R7 production ledger is untouched.

`SOURCE_MANIFEST.json` lists declared paths and BLAKE2b provenance. Its
`loaded_solver_modules=[]` is intentional for this zero-physics audit. Actual
paid jobs write a fresh `SOLVER_SOURCES.json` after loading the backend, including
the paths of the loaded A30/A24 solver modules. No new Maxwell call has run in
this adapter task.
