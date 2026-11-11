# Reusable dense field/current core

`field_current.py` ports the original Pro Section 4 projectors while retaining
the supplied physical backend. It has no truth input or truth-file access.

```python
from field_current import DenseFieldCurrentProjectors, PhysicalAudit, run_projection

projector = DenseFieldCurrentProjectors(
    G, E, S, Y, delta=global_whitened_noise_bound,
    scale=6., lower=0., upper=12., components=1,
    receiver_masks=[mask_for_source_0, mask_for_source_1],
    whiteners=[W_for_source_0, W_for_source_1],
)

def audit(m):
    before = existing_backend.costs.copy()
    objective, J, Q, metrics = existing_backend.evaluate_data_plus_prior(m)
    return PhysicalAudit(
        objective, J, Q, metrics=metrics,
        costs={key: existing_backend.costs[key] - before[key] for key in before},
    )

result = run_projection(
    "dr", projector, initial_material, audit,
    steps=3000, audit_every=50, omega=1.,
    cpu_budget_seconds=remaining_cpu_budget_after_setup,
    # checkpoint_callback receives a complete AuditRecord on each audit.
)
result.save(output_directory)
selected_material = result.material
```

The example callback is an adapter sketch: use the real existing backend API.
Its objective must match the common comparison's data term **and prior**;
selection never uses the auxiliary gap. Physical cost dictionaries are
increments for that callback invocation, not cumulative backend counters.

G has shape `(dofs,dofs)`, E has shape `(dofs,sources)`. S is either a common
receiver map or one map per source; Y is a receiver/source matrix or a list of
source vectors. Masks select rows of S and Y before whitening. The whitening
list supports scalar, diagonal-vector and matrix blocks, including rectangular
blocks. Cross-source correlated whitening is not implemented. All source blocks
use one delta and one Lagrange multiplier. Identical whitened port maps share
their H solves and port eigensystem. Null port eigenvalues define the reachable
data floor; impossible balls raise `UnreachableDataBall`, boundary balls use the
pseudoinverse limit rather than seeking an infinite finite root. Numerical port
rank uses `32*eps*max(dofs,receivers)*max_abs_eigenvalue`.

For scalar fields use `components=1`. For isotropic vector fields use
`components=3`: rows MUST be cell-major (`cell0-x,y,z,cell1-x,y,z,...`), and the
material length is `dofs/3`. The C projector aggregates all source AND component
entries in a cell before choosing one real material in its box. It uses stable
quadratic roots, evaluates both box endpoints and actual squared distances, and
deterministically chooses the first minimizer in tied/dark cases. Bounds may be
scalar or per-cell finite vectors.

AP and DR share projectors and physical audit selection. The initial material
is physically evaluated once. Each scheduled audit evaluates the C shadow's
material using the existing physical backend. The incumbent starts at the
initial material and changes only for a strictly smaller physical objective.
The final iteration is also audited if it is not a scheduled multiple. No
extra physical solve is forced after a budget exit; the last un-audited state
and the last audited incumbent are both retained. Set the same audit schedule
for both methods. DR requires `0 < omega < 2`; default omega is 1.

Each audit stores the pre-step and post-step trajectory state, both projected
shadows, physical candidate fields, candidate and incumbent materials, incumbent
physical fields, objective and metrics, CPU/wall timestamps and cumulative
physical/fixed-background cost ledgers. `save` writes `history.json` and
`states.npz` without pickle. A `checkpoint_callback(record)` can serialize each
audit immediately. The returned last state can be later than the last audit if
a budget expires between audits.

Setup CPU/wall includes forming H, one H Cholesky, setup H RHS and port
eigendecompositions. Runner CPU/wall begins before the initial audit and excludes
setup; the caller must subtract setup from its total proposal budget. Budget
checks occur at iteration boundaries: an in-progress iteration, audit or
checkpoint callback can overrun. Physical counters come only from the callback.
Fixed-background counters include H RHS, factor/eigen counts, G/G-adjoint and
receiver actions, A/C counts, and multiplier evaluations. Receiver-adjoint counts
refer to cached H-inverse receiver-adjoint actions, not new H solves. Result
`costs` reports runner increments; audit `costs` is the projector's cumulative
ledger including setup. Use a fresh projector per independently charged arm,
or explicitly account for shared setup once.

The constructor's optional `grouped_ports=True` selects v2. Default `False`
preserves the initially frozen v1 implementation. V2 groups identical port maps
into matrix RHS actions and caches aggregated eigenspace powers for multiplier
root finding. It does not alter the objective, metric, ball or mathematical
projection. Action RHS counts retain the same meaning; call counts decrease
when source columns share a matrix action. `benchmark_grouping.py` times one
existing Pro-toy background projection workload and checks the two outputs on
that workload. Its measured speedup is not a universal speed or recovery claim.

Small wiring check (no imaging study):

```sh
OPENBLAS_NUM_THREADS=1 python research/research_max_struct_r12/src/wiring_check.py
```
