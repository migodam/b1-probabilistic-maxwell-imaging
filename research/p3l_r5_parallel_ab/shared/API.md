# R5 shared physics API

Run a fresh process from repository root with `VECLIB_MAXIMUM_THREADS=1 OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1` and `../../.venv_nn/bin/python`.

```python
from research.p3l_r5_parallel_ab.shared.physics_adapter import Target
T = Target('A', 'declared_job', output_dir='research/p3l_r5_parallel_ab/stage_a_response_coordinates/physics/declared_job')
chi, aux = T.decoder(x66)                 # original P1; zero new Maxwell RHS
p = T.evaluate(x66, purpose='declared')    # 1 state +8 forward RHS
J2 = T.jacobian(p, V66x2)                 # white512x2;16 derivative RHS
J66 = T.jacobian(p)                      # white512x66;32 adjoint RHS on first transfer
jv = T.jvp(p, v66)                       # white512;8 derivative RHS
grad_phi = T.gradient(p, y_complex32x8)  #8 adjoint RHS
phi = T.phi(p, y_complex32x8)
y = T.load_observation(explicit_online_path)
print(T.book.totals())
```

`p.x` is Gaussianized66, `p.q` is native Gaussian64×Uniform2, `p.prediction` is complex32×8. All Jacobians are with respect to Gaussianized x, already whitened. Do not multiply by the probit derivative again. For an affine subspace x=x0+Vh the supplied V must contain full66-coordinate directions.

`retain_reference=True` on evaluate retains LU/current/field only for designated references. Otherwise compact endpoints persist immediately; returned Point holds one live LU for local derivative use. Caller must release unused Points. Adapter receiver transfer cache holds at most one Point.

`Target('B', ..., selection_frozen=explicit_file)` requires the concrete frozen selection/protocol path. This file is a necessary execution gate, not permission to change the fixed budget or historical anchor definition. New truth generation is not implemented or run by this adapter.

Accounting and qualification PASS are recorded in `audit/`. Qualification consumed3states and80RHS, included in A's fixed allocation. Historical strong contexts use96native+64probit optimizer limits, with median3376–3616RHS/parent; a40RHS supplied-x forward/J is a different context regime.

A5 inverse reference support (added after the initial qualification; no new physics run by the shared worker):

```python
R0 = T.inverse_action(p, np.eye(576, dtype=complex),
                      purpose='declared_A5_inverse_norm_probe',
                      persist_path='research/p3l_r5_parallel_ab/shared/A5_inverse_reference.npz')
```

This costs576 derivative RHS at the supplied live paid Point, with no new LU/state. Reacquiring a fresh Point costs one state/eight forward RHS. Two parent-authorized probes therefore reserve1152 derivative RHS plus up to16 forward RHS inside the unchanged A allocation. The optional stored inverse array is about5.3MB before compression. Any resulting Float64 norm/Neumann bound remains a numerical bound unless its rounding enclosure is established separately.
