# CPU vector adapter

```python
from vector_backend import TensorDDA, sphere_grid
import numpy as np
points, v = sphere_grid(radius=.12, h=.03)
model = TensorDDA(points, v, k=2*np.pi)
eps = np.broadcast_to(3*np.eye(3), (len(points),3,3)).copy()
state = model.solve(eps, [[0,0,1]], [[1,0,0]], receipt_path='solve.json')
near = state.near([[0,0,.5]])
far = state.far([[0,0,1]])
```

Arrays `state.E/H/j/J/p/E_exc/incident` have shape `(3N,P)` in cell-major,
component-minor rows. `E` is physical material internal E under CM+RR, while
`E_exc` excludes the cell self field. `j=(epsilon-I)E`, `p=v*j`, and normalized
physical polarization current `J=-ik*j`. `state.matrix` and `state.factor`
retain the single all-source CPU complex128 LU problem. Background identity
permittivity is legal; epsilon-I is never inverted.

`model.G` maps density j to internal scattered electric field, including
`(-1/3+i k^3 v/(6pi))*I` self. `model.Goff` maps dipole p to off-diagonal
electric field. `model.external_electric_operator(receivers)` and
`external_magnetic_operator` map density j to external scattered E/H.

`near()` returns E/H and E_scattered/H_scattered, each `(3R,P)`, and excludes
receiver points in occupied voxel interiors. `far()` returns transverse E/H
amplitudes `(3R,P)` with `Esc(r*n)=exp(ikr)/r * F_E(n)+O(r^-2)`.
Near and far contracts are distinct. Near receiver positions need not be on
an artificial sphere; far directions must be unit vectors.

Input epsilon is `(N,3,3)` and must be finite complex symmetric (transpose,
not conjugate transpose). The caller controls SPD bounds/passivity. Plane
wave directions are unit real `(P,3)`, polarizations are complex `(P,3)`
transverse to each direction. Units are vacuum epsilon=mu=Z0=1, omega=k,
time dependence exp(-i omega t).

Source reuse is standalone: `a9_green.py` is the unchanged selected A9 Green
kernel copied into this study. The complete read-only source snapshot is
under `audit/sources/a9_engine_source_snapshot.py`. See implementation audit
and M_QUALIFICATION for approximation and H(curl) limits.
