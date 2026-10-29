# R2 ordinary FEM adapter

`src/r2/physics.py` imports the unchanged GQ/BGGA/W-SLAM code. Put the new
`src`, old GQ `src`, old GQ `vendor`, and existing continuation `.deps` on the
Python path; the R2 job wrapper supplies those paths. No new solver is defined.

```python
from r2.physics import Eval, WallGrid, Target, S0_TX
wall = WallGrid(x_edges, y_edges, epsinf_yx, ell_yx)
target = Target((0., 2.), radius=.12, contrast=3., loss=.02,
                shape='ellipse', aspect=2/3, orientation=0.)
ev = Eval(book, 'E1', f=1., tx=S0_TX, wall_x=x_edges, wall_y=y_edges,
          wall_refine=1, h=.17, targetstep=.125)
state = ev.evaluate(wall, target, label='declared candidate')
```

Wall axes span exactly `[-1.5,1.5] x [.5,1.25]`. Arrays have shape
`(len(y)-1,len(x)-1)` and are copied into immutable storage. Cells are ambient
`(epsinf,ell)=(1,0)` or occupied with `epsinf in [1.5,6]`, `ell in [0,.25]`.
Loss enters literally as `ell/f`. Every wall axis is a constrained segment;
`wall_refine=2` subdivides each public interval without changing its material.
The mesh depends only on public axes, refinement, h, targetstep and PML.
It uses the inherited P3, physical/PML boxes, line sources and sixteen receivers.
`tx` must equal S0_TX (two sources) or S1_TX (four sources).

Target pose is in `[-1,1] x [1.6,2.6]`, maximum semiaxis in `[.08,.16]`,
contrast in `[0,8]` and loss in `[0,.25]`. Ellipse aspect is exactly `2/3`,
orientation is `0` or `pi/2`. Its affine inscribed 384-gon uses the pinned
triangle/polygon intersection integral; a nodal mask never assembles a target.

`state` exposes `data (ns,16)`, `u (nfree,ns)`, `matrix`, `rhs`, `factor`,
`target` (cut assembly with Vtarget/bX), `wall`, `world`/`target_spec` (Target),
`residual`, `timings`, `metadata` and `derivative_audit`. It does not retain
sampled material fields, global scattered DOFs, currents or adjoint fields.
One exact matching state is cached; invalid/new requests clear that slot.

The assembly and return-operator hooks are:

- `background(wall) -> (L_W,b_W)`, including the analytic incident weak forcing.
- `free_operator() -> (L0,zero_rhs)`; `ev.p._background_matrix` is the same L0.
- `assemble_target(target)` returns Vtarget and analytic cut-integral bX.
- `factor(matrix,label) -> (factor,immutable_matrix,factor_bytes)`.
- `ev.R` is the receiver probe restricted to free DOFs.
- `ev.ass.public_roi_free_indices` contains the common public target ROI ports.
- `ev.p._free_dofs`, `p.tx`, `p.rx`, `p.k`, `p.freq`, `p.basis` retain the pinned meanings.

`factor.solve(rhs,trans='N'|'T'|'H',label=...)` reserves its RHS automatically,
then checks every column against backward tolerance `1e-10`. Callers must not
precharge the same RHS. Factorization separately reserves one LU before
SuperLU. CPU/RSS accounting and process limits remain owned by the outer job
wrapper. The default DOF/factor guards are 80000/4e9 bytes, independently of
the stricter campaign RSS limit.

| Action | LU | Complex RHS |
|---|---:|---:|
| Background/target/differential assembly | 0 | 0 |
| New `evaluate` | 1 | ns |
| Exact state cache hit | 0 | 0 |
| `wall_jvp` or `target_jvp` | 0 | ns |
| `wall_jacobian` or `target_jacobian` | 0 | 16 |
| Joint `jacobians` | 0 | 16 total |
| `objective_gradient` | 0 | ns |

`wall_mass_direction(wall,deps,dell=0)` returns positive
`V'=k² integral(delta_epsilon phi phi)`; the system derivative is `A'=-V'`.
`iter_wall_mass_derivatives` streams epsinf then ell matrices.
`wall_differential_rhs` returns `V'u+b_W'`, including incident forcing.
All wall derivatives hold the occupied support fixed: air columns are zero.

`wall_jacobian(state,selected_cells=None,basis=None,include_loss=True)` supports
flat cell indices, `(y,x)` pairs, or arbitrary real `(ncell,ncolumn)` basis
columns. Rows flatten source then receiver; columns are supplied epsinf
directions then their ell directions. `jacobians` returns `(wall_J,target_J)`
using one shared receiver adjoint. Target order is
`('px','py','contrast','loss','radius')`; `target_loss` is accepted as a name alias.

`target_differentials(state,names=...)` returns a dictionary of source-resolved
`V_X'u+b_X'` arrays without a PDE action. Material derivatives use identical
cut quadrature and remain available at zero contrast/loss. Shape derivatives
reuse the inherited degree-4 polygon-edge rule, with fixed aspect/orientation.
No finite-difference PDE solve is hidden in these methods.

`objective_gradient(state,weights,...)` uses
`delta_phi=Re(vdot(weights,delta_data))`; it returns real `wall`, grid-shaped
`epsinf`/`ell`, `target`, names and audit. `wall_objective_gradient` returns only
the flat wall vector. Whitening, priors, nuisance gain and optimizer are caller
choices. These derivatives concern the same finite discrete equations;
continuum, grid-convergence, finite-difference and scientific gates remain
separate parent-run audits.

`tests/test_physics_algebra.py` forbids FEM and sparse LU, and checks geometry,
immutable coordinates, product-kernel/RHS duality, existing-basis interpolation,
and reservation/backward-error behavior. It is not a physical qualification.
