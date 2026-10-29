# R2 FEM adapter implementation receipt

Status: IMPLEMENTED; zero-PDE algebra validation PASS. Scientific qualification
and full derivative/grid audits are parent-owned and are not asserted here.

## Scope and ownership

This worker implemented the ordinary R2 adapter, its minimal API documentation,
and necessary algebra checks. It used the old
`research/research_w_gq/vendor/{wslam,bgga}` and GQ physics as read-only imports.
It did not edit old evidence, accounting, campaign configuration, protocols,
or publication; did not compute hashes, install dependencies, launch DeepSeek,
run a scientific PDE job, or change the default provider/model.

New files:

- `src/r2/physics.py`
- `docs/PHYSICS_API.md`
- `tests/test_physics_algebra.py`
- this receipt under `research/delegated/backend/`

The parent froze the contract and ran its finest-grid benchmark. That benchmark
is not this worker's evidence. This worker made no `evaluate` change after the
parent reported the benchmark complete.

## Implemented behavior

1. `WallGrid` stores arbitrary strictly increasing public x/y boundaries over
   the exact wall box, with immutable `(ny,nx)` epsinf/ell arrays. Ambient cells
   are exactly `(1,0)`; occupied cells require epsinf `[1.5,6]`, ell `[0,.25]`.
   Complex permittivity uses literal `epsinf+i*ell/f`, including endpoints.
2. `public_mesh` adapts the pinned Triangle PSLG recipe. All public wall axes
   are atomic constrained segments; optional refinement 2 subdivides each
   interval. A runtime triangle/interface crossing check rejects a misplaced
   public interface. Mesh input contains no material value, target position,
   target size, observation, occupancy, or truth-dependent quantity. The fixed
   target ROI, material/PML boxes, P3 basis, line sources and RX16 are inherited.
3. `Target` separates pose, maximum semiaxis, contrast, loss, disk/ellipse,
   aspect and orientation. R2 ellipse aspect is `2/3`; orientations are exact
   `0` or `pi/2`. Affine cut-polygon integration reuses the pinned target mass
   and incident forcing assembly. No nodal occupancy mask assembles a target.
4. `Eval.background(wall)` returns the actual `L_W,b_W`, including weak
   analytic-incident forcing. Full evaluation solves
   `(L_W-V_X)u=b_W+b_X`. Public free/target/factor hooks support the return
   worker: `free_operator`, `assemble_target`, `factor`, R, pinned free DOFs,
   and public ROI ports.
5. Every factorization calls `book.reserve(...factorizations=1,rhs=0)` before
   SuperLU. Every factor `solve`, including externally requested N/T/H actions,
   reserves its own column count before the numerical action. Callers must not
   precharge those RHS. All columns pass inherited backward checks at `1e-10`.
6. States retain free scattered DOFs, receiver data, sparse operator/forcing,
   cut assembly and metadata. They do not retain global scattered DOFs,
   material fields, adjoint fields, currents or target field copies. Wall
   product kernels interpolate the existing P3 basis in bounded element chunks.
7. Wall mass direction matrices, streamed arbitrary-basis/cell matrices and
   source-resolved wall differential RHS are exposed. The differential RHS is
   `V'u+b_W'`, including the forcing derivative. Support is held fixed, so
   ambient columns are zero.
8. Wall and target Jacobians use a shared sixteen-receiver adjoint. Joint
   `jacobians` computes both with sixteen RHS total. Cell selection supports
   flat indices, `(y,x)` pairs, or arbitrary real basis columns. Source-major
   rows preserve all source and receiver dimensions. Scalar objective gradients
   use ns weighted-adjoint RHS and the convention
   `delta_phi=Re(vdot(weights,delta_data))`.
9. Target differential RHS are assembled without PDE solves. Material terms
   use identical cut quadrature, including at zero target chi. Shape terms use
   the same 384-gon boundary with inherited degree-4 Gaussian edge integration.
   Target names are `(px,py,contrast,loss,radius)`; `target_loss` is a name alias.
   Aspect and orientation remain fixed hypotheses rather than derivative axes.

## Validation performed

Authorized wrapper command:

```text
<existing A25 .venv_a25/bin/python> tools/run_job.py E1 backend-zero-algebra \
    <existing A25 .venv_a25/bin/python> tests/test_physics_algebra.py
```

The wrapper supplies new/old src, read-only vendor, existing scikit-fem/Triangle
dependencies, single-thread environment and campaign CPU/RSS monitoring.
System/Homebrew Python import probes initially failed because their default
paths have no NumPy; no installation or environment mutation followed. The
parent supplied the established A25 environment, used successfully above.

Six tests PASS:

- actual public-cell preservation, immutable ownership, legal material domain,
  and exact loss/f values at all declared frequencies;
- every original/refined wall segment present in the target-independent PSLG;
- affine target polygon area partition for both ellipse orientations, with all
  supplied mesh vertices outside the target;
- product-kernel Jacobian equals an independently formed FE-quadrature
  `q^H(V'u+b')` for multiple complex sources/receivers, ambient cells and signed
  arbitrary basis columns;
- bounded chunk interpolation equals the existing P3 `Basis.interpolate`,
  including zero outer-boundary values;
- factor wrapper reserves before solve for N/T/H, tracks column counts and
  rejects a defect confined to a tiny RHS that an aggregate norm could conceal.

FEM construction and sparse LU are patched to raise in this suite. The tiny
P3 Basis and explicit algebra fixture do not assemble or solve Helmholtz.
The fake factor uses an identity action and a deliberately wrong action only.

Recorded job receipt: `results/jobs/backend-zero-algebra/RECEIPT.json`.

| Quantity | Recorded value |
|---|---:|
| Exit status | 0 |
| Tests | 6 PASS |
| Physical LU | 0 |
| Physical RHS | 0 |
| Wrapper process CPU seconds | 0.336211 |
| Wrapper wall seconds | 0.5156872501 |
| Polled peak RSS bytes | 58163200 |
| Memory killed | false |

## Remaining parent audits and limitations

No worker scientific result is produced. Full-PDE primal/adjoint duality,
shape finite differences, independent-grid checks, PML/cut-forcing convergence,
optimization correctness, rank interpretation and scientific claims remain
NOT_RUN by this worker. Polygon/FE derivatives concern these finite discrete
equations and the inherited boundary integration rule; they are not a
continuum or uniformly certified derivative claim. There is no nuisance gain,
prior, whitening, optimization step or scientific threshold in the adapter.

Repeated calls to separate wall/target Jacobian methods each legitimately pay
sixteen RHS. Use joint `jacobians` when both are required. Cached evaluation
returns the identical immutable state with zero new LU/RHS; an invalid/new
request clears the one-slot cache. CPU and RSS are caller-owned, and the
campaign's 4-GiB RSS limit is independent of the backend's configurable 80000
DOF / 4e9 factor-byte guards.
