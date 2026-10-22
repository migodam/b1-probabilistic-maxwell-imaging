# Independent scalar-TM FEM/PML route

Implementation: `src/wq/fem.py`. This route imports no P0 VIE code, integrated Green matrix, or cell averaging routine. It requires scikit-fem in the project dependency path and the existing NumPy/SciPy runtime.

## Equation and source convention

The time convention is `exp(-i omega t)`. At normalized frequency `f`, `k=2*pi*f`, and passive material is `epsilon=epsilon_infinity+i*sigma/f`. The unit line-source incident field is `i*H0^(1)(k*|x-s|)/4`, including its phase and absolute normalization. For plane waves it is `exp(i*k*d.(x-origin))`, evaluated exactly at each material quadrature point.

Writing `u=u_inc+u_s` gives `(-Delta-k^2*epsilon)u_s=k^2*(epsilon-1)*u_inc`. Only `u_s` is a finite-element unknown. The line source is outside the material region; the forcing therefore has no line-source singularity. The free-space incident solution is never approximated by a mesh delta or a center-sampled source.

The material support box is `[-2,2] x [0,3]`. The coefficient callback receives an `(n,2)` array of material triangle quadrature points and returns a complex vector `(n,)` or one scalar. The coefficient is set exactly to one outside that box. The same frozen physical bounds as P0 are enforced: `1<=Re(epsilon)<=9` and `0<=f*Im(epsilon)<=0.15`.

## PML and mesh

For each coordinate the stretch is `s=1+i*alpha*(distance/thickness)^power` inside its PML, and one inside the physical rectangle. The transformed weak operator uses `A=diag(s_y/s_x,s_x/s_y)` and `B=s_x*s_y`: `integral(A grad(u_s).grad(v)-k^2*B*epsilon*u_s*v)`. The outer computational rectangle has homogeneous Dirichlet data. With outgoing `exp(+i*k*r)`, the positive imaginary coordinate stretch supplies damping.

The default physical rectangle is `[-2.5,2.5] x [-0.8,3.5]`, the default thickness is one on each side, the maximum imaginary stretch is four, and the polynomial power is three. These settings are an exposed starting point, not an error certificate. `PMLConfig` accepts separate left/right/bottom/top thicknesses and strengths; changing it requires a new backend and new counted solve.

The Cartesian tensor triangulation contains the material-box edges and all PML interfaces exactly. Each Cartesian interval is at most requested `h`; triangular diameter is at most `sqrt(2)*h`. P2 and P3 continuous Lagrange elements are available. With the uniform mesh, curved walls, cavities, and disk boundaries remain pointwise coefficient discontinuities inside triangles. Their resolution must be checked by decreasing `h` and increasing triangle quadrature order. No area averaging or reuse of the VIE square grid is involved. The default integration order is `max(8,2*order+2)`.

The constructor additionally accepts `mesh=MeshTri(...)`. A supplied mesh must fill exactly the configured outer rectangle, have no interior boundary holes, and conform to the four edges of the rectangular material box. These checks run before `Basis` allocation. Internal material interfaces can be supplied through that mesh by the parent's separate mesher; the backend still evaluates the original coefficient callback at actual quadrature points. For a custom mesh, `h` is a label and the actual maximum triangular edge length is reported separately.

## API

```python
from wq.fem import FEM, PMLConfig

p = FEM(freq, tx, rx, book, h=0.1, order=2,
        quadrature_order=8, pml=PMLConfig(), phase='validation')
s = p.solve(epsilon_callable)                  # all constructor tx line sources
t = p.solve(epsilon_callable, angles=angles,   # all requested plane waves
            plane_wave_origin=(0.0, 1.5))
F = p.far_field(s, observation_angles)          # postprocess only; no PDE solve
values = p.evaluate(s, observation_points)      # source-major total fields
```

This example is illustrative and has not been executed by the implementation worker. The parent owns physical jobs and the external budget wrapper.

Receiver arrays `data`, `scattered`, `direct`, and `incident_at_receivers` have shape `(source,receiver)`. Line-source `data=direct+scattered`. For plane-wave cylinder qualification, `direct` is exactly zero and `data=scattered`; the nonzero analytic incident receiver field is available as `incident_at_receivers`. `evaluate(...,total=True)` returns the actual total field for either source mode.

State `epsilon`, `incident`, `field` (total), and `current=(epsilon-1)*field` are sampled at `material_points`; source-bearing arrays have shape `(material_quadrature_point,source)`. `material_weights` supply the triangle integration weights in the same flattened order. `scattered_dofs` are `(global_FEM_dof,source)` and include exactly zero outer boundary values. No full Jacobian matrix is formed.

`jvp(state,delta_epsilon_callable)` solves `L*delta_u=k^2*delta_epsilon*u_total` with the state's existing LU and returns source-major complex receiver data. `vjp(state,data_cotangent)` solves the conjugate-transpose system and returns a complex gradient `g` at material quadrature points, including integration weights. Its convention is `Re(vdot(y,J*delta_epsilon))=Re(vdot(g,delta_epsilon(points)))`. Thus real/imaginary sample-coordinate gradients are `g.real` and `g.imag`; a material parameterization still needs its own chain rule. Both actions reserve exactly one RHS per source before forming or solving that group, and retain their state's original matrix/factor even if a later material has replaced the backend cache.

The far-field coefficient is `F(theta)=k^2*integral_D chi*u_total(x)*exp(-i*k*d(theta).x) dx`. Its asymptotic convention is `u_s(r*d) ~ exp(i*(k*r+pi/4))/sqrt(8*pi*k*r)*F(theta)`. It uses the finite-element total field and material quadrature, so it is postprocessing rather than an independent scattering reference.

## Accounting and failure behavior

Every new material/frequency sparse LU first calls `book.reserve(phase,'fem_material_frequency_LU',rhs=0,factorizations=1,...)`. Every group of incident RHS first calls `book.reserve(phase,'fem_forward_sources',rhs=n_sources,factorizations=0,...)`, before assembly or numerical solution of that group. Repeating a coefficient with exactly identical material quadrature values reuses the same LU; it still reserves all new RHS. There is one cached material factor, with no tolerance, rounding, or implicit equivalence criterion.

The required `book` is external and must expose `reserve`; optional `event` receives timings with `charged_cpu_seconds=0`. There are no local accounting caps, CPU charges, subprocess jobs, convergence loops, or hidden budget. CPU charging remains entirely in the parent's wrapper.

The solver uses complex sparse SuperLU without diagonal jitter, pseudoinverse, least-squares fallback, coefficient clipping, artificial loss, or mesh-dependent source normalization. LU failure propagates. Nonfinite coefficients, matrix entries, RHS, solutions, or excessive RHS-relative linear residual stop the solve. The reported residual is the maximum over sources of `||A*u-b||_2/max(||b||_2,tiny)`. It diagnoses only the algebraic solve, not discretization or radiation error.

Before constructing `Basis`, exact scalar P2/P3 DOF counts are compared with `max_dofs` (default 60000). An exposed conservative 2D fill heuristic estimates factor storage as `64*N*ceil(log2(N))*24` bytes and rejects values exceeding `max_factor_bytes` (default 1.5 billion). Actual CSC matrix, SuperLU L/U/permutation bytes, and factor nonzero counts are reported after factorization; an actual storage overrun also stops the job. The estimate and measured sparse storage are not rigorous process-RSS bounds: opaque SuperLU workspaces and live Python/basis arrays also use memory, so the parent's process monitor remains necessary.

## Qualification required from the parent

Worker status: **NO PHYSICAL SOLVES RUN**. Import and compilation passed using scikit-fem 12.0.2 from project `.deps`; no FEM constructor was executed by the worker. No FEM error, derivative identity, or PML convergence has been established by implementation.

The parent should retain separate comparisons for (a) mesh/order refinement with fixed PML, (b) quadrature refinement for internal discontinuities, and (c) PML thickness/strength/physical-box variation at an already resolved mesh. A change in all three at once does not isolate their errors. The independent cylindrical plane-wave expansion can test sign, amplitude, material integration, and full receiver scattering data. A line-source comparison/reciprocity test then covers the actual four-source geometry. Vacuum scattering is an exact-zero-forcing case and cannot on its own validate PML performance or source normalization.

Comparison values, tolerances, selected settings, failures, and consumed physical budget belong in the parent's qualification artifacts. This module never marks itself converged.

API references checked during implementation: [scikit-fem official API](https://scikit-fem.readthedocs.io/en/stable/api.html) and [official assembly guide](https://scikit-fem.readthedocs.io/en/stable/advanced.html).
