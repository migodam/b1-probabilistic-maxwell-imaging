# Fixed public mesh and target cut integration

Status: implementation delivered; geometry/mass tests **NOT_RUN** by this worker. Source AST parsing passed. No FEM backend, PDE solve, LU, or physical experiment was executed. Numerical error, runtime, RSS, and scientific qualification remain **NOT_RUN**.

## Interface

- `public_fixed_mesh(h=.17, grid_refinement=2, target_step=.125, pml=None) -> MeshTri`. Imports `wslam.fem_backend` from the campaign's frozen `vendor/wslam` copy. Preserves all public wall grid lines. Adds an axis-aligned target grid covering `[-1.2,1.2] x [1.4,2.8]`. Each grid emits only atomic segments, without a second overlapping long rectangle outline. There are no target circle edges.
- `public_mesh_key(...) -> tuple`, also attached to `mesh.bgga_public_key`. Only public mesh/PML settings enter the key; there is no pose, radius, target contrast, wall occupancy, or observation argument.
- `assembler = TargetAssembler(basis, free_dofs, degree=8, chunk_points=4096, max_points=250_000, max_work_bytes=128_000_000)`. Requires scalar affine-triangle P3. Its one-time triangle bounding-box pass builds a public-ROI spatial index.
- `assembler.assemble(pose, radius, k, chi, incident, source_count=None) -> TargetAssembly`. Here `chi` is the complete complex target contrast (caller supplies frequency-dependent loss), and `incident(points_Nx2)` directly evaluates analytic fields with output `(N,s)` or `(N,)`. Allowed centers are `[-1,1] x [1.6,2.6]`, radii `[.06,.16]`.
- Convenience: `assemble_target(basis, free_dofs, pose, radius, k, chi, incident, assembler=...)`; `public_roi_ports(basis, free_dofs)` supplies the full frozen target-ROI free-space index array for `FixedSpaceCache.U`.

`TargetAssembly` returns `Vtarget` (free x free CSC), `bX` (free x source), `target_area`, actual `touched_global_dofs` / `touched_free_indices`, complete `public_roi_global_dofs` / `public_roi_free_indices`, `touched_rank`, `Vtarget_rank`, `quadrature`, and a metadata dictionary. Every free index addresses the order in the supplied `free_dofs`; it is not a global FEM index. The public ROI support is frozen before any target is supplied.

## Integrals and resource boundaries

The target is the inscribed regular 384-gon, intersected with each candidate FEM triangle. The resulting convex polygon is split into triangles. An independent positive Duffy Gauss rule supplies the retained physical quadrature coordinates/weights and parent-element indices. Reference coordinates and weights are returned in `quadrature.reference_points` / `reference_weights`; reference weights sum to `1/2`. The rule maps `x=u, y=(1-u)v`, with Jacobian `1-u`; degree 8 uses 5 x 5 Gauss points and integrates polynomial degree 8. P3 mass products have degree 6.

Assembly uses sparse point evaluation `P`: `M = P.T @ diag(w) @ P`, `F = P.T @ (w * incident(points))`, and returns `Vtarget = k**2 * chi * M`, `bX = k**2 * chi * F`. The forcing is directly integrated; no nodal incident interpolant or `Vtarget @ incident_at_dofs` substitution is used.

Candidate lookup uses a cKDTree of ROI triangle centers, conservative triangle radii, and local bounding boxes. Pinned scikit-fem `Basis.probes` has an all-elements locator fallback. This module uses its same `mapping.invF` / `elem.gbasis` sparse-probe formula with the parent cells already known from clipping, avoiding that fallback. Point counts and owned-array estimates are checked, and sparse probes/analytic fields are evaluated in bounded chunks. These limits are not a process-RSS certificate; opaque dependency allocations, Python object overhead, the caller's Basis, and prior process memory are excluded.

## Limits and tests

- The polygon area is `384/2 * r**2 * sin(2*pi/384)`, not `pi*r**2`. The analytic area deficit and mesh-partition area difference are retained in metadata.
- `touched_rank` is the exact algebraic mass rank on the actual positive-area intersected P3 support. Restricting a P3 polynomial to any open part of a triangle preserves local unisolvence. `Vtarget_rank` is this count for nonzero contrast, zero for zero contrast. Neither field is a numerical SVD, conditioning measurement, or small-cut stability certificate; numerical rank is explicitly **NOT_RUN**.
- Tests cover public PSLG non-overlap/coverage, independent-rule monomials, polygon area/centroids at extreme allowed poses/radii, mass symmetry/positivity/partition of unity, constant/linear analytic forcing, an independent quartic polygon moment that differs from P3 nodal interpolation, frozen ROI ports, zero contrast, and point-limit rejection. They construct integration meshes/Bases only and never a PDE operator or factorization.
- Worker-edited files: `src/bgga/fixed_mesh.py`, `tests/test_fixed_mesh.py`, and this summary. No old source was changed. Main-thread test execution and accounting remain outstanding.
