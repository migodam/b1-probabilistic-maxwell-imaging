# Campaign FEM adapter

`src/wslam/fem_backend.py` was copied directly using local git show from `483d260c6b5b111db9fb814b096685fccb362448:research/research_w_continuation/src/wq/fem.py`.

Only physical contract adaptations: material box is (-2,2,.25,3), default PML physical box is (-2.75,2.75,-.8,3.5), passive conductivity upper bound is .25. Mesh conformity cuts and automatically constructed Cartesian axes use MATERIAL_BOX including .25. Thickness1, strength4, power3 remain defaults; caller may provide PMLConfig values. All inherited Helmholtz/PML formulae, complex128 assembly, residual checks, LU/RHS pre-reservations, memory caps, exact factor reuse and no-fallback behavior are retained. Default residual tolerance is1e-9. New geometry/material scope is **not** covered automatically by historical qualification.

## API

- FEM(freq,tx,rx,book,*,h,order=2,quadrature_order=None,pml=None,...,receiver_mode='probe'); use receiver_mode='volume' for efficient receiver-transpose Jacobians. Caller supplies a NEW ledger and mesh, not an old Book or launcher.
- solve(epsilon) -> State, with source-major data and live factor. epsilon(points Nx2) -> complex N.
- physics.acquisition(outer=False) -> TX(-1.5,0),(1.5,0), optionally appended(-2.25,0),(2.25,0);16 RX evenly spread x[-2,2],y0.
- physics.pack(data,sigma=.01) flattens frequency/source/receiver C-order, concatenates real then imaginary, divides sigma. packj(complex_matrix,sigma=.01) stacks real/imag rows.
- coefficient_jacobian(p,state,derivative) -> complex(nsource*nreceiver,nparameter), exact volume-output chain rule,16 transpose RHS for declared RX. derivative(points,frequency) returns complex(nquad,nparameter). No full fields exported.
- material_jvp_jacobian(p,state,derivative) chooses min(nparam*nsource,nreceiver) RHS for volume receivers (ties choose receiver transpose); probe receivers use nparam source-JVPs. No new factorization. Every RHS reserved before solve. This is a coefficient derivative on fixed quadrature, not a shape/antenna derivative.
- apply_target_mass(p,wall_state,full_state,total_field=None) assembles free-DOF forcing k² M_delta E; no physical solve.
- paired_feedback_residual verifies Lwall (us_full-us_wall)=k² M_delta Efull on the paired identical backend/source set; returns residuals/pass rather than silently accepting. No extra solve.
- target_feedback_terms returns exact first-interaction and repeated-feedback RHS split; algebra only, no standalone ROM or four-block scientific claim.
- paired_data_decomposition returns Fwall,delta,Ffull from already-paid paired states.

Both states for feedback helpers must be solved on **one identical backend/mesh**; retained state.factor and state._matrix keep their own operator after a subsequent material solve. No serialized LU reconstruction, approximate material-key cache, truth loader, geometry builder, run loop or old-root output writer is included. Caller owns frequency stacking and output serialization.

Validation performed: AST parsing of both files and direct inspection of changed constants/cuts. No imports, PDE, LU, RHS, physical tests, networking or hash checks were executed. Parent must perform authorized E0 checks before treating new-domain derivatives or feedback identities as qualified.

## Point-output extension

The main contract now uses receiver_mode='probe'; volume extraction remains an independent check. `point_receiver_adjoint(p,s)` reserves nreceiver transpose RHS and returns (global-DOF×receiver complex fields, maximum relative residual). It uses the bilinear transpose T, checks against the backend residual_tolerance (default1e-9), and emits the paid solve event.

`point_coefficient_jacobian(p,s,derivative,qnodes=None)` returns complex source-major receiver×parameter rows. Pass qnodes from one paid point_receiver_adjoint to reuse the same fields across all wall/material/shape parameters; omission pays for a receiver transpose set. This leaves the earlier volume helper unchanged.

`surface_shape_jacobian(p,s,qnodes,boundary_points,normals,ds,radius_weights=None,*,contrast)` returns columns x translation,y translation,then supplied normal-velocity radius columns. contrast is epsilon_inside-epsilon_outside scalar or boundary vector. E is analytic incident line Green times actual state amplitudes plus sampled scattered FE field; q is sampled receiver transpose field. The expression is k² sum(ds*contrast*E*q*normal_velocity), with no complex conjugation. No additional PDE actions. The caller owns scalar-TM shape-formula validity, polygon-consistent boundary quadrature and substantive FD qualification; this implementation is not a geometry-independent proof. Both imports and physical execution remain NOT_RUN; AST parsing passed.
