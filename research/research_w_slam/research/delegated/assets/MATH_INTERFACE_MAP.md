# Math interface map

Paths prefixed C are under `Gaussian/Research_W/research-w-continuation`; paths prefixed R are under `Gaussian/B1/b1-w-p0r-to-p1/research/research_w_p0r_to_p1`.

| Need | Existing interface and lines | Limits / adaptation |
|---|---|---|
| Scalar-TM scattered operator | C/src/wq/fem.py:1–24,137–150 | exp(-i wt), k=2πf, epsilon=epsilon∞+i sigma/f; PML, outgoing line source |
| P2/P3 high-order elements | C/src/wq/fem.py:274–327 | explicit h, quadrature, PML, DOF/memory caps; not automatically qualified |
| Material map | C/src/wq/fem.py:441–454 | callable Nx2→complex N; 1≤Re eps≤9, 0≤f Im eps≤.15 |
| LU reuse | C/src/wq/fem.py:456–501 | exact quadrature coefficient equality only; one live material factor |
| Forward state | C/src/wq/fem.py:536–619 | State holds epsilon, total E, induced current, DOFs, data, residuals, factor and matrix |
| Coefficient JVP | C/src/wq/fem.py:652–692 | delta eps callable; solves L delta u=k² delta eps E with four source RHS; returns source-major receiver derivative |
| Coefficient VJP | C/src/wq/fem.py:694–743 | real-adjoint gradient includes weights; must apply low-dimensional decoder chain rule |
| Full decoder Jacobian | R/src/wr/physics.py:19–37 | volume receiver extraction, transpose T (not H) complex bilinear chain; 12 RHS/frequency |
| Differentiable wall chart | R/src/wr/native_world.py:14–104 | 108 tensor-hat coefficients; restrict via explicit fixed linear map q=q_base+B a for 2–8 parameters; positivity/range must be checked; known support is conditional |
| Simple material/target override | R/src/wr/physical_witness.py:18–53 | fixed target disk; material derivative excludes target interior. Sharp target pose derivative NOT implemented |
| Interface fitted meshes | C/src/wq/meshes.py:20–55; C/src/wq/target_mesh.py:18–32 | triangle PSLG wall/circle interfaces; fit_probe_mesh preserves base facets. Remeshing is new geometry and new LU |
| Source/receiver pose | C/src/wq/fem.py:127–134,295–299,536–542 | coordinates passed at constructor; no source/receiver coordinate derivative interface |
| Packing/whitening | R/src/wr/physics.py:15–17 | frequency-major then source-major complex flatten, concatenate all real then all imag, divide .01, no count normalization |
| Restricted posterior math | R/src/wr/posterior_math.py:87,141,177 | log_marginal_likelihood shared nuisance gamma; finite_grid_posterior, finite_grid_summary pure numerical math; adapt observation_shape explicitly for SLAM |
| Historical two-branch family | R/src/wr/probability_family.py:18–27,145–159,249–284 | fixed two walls, x-offset/contrast quadrature; point() is old-root bound and triggers new physical solves |

For four occupancy states a new decoder must represent 00,01,10,11 including joint multiple scattering in 11. Compute F_h(theta) separately for each declared state. `log_marginal_likelihood` can be reused only after matching the likelihood/observation contract; positive prior quadrature weights must sum to each state's prior mass and posterior normalization must span all states. No cached historical F implies the four-state model.

For E0 local J use explicit coefficient derivatives and the existing JVP/VJP on a fixed discretization. Moving discontinuous interfaces or antenna probes are not ordinary coefficient derivatives; use a declared fitted-remesh finite difference or derive an independently checked shape/pose method. For E1 local information, pack J with sigma consistently and remove/gauge-fix pose modes according to the new contract; solver rank does not by itself establish global identifiability.
