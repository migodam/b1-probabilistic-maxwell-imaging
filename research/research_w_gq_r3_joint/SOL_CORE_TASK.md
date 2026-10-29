# R3 Sol core implementation contract

All previous science and R2 files are read-only. This worker implemented src/r3/{physics,inverse,charts,accounting,io}.py, tools/run_r3.py and configs/INVERSE_FROZEN.json. Observation generation, private truth and scoring are independent and were not read by core.

## Model

Frozen public wall: x=[-1.5,-.5,.5,1.5], y=[.5,.875,1.25], all six cells active unknown, epsilon real [1,6], ell [0,.25]. z order: epsilon6 y-major, ell6 y-major, target x,y,rx,ry,contrast,loss,angle. Target epsilon=1+contrast+i loss/f; wall epsilon=eps+i ell/f. Three frequencies .85,1,1.15; 4TX16RX, real/imag independent sigma .01. Target bounded x[-.9,.9],y[1.7,2.5], independent radii[.07,.18], contrast[.5,8],loss[0,.25],angle[-pi/2,pi/2]. A circle is nested continuously.

Exact 384-gon ellipse/triangle cut integrals, analytic incident weak forcing, unchanged scalar-TM P3 FEM and PML. All target shape derivatives use the polygon boundary with affine radius/rotation velocities, and all material derivatives include incident forcing. No nodal target occupancy in the PDE. Raster images are actual inferred materials evaluated on an independent pixel grid; no post-hoc clipping or truth selection.

JointModel(config,book).forward(z,jacobian=False) returns Prediction(data (nf,4,16),jacobian (64nf,19),diagnostics). material_maps(z) supplies full material pixel maps. default_config and parameter_spec define standalone public schema. Model initialization and physical actions hold the shared R3/.fem.lock; new heavy jobs after 18:15:21UTC and all physical actions after18:35:21UTC are forbidden. Every LU/RHS is reserved in the fresh R3 ledger before action.

## Inversion

Full-GN baseline: measurement-only wall calibration with target contrast0 (initializer only), 18 public continuous x/y/geometry seeds, best4 spatial-diverse starts selected only by full objective. Each joint start allows120 exact objective evaluations and100 accepted steps. Bound-constrained scaled linearized least-squares subproblems with adaptive trust radius/damping; full nonlinear objective decides all accepted steps. Priors and all optimization options frozen in configs/INVERSE_FROZEN.json before blind scoring. Proper independent truncated-Gaussian wall priors and proper uniform target priors, no finite12 target templates or truth input. Stops distinguish projected gradient, objective tolerance, maxbudget and trust failure.

Artifacts: TRACE.jsonl plus candidate_*.npz preserve every evaluated full objective; STEPS.jsonl preserves all accept/reject decisions and failures; START_*_ENDPOINT.npz preserves all converged/stopped branches. endpoint.npz and ONLINE_SEALED.json select minimum full objective including prior. CLI: Python tools/run_r3.py run --measurement observation.npz --config config.json --output destination --method FULL_GN. Eval guarded_online.py should wrap all blind runs and log file access.

## Forcing chart, not native SOM

Scientific label: SOM-inspired material-forcing chart. Fixed free-reference Green/source operators need not change just because wall material is unknown; this implementation does not claim otherwise. It projects differentiated weak loads q_p=b_p+V_p u_sc, never material coordinates. D is **the diagonal of the consistent P3 FE mass integrated on the fixed public WALL_BOX=(-1.5,1.5,.5,1.25) union TARGET_BOX=(-1.2,1.2,1.4,2.8) using the inherited material basis and its quadrature**, restricted to free DOFs. D_i is positive on DOFs intersecting that public union. D^-1/2 is zero off this support. This is a positive diagonal FE metric, not a row-sum lumping formula (P3 row sums need not be positive).

At public homogeneous wall eps3,ell.06 with zero target contrast, use receiver adjoint to form S_ref=R L_ref^-1 D^1/2. Its thin SVD gives right forcing modes U_r. At each current state, solve L(current) with r loads D^1/2 U_r. Approximate J_p=R L(current)^-1 D^1/2 U_r U_r^* D^-1/2 q_p. Thus current nonlinear material interaction stays in L(current). Primary fixed r16 retains all16 receiver-observable directions of the public reference. Rank8 is retained only as a limited development ablation after a preblind parent amendment. Adaptive checks each proposed step with exact4-source differentiated RHS/frequency and relative defect .1; a failure triggers full current16-adjoint Jacobian, current rank16 receiver chart rebuild (the earlier rank8 ablation code also permits residual-forcing enrichment up to rank16). All actions, failures, and basis events are charged and retained. No physical acceleration or novelty claimed before measured task/cost comparison.

## Necessary derivative audit

results/benchmark/BENCHMARK.json: random scaled full19D analytic directional derivative relative discrepancy1.73e-8, h1e-5. Three frequencies,24,736DOFs each; audit9FE84RHS,5.08CPU sec,~919MBpeak. This verifies derivative implementation at one public point, not global inverse correctness or physical identification.

## Pre-score correction

Independent Eval found a target-loss sampler mismatch before scoring. Core immediately stopped its initial dev run (exit130), retained complete partial outputs and receipt under results/online/dev_r0/INVALID_SAMPLER_FULL_GN, and restarted only after corrected-ready. All spent physics remains in the R3 ledger. This was a contract correction, not scientific tuning after truth inspection.

Final exact full-J stationarity diagnostics are added uniformly to every arm, outside optimizer step selection. They do not change branches or iterations. A chart stop alone is not called converged if its final exact projected gradient exceeds .005. The first chart audit included the inherited full material box instead of the requested wall/target union; this support mismatch was repaired before chart-arm scoring and the attempted audit/cost preserved.

Standalone definitive run aliases also include target_material.npy (y,x,channels=[real epsilon,ell] on public target ROI), predicted_measurement.npz, reconstruction.json, optimization_trace.csv, runtime_cost.json, reconstruction.png and target.png. Every final start saves full complex Jacobian, residual, prior residual/Jacobian and parameter scales using already-paid final stationarity actions. Same-endpoint approximate Jacobian is saved only when already in memory; otherwise the unavailable status is explicit. Reference and adaptive fallback approximate/full matrix pairs are captured from already-paid adjoints for offline Schur/step analysis, with no new PDE for artifact requests.

After the research deadline, the default campaign entry point correctly refuses new physics. Explicit user reconstruction uses tools/deploy_r3.py --deployment-budget-seconds N (1..10800), a new empty output directory and frozen scientific source. It creates output-local deployment_accounting/RUN_CLOCK.json and PHYSICS_LEDGER.jsonl, shares global R3 physical lock, and never extends the original campaign clock or ledger. Interface only was added; no extra deployment PDE was run for validation.
