# Research F source and asset audit

Status: completed source/metadata audit, zero new project physics. Read Research_F package report (all 1804 lines through chunked reads), prompt, protocol, source audit, README and algebra script; inspected original finite_feedback/common.py, e0.py, diagnose.py interfaces, local2d.py MaterialPrior and prior_constrained/prior.py, original A25 dataset generator/rasterizer and background archive keys. This audit makes no scientific novelty/performance judgment.

## Frozen exposed DEVELOPMENT cases

Deterministic eligibility: existing `validation.npz` row order, split=validation, family=smooth, ood_type=ID, amplitude=1.2 (`a1`), noise_db=20, status=OK. First two are:

| Row | Case | Background | Fixed circular-complex sigma |
|---:|---|---|---:|
| 67 | p250073_ring_a1_snr20 | ring | 0.0007590054418520374 |
| 71 | p250073_weak_a1_snr20 | weak | 0.00119211798514519 |

Selection used metadata only, no truth material or outcome/error. These are the SAME smooth parent 250073 in two backgrounds, not two independent objects. Both are exposed DEVELOPMENT conditions, not blind validation. Do not substitute TRAIN anchors 250033_weak/250034_ring: their old starts equal their own known truth and establish calibration only. physics_R_base is a single ROM update, not an independently converged physics reconstruction.

## Conditional metadata provenance and legitimate scope

Original public manifest intentionally omits parent geometry. No local persisted unredacted generator manifest was found. Exact original `src/a25/dataset.py:generate_manifest` was isolated using AST, with only NumPy and constants; train/validation metadata were regenerated preserving split indices and parent numbering. No material rasterizer, label solver, current solve, or locked split was invoked. Parent seed is [251008,250,1,2,0]. This is deterministic metadata replay, not direct reading of a persisted full original row.

Known conditioned parameters for this experiment: family=smooth, ID law, background, amplitude=1.2, global center=[0.012611600067911696,0.07317843756864494], rotation=-2.235014898736308, scale=0.9542157111194681, axis_ratio=1.1343613140435578. These conditions are revealed from generator metadata and MUST be disclosed; this is not general unknown-geometry inversion. Shell/disk/rectangle/subtype/segment and third Gaussian parameters are unused by ID smooth rasterization and can be retained as inert known metadata. Original unknown first-two Gaussian coordinates are placed only in offline_only and never in the online loader.

## Exact decoder interface and prior-only initialization

Reuse `src/lxz/local2d.py:MaterialPrior(row,bg)` and `prior_constrained.prior.decode`. Coordinates in order: (g0_x,g0_y,g1_x,g1_y,g0_width,g1_width,g0_weight,g1_weight). Uniform bounds: centers [-.13,.13], widths [.05,.095], weights [.4,.9]. Mapping is `a=lo+(hi-lo)*ndtr(u)` with independent standard normal u for the original uniform parameter law. `decode(u)` returns (complex chi, updated generator row). Keep actual raster maximum normalization, phase (1+.2i)/sqrt(1.04), and amplitude1.2.

The literal support midpoint u=0 yields two coincident Gaussians (centers0, widths.0725, weights.65): additional singular degeneracy, not a generic rank7 point. Parent therefore froze a deterministic prior-only noncoincident initialization: centers(-.052,0),(.052,0), widths(.0725,.0725), weights(.65,.65). Online JSON replaces ALL eight unknown values before serialization, and redacts original parent_seed to null to prevent accidental regeneration of hidden coordinates. This start uses no target parameters/GT/N3/measurement fitting. Common weight-scale gauge remains; do not assume full coordinate rank or infer uniform seven-dimensional quotient measure. Geometry may remove the gauge; posterior density stays the original eight-coordinate pushforward unless correctly derived otherwise.

`coordinates()` reads the row, so pass ONLY sanitized online generator, never offline original generator. For initialization, calling it on sanitized row is allowed. Analytic derivative includes max normalization and fixed-max-cell assumption; monitor active max and rank along any future path. Finite differences decoder-only are not Maxwell solves.

## Online/offline asset contract

`inputs/online/<case>.npz` has exactly `observation` (32,8) complex and scalar `sigma`; no chi, clean truth prediction, initial_mass, or original unknown parameters. Adjacent JSON supplies sanitized generator and known conditions. Offline evaluator assets are `inputs/offline_only/<case>.npz` (truth, truth_prediction) and JSON (original_generator). Selection preceded offline packaging. No GT errors were computed during selection.

Original acquisition: 24x24 grid, 8 sources,32 receivers, gain1; material mass coordinates `[Re chi,Im chi]/24`. Pack data source-major `[Re32,Im32]` per source; whiten sqrt(2)/sigma once. Background `weak.npz`/`ring.npz` are present under A30 `inputs/a29_release/inputs/a25/background/`, with chi_b, eb,yb,Db,G,S,incident,L_b,volume,points/source_points/receiver_points, and old J/step. G(576,576), S(32,576), incident(576,8), L_b(576,576), old background J(512,1152). Reading archive keys is not a new factorization/solve. Their availability does NOT authorize invoking R0/LU solves in this zero-physics stage.

## TRAIN scale and structure gate

Read only public TRAIN `chi` and background IDs for all512 rows and form physical complex increment RMS=sqrt(mean(abs(chi-chi_b)^2)). Median=0.10551092559851608. Frozen 5% scale=0.005275546279925804. Each material occurs equally over noise variants, so retaining all rows leaves median unchanged from equal multiplicity deduplication. Both amplitudes/backgrounds retained; no DEVELOPMENT material used to set scale. Details: TRAIN_INCREMENT_SCALE.json.

For structures reuse prior.structural_distance (24x24, smoothing scales0/1/2 cells, thresholds .1/.25/.5/.75/.9). Its existing `persistent_thresholds_ge_005` means mean occupancy disagreement >=.05 at a threshold after averaging smoothing scales, and reports how many thresholds meet it; it is not probability mass and does not establish a new topology. Report continuous multiscale RMS, boundary/occupancy differences and component/hole differences separately. Any future requirement such as persistent_thresholds_ge_005>=2 must be frozen by parent before physics, not selected after results. Smooth local displacement may be meaningful without a topology change.

## Missing caches and evidence boundaries

Local finite_feedback inputs contain N3, physics_R_base and TRAIN anchors, but no results/*JACOBIAN.npz, *CHART.npz, results/endpoints/*.npz. Old JSON metrics and old source cannot supply missing curvature vectors. Release retrieval is needed to execute the original full-cache three-term audit; never reconstruct missing arrays through LU/RHS and call it zero physics. New independent anchors are not present yet. No claim of independent converged reconstruction is possible from this audit alone.

Copied package algebra_checks.py to evidence/package_checks and ran there with existing Gaussian/.venv_nn. Assertions PASS. It uses synthetic six-dimensional/two-cell matrix algebra only; those solves are not project Maxwell states or RHS. The finite-resolvent relative check was 1.235333123192667e-15 in this runtime; minor floating differences from package recorded result are expected. No network training, Mie, mesh physics or historical trajectory replay ran. Original package and closed ledgers untouched. New project full states=0, complex RHS=0, Mie=0,GPUtraining=0. Decoder evaluations in this worker=0; only source/metadata/array reads and TRAIN array arithmetic. No SHA256 computations or checks.
