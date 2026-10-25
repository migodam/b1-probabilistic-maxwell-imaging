# Existing physics and cache inventory

Read-only inspection; no PDE/LU/RHS, network, hash checks or old-input edits. NPZ inspection read ZIP metadata and NPY headers only. One header-reader compatibility error was repaired using NumPy's public read_array_header_1_0 API.

## Recommendation

Use the pinned continuation `src/wq/fem.py` as the small numerical engine, with a NEW campaign-owned coefficient/derivative model, ledger and output root. Its callable coefficient accepts arbitrary low-dimensional (2–8 wall-material parameters plus target parameters) maps. There is no ready-made SLAM decoder, analytic receiver-position derivative, or sharp-target shape derivative. Those remain new implementation and verification work. Enumerating four occupancy states is an outer model loop; this is not provided by the two-branch historical family.

`FEM.solve` gives total data, quadrature fields/current and exact live SuperLU factor. `jvp` and `vjp` reuse that live factor. `wr.physics.coefficient_jacobian` is an efficient full parameter-Jacobian pattern: 12 receiver-transpose RHS per frequency rather than four source RHS per parameter. It is **volume-output-specific** and applies the decoder chain rule. Copy/adapt the pure routine into the NEW campaign, rather than importing `wr.physics.Forward`: that class is bound to its old ROOT and writes there.

## Exact historical sources and usable files

Git object `958a478739d9553269316ded8125fda9df225c0c:research/research_w_p0r_to_p1` is available in local `Gaussian/Research_W/publish-w-p0/.git`, although that stage is absent from its working-tree HEAD. Use `git show` for exact historical source text. The stage and caches physically exist at `Gaussian/B1/b1-w-p0r-to-p1/research/research_w_p0r_to_p1`. This inventory did not assert byte identity of that working copy to the pinned commit.

The continuation source is at `Gaussian/Research_W/research-w-continuation` and is also represented by the pinned `483d260c6b5b111db9fb814b096685fccb362448` tree. Source lines below refer to the inspected local files; retrieve pinned blobs for the new source snapshot.

Exact saved joint witness: `results/B2/validation/joint_seed2_narrow_fullj/RESULTS.json` reports d_joint=0.7987161019689071, d_cal=3.5605119112114565, target_shift=.1, numerical_pass=true, cal_physics_threshold=false, joint_physics_threshold=true, strict_original_prior_GO=false, posterior_probability=null. It is conditional fixed-support evidence, not calibration equivalence or a broad-prior posterior statement.

Paired means are already stored in `joint8_PAIR.npz` and `cal_PAIR.npz`: coarse0/coarse1/fine0/fine1 each complex128 (5,4,12). Full endpoint fields exist: sample `fields/joint8_world0_fine/endpoint_0000/f0.npz` has 194656 material quadrature points, field/current (194656,4), scattered DOFs (81685,4), mesh vertices (2,9155) and triangles (3,18074). No LU factor is serialized in these files. Saved fields cannot provide free new-material derivatives without reconstructing/factoring the exact operator.

`inputs/W_FAMILY_WALLS.npz` contains q0/q1 float64 (108,). `results/P1A/rules/full_5.npz` contains F (50,5,4,12), branch/u/contrast/x/weights (50,). These are an existing two-wall/one-target restricted family, not a 2–8 parameter wall posterior and not four occupancy states. Use only as explicitly declared historical checks, not as new experiment data.

Continuation caches include `results/wall-qualification-resume/WORLD_2026102202_h0.145.npz`: nine frequency cal/task (4,12) means plus single saved cal/task fields (128176,4), DOFs (57985,4), geometry and epsilon. The 22-byte `results/fem-derivatives/ARRAYS.npz` is empty; derivative numerical evidence is in RESULTS.json (adjoint defect 3.15e-14; best FD relative error 1.49e-9), not cached derivative fields.

## Runtime and qualification boundary

Existing interpreter: `Gaussian/.venv_nn/bin/python` (symlink exists and successfully read NPZ headers). Existing FEM packages: continuation `.deps/skfem` and `.deps/triangle`. Minimal paths: NEW campaign src, pinned continuation src, continuation vendor (contains p0), continuation .deps. Set all BLAS thread variables to 1 and PYTHONDONTWRITEBYTECODE=1. Recorded environment is Python3.13, NumPy2.5.3, SciPy1.18.1, skfem12.0.2, Triangle20250106; imports/versions were not rerun.

Do not invoke old tools/run_job.py: it writes old ledgers and consumes the closed old budget. NEW wrapper must charge child CPU once and reserve every LU/RHS before action.

The historical Q result is PASS / QUALIFIED_DISCRETE_FULLWAVE only, for the eight frozen data contexts, not continuum certification or arbitrary material/geometry/pose paths. New candidate geometries, four occupancy combinations, and pose derivatives need their own substantive E0 checks within the new authorization.
