# R12 inherited evidence and callable interfaces

Read-only inventory, 2026-11-11. No inversion, forward evaluation, GT access, hashes, publication, or process control. `INVENTORY.json` records exact absolute paths, array keys/shapes, prediction metadata, task and frozen configuration copies. Sources actually inspected: R11 REPRODUCE, endpoint-path/face-curvature scripts and summaries; R10 adapter/runner; R8 warm-network/features, Fresnel adapter/features/solver, vector classes; R1 TensorDDA and R7 vector inverse. This is targeted inspection, not a claim of reading all R4–R11 or the Atlas.

Use bases below (all exist):

- `R12=/Volumes/migodam's-external-brain/Research/Inv_SLAM/research/delegated/max_struct_r1/r12_worktree`
- `R10=/Volumes/migodam's-external-brain/Research/Inv_SLAM/research/delegated/max_struct_r1/r10_worktree/research/research_max_struct_r10`
- `R8=/Volumes/migodam's-external-brain/Research/Inv_SLAM/research/delegated/max_struct_r1/r8_worktree/research/research_max_struct_r8`
- `PY=/Volumes/migodam's-external-brain/Research/Inv_SLAM/Gaussian/.venv_nn/bin/python`

## Exact three inherited wrong endpoints

| Alias | Case | Original saved endpoint relative to R10 | Saved original Phi | Relative data residual |
|---|---|---|---:|---:|
| ring CIE0001 | r10_twod_0001 | RESULTS/twod/r10_twod_0001__cie_single__noise_ball_99/prediction.npz | .009385376006836756 | .13684356879867443 |
| shell G190004 | r10_twod_0004 | RESULTS/twod/r10_twod_0004__matrix_green19__noise_ball_99/prediction.npz | .17393421915111248 | .5897984853416794 |
| bridge G190022 | r10_twod_0022 | RESULTS/twod/r10_twod_0022__matrix_green19__noise_ball_99/prediction.npz | .014208880543995337 | .16849809526410067 |

Each original NPZ contains `m(576 float64), j/Q(576,12 complex128), Y_prediction(48,12 complex128), gradient(576 float64)`; sibling `prediction.json`, `PUBLIC_TASK.json`, `FROZEN_CONFIG.json`, `NUMERICAL_SCHEDULE.json`, receipts retain provenance. The tasks record remote Windows original paths; use the verified local R10 paths above, not those remote strings.

Already portable inside R12: `research/research_max_struct_r11/RESULTS/exposed_face_curvature/{case_id}.npz`, key `m`, is the wrong endpoint; `RESULTS/exposed_endpoint_paths/{case_id}.npz`, key `m[0]`, is also the wrong endpoint. That path array has nine sampled materials and **must not be treated as an optimizer trajectory**. The path's `m[-1]` is a separately paid alternative endpoint. Exact mapping is independently recorded in R11 `RESULTS/exposed_endpoint_paths/SUMMARY.json` and its creating script. Endpoint curvature is positive only for ring/shell; bridge has negative free-face eigenvalues and non-negligible PG. Do not label all three equally converged minima.

## Acquisition, noise, prior, initialization

Load `R10/data_public/twod/operators.npz` then overlay `R10/data_public/twod/{case_id}.npz`. Common operators: `G(576,576), E0(576,12), S(48,576)` complex128; `centers/coords(576,2), sources(12,2), receivers(48,2), grid_shape(2), k(), radius()` and aliases are retained. Public observations are `Y(48,12)` with scalar sigma and noise radius, bounds, grid shape; no material truth appears there.

Noise is inherited raw iid circular complex with E|N|²=sigma², no whitening/source normalization. Sigma is `.0004179219856601356`, `.00022659487662092322`, `.000377920123373535` respectively; all public noise_radius values are `25.16542219663291`. Do not replace sigma with this radius. Noise99 adapter implements one full 99% noise-ball event and supported receiver-prefix selection. Frozen DPCS ranks are 13,13,12.

The original starts were uniform m=.1, real bounds [.02,20]. Frozen per-task continuation is mf12/maxiter1500, corrector maxiter3000. The class `research.research_fbg_r3_experiments.src.physics.MaxwellModel(public, shape=None)` reuses the R4 disk-pulse operator. Useful calls: `state(m)` returns `Y,j,Q,L,lu`; `forward(m)`, `data_value_gradient(m,Y,beta=None,reference=None)`, `physical_score(m,Y)`, `tangent(m,h)` (returns current tangents, not measurement tangents). Its exact original objective is `.5*||Ypred-Y||²/||Y||²+R(m)`. `model.reg` is `research_z_alg_r4_diversity/src/physics.py:Regularizer(grid_shape,smooth=1e-4,tv=1e-4,epsilon=.01)`: mean squared neighbor differences plus smoothed TV; retain its exact curvature when needed.

Historical reconstruction command form, **not run here**: `$PY research/research_max_struct_r10/experiments/twod_run.py --manifest <portable manifest> --task-id <one of above task IDs> --config <matching frozen config> --output-dir <NEW directory> --r8-repo "$R12" --dpcs-rule noise_ball_99`. Existing output paths are immutable. G19 material hook is `research/research_max_struct_r10/src/lowrank_denominator.py`, prepare `prepare_green_strong`, callable `matrix_material_readout`, label matrix_green19; its public task arm is a_plus. Prefer loading the endpoint over rerunning history.

## Frozen R8 neural initializer

Weights are present **inside R12** at `research/research_max_struct_r8/RESULTS/training/{physics,born,domain_augmented}/weights.pt`. `src/warm_network.py:load_predictor(path,device='cpu') -> (predict,checkpoint)`, `predict(features)` accepts `(10,H,W)` or `(B,10,H,W)` and returns `(B,H,W)`. Checkpoint defines architecture, state_dict, mean/std and feature_mode. TinyWarmStart has 24,737 parameters, 10 channels, width16, and bounded output `.02+19.98*sigmoid`. Physics training receipt says 1536 training/192 validation objects, 120 epochs, best epoch86, CPU, no test objects read. This is an initialization prior, not an unsupervised constitutive update.

`src/warm_features.py:FeatureOperator(public)` followed by `op(Y) -> features, {born,dpcs}, cost`. Ten frozen channels use Born, DPCS readout, current/field energy, source aperture/count and coordinates. This API assumes one shared S and must not be applied blindly to a masked Fresnel acquisition. Missing local R12 training arrays are available at `R8/data_training/train.npz` and `validation.npz`; manifests are present in R12. Train script expects keys X,m; no labels were loaded here. Existing frozen weights require no retraining.

## Fresnel masked acquisition

Raw original measured files are **absent in R12** but verified present at `R8/data_local/fresnel/{dielTM_dec8f,twodielTM_8f,FoamTwinDielTM}.exp`. Do not publish them as derived arrays. `research/research_max_struct_r8/src/fresnel_adapter.py:prepare(path,frequency=2.,side=24,width=.16,illumination='multipole',matrix_free=False)` returns model/public/Y/I/Y_train/train_sources/train_mask/withheld_mask/held_sources/centers/calibration/incident_evaluator/setup cost. Source fit uses empty-target fields; scattered heldout channels are not fitted.

`SourceReceiverModel(public,source_indices,receiver_indices)` has `receive(J)`, `receive_adjoint(R)`, `state(m)`, `data_value_gradient(m,Y,beta=None,reference=None)`; Ss shape is `(gathered receivers,n,selected sources)`. This is the exact source-specific mask; `model.S@J` is not the correct measured forward. Matrix-free counterpart has the same receive calls and bounded GMRES then retained dense-LU fallback. `src/fresnel_features.py:features(model,Y,centers,halfwidth=.08)` supplies the frozen 10-channel masked NN transform; requires dense model.G. `src/fresnel_solver.py:run(model,Y,m0,arm,config)` is the existing baseline interface. Foam historical entry uses frequency8, side60,width.2,matrix_free=True. Reusing R8's synthetic NN on Fresnel is explicitly out of distribution; no pixel GT.

## 3D TensorDDA and baseline APIs

`research/research_max_struct_r1/src/vector_backend.py:TensorDDA(points,volume,k=2*pi)` constructs full dyadic G with self `-1/3+i*k³*V/(6*pi)` and `G=V*Goff+self*I`. `external_electric_operator(receivers)` returns `(3Nr,3N)`. `solve(epsilon,directions,polarizations,receipt_path=None)` requires epsilon `(N,3,3)` complex symmetric, returns VectorState carrying E,H,j,p,inc,L,factor,receipt; forward epsilon is absolute relative epsilon, not chi. One LU shares all source RHS. `plane_waves(points,k,directions,polarizations)` supplies incident fields.

Reusable small inverse: `research/research_max_struct_r7/src/vector_inverse.py:VectorInverse(nside=5,h=.08)`: N125 real isotropic chi, six plane-wave sources, 30 vector receivers, bounds [0,4], prior lam2e-4/l2=1e-6; `physical(x,Y=None,gradient=True)` returns `(pred,J,E,lin)` without data, or `(Phi,g)` with gradient, or `(Phi,data_residual_squared,J,E,lin)` without gradient. `optimize(Y,kind,outer=15,inner=12,corrector=40)` supports direct_physical, ordinary_joint, dual_like_joint, A_material_joint, includes common corrector and paid selection; starts .25.

Larger R8 adaptation: `src/vector_cie3d.py:VectorCIE3D(ns=6,width=.6,lossy=False)` builds 60 dipole sources on three circles and 60 vector receivers, `physical(z,Y=None,gradient=True)`, `prepare(Y)`, `material(z,J,kind,beta)`, `current(z,Y,kind,beta,mf)`. It uses per-source data normalization and a different baseline contract from R7. R10 `src/vector_grid.py:VectorGrid(side,width,acquisition,gmres_rtol=2e-6,gmres_restart=60,gmres_cycles=20)` retains dyadic FFT operator for even grids; full unknown support. These operators are physically distinct from scalar-TM and must not be interchanged.

## Blockers and practical limits

No blocker for starting endpoint experiments or frozen NN inference. Historical exact task scripts contain remote paths and overwrite protection; new manifests/output paths are required for repeats. R12 lacks raw Fresnel and training NPZ files but sibling R8 originals are accessible. Numerical APIs were inspected, not executed or newly qualified. No claim that older derivative/physical residual checks imply new experiment correctness. Final direction, claim support, novelty and experiment authorization remain parent Codex decisions.
